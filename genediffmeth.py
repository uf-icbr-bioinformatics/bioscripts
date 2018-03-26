#!/usr/bin/env python

import sys
import Script
import GeneList

# TODO: option to output results with sites in separate lines instead of combined

# Utils

REGIONS = [['p', 'promoter'],
           ['b', 'gene body'],
           ['e', 'exons'],
           ['i', 'introns'],
           ['d', 'downstream']]

MODES = {'avg': 'AvgDiffMeth',
         'abs': 'AvgAbsMeth',
         'max': 'MaxDiffMeth',
         'min': 'MinDiffMeth',
         'bal': 'DiffMethBalance',
         'none': 'DiffMeth'}

def decodeRegions(regs):
    res = []
    for r in REGIONS:
        if r[0] in regs:
            res.append(r[1])
    if len(res) == 0:
        return "no regions. I can predict your results will not be very interesting."
    else:
        return ", ".join(res)

def outColName(mode):
    return MODES[mode]

class Region():
    regwanted = ""
    updistance = 2000
    dndistance = 2000

    def __init__(self, reg, up, dn):
        self.regwanted = reg
        self.updistance = up
        self.dndistance = dn

# Produce table of differential methylation by gene

class Params(Script.Script):
    methfile = None
    genesfile = None
    outfile = None
    chromcol = 0                # Column for chromosome
    poscol = 1                  # Column for position
    datacol = 4                 # Column for data
    regions = "b"               # b = gene body, B = gene body + distance, p = promoter, u = upstream, d = downstream
    updistance = 2000           # number of nucleotides upstream of gene
    dndistance = 2000           # number of nucleotides downstream of gene
    mode = 'avg'                # how to combine methyl sites in score for gene. Other options: max, min, bal, abs, none
    minsites = 1                # minimum number of sites
    classify = False            # If true, run in classify mode (all other options are ignored)

    # Runtime
    genelist = []               # All genes
    diffsorted = False          # If true, output will be in order of decreasing differential methylation
    results = []

    def init(self):
        self.genelist = []
        self.results = []

    def parseArgs(self, args):
        next = ""
        self.standardOpts(args)
        for a in args:
            if next == '-r':
                self.regions = a
                next = ""
            elif next == '-d':
                self.updistance = int(a)
                self.dndistance = int(a)
                next = ""
            elif next == '-dup':
                self.updistance = int(a)
                next = ""
            elif next == '-ddn':
                self.dndistance = int(a)
                next = ""
            elif next == '-m':
                self.mode = a
                next = ""
            elif next == '-l':
                self.minsites = int(a)
                next = ""
            elif next == '-x':
                self.datacol = int(a) - 1
                next = ""
            elif a in ['-r', '-d', '-dup', '-ddn', '-m', '-l', '-x']:
                next = a
            elif a == '-s':
                self.diffsorted = True
            elif a == '-c':
                self.classify = True
            elif self.methfile == None:
                self.methfile = a
            elif self.genesfile == None:
                self.genesfile = a
            else:
                self.outfile = a
    
    def describe(self):
        sys.stderr.write("""Input file: {}
Genes file: {}
Output file: {}
Regions: {}
Upstream/downstream: {} / {} bp
Mode: {}
Min sites: {}
""".format(self.methfile, self.genesfile, self.outfile or "standard output", decodeRegions(self.regions), self.updistance, self.dndistance, self.mode, self.minsites))

    def genesOnChrom(self, chrom):
        return self.genelist.genesOnChrom(chrom)

    # def loadGenes(self, filename):
    #     # Autodetect file type
    #     with open(filename, "r") as f:
    #         line = f.readline().rstrip("\r\n")
    #         if line[0:5] == 'LOCUS':    # is it Genbank?
    #             self.genelist = genes.readGenesFromGenbank(filename)
    #         elif line.count("\t") == 10: # is it refFlat?
    #             self.genelist = genes.readGenesFromRefFlat(filename)
    #         elif filename.endswith("gff3") or filename.endswith("gff"):
    #             self.genelist = genes.readGenesFromGFF3(filename)
    #         elif filename.endswith("gtf"):
    #             self.genelist = genes.readGenesFromGTF(filename)
    #         else:
    #             sys.stderr.write("Format of genes file not recognized.")
    #             exit(-1)
    #     return self.genelist

    def main(self):
        genelist = GeneList.loadGenes(self.genesfile, preload=True)

        reg = Region(self.regions, self.updistance, self.dndistance)

        nwritten = 0
        chrcol = self.chromcol
        poscol = self.poscol
        valcol = self.datacol
        thischrom = ""
        chromdata = []
        if self.outfile:
            out = open(self.outfile, "w")
        else:
            out = sys.stdout
        try:
            out.write("# Gene\tChrom\tStart\tEnd\tStrand\tSites\t{}\n".format(outColName(self.mode)))
            with open(self.methfile, "r") as f:
                f.readline()            # skip header
                for line in f:
                    parsed = line.rstrip("\r\n").split("\t")
                    chrom = parsed[chrcol]
                    if chrom != thischrom:
                        if thischrom != "":
                            nwritten += self.writeChromGenes(out, reg, chromdata, genelist.genesOnChrom(thischrom), thischrom)
                        thischrom = chrom
                        # sys.stderr.write("Found chrom {}\n".format(chrom))
                        v = float(parsed[valcol])
                        chromdata = [ (int(parsed[poscol]), v) ]
                    else:
                        v = float(parsed[valcol])
                        chromdata.append((int(parsed[poscol]), v))
            nwritten += self.writeChromGenes(out, reg, chromdata, genelist.genesOnChrom(thischrom), thischrom)
            
            # If we wanted sorted output, nothing has been written so far
            if self.diffsorted:
                self.writeSortedGenes(out)
        finally:
            if self.outfile:
                out.close()
        sys.stderr.write("{} genes written.\n".format(nwritten))

    def writeChromGenes(self, out, reg, data, genes, chrom):
        sys.stderr.write("Writing values for {} genes on {}.\n".format(len(genes), chrom))
        #print data[:100]
        mode = self.mode
        nwritten = 0

        for g in genes:
            grange = g.getRegion(reg)
            if not grange:
                continue
            #    print g.name
            #    print reg.regwanted
            #    raw_input()
            nsites = 0
            totdmc = 0.0
            totabs = 0.0
            maxdmc = 0.0
            mindmc = 0.0
            balance = 0
            datarow = [g.name, g.chrom, grange[0], grange[1], g.strand, 0, 0]
            genesites = []
            for x in data:
                pos = x[0]
                if pos > grange[1]:
                    break            # Too far...
                if pos >= grange[0]: # we're inside
                    if mode == 'none':
                        datarow[2] = pos
                        datarow[3] = pos+1
                        datarow[5] = 1
                        datarow[6] = x[1]
                        out.write("\t".join([str(x) for x in datarow]) + "\n")
                    else:
                        genesites.append(pos)
                        v = x[1]
                        nsites += 1
                        totdmc += v
                        totabs += abs(v)
                        if v > maxdmc:
                            maxdmc = v
                        elif v < mindmc:
                            mindmc = v
                        if v > 0:
                            balance += 1
                        else:
                            balance += -1
            # if nsites > 0:
            #     sys.stderr.write("{}: {} sites\n".format(g.name, nsites))
            #     print datarow
            #     print genesites
            #     print totdmc
            #     raw_input()
            if nsites >= self.minsites:
                if mode == 'avg':
                    score = totdmc / nsites
                elif mode == 'max':
                    score = maxdmc
                elif mode == 'min':
                    score = mindmc
                elif mode == 'bal':
                    score = balance
                elif mode == 'abs':
                    score = totabs / nsites
                datarow[5] = nsites
                datarow[6] = score
                if self.diffsorted:
                    self.results.append(datarow + genesites)
                else:
                    out.write("\t".join([str(x) for x in datarow + genesites]) + "\n")
                nwritten += 1
        return nwritten

    def writeSortedGenes(self, out):
        self.results.sort(key=lambda x: x[6], reverse=True)
        for datarow in self.results:
            out.write("\t".join([str(x) for x in datarow]) + "\n")

### Classify mode

    def findAllClasses(self, genelist, chrom, pos, distance):
        """Find all possible classifications for a site at position `pos' of chromosome `chrom'
using gene database `genelist'."""
        res = []
        names = []
        genes = genelist.allIntersecting(chrom, pos-distance, pos+distance)
        for g in genes:
            names.append(g[0].name)
            c = g[0].classifyPosition(pos, distance)
            if c not in res:
                res.append(c)
        return (res, names)

    def classifyMain(self):
        counters = {'p': 0, 'd': 0, 'E': 0, 'e': 0, 'i': 0, 'n': 0}
        nseen  = 0
        self.datacol = 3        # Assume BEDgraph file

        genelist = self.loadGenes(self.genesfile)
        sys.stderr.write("{} genes loaded from file {}\n".format(genelist.ngenes, self.genesfile))

        if self.outfile:
            out = open(self.outfile, "w")
        else:
            out = sys.stdout
        try:
            out.write("#Chrom\tStart\tEnd\tValue\tClassification\tGenes\n")
            with open(self.methfile, "r") as f:
                f.readline()            # skip header
                for line in f:
                    parsed = line.rstrip("\r\n").split("\t")
                    chrom = parsed[self.chromcol]
                    pos = int(parsed[self.poscol])
                    (cls, names) = self.findAllClasses(genelist, chrom, pos, max(self.updistance, self.dndistance))
                    out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, pos, parsed[self.poscol+1], parsed[self.datacol], "".join(cls), ",".join(names)))
                    if cls == []:
                        nseen += 1
                        counters['n'] += 1
                    else:
                        for c in cls:
                            nseen += 1
                            counters[c] += 1
        finally:
            if self.outfile:
                out.close()
        for c in ['p', 'e', 'E', 'i', 'd', 'n']:
            sys.stderr.write("{}\t{}\t{}%\n".format(c, counters[c], genes.f2dd(100.0*counters[c]/nseen)))

def usage(what=None):
    if what and what.startswith('reg'):
        sys.stderr.write("""This page shows the regions analyzed around a transcript when using
each one of the possible values for -r (b, B, p, u, d). The first example transcript is on the
top strand, the second one on the bottom strand. -dup is set to the equivalent of four characters
and -dn to the equivalent of two characters.

        ------[>>>>>>>>>>>>>>>>>>>>>>>>>>]------
              |bbbbbbbbbbbbbbbbbbbbbbbbbb|
          BBBB|BBBBBBBBBBBBBBBBBBBBBBBBBB|BB
          pppp|                          |
          uuuu|uu                        |
              |                      dddd|dd

        ------[<<<<<<<<<<<<<<<<<<<<<<<<<<]------
              |bbbbbbbbbbbbbbbbbbbbbbbbbb|
            BB|BBBBBBBBBBBBBBBBBBBBBBBBBB|BBBB
              |                          |pppp
              |                        uu|uuuu
          dddd|dd                        |

""")
    else:
        sys.stderr.write("""Usage: genediffmeth [options] diffmeth genesfile [outfile]

Compute summary of differential methylation for the genes in `genesfile' using the
data in file `diffmeth'. This file is expected to have chromosome names in the column 
1, positions in column 2, and differential methylation in column 5 (change with the -x
option). Results are written to `outfile', if provided, or to standard output.

The genes file should be either in Genbank format, in the UCSC Genome Browser 'refFlat'
format, in GTF or GFF3 format.

Options:
  -r region | Limit analysis to sites in the specified gene region. The value should
              be one of: b (gene body, ie from start of transcript to end of transcript),
              B (gene body plus D basepairs upstream and downstream), p (promoter, D
              basepairs upstream of transcript start), u (upstream, ie D basepairs around
              transcript start, d (downstream, ie D basepairs around transcript end).
              Use the -d option to set the distance D. Use "-h regions" to show a help page
              describing the regions produced by the different options.

  -dup D    | Number of nucleotides upstream of the transcript, when using the 'B', 'p', 
              'u' or 'd' regions.

  -ddn D    | Number of nucleotides downstream of the transcript, when using the 'B', 'p', 
              'u' or 'd' regions.

  -d D      | Number of nucleotides upstream or downstream of the transcript, when using
              the 'B', 'p', 'u' or 'd' regions (equivalent to setting -dup and -dn to the
              same value).

  -m mode   | One of 'avg', 'abs', 'max', 'min', 'bal', or 'none'. This option determines
              how the differential methylation values for the sites in a gene are combined 
              to produce the final score. 'avg' (the default) computes the average. 'abs'
              computes the average of the absolute values (measures amount of change 
              irrespective of direction). 'max' and 'min' return the largest positive 
              or negative value, respectively. 'bal' returns the difference between 
              the number of positive sites and the number of negative sites. If 'none',
              then values will NOT be combined, and the output file will contain one line
              for each site.

  -x col    | Specify column containing differential methylation values.

  -l nsites | Limit analysis to regions containing at least `nsites' sites. Default is 1.

  -s        | If specified, output will be sorted by differential methylation (high to
              low) instead of the default order (chromosome, gene position).

  -c        | If specified, switch to 'classify' mode. All other options are ignored except
              for -o and -d. 

""")

if __name__ == "__main__":

    par = Params("genediffmeth.py", version="1.0", usage=usage)
    par.parseArgs(sys.argv[1:])
    if par.genesfile == None:
        usage()
        exit(-1)
    if par.classify:
        par.classifyMain()
    else:
        par.describe()
        par.main()
