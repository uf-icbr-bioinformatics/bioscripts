#!/usr/bin/env python

import sys
import os.path
import Utils
import Script
import GeneList

def usage(what=None):
    if what == None:
        sys.stderr.write("""genes.py - Create and query databases of genes and transcripts

Usage: genes.py [options] command command-arguments...

where `command' can be one of: makedb, classify, region, transcripts, overlap.

Common options:

  -db D  | Load gene database from database file D. D can be in one of the following formats:
           - gtf (extensions .gtf, .GTF)
           - gff3 (extensions .gff, .gff3, .GFF, .GFF3)
           - genbank (extension .gb)
           - UCSC refFlat (extension .csv)
           - sqlite3 (extension .db)

Use genes.py -h <command> to get help on <command> and its specific options.

""")

    elif what == 'makedb':
        sys.stderr.write("""genes.py - Create and query databases of genes and transcripts

Usage: genes.py classify [options] dbfile

Convert a gene database in gtf/gff/genbank/refFlat format to sqlite3 format.

Options

  -db D | Load genes from file D.

""")

    elif what == 'classify':
        sys.stderr.write("""genes.py - Create and query databases of genes and transcripts

Usage: genes.py classify [options] chrom:position ... 

Classify the specified chromosome positions according to the transcripts they fall in.

Options:

  -d N   | Upstream/downstream distance from gene; if specified, sets -dup and -ddn 
           to the same value (default: {}).
  -dup N | Upstream distance from gene (default: {}).
  -ddn N | Downstream distance from gene (default: {}).
  -t     | List individual transcripts.
  -s     | Summary classification only (number and percentage of region classes).

""".format(Prog.updistance, Prog.updistance, Prog.dndistance))

    elif what == 'region':
        sys.stderr.write("""genes.py - Create and query databases of genes and transcripts

Usage: genes.py region [options] spec ... 

Display the genomic region for the genes (or transcripts, if -t is specified) identified by 
`spec'. If `spec' is the string `all', outputs all genes or transcripts in the database. If a
spec has the form @filename, gene or transcript ids are read from the first column file `filename', 
one per line (a different column can be specified using @filename:col). Otherwise, each `spec' is
treated as a gene or transcript identifier.

Options:

  -r R   | Desired gene region: either `b' (gene body), `u' (upstream), or `d' (downstream). 
           Default: {}.
  -f F   | Output format: either `c' (for coordinates, e.g. chr1:1000-1500) or `t' (tab-delimited).
           Default: {}.
  -d N   | Upstream/downstream distance from gene; if specified, sets -dup and -ddn 
           to the same value (default: {}).
  -dup N | Upstream distance from gene (default: {}).
  -ddn N | Downstream distance from gene (default: {}).
  -t     | Show transcripts instead of genes.

""".format(Prog.regwanted, Prog.oformat, Prog.updistance, Prog.updistance, Prog.dndistance))

    elif what == 'overlap':
        sys.stderr.write("""genes.py - Create and query databases of genes and transcripts

Usage: genes.py overlap [options] bedfile spec...

This command reads regions from BED file `bedfile', and writes to standard output those that 
overlap a region associated with one or more of the transcripts in the current database. The transcript
region is specified using the same options as for the `region' command. If `spec' is the string 
`all', considers all  genes or transcripts in the database. If a spec has the form @filename,
transcript ids are read from the first column of file `filename', one per line (a different 
column can be specified using @filename:col). Otherwise, each `spec' is treated as a gene or 
transcript identifier.

A region from the BED file is written to output if its overlap with at least one of the transcript
regions is larger than the number of bases specified with the 

Options:

  -r R   | Desired gene region: either `b' (gene body), `u' (upstream), or `d' (downstream). 
           Default: {}.
  -d N   | Upstream/downstream distance from gene; if specified, sets -dup and -ddn 
           to the same value (default: {}).
  -dup N | Upstream distance from gene (default: {}).
  -ddn N | Downstream distance from gene (default: {}).
  -b B   | Minimum number of bases in minimum overlap (default: {}).

""".format(Prog.regwanted, Prog.updistance, Prog.updistance, Prog.dndistance, Prog.ovbases))

    elif what == 'transcripts':
        sys.stderr.write("""genes.py - Create and query databases of genes and transcripts

Usage: genes.py transcripts spec ... 

Display the transcripts for the gene identified by the supplied `spec's. If `spec' is the string 
`all', outputs all transcripts in the database. If a spec has the form @filename, gene ids are 
read from the first column file `filename', one per line (a different column can be specified 
using @filename:col). and transcripts for those genes are printed. Otherwise, each `spec' is
treated as a gene identifier, and its transcripts are printed.

""")

### Program object

class Prog(Script.Script):
    source = None
    sourcetype = ""
    args = []
    updistance = 2000
    dndistance = 2000
    oformat = "c"
    regwanted = "b"             # Gene body by default
    gl = None                   # Gene list
    mode = "g"                  # used to be cltrans - also "t" for transcript-level, "s" for summary
    ovbases = 100               # Number of overlap bases
    excel = False

    def parseArgs(self, args):
        cmd = None
        next = ""
        
        self.standardOpts(args)
        for a in args:
            if next == '-db':
                self.source = P.isFile(a)
                self.sourcetype = 'DB'
                next = ""
            elif next == '-d':
                self.updistance = P.toInt(a)
                self.dndistance = self.updistance
                next = ""
            elif next == '-dup':
                self.updistance = P.toInt(a)
                next = ""
            elif next == '-ddn':
                self.dndistance = P.toInt(a)
                next = ""
            elif next == '-f':
                self.oformat = a
                next = ""
            elif next == '-x':
                self.excel = a
                next = ""
            elif next == '-b':
                self.ovbases = P.toInt(a)
                next = ""
            elif next == '-r':
                if a in ['b', 'u', 'd']:
                    self.regwanted = a
                    next = ""
                else:
                    P.errmsg('BADREGION')
            elif a in ["-db", "-d", "-dup", "-ddn", "-f", "-r", "-x", "-b"]:
                next = a
            elif a == '-t':
                self.mode = "t"
            elif a == '-s':
                self.mode = "s"
            elif cmd:
                self.args.append(a)
            else:
                cmd = a
        if cmd not in ['region', 'transcripts', 'classify', 'split', 'makedb', 'overlap']:
            P.errmsg(P.NOCMD)
        return cmd

P = Prog("genes.py", version="1.0", usage=usage, 
         errors=[('BADSRC', 'Missing gene database'),
                 ('BADREGION', 'Bad gene region', "Region should be one of b, u, d."),
                 ('NOOUTDB', 'Missing database filename', "Please specify the name of the output database file."),
                 ('NOCMD', 'Missing command', "Please specify a command (one of 'makedb', 'region', 'transcripts', 'classify', 'split').") ])

# Utils

def parseRegion(c):
    """Parse a region specification of the form chr:start-end and return a
tuple containing the three elements. If -end is omitted, the third element
will be equal to the second one."""
    dc = c.find(":")
    if dc > 0:
        chrom = c[:dc]
        dd = c.find("-")
        if dd > 0:
            return (chrom, int(c[dc+1:dd]), int(c[dd+1:]))
        else:
            pos = int(c[dc+1:])
            return (chrom, pos, pos)
    else:
        return None

def readRegions(atfile):
    """Read regions from a file indicated as @filename:col."""
    regs = []
    afr = Utils.AtFileReader(atfile)
    while True:
        line = afr.stream.readline()
        if line == '':
            afr.stream.close()
            return regs
        if line[0] != afr.ignorechar:
            line = line.rstrip("\r\n").split("\t")
            regs.append( (line[afr.col], int(line[afr.col+1]), int(line[afr.col+2])) )

### Classifier    

class Classifier():
    regnames = [ ('Upstream', 'u'), ('Exon', 'e'), ('CodingExon', 'E'), ('Intron', 'i'), ('Downstream', 'd'), ('Intergenic', 'o') ]
    totals = {'n': 0, 'u': 0, 'd': 0, 'E': 0, 'e': 0, 'i': 0, 'o': 0}
    streams = {'u': None, 'd': None, 'E': None, 'e': None, 'i': None, 'o': None, 's': None}

    def add(self, key):
        self.totals['n'] += 1
        self.totals[key] += 1

    def report(self, stream):
        if stream == 's':
            stream = self.streams['s']
        stream.write("Class\tNumber\tPercentage\n")
        tot = self.totals['n']
        for regs in self.regnames:
            stream.write("{}\t{}\t{:.2f}%\n".format(regs[0], self.totals[regs[1]], 100.0 * self.totals[regs[1]] / tot))
        stream.write("Total\t{}\t\n".format(tot))

    def initStreams(self, filename, header=None, stream=None):
        base = os.path.splitext(filename)[0]
        filename = base + "-Summary.csv"
        if stream:
            if P.excel:
                stream.write("{} -firstrowhdr -name {}-Summary \n".format(filename, P.excel))
            else:
                stream.write(filename + "\n")
        self.streams['s'] = open(filename, "w")
        for rn in self.regnames:
            filename = base + "-" + rn[0] + ".csv"
            if stream:
                if P.excel:
                    stream.write("{} -firstrowhdr -name {}-{} \n".format(filename, P.excel, rn[0]))
                else:
                    stream.write(filename + "\n")
            self.streams[rn[1]] = out = open(filename, "w")
            if header:
                out.write(header + "\n")

    def closeStreams(self):
        for rn in self.regnames:
            self.streams[rn[1]].close()

    def writeLine(self, key, line):
        self.add(key)
        stream = self.streams[key]
        if stream:
            stream.write("\t".join(line) + "\n")

### Main

def loadGenes(filename):
    loaderClass = GeneList.getLoader(filename)
    if not loaderClass:
        P.errmsg(P.BADSRC)
    loader = loaderClass(filename)
    sys.stderr.write("Loading genes database from {}... ".format(filename))
    gl = loader.load(preload=False)
    sys.stderr.write("{} genes loaded.\n".format(gl.ngenes))
    return gl

def doMakeDB():
    if len(P.args) == 0:
        P.errmsg(P.NOOUTDB)
    dbfile = P.args[0]
    sys.stderr.write("Saving gene database to {}...\n".format(dbfile))
    GeneList.initializeDB(dbfile)
    ng = P.gl.saveAllToDB(dbfile)
    sys.stderr.write("done, {} genes written.\n".format(ng))

def doRegion():
    if len(P.args) == 0:
        P.errmsg(P.BADSRC)

    if P.args[0][0] == '@':
        source = Utils.AtFileReader(P.args[0])
    elif P.args[0] == "all":
        if P.mode == "t":
            source = P.gl.allTranscriptNames()
        else:
            source = P.gl.allGeneNames()
    else:
        source = P.args

    for name in source:
        if P.mode == "t":
            tx = P.gl.findTranscript(name)
            if tx:
                reg = tx.getRegion(P)
                if not reg:
                    continue
                (start, end) = reg
                if P.oformat == "c":
                    sys.stdout.write("{}\t{}:{}-{}\t{}\n".format(name, tx.chrom, start, end, "+" if tx.strand == 1 else "-"))
                elif P.oformat == "t":
                    sys.stdout.write("{}\t{}\t{}\t{}\t{}\n".format(tx.chrom, start, end, "+" if tx.strand == 1 else "-", name))
            else:
                sys.stderr.write("No transcript `{}'.\n".format(name))
        else:
            gene = P.gl.findGene(name)
            if gene:
                reg = gene.getRegion(P)
                if not reg:
                    continue
                (start, end) = reg
                if P.oformat == "c":
                    sys.stdout.write("{}\t{}:{}-{}\t{}\n".format(name, gene.chrom, start, end, "+" if gene.strand == 1 else "-"))
                elif P.oformat == "t":
                    sys.stdout.write("{}\t{}\t{}\t{}\t{}\n".format(gene.chrom, start, end, "+" if gene.strand == 1 else "-", name))
            else:
                sys.stderr.write("No gene `{}'.\n".format(name))

def findOverlapping(start, end, regions, minover):
    result = []
    for reg in regions:
        if reg[0] > end:
            break
        if (start <= reg[0] <= end) or (start <= reg[1] <= end):
            overlap = min(end, reg[1]) - max(start, reg[0])
            if overlap >= minover:
                result.append((overlap, reg))
    return result

def doOverlap():
    if len(P.args) < 2:
        P.errmsg()
    bedfile = P.args[0]
    if P.args[1] == "all":
        source = P.gl.allTranscriptNames()
    else:
        source = P.args[1:]

    regions = {}
    sys.stderr.write("Computing regions... ")
    for name in source:
        tx = P.gl.findTranscript(name)
        if tx:
            reg = tx.getRegion(P)
            if not reg:
                continue
            (start, end) = reg
            chrom = tx.chrom
            if chrom not in regions:
                regions[chrom] = []
            regions[chrom].append((start, end, name, tx.strand))
    for chrom in regions.keys():
        regions[chrom].sort()
    sys.stderr.write("done.\n")

    with open(bedfile, "r") as f:
        reader = Utils.CSVreader(f)
        for line in reader:
            chrom = line[0]
            start = int(line[1])
            end   = int(line[2])
            txregions = regions[chrom]
            overs = findOverlapping(start, end, txregions, P.ovbases)
            if overs:
                overs.sort(key=lambda n: n[0])
                for ovr in overs:
                    ov = ovr[1]
                    sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, ov[0], ov[1], ov[2], "+" if ov[3] == 1 else "-", ovr[0]))

def doTranscripts():
    if len(P.args) == 0:
        P.errmsg()
    sys.stdout.write("Gene\tID\tName\tAccession\tChrom\tTXstart\tTXend\tExons\n")
    if P.args[0] == "all":
        for tx in P.gl.getAllTranscripts():
            sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(tx.gene, tx.ID, tx.name, tx.accession, tx.chrom, tx.txstart, tx.txend, ",".join(["{}-{}".format(e[0], e[1]) for e in tx.exons])))

    for name in P.args:
        gene = P.gl.findGene(name)
        if gene:
            for tx in gene.transcripts:
                sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(name, tx.ID, tx.name, tx.accession, tx.chrom, tx.txstart, tx.txend, ",".join(["{}-{}".format(e[0], e[1]) for e in tx.exons])))
        else:
            sys.stderr.write("No gene `{}'.\n".format(name))

def doClassify():
    C = Classifier()
    maxd = max(P.updistance, P.dndistance)

    regions = []
    for a in P.args:
        if a[0] == '@':
            regs = readRegions(a)
            for r in regs:
                regions.append(r)
        else:
            regions.append(a)

    sys.stderr.write("Classifying {} regions.\n".format(len(regions)))

    if P.mode == "t":
        sys.stdout.write("Gene\tID\tAccession\tClass\n")
    elif P.mode == "g":
        sys.stdout.write("Gene\tID\tClass\n")
    
    for reg in regions:
        if type(reg).__name__ == 'str':
            reg = parseRegion(reg)
        pos = reg[1]
        genes = P.gl.allIntersecting(reg[0], pos - maxd, pos + maxd)
        if P.mode == "t":
            for g in genes:
                for tr in g.transcripts:
                    c = tr.classifyPosition(pos, P.updistance, P.dndistance)
                    sys.stdout.write("{}\t{}\t{}\t{}\n".format(g.name, tr.ID, tr.accession, c))
        elif P.mode == "g":
            for g in genes:
                name = g.name
                c = g.classifyPosition(pos, P.updistance, P.dndistance)
                sys.stdout.write("{}\t{}\t{}\n".format(g.name, g.ID, c))
        elif P.mode == "s":
            if len(genes) == 0:
                C.add('o')
            else:
                for g in genes:
                    c = g.classifyPosition(pos, P.updistance, P.dndistance)
                    for x in c:
                        C.add(x)
    if P.mode == "s":
        C.report(sys.stdout)
    
def doSplit():
    C = Classifier()
    maxd = max(P.updistance, P.dndistance)
    for a in P.args:
        if a[0] != '@':
            a = '@' + a
        afr = Utils.AtFileReader(a)
        header = afr.stream.readline().rstrip("\r\n")
        C.initStreams(afr.filename, header=header, stream=sys.stdout)
        while True:
            line = afr.readline()
            if line == '':
                afr.close()
                break
            reg = (line[afr.col], int(line[afr.col+1]), int(line[afr.col+2]))
            pos = reg[1]
            genes = P.gl.allIntersecting(reg[0], pos - maxd, pos + maxd)
            if len(genes) == 0:
                C.writeLine('o', line)
            else:
                for g in genes:
                    c = g.classifyPosition(pos, P.updistance, P.dndistance)
                    for x in c:
                        C.writeLine(x, line)
        C.closeStreams()
        C.report('s')

def main(args):
    cmd = P.parseArgs(args)
    P.gl = loadGenes(P.source)
    try:
        if cmd == 'region':         # Print the region for the genes passed as arguments.
            doRegion()
        elif cmd == 'transcripts':  # Display the transcripts for the genes passed as arguments.
            doTranscripts()
        elif cmd == 'classify':  # Classify the given position (chr:pos) according to the transcripts it falls in
            doClassify()
        elif cmd == 'split':
            doSplit()
        elif cmd == 'makedb':
            doMakeDB()
        elif cmd == 'overlap':
            doOverlap()
    except IOError:
        return

if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) > 0:
        main(args)
    else:
        P.errmsg(P.NOCMD)
