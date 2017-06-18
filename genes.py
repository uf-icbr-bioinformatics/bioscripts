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

where `command' can be one of: makedb, classify, region, transcripts.

Common options:

  -db D  | Load gene database from sqlite3 file D.
  -gff G | Load gene database from GFF file G.

Use genes.py -h <command> to get help on <command> and its specific options.

""")

    elif what == 'classify':
        sys.stderr.write("""genes.py - Create and query databases of genes and transcripts

Usage: genes.py classify [options] chrom:position ... 

Classify the specified chromosome positions according to the transcripts they fall in.

Options:

  -d N   | Upstream/downstream distance from gene; if specified, sets -dup and -ddn 
           to the same value (default: {}).
  -dup N | Upstream distance from gene (default: {}).
  -ddn N | Downstream distance from gene for classify or region (default: {}).
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

  -r R   | Desired gene region for `region' command; either `b' (gene body),
           `u' (upstream), or `d' (downstream). Default: {}.
  -f F   | Use format F in the region command, either `c' (for coordinates, e.g.
           chr1:1000-1500) or `t' (tab-delimited). Default: {}.
  -d N   | Upstream/downstream distance from gene; if specified, sets -dup and -ddn 
           to the same value (default: {}).
  -dup N | Upstream distance from gene (default: {}).
  -ddn N | Downstream distance from gene for classify or region (default: {}).
  -t     | Show transcripts instead of genes.

""".format(Prog.updistance, Prog.updistance, Prog.dndistance, Prog.regwanted, Prog.oformat))

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
    excel = False

    def parseArgs(self, args):
        cmd = None
        next = ""
        
        self.standardOpts(args)
        for a in args:
            if next == '-gff':
                self.source = P.isFile(a)
                self.sourcetype = 'GFF'
                next = ""
            elif next == '-db':
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
            elif next == '-r':
                if a in ['b', 'u', 'd']:
                    self.regwanted = a
                    next = ""
                else:
                    P.errmsg('BADREGION')
            elif a in ["-gff", "-db", "-d", "-dup", "-ddn", "-f", "-r", "-x"]:
                next = a
            elif a == '-t':
                self.mode = "t"
            elif a == '-s':
                self.mode = "s"
            elif cmd:
                self.args.append(a)
            else:
                cmd = a
        if cmd not in ['region', 'transcripts', 'classify', 'split']:
            P.errmsg('NOCMD')
        return cmd

P = Prog("genes.py", version="1.0", usage=usage, 
         errors=[('BADSRC', 'Missing gene database'),
                 ('BADREGION', 'Bad gene region', "Region should be one of b, u, d."),
                 ('NOCMD', 'Missing command', "Please specify a command (either 'region', 'transcripts', or 'classify').") ])

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

def loadGenes(source, format):
    if format == 'GFF':
        l = GeneList.GFFloader(source)
    elif format == 'DB':
        l = GeneList.DBloader(source)
    else:
        P.errmsg(P.BADSRC)
    sys.stderr.write("Loading genes from {} database {}... ".format(format, source))
    gl = l.load()
    sys.stderr.write("{} genes loaded.\n".format(gl.ngenes))
    return gl

def doRegion():
    if len(P.args) == 0:
        P.errmsg("BADSRC")

    if P.args[0][0] == '@':
        source = Utils.AtFileReader(P.args[0])
    elif P.args[0] == "all":
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
    P.gl = loadGenes(P.source, P.sourcetype)
    try:
        if cmd == 'region':         # Print the region for the genes passed as arguments.
            doRegion()
        elif cmd == 'transcripts':  # Display the transcripts for the genes passed as arguments.
            doTranscripts()
        elif cmd == 'classify':  # Classify the given position (chr:pos) according to the transcripts it falls in
            doClassify()
        elif cmd == 'split':
            doSplit()
    except IOError:
        return

if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) > 0:
        main(args)
    else:
        P.errmsg(P.NOCMD)
