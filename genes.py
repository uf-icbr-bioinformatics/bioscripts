#!/usr/bin/env python

import sys
import csv
import os.path
import Utils
import Script
import GeneList

# Utils

def loadGenes(filename):
    loaderClass = GeneList.getLoader(filename)
    if not loaderClass:
        P.errmsg(P.BADSRC)
    loader = loaderClass(filename)
    sys.stderr.write("[Loading genes database from {}... ".format(filename))
    gl = loader.load(preload=False)
    sys.stderr.write("{} genes loaded.]\n".format(gl.ngenes))
    return gl

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
    """Read regions from a file indicated as @filename:col. Lines that don't have
the required number of columns or that contains start and end positions that are
not numbers are silently ignored."""
    regs = []
    afr = Utils.AtFileReader(atfile)
    while True:
        line = afr.stream.readline()
        if not line:
            afr.stream.close()
            return regs
        if line[0] != afr.ignorechar:
            line = line.rstrip("\r\n").split("\t")
            try:
                regs.append( (line[afr.col], int(line[afr.col+1]), int(line[afr.col+2])) )
            except ValueError:
                pass
            except IndexError:
                pass

def findOverlapping(start, end, regions, minover):
    """Return all itntervals from the list `regions' that overlap the interval `start'-`end'
by at least `minover'."""
    
    result = []
    for reg in regions:
        if reg[0] > end:
            break
        if (start <= reg[0] <= end) or (start <= reg[1] <= end) or (reg[0] <= start <= end <= reg[1]):
            overlap = min(end, reg[1]) - max(start, reg[0])
            if overlap >= minover:
                result.append((overlap, reg))
    return result

### Command classes

class MakeDB(Script.Command):
    _cmd = "makedb"

    def usage(self, P):
        sys.stderr.write("""Usage: genes.py makedb [options] dbfile

Convert a gene database in gtf/gff/genbank/refFlat format to sqlite3 format. 
The -db option specifies the input database, the -o option specifies the output 
database. Both are required.

""")

    def run(self, P):
        if len(P.args) == 0:
            P.errmsg(P.NOOUTDB)
        dbfile = P.args[0]
        sys.stderr.write("Saving gene database to {}...\n".format(dbfile))
        GeneList.initializeDB(dbfile)
        ng = P.gl.saveAllToDB(dbfile)
        sys.stderr.write("done, {} genes written.\n".format(ng))

class Classify(Script.Command):
    _cmd = "classify"
    regnames = [ ('Upstream', 'u'), ('Exon', 'e'), ('CodingExon', 'E'), ('Intron', 'i'), ('Downstream', 'd'), ('Intergenic', 'o') ]
    totals = {}
    streams = {}
    
    def __init__(self):
        self.totals = {'n': 0, 'u': 0, 'd': 0, 'E': 0, 'e': 0, 'i': 0, 'o': 0}
        self.streams = {'u': None, 'd': None, 'E': None, 'e': None, 'i': None, 'o': None, 's': None}

    def add(self, key):
        self.totals['n'] += 1
        self.totals[key] += 1

    def report(self, stream):
        if stream == 's':
            stream = self.streams['s']
        stream.write("Class\tNumber\tPercentage\n")
        tot = self.totals['n']
        for regs in self.regnames:
            if tot > 0:
                pct = 100.0 * self.totals[regs[1]] / tot
            else:
                pct = 0
            stream.write("{}\t{}\t{:.2f}%\n".format(regs[0], self.totals[regs[1]], pct))
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

    def usage(self, P):
        sys.stderr.write("""Usage: genes.py classify [options] chrom:position ... 

Classify the specified chromosome position(s) according to the transcripts they fall in. If
an argument has the form @filename:col, read regions from column `col' of file `filename'
(col defaults to 1).

Options:

  -d N   | Upstream/downstream distance from gene; if specified, sets -dup and -ddn 
           to the same value (default: {}).
  -dup N | Upstream distance from gene (default: {}).
  -ddn N | Downstream distance from gene (default: {}).
  -a     | Preserve original contents of input file.
  -t     | List individual transcripts.
  -s     | Summary classification only (number and percentage of region classes).

""".format(P.updistance, P.updistance, P.dndistance))
    
    def run(self, P):
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

        with Utils.Output(P.outfile) as out:
            if P.mode == "t":
                out.write("Gene\tID\tAccession\tClass\n")
            elif P.mode == "g":
                out.write("Gene\tID\tClass\n")

            for reg in regions:
                if type(reg).__name__ == 'str':
                    reg = parseRegion(reg)
                start = reg[1]
                end   = reg[2]
                pos   = (start + end) / 2
                genes = P.gl.allIntersecting(reg[0], start - maxd, end + maxd)
                if P.mode == "t":
                    for g in genes:
                        for tr in g.transcripts:
                            c = tr.classifyPosition(pos, P.updistance, P.dndistance)
                            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(reg[0], start, end, g.name, tr.ID, tr.accession, c))
                elif P.mode == "g":
                    for g in genes:
                        name = g.name
                        c = g.classifyPosition(pos, P.updistance, P.dndistance)
                        out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(reg[0], start, end, g.name, g.ID, c))
                elif P.mode == "s":
                    if len(genes) == 0:
                        self.add('o')
                    else:
                        for g in genes:
                            c = g.classifyPosition(pos, P.updistance, P.dndistance)
                            for x in c:
                                self.add(x)
                elif P.mode == "st":
                    if len(genes) == 0:
                        self.add('o')
                    else:
                        for g in genes:
                            for tr in g.transcripts:
                                c = tr.classifyPosition(pos, P.updistance, P.dndistance)
                                if c:
                                    for x in c:
                                        self.add(x)
                                else:
                                    self.add('o')

            if "s" in P.mode:
                self.report(out)

class Region(Script.Command):
    _cmd = "region"

    def usage(self, P):
        sys.stderr.write("""Usage: genes.py region [options] spec ... 

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

""".format(P.regwanted, P.oformat, P.updistance, P.updistance, P.dndistance))

    def run(self, P):
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

class Transcripts(Script.Command):
    _cmd = "transcripts"

    def usage(self, P):
        sys.stderr.write("""Usage: genes.py transcripts spec ... 

Display the transcripts for the gene identified by the supplied `spec's. If `spec' is the string 
`all', outputs all transcripts in the database. If a spec has the form @filename, gene ids are 
read from the first column file `filename', one per line (a different column can be specified 
using @filename:col). and transcripts for those genes are printed. Otherwise, each `spec' is
treated as a gene identifier, and its transcripts are printed.

""")
        
    def run(self, P):
        if len(P.args) == 0:
            P.errmsg(P.NOSPECS)
        sys.stdout.write("Gene\tID\tName\tAccession\tChrom\tTXstart\tTXend\tExons\n")
        if P.args[0] == "all":
            for tx in P.gl.getAllTranscripts():
                sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(tx.gene, tx.ID, tx.name, tx.accession, tx.chrom, tx.txstart, tx.txend, ",".join(["{}-{}".format(e[0], e[1]) for e in tx.exons])))
        else:
            for name in P.args:
                gene = P.gl.findGene(name)
                if gene:
                    for tx in gene.transcripts:
                        sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(name, tx.ID, tx.name, tx.accession, tx.chrom, tx.txstart, tx.txend, ",".join(["{}-{}".format(e[0], e[1]) for e in tx.exons])))
                else:
                    sys.stderr.write("No gene `{}'.\n".format(name))

class Overlap(Script.Command):
    _cmd = "overlap"

    def usage(self, P):
        sys.stderr.write("""Usage: genes.py overlap [options] bedfile spec...

This command reads regions from BED file `bedfile', and writes to standard output those that 
overlap a region associated with one or more of the genes (or transcripts if -t is specified)
in the current database. The transcript region is specified using the same options as for the 
`region' command. If `spec' is the string `all', considers all genes or transcripts in the 
database. If a spec has the form @filename, gene or transcript ids are read from the first 
column of file `filename', one per line (a different column can be specified using @filename:col). 
Otherwise, each `spec' is treated as a gene or transcript identifier.

A region from the BED file is written to output if its overlap with at least one of the
regions is larger than the number of bases specified with the -b option.

Options:

  -r R   | Desired gene region: either `b' (gene body), `u' (upstream), or `d' (downstream). 
           Default: {}.
  -d N   | Upstream/downstream distance from gene; if specified, sets -dup and -ddn 
           to the same value (default: {}).
  -dup N | Upstream distance from gene (default: {}).
  -ddn N | Downstream distance from gene (default: {}).
  -t     | Show transcripts instead of genes.
  -a     | Preserve original contents of input file.
  -b B   | Minimum number of bases in minimum overlap (default: {}).

""".format(P.regwanted, P.updistance, P.updistance, P.dndistance, P.ovbases))



    def run(self, P):
        if len(P.args) < 2:
            P.errmsg(P.NOFILE)
        bedfile = P.args[0]
        if P.args[1] == "all":
            if P.mode == "t":
                source = P.gl.allTranscriptNames()
            else:
                source = P.gl.allGeneNames()
        else:
            source = P.args[1:]

        regions = {}
        sys.stderr.write("Computing regions... ")
        for name in source:
            if P.mode == "t":
                tx = P.gl.findTranscript(name)
            else:
                tx = P.gl.findGene(name)
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
                start = Utils.safeInt(line[1])
                end   = Utils.safeInt(line[2])
                if start is None or end is None:
                    continue
                if chrom not in regions:
                    continue
                txregions = regions[chrom]
                overs = findOverlapping(start, end, txregions, P.ovbases)
                if overs:
                    overs.sort(key=lambda n: n[0])
                    for ovr in overs:
                        ov = ovr[1]
                        sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}".format(chrom, ov[0], ov[1], ov[2], "+" if ov[3] == 1 else "-", ovr[0]))
                        if P.addBedFields:
                            sys.stdout.write("\t" + "\t".join(line[1:]))
                        sys.stdout.write("\n")

class Split(Script.Command):
    _cmd = "split"

    def usage(self, P):
        sys.stderr.write("""Usage: genes.py [options] split filename...

Separate regions contained in the input files into multiple files depending on their
classification. For each input file, six output files are created concatenating the
filename with each of following labels:

{}

Each row of input file is interpreted as a region (chrom, start, end), classified as 
per the `classify' command, and written to the appropriate output file.

""".format(", ".join([w[0] for w in Classify.regnames])))
    
    def run(self, P):
        C = Classify()
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

class Closest(Script.Command):
    _cmd = "closest"

    def usage(self, P):
        sys.stderr.write("""Usage: genes.py closest [options] spec...

Find the closest gene (or transcript, if -t is specified) to each one of the 
regions specified on the command line. Each `spec' can be a position (chrom:pos),
a region (chrom:start-end), or a file specification of the form @filename:col, 
in which case regions are read from filename, assuming that the chromosome is in 
column `col' (defaulting to 1) and start and end positions are in the two next columns.

If a region was read from the command-line, output consists of the following 
tab-separated fields:

        chrom, start, end, distance from nearest gene, gene ID, name of gene

If regions are read from a file, output consist of each line of the input
file followed by distance from nearest gene, gene ID, name of gene.

Output is written to standard ouptut, unless an output file is specified with -o.
""")
    
    def run(self, P):
        regions = []
        transcript = "t" in P.mode
        biotype = "protein_coding" if P.codingOnly else None

        with Utils.Output(P.outfile) as out:
            nin = 0
            nout = 0

            with P.gl:
                for a in P.args:
                    if a[0] == '@':
                        (nin, nout) = self.runFile(a[1:], out, transcripts=transcript, biotype=biotype, canonical=P.canonical)
                    else:
                        regions.append(parseRegion(a))

                for reg in regions: # Does nothing if we're in @ mode
                    nin += 1
                    (geneID, dist) = P.gl.findClosestGene(reg[0], reg[1], reg[2], transcripts=transcript, biotype=biotype, canonical=P.canonical)
                    if geneID:
                        nout += 1
                        gene = P.gl.findGene(geneID)
                        out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(reg[0], reg[1], reg[2], dist, gene.ID, gene.name))
        sys.stderr.write("Classified {}/{} regions.\n".format(nin, nout))

    def runFile(self, infile, out, transcripts=False, biotype=None, canonical=False):
        if not os.path.isfile(infile):
            P.errmsg(P.NOFILE)
        nin = 0
        nout = 0
        with open(infile, "r") as f:
            c = csv.reader(f, delimiter='\t')
            for line in c:
                nin += 1
                (ID, dist) = P.gl.findClosestGene(line[0], int(line[1]), int(line[2]), transcripts=transcripts, biotype=biotype, canonical=canonical)
                if ID:
                    nout += 1
                    if transcripts:
                        tx = P.gl.findTranscript(ID)
                        outrow = line + [str(dist), ID, tx.name]
                    else:
                        gene = P.gl.findGene(ID)
                        outrow = line + [str(dist), gene.ID, gene.name]
                    out.write("\t".join(outrow) + "\n")
        return (nin, nout)

### Annotate

class Annotate(Script.Command):
    _cmd = "annotate"

    def usage(self, P):
        sys.stderr.write("""Usage: genes.py annotate [options] infile

Rewrite input file `infile' adding gene annotations. Options:

  -i I | Gene identifiers are in this column (default: 1).
  -w W | Comma-separated list of wanted fields from db (default: 'name').

""")

    def run(self, P):
        if len(P.args) == 0:
            P.errmsg(P.NOFILE)
        infile = P.args[0]
        idcol = P.idcol
        inscol = idcol + 1
        fieldnames = [ s.upper() for s in P.wanted ]
        missing = [ "???" for s in P.wanted ]
        query = "SELECT {} FROM Genes WHERE ID=?".format(",".join(P.wanted))
        with Utils.Output(P.outfile) as out:
            with P.gl:
                with open(infile, "r") as f:
                    c = csv.reader(f, delimiter='\t')
                    hdr = c.next()
                    hdr[inscol:inscol] = fieldnames
                    out.write("\t".join(hdr) + "\n")
                    for row in c:
                        gid = row[P.idcol].strip('"')
                        data = P.gl.getGeneInfo(gid, query)
                        if data:
                            row[inscol:inscol] = data
                        else:
                            row[inscol:inscol] = missing
                        out.write("\t".join(row) + "\n")

### Main

def usage(what=None):
    if what == None:
        sys.stderr.write("""genes.py - Create and query databases of genes and transcripts

Usage: genes.py [common-options] [options] command command-arguments...

where `command' can be one of: {}.

Common options (applicable to all commands):

  -o  O | Write output to file O (default: stdout)
  -db D | Load gene database from database file D. D can be in one of the following formats:
          - gtf (extensions .gtf, .GTF)
          - gff3 (extensions .gff, .gff3, .GFF, .GFF3)
          - genbank (extension .gb)
          - UCSC refFlat (extension .csv)
          - sqlite3 (extension .db)

Use genes.py -h <command> to get help on <command> and its specific options.

""".format(", ".join(P._commandNames)))

    else:
        c = P.findCommand(what)
        if c:
            c().usage(P)
        
### Program object

class Prog(Script.Script):
    outfile = None
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
    addBedFields = False        # If True, add contents of original BED file
    excel = False
    idcol = 0                   # Column containing gene id for Annotate
    wanted = ['name']           # Wanted field(s) from db for Annotate
    codingOnly = False          # If True (-pc) only look at protein coding genes in Closest
    canonical = False           # If True (-ca) only look at canonical transcript for each gene in Closest

    def parseArgs(self, args):
        cmd = None
        next = ""
        
        self.standardOpts(args)
        for a in args:
            if next == '-db':
                self.source = self.isFile(a)
                self.sourcetype = 'DB'
                next = ""
            elif next == '-o':
                self.outfile = a
                next = ""
            elif next == '-d':
                self.updistance = self.toInt(a)
                self.dndistance = self.updistance
                next = ""
            elif next == '-dup':
                self.updistance = self.toInt(a)
                next = ""
            elif next == '-ddn':
                self.dndistance = self.toInt(a)
                next = ""
            elif next == '-f':
                self.oformat = a
                next = ""
            elif next == '-x':
                self.excel = a
                next = ""
            elif next == '-b':
                self.ovbases = self.toInt(a)
                next = ""
            elif next == '-r':
                if a in ['b', 'u', 'd']:
                    self.regwanted = a
                    next = ""
                else:
                    self.errmsg('BADREGION')
            elif next == "-c":
                self.idcol = self.toInt(a) - 1
                next = ""
            elif next == "-w":
                self.wanted = a.split(",")
                next = ""
            elif a in ["-db", "-d", "-dup", "-ddn", "-f", "-r", "-x", "-b", "-o", "-c", "-w"]:
                next = a
            elif a == '-t':
                if self.mode == "s":
                    self.mode = "st"
                else:
                    self.mode = "t"
            elif a == '-s':
                if self.mode == "t":
                    self.mode = "st"
                else:
                    self.mode = "s"
            elif a == "-a":
                self.addBedFields = True
            elif a == "-pc":
                self.codingOnly = True
            elif a == "-ca":
                self.canonical = True
            elif cmd:
                self.args.append(a)
            else:
                cmd = a
        if cmd not in self._commandNames:
            P.errmsg(self.NOCMD)
        if not self.source:
            P.errmsg(self.BADSRC)
        return cmd

P = Prog("genes.py", version="1.0", usage=usage, 
         errors=[('BADSRC', 'Missing gene database'),
                 ('NOFILE', 'The input file does not exist.'),
                 ('BADREGION', 'Bad gene region', "Region should be one of b, u, d."),
                 ('NOOUTDB', 'Missing output database filename', "Please specify the name of the output database file."),
                 ('NOCMD', 'Missing command', "Please specify a command."),
                 ('NOSPECS', 'Missing specs', "Please provide at least one gene or region spec.") ])
P.addCommand(MakeDB)
P.addCommand(Classify)
P.addCommand(Region)
P.addCommand(Transcripts)
P.addCommand(Overlap)
P.addCommand(Split)
P.addCommand(Closest)
P.addCommand(Annotate)

### Main

def main(args):
    cmd = P.parseArgs(args)
    P.gl = loadGenes(P.source)
    try:
        c = P.findCommand(cmd)
        if c:
            c().run(P)
        else:
            usage(P)
    except IOError:
        return

if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) > 0:
        main(args)
    else:
        P.errmsg(P.NOCMD)
