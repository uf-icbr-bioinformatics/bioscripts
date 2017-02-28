#!/usr/bin/env python

# (c) 2016, A. Riva, DiBiG, ICBR Bioinformatics, University of Florida

__doc__ = """Reimplementation of RSEM's rsem-generate-data-matrix that also handles ERCC normalization."""

# Functions:
# 1. generate-data-matrix
#    - Read two or more .genes.results (or .isoforms.results) files, create single matrix
#    - If ERCC normalization was requested, do linear regression normalization
#      - Needs to specify which mix is used by which sample
#      - Default ERCC data or loaded from file
#    - Output normalized dataset
#
#    Example cmdline: generate-data-matrix -e -ercc <path> -mix 1,1,2,2 file1 file2 file3 file4
#
# 2. merge multiple Diff files into matrix for heatmap

import sys
import math
import ercc
import os.path
import Utils
import Script

def usage(what=None):
    progname = os.path.basename(sys.argv[0])
    if what == 'matrix':
        sys.stderr.write("""Usage: {} matrix [-e] [-ercc filename] [-mix a,b,c...] files...

Combine multiple .genes.results or .isoforms.results from RSEM into a single data matrix. 
This command is equivalent to rsem-generate-data-matrix, but handles ERCC-based normalization
as well. Options:

 -c name        | Output column 'name' (default: expected_count).
 -e             | Enable ERCC normalization.
 -ercc filename | Load ERCC data from `filename', instead of using defaults.
 -mix a,b,c...  | Specify ERCC mix used by each sample. There should be one
                  value for each sample, either '1' or '2'.

""".format(progname))

    elif what == 'merge':
        sys.stderr.write("""Usage: {} merge labels...

Combine the log2(FC) values for the differentially expressed genes and isoforms for the specified 
labels into two single files. Options:

 -gout filename | Name of output file for genes (default: genes.merged.csv)
 -iout filename | Name of output file for isoforms (default: isoforms.merged.csv)

""".format(progname))

    else:
        sys.stderr.write("""Usage: {} command [args...]

where command is one of:

matrix - combine multiple .genes.results or .isoforms.results into a single data matrix
merge  - merge multiple .gfdr.csv or .ifdr.csv files into a single matrix
allexp - merge multiple .gmatrix.csv files into a single matrix

Call with command and no arguments to get help about a specific command.

""".format(progname))

S = Script.Script("dmaptools.py", version="1.0", usage=usage,
                  errors=[('NOCMD', 'Missing command', 'The first argument should be one of: merge, matrix.')])

# Utils

def ew(string, *args):
    sys.stderr.write(string.format(*args))

def square(a):
    return a * a

def linreg(xs, ys):
    n = len(xs)
    xbar = sum(xs) / n
    ybar = sum(ys) / n
    Lxx = sum([square(xi - xbar) for xi in xs])
    Lxy = sum([(xi - xbar) * (yi - ybar) for (xi, yi) in zip(xs, ys)])
    slope = Lxy / Lxx
    intcp = ybar - (xbar * slope)
    return (slope, intcp)

# Main class

class MatrixGenerator():
    infiles = []                # List of input files
    ncols = 0                   # Number of input files
    rows = []                   # List of rows being assembled
    nrows = 0
    erccrows = []
    nercc = 0
    ngroup1 = 0                 # Number of inputs in group 1 (for normalization)
    ngroup2 = 0                 # Number of inputs in group 2 (for normalization)
    doERCC = False
    ERCCdb = None
    mixes = []
    column = ""
    colidx = 0

    def __init__(self):
        self.rows = []
        self.nrows = 0
        self.erccrows = []
        self.nercc = 0
        self.column = "expected_count"
        self.colidx = 4

    def setFiles(self, infiles):
        self.infiles = infiles
        self.ncols = len(infiles)

    def initERCC(self, mixes=[], ERCCfile=None):
        """This method is called if we want ERCC normalization done before writing the
output matrix. Sets the doERCC field to True, and stores the `mixes' list in the mixes
field. Initializes the ERCCdb to the default list of ERCC controls, unless `ERCCfile' is 
specified, in which case it gets loaded from that file."""
        self.doERCC = True
        if mixes == []:
            self.mixes = ['1'] * self.ncols
        else:
            self.mixes = [ m.strip(" ") for m in mixes.split(",")]
            if len(self.mixes) <> len(self.infiles):
                ew("Error: the number of mixes ({}) should be equal to the number of input files ({}).\n", len(self.mixes), len(self.infiles))
                exit(-1)
            for m in self.mixes:
                if not (m == '1' or m == '2'):
                    ew("Error: mixes can be only '1' or '2', not '{}'.\n",m)
                    exit(-1)
        self.ERCCdb = ercc.ERCCdb()
        if ERCCfile:
            self.ERCCdb.init(ERCCfile)

    def isERCC(self, fields):
        if fields[0].startswith("ERCC-"):
            return fields[0]
        elif fields[1].startswith("ERCC-"):
            return fields[1]
        else:
            return False

    def loadFiles(self):
        """Load the contents of the files listed in the infiles slot into this object. Values
for genes/transcripts and ERCC controls are stored separately."""
        first = self.infiles[0]
        #sys.stderr.write("Reading file {}...\n".format(first))
        ew("Reading file {}...\n", first)
        with open(first, "r") as f:
            hdr = Utils.parseLine(f.readline())
            if self.column in hdr:
                self.colidx = hdr.index(self.column)
                sys.stderr.write("Reading column {} ({})\n".format(self.colidx, self.column))
            else:
                sys.stderr.write("Warning: column {} not found.\n".format(self.column))

            for line in f:
                fields = Utils.parseLine(line)
                isE = self.isERCC(fields)
                if isE:
                    self.erccrows.append([isE, float(fields[self.colidx])] + [0]*(self.ncols-1))
                    self.nercc += 1
                else:
                    self.rows.append([fields[0], float(fields[self.colidx])] + [0]*(self.ncols-1))
                    self.nrows += 1
        idx = 2
        for infile in self.infiles[1:]:
            ew("Reading file {}...\n", infile)
            with open(infile, "r") as f:
                f.readline()        # skip header line
                for i in range(self.nrows):
                    fields = Utils.parseLine(f.readline())
                    row = self.rows[i]
                    if row[0] != fields[0]:
                        ew("ERROR: wrong row order in file {}, expected {}, found {}.", infile, row[0], fields[0])
                        sys.exit(-1)
                    row[idx] = float(fields[self.colidx])
                for i in range(self.nercc):
                    fields = Utils.parseLine(f.readline())
                    row = self.erccrows[i]
                    E = self.isERCC(fields)
                    if row[0] != E:
                        ew("ERROR: wrong row order in file {}, expected {}, found {}.", infile, row[0], E)
                        sys.exit(-1)
                    row[idx] = float(fields[self.colidx])
            idx += 1

    def writeDataMatrix(self):
        """Write the full (normalized) data matrix to stdout."""
        sys.stdout.write("\t" + "\t".join([ '"' + f + '"' for f in self.infiles]) + "\n")
        for r in self.rows:
            sys.stdout.write('"' + r[0] + '"\t' + "\t".join([str(x) for x in r[1:]]) + "\n")

    def ERCCnormalize(self):
        ew("Performing ERCC normalization. Mixes:\n")
        for i in range(self.ncols):
            ew("  {}: mix{}\n", self.infiles[i], self.mixes[i])
        for idx in range(2, self.ncols+1):
            ew("Normalizing column {}:\n", idx)
            ereg = self.ERCCregression(idx)
            if ereg:
                ew("  ERCC regression: {} {}\n", ereg[0], ereg[1])
                self.lnorm(self.rows, idx, ereg[0], ereg[1])

    def regression(self, data, idx):
        """Compute the linear regression between the elements of columns 1 and `idx' in `data'.
Rows in which both values are 0 are ignored."""
        base = [ [r[1], r[idx]] for r in data]
        fx = []
        fy = []
        for b in base:
            if b[0] <> 0 and b[1] <> 0:
                fx.append(b[0])
                fy.append(b[1])
        reg = linreg(fx, fy)
        return reg

    def ERCCregression(self, idx):
        """Compute the linear regression between the elements of columns 1 and `idx' in `data',
containing expression values for ERCC controls. Rows in which both values are 0 are ignored. 
Also handles different concentrations due to the use of different mixes."""        
        base = [ [r[0], r[1], r[idx]] for r in self.erccrows]
        fx = []
        fy = []
        for b in base:
            if b[1] <> 0 and b[2] <> 0:
                e = self.ERCCdb.find(b[0])
                if e:
                    fact = e.factor(self.mixes[0], self.mixes[idx-1])
                    # ew("ERCC: {}, {}, {} => {}\n", b[0], self.mixes[0], self.mixes[idx-1], fact)
                else:
                    fact = 1
                fx.append(b[1])
                fy.append(b[2] * fact)
        if len(fx) > 5:
            reg = linreg(fx, fy)
            return reg
        else:
            ew("No ERCC controls found, unable to normalize.\n")
            return False

    def ERCCregressions(self):
        for idx in range(2, self.ncols+1):
            sys.stderr.write("Columns 1-{}:\n".format(idx))
            regbefore = self.regression(self.rows, idx)
            sys.stderr.write("  Regression before: {} {}\n".format(regbefore[0], regbefore[1]))
            ereg = self.regression(self.erccrows, idx)
            sys.stderr.write("  ERCC regression: {} {}\n".format(ereg[0], ereg[1]))
            self.lnorm(self.rows, idx, ereg[0], ereg[1])
            self.lnorm(self.erccrows, idx, ereg[0], ereg[1])
            eregafter = self.regression(self.erccrows, idx)
            sys.stderr.write("  ERCC regression after: {} {}\n".format(eregafter[0], eregafter[1]))
            regafter = self.regression(self.rows, idx)
            sys.stderr.write("  Regression after: {} {}\n".format(regafter[0], regafter[1]))
            
    def lnorm(self, rows, idx, slope, intercept):
        # sys.stderr.write("Normalizing column {}\n".format(idx))
        for r in rows:
            if r[1] <> 0 and r[idx] <> 0:
                r[idx] = max(0.0, (r[idx] - intercept) / slope) # RSEM doesn't like negative counts... ;)

    def run(self):
        self.loadFiles()
        if self.doERCC:
            self.ERCCnormalize()
        self.writeDataMatrix()

class DiffMerger():
    """Files are either .gdiff.csv or .idiff.csv"""
    labels = []
    ncols = 0
    rows = []
    nrows = 0
    # suffixes = {'gf': ".gfdr.csv",      # differentially expressed genes
    #             'gd': ".gdiff.csv",     # gene fold changes
    #             'if': ".ifdr.csv",      # differentially expressed isoforms
    #             'id': ".idiff.csv"}     # isoform fold changes
    suffixes = {'gf': ".geneDiff.csv",       # differentially expressed genes
                'gd': ".gdiff.csv",          # gene fold changes
                'cf': ".codinggeneDiff.csv", # differentially expressed coding genes
                'cd': ".gdiff.csv",          # gene fold changes
                'if': ".isoDiff.csv",        # differentially expressed isoforms
                'id': ".idiff.csv"}          # isoform fold changes

    gmerged = "merged.geneDiff.csv"
    cmerged = "merged.codinggeneDiff.csv"
    imerged = "merged.isoDiff.csv"

    def __init__(self):
        self.rows = []
        self.nrows = 0

    def setLabels(self, labels):
        self.labels = labels
        self.ncols = len(labels)

    def labelFilename(self, label, key):
        return label + self.suffixes[key]

    def addIdsFromFile(self, filename, ids):
        """Read the identifiers in the first column of `filename' and
add them to set `ids'. Returns the number of identifiers seen."""
        ng = 0
        with open(filename, "r") as f:
            f.readline()
            for line in f:
                parsed = Utils.parseLine(line)
                ids.add(parsed[0].strip('"'))
                ng += 1
        return ng

    def fillMatrixColumn(self, filename, col, ids, matrix, fccol=3, dolog=True):
        """Read column `fccol' from `filename', containing fold changes, 
and add it to the `col'th element of the appropriate row in `matrix' if 
the identifier in the first column is in the set `ids'."""
        with open(filename, "r") as f:
            f.readline()
            for line in f:
                parsed = Utils.parseLine(line)
                g = parsed[0].strip('"')
                if g in ids:
                    row = matrix[g]
                    if dolog:
                        row[col] = math.log(float(parsed[fccol]), 2)
                    else:
                        row[col] = float(parsed[fccol])

    def mergeOne(self, desc, key1, key2, fccol, pcol, outfile):
        wanted = set()

        # Read files with significant items
        for l in self.labels:
            filename = self.labelFilename(l, key1)
            ew("Reading {} from file '{}'... ", desc, filename)
            ng = self.addIdsFromFile(filename, wanted)
            ew("{} {}.\n", ng, desc)
        ew("{} total {} found.\n", len(wanted), desc)

        # Build and initialize matrices
        matrix = {}
        pmatrix = {}
        for g in wanted:
            matrix[g] = [0.0]*self.ncols
            pmatrix[g] = [0.0]*self.ncols

        # Fill matrix
        col = 0
        for l in self.labels:
            filename = self.labelFilename(l, key2)
            ew("Reading fold changes and P-values from file '{}'.\n", filename)
            self.fillMatrixColumn(filename, col, wanted, matrix, fccol=fccol, dolog=True)
            self.fillMatrixColumn(filename, col, wanted, pmatrix, fccol=pcol, dolog=False)
            col += 1

#        if 'ENSG00000167741' in matrix:
#            print matrix['ENSG00000167741']

        # Write output file
        ew("Writing fold changes for {} to file '{}'.\n", desc, outfile)
        with open(outfile, "w") as out:
            # Write header
            out.write("#ID")
            for l in self.labels:
                out.write("\t" + l + "\tp(" + l + ")")
            out.write("\n")

            for g in matrix:
                fvals = matrix[g]
                pvals = pmatrix[g]
                out.write(g)
                for (fv, pv) in zip(fvals, pvals):
                    out.write("\t" + str(fv) + "\t" + str(pv))
                out.write("\n")

    def merge(self):
        self.mergeOne('genes', 'gf', 'gd', 3, 1, self.gmerged)
        self.mergeOne('coding genes', 'cf', 'cd', 3, 1, self.cmerged)
        self.mergeOne('isoforms', 'if', 'id', 3, 1, self.imerged)

    def run(self):
        self.merge()

### ExpMerger

class ExpMerger():
    filenames = []
    labels = []
    table = {}
    outfile = None

    def __init__(self):
        self.filenames = []
        self.labels = []
        self.table = {}

    def readOneFile(self, filename):
        sys.stderr.write("Reading {}... ".format(filename))
        with open(filename, "r") as f:
            hdr = Utils.parseLine(f.readline())
            hdr = [ h.strip('"') for h in hdr ]
            ncols = len(hdr)
            for h in hdr[1:]:
                if h not in self.labels:
                    self.labels.append(h)
            for line in f:
                fields = Utils.parseLine(line)
                g = fields[0]
                if g in self.table:
                    gdata = self.table[g]
                else:
                    gdata = {}
                    self.table[g] = gdata
                for i in range(1, ncols):
                    gdata[hdr[i]] = fields[i]
        sys.stderr.write("done.\n")

    def writeAllExp(self, out):
        out.write("#Gene\t" + "\t".join(self.labels) + "\n")
        for g, gdata in self.table.iteritems():
            out.write(g)
            for l in self.labels:
                out.write("\t")
                out.write(gdata[l])
            out.write("\n")

    def run(self):
        for f in self.filenames:
            self.readOneFile(f)
        sys.stderr.write("{} genes, {} labels.\n".format(len(self.table), len(self.labels)))
        if self.outfile:
            sys.stderr.write("Writing results to {}\n".format(self.outfile))
            with open(self.outfile, "w") as out:
                self.writeAllExp(out)
        else:
            self.writeAllExp(sys.stdout)

def parseArgs(cmd, args):
    if cmd == 'matrix':
        P = MatrixGenerator()
        files = []
        e = False
        ercc = None
        mix = []
        next = ""
        for a in args:
            if a == '-e':
                e = True
            elif next == '-c':
                P.column = a
                next = ""
            elif next == '-ercc':
                ercc = a
                next = ""
            elif next == '-mix':
                mix = a
                next = ""
            elif a in ['-ercc', '-mix', '-c']:
                next = a
            else:
                files.append(a)
        if files == []:
            return usage(cmd)
        P.setFiles(files)
        if e:
            P.initERCC(mix, ERCCfile=ercc)
        return P

    elif cmd == 'merge':
        P = DiffMerger()
        labels = []
        next = ""
        for a in args:
            if next == '-gout':
                P.gmerged = a
                next = ""
            elif next == '-cout':
                P.cmerged = a
                next = ""
            elif next == '-iout':
                P.imerged = a
                next = ""
            elif a in ['-gout', '-cout', '-iout']:
                next = a
            else:
                labels.append(a)
        if labels == []:
            return usage(cmd)
        else:
            P.setLabels(labels)
            return P
    elif cmd == 'allexp':
        P = ExpMerger()
        filenames = []
        next = ""
        for a in args:
            if next == '-o':
                P.outfile = a
                next = ""
            elif a == '-o':
                next = a
            else:
                filenames.append(S.isFile(a))
        P.filenames = filenames
        return P

    else:
        return usage()

if __name__ == "__main__":
    args = sys.argv[1:]
    nargs = len(args)
    if nargs == 0:
        S.errmsg(S.NOCMD)

    S.standardOpts(args)

    cmd = args[0]
    arglist = args[1:]
    P = parseArgs(cmd, arglist)
    if P:
        P.run()

