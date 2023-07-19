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
import csv
import math
import ercc
import os.path
import scipy.stats
import Utils
import Script

def usage(what=None):


    progname = os.path.basename(sys.argv[0])
    if what == 'matrix':
        sys.stdout.write("""Usage: {} matrix [options] files...

Combine columns from multiple files into a single data matrix. This command is equivalent to 
rsem-generate-data-matrix, but handles ERCC-based normalization as well. Options:

 -c name        | Output column 'name' (default: expected_count).
 -l A,B,C...    | Use supplied labels for columns instead of filenames.
 -e             | Enable ERCC normalization.
 -ercc filename | Load ERCC data from `filename', instead of using defaults.
 -mix a,b,c...  | Specify ERCC mix used by each sample. There should be one
                  value for each sample, either '1' or '2'.
 -f F           | Omit genes with intra-condition variability higher than F (in log2 FC).
                  Requires the -n option.
 -u             | Invert the meaning of -F: output only the high-variability genes.
 -n X,Y,Z,...   | Specifies the number of samples for each condition. Eg:
                  -n 2,3 = samples 1-2 from condition A, samples 3-5 from condition B.
 -r             | Round values to closest integer before writing matrix.

""".format(progname))

    elif what == 'merge':
        sys.stdout.write("""Usage: {} merge labels...

Combine the log2(FC) values for the differentially expressed genes and isoforms for the specified 
labels into two single files. Options:

 -gout filename | Name of output file for genes (default: genes.merged.csv)
 -iout filename | Name of output file for isoforms (default: isoforms.merged.csv)

""".format(progname))

    elif what == 'sort':
        sys.stdout.write("""Usage: {} sort [options] filename

Sort items from `filename' according to the average value(s) in one or
more specified columns. Options:

  -o F          | Write output to file F (default: stdout).
  -i I          | Identifiers are in column I (default: {}).
  -c C1,C2, ... | Values are in columns C1, C2, etc.
  -t T          | Output only the top T % of genes (default: {}).

""".format(progname, GeneSorter.idcol + 1, GeneSorter.topperc))
    elif what == 'extract':
        sys.stdour.write("""Usage: {} extract [options] datafile filenames...

Extract rows from `datafile' corresponding to all IDs found in `filenames' (removing duplicates)."

Options:

  -o O | Write output to file O (default: standard output).

""".format(progname))

    else:
        sys.stdout.write("""Usage: {} command [args...]

where command is one of:

matrix - combine multiple .genes.results or .isoforms.results into a single data matrix
merge  - merge multiple .gfdr.csv or .ifdr.csv files into a single matrix
allexp - merge multiple .gmatrix.csv files into a single matrix
sort   - sort genes from input file according to values in specified columns

Call with command and no arguments to get help about a specific command.

""".format(progname))

S = Script.Script("rnaseqtools.py", version="1.0", usage=usage,
                  errors=[('NOCMD', 'Missing command', 'The first argument should be one of: merge, matrix, allexp, extract, sort.')])

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
    labels = []                 # Labels for data columns (if not specified, uses filename)
    conditions = []             # Condition that each sample belongs to
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
    maxfc = None
    minval = False
    condsmpls = []
    filtover = True             # Filter genes over maxfc threshold?
    roundvals = False           # Round values to closest integer?

    def __init__(self):
        self.labels = []
        self.conditions = []
        self.rows = []
        self.nrows = 0
        self.erccrows = []
        self.nercc = 0
        self.column = "expected_count"
        self.colidx = 4
        self.condsmpls = []

    def parseArgs(self, args):
        files = []
        e = False
        ercc = None
        mix = []
        prev = ""
        condsmpls = None

        for a in args:
            if a == '-e':
                e = True
            elif prev == '-c':
                self.column = a
                prev = ""
            elif prev == '-ercc':
                ercc = a
                prev = ""
            elif prev == '-mix':
                mix = a
                prev = ""
            elif prev == '-f':
                self.maxfc = float(a)
                prev = ""
            elif prev == '-n':
                condsmpls = [ int(s) for s in a.split(",") ]
                prev = ""
            elif prev == '-m':
                self.minval = float(a)
                prev = ""
            elif prev == '-l':
                self.labels = a.split(",")
                prev = ""
            elif prev == '-g':
                self.conditions = a.split(",")
                prev = ""
            elif a in ['-ercc', '-mix', '-c', '-n', '-f', '-m', '-l', '-g']:
                prev = a
            elif a == "-u":
                self.filtover = False
            elif a == "-r":
                self.roundvals = True
            else:
                files.append(S.isFile(a))
        if files == []:
            S.usage(cmd)
        self.infiles = files
        self.ncols = len(files)
        if condsmpls:
            self.condsmpls = self.condIndexes(condsmpls)
        if e:
            self.initERCC(mix, ERCCfile=ercc)

    def condIndexes(self, condsmpls):
        result = []
        idx = 1
        for cs in condsmpls:
            r = []
            for i in range(cs):
                r.append(idx)
                idx += 1
            result.append(r)
        return result

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
            if len(self.mixes) != len(self.infiles):
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
                        ew("ERROR: wrong row order in file {}, expected {}, found {}.\n", infile, row[0], fields[0])
                        sys.exit(-1)
                    row[idx] = float(fields[self.colidx])
                for i in range(self.nercc):
                    fields = Utils.parseLine(f.readline())
                    row = self.erccrows[i]
                    E = self.isERCC(fields)
                    if row[0] != E:
                        ew("ERROR: wrong row order in file {}, expected {}, found {}.\n", infile, row[0], E)
                        sys.exit(-1)
                    row[idx] = float(fields[self.colidx])
            idx += 1

    def rowMaxFC(self, row):
        """Determine the highest log2(FC) between samples of the same condition."""
        maxfc = 0
        for idxs in self.condsmpls:
            vmin = row[idxs[0]]
            vmax = row[idxs[0]]
            for idx in idxs[1:]:
                x = row[idx]
                if x < vmin:
                    vmin = x
                if x > vmax:
                    vmax = x
            if vmin > 0:
                fc = math.log(vmax/vmin, 2)
                if fc > maxfc:
                    maxfc = fc
        return maxfc

    def rowMinAvg(self, row):
        """Determine the smallest average between samples of the same condition."""
        minavg = 100000000
        for idxs in self.condsmpls:
            s = 0
            for i in idxs:
                s += row[i]
            a = s / len(idxs)
            if a < minavg:
                minavg = a
        return minavg

    def writeDataMatrix(self):
        """Write the full (normalized) data matrix to stdout."""
        nfiltered1 = 0
        nfiltered2 = 0

        if self.labels:
            if self.conditions:
                sys.stdout.write("#NAME\t" + "\t".join(self.labels) + "\n")
            else:
                sys.stdout.write("\t" + "\t".join(self.labels) + "\n")
        else:
            sys.stdout.write("\t" + "\t".join([ '"' + f + '"' for f in self.infiles]) + "\n")
        if self.conditions:
            sys.stdout.write("#CLASS:Condition\t" + "\t".join(self.conditions) + "\n")
        for r in self.rows:
            if self.minval:
                #if min(r[1:]) < self.minval:
                if self.rowMinAvg(r) < self.minval:
                    nfiltered1 += 1
                    continue
            if self.maxfc:
                over = (self.rowMaxFC(r) > self.maxfc)
                if over == self.filtover:
                    nfiltered += 1
                    continue
            if self.roundvals:
                sys.stdout.write('"' + r[0] + '"\t' + "\t".join([str(int(round(x))) for x in r[1:]]) + "\n")
            else:
                sys.stdout.write('"' + r[0] + '"\t' + "\t".join([str(x) for x in r[1:]]) + "\n")
        if nfiltered1 > 0:
            ew("{} genes filtered because value < {}.\n", nfiltered1, self.minval)
        if nfiltered2 > 0:
            ew("{} genes filtered because variability > {}.\n", nfiltered2, self.maxfc)
#        for r in self.erccrows:
#            sys.stdout.write('"' + r[0] + '"\t' + "\t".join([str(x) for x in r[1:]]) + "\n")

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
            if b[0] != 0 and b[1] != 0:
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
            if b[1] != 0 and b[2] != 0:
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
            if r[1] != 0 and r[idx] != 0:
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
    pval = 1
    log2fc = 0
    # suffixes = {'gf': ".gfdr.csv",      # differentially expressed genes
    #             'gd': ".gdiff.csv",     # gene fold changes
    #             'if': ".ifdr.csv",      # differentially expressed isoforms
    #             'id': ".idiff.csv"}     # isoform fold changes
    suffixes = {'gf': ".geneDiff.csv",       # differentially expressed genes
                'gd': ".gdiff.csv",          # gene fold changes
                'cf': ".codinggeneDiff.csv", # differentially expressed coding genes
                'cd': ".gdiff.csv",          # coding gene fold changes
                'if': ".isoDiff.csv",        # differentially expressed isoforms
                'id': ".idiff.csv"}          # isoform fold changes

    gmerged = "merged.geneDiff.csv"
    cmerged = "merged.codinggeneDiff.csv"
    imerged = "merged.isoDiff.csv"

    def __init__(self):
        self.rows = []
        self.nrows = 0

    def parseArgs(self, args):
        labels = []
        prev = ""
        for a in args:
            if prev == '-gout':
                self.gmerged = a
                prev = ""
            elif prev == '-cout':
                self.cmerged = a
                prev = ""
            elif prev == '-iout':
                self.imerged = a
                prev = ""
            elif prev == "-p":
                self.pval = float(a)
                prev = ""
            elif prev == "-f":
                self.log2fc = float(a)
                prev = ""
            elif a in ['-gout', '-cout', '-iout', "-p", "-f"]:
                prev = a
            else:
                labels.append(a)
        if labels == []:
            S.usage(cmd)
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
                        row[col] = parsed[fccol]

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
            ew("Reading fold changes (col {}) and P-values (col {}) from file '{}'.\n", fccol, pcol, filename)
            self.fillMatrixColumn(filename, col, wanted, matrix, fccol=fccol, dolog=False)
            self.fillMatrixColumn(filename, col, wanted, pmatrix, fccol=pcol, dolog=False)
            col += 1

#        if 'ENSG00000167741' in matrix:
#            print(matrix['ENSG00000167741'])

        # Write output file
        ew("Writing fold changes for {} to file '{}'.\n", desc, outfile)
        with open(outfile, "w") as out:
            # Write header
            out.write("#ID")
            for l in self.labels:
                out.write("\t" + l + "\tp(" + l + ")")
                #out.write("\t" + l)
            out.write("\n")

            for g in matrix:
                fvals = matrix[g]
                pvals = pmatrix[g]
                out.write(g)
                for (fv, pv) in zip(fvals, pvals):
                    if pv == "NA":
                        out.write("\t\t")
                    elif abs(float(fv)) >= self.log2fc and float(pv) <= self.pval:
                        out.write("\t" + str(fv) + "\t" + pv)
                    else:
                        out.write("\t\t")
                #out.write("\t".join([str(fv) for fv in fvals]))
                out.write("\n")

    def merge(self):
        self.mergeOne('genes', 'gf', 'gd', 2, 6, self.gmerged)
        self.mergeOne('coding genes', 'cf', 'cd', 2, 6, self.cmerged)
        self.mergeOne('isoforms', 'if', 'id', 2, 6, self.imerged)

    def run(self):
        self.merge()

### ExtractRows

class ExtractRows():
    datafile = ""
    wantedfiles = []
    outfile = "/dev/stdout"

    def __init__(self):
        self.wantedfiles = []

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-o":
                self.outfile = a
                prev = ""
            elif a in ["-o"]:
                prev = a
            elif self.datafile:
                self.wantedfiles.append(a)
            else:
                self.datafile = a
        if not self.wantedfiles:
            S.usage("extract")

    def run(self):
        wanted = []
        for wf in self.wantedfiles:
            with open(wf, "r") as f:
                for line in Utils.CSVreader(f):
                    gid = line[0]
                    if gid not in wanted:
                        wanted.append(gid)
        sys.stderr.write("{} genes to extract.\n".format(len(wanted)))
        with open(self.outfile, "w") as out, open(self.datafile, "r") as f:
            out.write(f.readline()) # Copy header (*** make this an option)
            for line in Utils.CSVreader(f):
                if line[0] in wanted:
                    out.write("\t".join(line) + "\n")

### ExpMerger

class ExpMerger():
    filenames = []
    labels = []
    detectLabels = True
    table = {}
    outfile = None

    def __init__(self):
        self.filenames = []
        self.labels = []
        self.table = {}

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == '-o':
                self.outfile = a
                prev = ""
            elif prev == '-l':
                self.labels = a.split(",")
                self.detectLabels = False
                prev = ""
            elif a in ['-o', '-l']:
                prev = a
            else:
                self.filenames.append(S.isFile(a))
        if self.filenames == []:
            S.usage("merge")

    def readOneFile(self, filename):
        sys.stderr.write("Reading {}... ".format(filename))
        with open(filename, "r") as f:
            c = csv.reader(f, delimiter='\t')
            hdr = next(c)
            hdr = [ h.strip('"') for h in hdr ]
            ncols = len(hdr)
            if self.detectLabels:
                for h in hdr[1:]:
                    if h not in self.labels:
                        self.labels.append(h)
            for fields in c:
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
        for g, gdata in Utils.get_iterator(self.table):
            out.write(g)
            for l in self.labels:
                out.write("\t")
                if l in gdata:
                    out.write(gdata[l])
                else:
                    out.write("0.0")
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

### GeneSorter

class GeneSorter():
    filename = None
    outfile  = None              # -o option
    idcol    = 0                 # -i option
    valcols  = []                # -c option
    topperc  = 100.0             # -t option
    nd       = 0                 # number of values read

    def __init__(self):
        self.valcols = []

    def parseArgs(self, args):
        next = ""

        for a in args:
            if next == "-o":
                self.outfile = a
                next = ""
            elif next == "-i":
                self.idcol = S.toInt(a) - 1
                next = ""
            elif next == "-c":
                self.valcols = Utils.parseColspec(a)
                next = ""
            elif next == "-t":
                self.topperc = S.toFloat(a)
                next = ""
            elif a in ["-o", "-i", "-c", "-t"]:
                next = a
            else:
                self.filename = S.isFile(a)
        self.nvals = len(self.valcols)
        if not self.filename:
            S.usage("sort")

    def getAvgValue(self, row):
        """Returns the average of the values indexed by valcols in `row',
ignoring out of bounds errors and fields that don't contain numbers."""
        v = 0
        n = 0
        for c in self.valcols:
            try:
                v += float(row[c])
                n += 1
            except (ValueError, IndexError):
                pass
        if n == 0:
            return 0
        else:
            return v/n

    def doSort(self, out):
        data = []
        self.nd = 0
        f = None
        try:
            if self.filename == "-":
                f = sys.stdin
            else:
                f = open(self.filename, "r")
            for row in csv.reader(f, delimiter="\t"):
                if len(row) > 0 and len(row[0]) > 0 and row[0][0] != "#":
                    a = self.getAvgValue(row)
                    data.append( (a, row[self.idcol]) )
                    self.nd += 1
        finally:
            if f:
                f.close()
        data.sort(key=lambda r: r[0], reverse=True)

        toWrite = round(self.nd * self.topperc / 100.0)
        for d in data:
            if toWrite == 0:
                break
            out.write("{}\t{}\n".format(d[1], d[0]))
            toWrite -= 1

    def run(self):
        sys.stderr.write("Averaging columns {} from file {}\n".format(", ".join([str(c) for c in self.valcols]), self.filename))
        sys.stderr.write("Gene IDs in column {}\n".format(self.idcol))
        if self.outfile:
            with open(self.outfile, "w") as out:
                self.doSort(out)
        else:
            self.doSort(sys.stdout)

class Ttester(object):
    filename = ""
    range1 = None
    range2 = None
    outfile = "/dev/stdout"

    def parseArgs(self, args):
        self.filename = args[0]
        self.range1 = self.parseRange(args[1])
        self.range2 = self.parseRange(args[2])
        return True

    def parseRange(self, r):
        if "+" in r:
            parts = r.split("+")
            a = int(parts[0])
            b = int(parts[1])
            return range(a, a+b)
        elif "-" in r:
            parts = r.split("+")
            a = int(parts[0])
            b = int(parts[1])
            return range(a, a+b)
        else:
            return [int(r)]

    def run(self):
        with open(self.outfile, "w") as out, open(self.filename, "r") as f:
            out.write("Gene\tAvg1\tAvg2\tLog2FC\tP-value\n")
            f.readline()        # skip header
            c = csv.reader(f, delimiter='\t')
            for line in c:
                gid = line[0]
                vals1 = [ float(line[i]) for i in self.range1 ]
                vals2 = [ float(line[i]) for i in self.range2 ]
                avg1 = sum(vals1) / len(vals1)
                avg2 = sum(vals2) / len(vals2)
                if avg1 != 0.0 and avg2 != 0.0:
                    log2fc = math.log(avg1 / avg2, 2.0)
                    (ttest, pval) = scipy.stats.ttest_ind(vals1, vals2)
                    out.write("{}\t{}\t{}\t{}\t{}\n".format(gid, avg1, avg2, log2fc, pval))

### Main

def parseArgs(cmd, args):
    if cmd == 'matrix':
        P = MatrixGenerator()
    elif cmd == 'merge':
        P = DiffMerger()
    elif cmd == 'allexp':
        P = ExpMerger()
    elif cmd == 'sort':
        P = GeneSorter()
    elif cmd == 'extract':
        P = ExtractRows()
    elif cmd == 'ttest':
        P = Ttester()
    else:
        return S.usage()
    P.parseArgs(args)
    return P

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

