#!/usr/bin/env python

import sys
import csv
import numpy as np
import scipy.stats
from collections import defaultdict

from Utils import BEDreader, DualBEDreader, MATreader, METHreader, REGreader, readDelim, filenameNoExt, safeInt, stringp

import Script
from Regions import BEDdict

# COMMANDS = "merge, avgmeth, histmeth, dmr, dmr2, winavg, winmat, cmerge, regavg, corr, dodmeth"

# def usage(what=None):
#     if what == None:
#         allcmd = allCommands()
#         sys.stderr.write("""dmaptools.py - Operate on methylation data.

# Usage: dmaptools.py command arguments...

# where command is one of: {}

# """.format(", ".join(allcmd)))
#         for cmd in allcmd:
#             cl = CLASSES[cmd]
#             sys.stderr.write("{} - {}.\n".format(cmd, cl.__doc__))

#         sys.stderr.write("""
# Use `dmaptools.py -h command' to display usage for `command'.

# """)
#     elif what == 'regavg':
#         sys.stderr.write("""dmaptools.py - Operate on methylation data.

#""")
#    elif what in CLASSES:
#        cl = CLASSES[what]
#        cl.usage()

#def dummy(what):
#    if what == 'regavg':
#        sys.stderr.write("""dmaptools.py - Operate on methylation data.

# Usage: dmaptools.py regavg [options] bedfile regionsfile

# Compute average methylation over a set of regions (e.g. gene transcripts) mapping
# site positions to a fixed-size vector. File `bedfile' should have at least four
# columns: chromosome, site start, site end, methylation. File `regionsfile' should 
# have at least four columns: chromsome, region start, region end, strand (+ or -). 

# Output is tab-delimited with three columns: position, average methylation at that 
# position, moving average at that position (computed over a window extending S 
# positions in both directions, where S is specified with the -s option).

# Options:

#   -o O | Write results to file O (default: stdout).
#   -w W | Map regions to a vector of size W (default: {}).
#   -m M | Map up/downstream of region to a vector of size M (default: {}).
#   -s S | Use smooting window of S positions (default: {}).
#   -u U | Size of up/downstream regions in bp in non-scaled mode (default: {}).
#   -f   | Do not scale up/downstream regions (in this case, -u is used to specify
#          size of up/downstream regions in bp).

# """.format(REGAVG.vectsize, REGAVG.margsize, REGAVG.winsize, REGAVG.upsizebp))

#     elif what == 'merge':
#         sys.stderr.write("""dmaptools.py - Operate on methylation data.

# Usage: dmaptools.py merge mcompfile matfile1 matfile2 [outfile]

# Merge methylation data from two "mat" files `matfile1' and `matfile2' at the sites
# listed in `mcompfile', containing differentially-methylated C positions. Write
# the results to standard output or to `outfile' if specified.

# """)
#     elif what == 'avgmeth':
#         sys.stderr.write("""dmaptools.py - Operate on methylation data.

# Usage: dmaptools.py avgmeth matfile1 matfile2 [outfile]

# Compute the average global methylation rates for all replicates in `matfile1' and `matfile2'
# and report them, then test the significance of the difference between the two groups using
# a two-sample t-test. The output (written to standard output or to `outfile' if specified)
# has the following format:

# File1: (name of matfile1)
# Replicates1: (number of replicates in matfile1)
# MethRates1: (comma-separated average methylation value for each replicate in matfile1)
# File2: (name of matfile2)
# Replicates2: (number of replicates in matfile2)
# MethRates2: (comma-separated average methylation value for each replicate in matfile2)
# Tstatistics: (T statistics from t-test)
# P-value: (P-value from t-test)
# Significant: (Y or N indicating if P-value < 0.01)

# """)
#     elif what == 'histmeth':
#         sys.stderr.write("""dmaptools.py - Operate on methylation data.

# Usage: dmaptools.py histmeth matfile1 matfile2 [outfile]

# (more to come)

# """)
#     elif what == 'dmr':
#         sys.stderr.write("""dmaptools.py - Operate on methylation data.

# Usage: dmaptools.py dmr [options] testbed ctrlbed

# This command compares two BED files containing methylation rates for two different conditions, 
# (test and control respectively) and detects regions of differential methylation (DMRs). The 
# genome is divided into consecutive windows of `winsize' nucleotides. A window is identified as 
# a DMR if: it contains at least `minsites1' sites in the test dataset and `minsites2' sites in
# the control dataset having coverage higher than `mincov'; if the difference between the methylation
# rates in test and control is higher than `methdiff'; and if the P-value of this difference, 
# computed using Fisher's exact test, is smaller than `pval'.

# Consecutive DMRs will be joined if the are separated by not more than `gap' windows. Only DMRs 
# with a differential methylation in the same direction (ie, both positive or both negative) will
# be joined, unless the -a option is specified. When regions are joined, the differential methylation
# value is the average of those of the joined regions (or the average of their absolute values if
# -a was specified) and the P-value is the highest of those in the joined regions.

# Options:

#  -o outfile   | Write output to `outfile' instead of standard output.
#  -w winsize   | Set window size (default: {}).
#  -t minsites1 | Minimum number of sites from test in DMR (default: {}).
#  -s minsites2 | Minimum number of sites from control in DMR (default: {}).
#  -c mincov    | Minimum coverage of sites for -t and -s (default: {}).
#  -d methdiff  | Minimum difference of methylation rates (default: {}).
#  -p pval      | P-value threshold (default: {}).
#  -g gap       | Maximum gap for DMR joining (default: {}).
#  -a           | Allow joining of DMRs in different directions.

# """.format(DMR.winsize, DMR.minsites1, DMR.minsites2, DMR.mincov, DMR.methdiff, DMR.pval, DMR.gap))

#     elif what == 'dmr2':
#         sys.stderr.write("""dmaptools.py - Operate on methylation data.

# Usage: dmaptools.py dmr2 [options] testbed ctrlbed

# This command compares two BED files containing methylation rates for two different conditions, 
# (test and control respectively) and detects regions of differential methylation (DMRs). A DMR
# is defined as a sequence of at least `minsites' consecutive sites having coverage higher than
# `mincov', and all showing differential methylation in the same direction, over `methdiff'.

# The output file contains chromosome, start, and end position of each DMR, the average diffmeth
# of all sites in the DMR, and the number of sites.

# Options:

#  -o outfile   | Write output to `outfile' instead of standard output.
#  -t minsites  | Minimum number of sites from test in DMR (default: {}).
#  -c mincov    | Minimum coverage of sites for -t and -s (default: {}).
#  -d methdiff  | Minimum difference of methylation rates (default: {}).
#  -s maxdist   | Maximum distance between sites in a window (default: {}).
#  -w minsize   | Minimum size of a DMR (default: {}).
#  --insig num  | Maximum number of sites below methdiff allowed in DMR (default: {}).

# """.format(DMR2writer.minsites, DMR2.mincov, DMR2.methdiff, DMR2writer.maxsitedist, DMR2writer.mindmrsize, DMR2writer.insig_max))

#     elif what == 'winavg':
#         sys.stderr.write("""dmaptools.py - Operate on methylation data.

# Usage: dmaptools.py winavg [options] bedfile

# This command generates a BED file containing average metylation in consecutive windows from a BED
# file containing methylation rates. The input file `bedfile' should be in the format produced by
# mcall or cscall; in particular the first four columns should contain chromosome, start of site,
# end of site, site % methylation respectively. The output file will also have four columns: 
# chromosome, start of window, end of window, average % methylation in window.

# Options:
 
#  -o outfile   | Write output to `outfile' instead of standard output.
#  -w winsize   | Set window size (default: {}).
#  -t minsites  | Minimum number of sites in window (default: {}).
#  -c mincov    | Minimum coverage of sites for -t (default: {}).

# """.format(WINAVG.winsize, WINAVG.minsites, WINAVG.mincov))

#     elif what == 'winmat':
#         sys.stderr.write("""dmaptools.py - Operate on methylation data.

# Usage: dmaptools.py winmat [options] matfile

# This command generates a BED file containing average metylation in consecutive windows from a -mat
# file.

# Options:
 
#  -o outfile   | Write output to `outfile' instead of standard output.
#  -w winsize   | Set window size (default: {}).
#  -t minsites  | Minimum number of sites in window (default: {}).

# """.format(WINAVG.winsize, WINAVG.minsites, WINAVG.mincov))
    
#     elif what == 'cmerge':
#         sys.stderr.write("""dmaptools.py - Operate on methylation data.

# Usage: dmaptools.py cmerge [options] files...

# This command merges columns from multiple input files into a single output file. The input files
# should have chromosome, start, and end position in the first three columns, and a value associated
# with each region in an additional column (the fourth one by default). For each region in the input
# files, the output file will contain chromosome, start, end, and the target columns from all input
# files in order.

# Options:
 
#  -o outfile   | Write output to `outfile' instead of standard output.
#  -c targetcol | Specify column containing values to be merged (default: {}).
#  -x missing   | Value to use for missing elements (default: {}).
#  -s           | Skip first row in each input file (default: False).

# """.format(CMERGE.targetcol, CMERGE.missing))

#     else:
#         P.usage()

# Merger

class Merger(Script.Command):
    """report per-replicate methylation values at differentially methylated sites"""
    _cmd = "merge"
    mcompfile = None
    matfile1 = None
    matfile2 = None
    outfile = None

    def parseArgs(self, args):
        nfiles = 0
        for a in args:
            if nfiles == 0:
                self.mcompfile = P.isFile(a)
                nfiles += 1
            elif nfiles == 1:
                self.matfile1 = P.isFile(a)
                nfiles += 1
            elif nfiles == 2:
                self.matfile2 = P.isFile(a)
                nfiles += 1
            elif nfiles == 3:
                self.outfile = a
        if nfiles < 3:
            P.errmsg(P.NOFILE)

    def usage(self, parent, out=sys.stdout):
        out.write("""dmaptools.py - Operate on methylation data.

Usage: dmaptools.py merge mcompfile matfile1 matfile2 [outfile]

Merge methylation data from two "mat" files `matfile1' and `matfile2' at the sites
listed in `mcompfile', containing differentially-methylated C positions. Write
the results to standard output or to `outfile' if specified.

""")

    def makeMatMap(self, matfile):
        hdr = []
        mmap = {}
        with open(matfile, "r") as f:
            c = csv.reader(f, delimiter='\t')
            line = next(c)
            hdr = line[4:]
            p = f.tell()
            for line in c:
                key = line[0] + ":" + line[1]
                # print("{} -> {}".format(key, p))
                # raw_input()
                mmap[key] = p
                p = f.tell()
        return (hdr, mmap)

    def mergeMatFiles(self, out):
        """Merge the contents of `matfile1' and `matfile2' at all locations contained in `mcompfile'.
    Write results to `outfile'."""
        hdr1 = []
        hdr2 = []
        nhdr1 = 0
        nhdr2 = 0
        map1 = {}
        map2 = {}

        sys.stderr.write("Indexing sites in `{}'... ".format(self.matfile1))
        (hdr1, map1) = self.makeMatMap(self.matfile1)
        sys.stderr.write("done, {} sites found.\n".format(len(map1)))
        sys.stderr.write("Indexing sites in `{}'... ".format(self.matfile2))
        (hdr2, map2) = self.makeMatMap(self.matfile2)
        sys.stderr.write("done, {} sites found.\n".format(len(map2)))

        nhdr1 = len(hdr1)
        nhdr2 = len(hdr2)

        # print("Header lengths: {}, {}".format(nhdr1, nhdr2))
        out.write("Chrom\tPos\tC1:Avg\tC1:Stdev")
        for h in hdr1:
            out.write("\tC1:" + h)
        out.write("\tC2:Avg\tC2:Stdev")
        for h in hdr2:
            out.write("\tC2:" + h)
        out.write("\n")

        nwritten = 0
        sys.stderr.write("Writing data for sites in `{}'... ".format(self.mcompfile))
        with open(self.matfile1, "r") as m1:
            with open(self.matfile2, "r") as m2:
                with open(self.mcompfile, "r") as f:
                    f.readline()
                    for line in f:
                        line = line.split("\t")
                        key = line[0] + ":" + line[1]
                        if key in map1 and key in map2:
                            fp1 = map1[key]
                            m1.seek(fp1)
                            dl1 = m1.readline().rstrip("\r\n")
                            fp2 = map2[key]
                            m2.seek(fp2)
                            dl2 = readDelim(m2)
                            out.write(dl1 + "\t" + "\t".join(dl2[2:]) + "\n")
                            nwritten += 1
        sys.stderr.write("done, {} sites in output.\n".format(nwritten))
        return nwritten
    
    def run(self):
        if self.outfile:
            with open(self.outfile, "w") as out:
                self.mergeMatFiles(out)
        else:
            self.mergeMatFiles(sys.stdout)

### ColMerger

class ColMergerRecord():
    chrom = ""
    start = 0
    end = 0
    data = []

    def __init__(self, chrom, start, end, nsources, missing="NA"):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.data = [missing]*nsources

    def addValue(self, idx, value):
        self.data[idx] = value

class ColMerger():
    filenames = []
    nsources = 0                # Number of files we're reading from
    target = 0                  # Column containing value to be merged
    chroms = []                 # List of chromosomes seen
    reclists = {}               # Dictionary of records seen (by chromosome)
    missing = "NA"              # Value to use for missing data
    skiphdr = False

    def __init__(self, filenames, target, skiphdr=False):
        self.filenames = filenames
        self.colnames = [ self.cleanFilename(f) for f in filenames ]
        self.nsources = len(filenames)
        self.target = target
        self.chroms = []
        self.recids = []
        self.records = {}
        self.skiphdr = skiphdr

    def cleanFilename(self, f):
        """Move this somewhere else..."""
        p = f.rfind("/")
        if p > 0:
            f = f[p+1:]
        if f[-4:] == ".csv":
            f = f[:-4]
        return f

    def parseSource(self, idx):
        """Parse source number `idx'."""
        sys.stderr.write("Parsing `{}'... ".format(self.filenames[idx]))
        with open(self.filenames[idx], "r") as f:
            if self.skiphdr:
                f.readline()
            for line in f:
                if len(line) > 0 and line[0] != '#':
                    parsed = line.rstrip("\r\n").split("\t")
                    chrom = parsed[0]
                    start = int(parsed[1])
                    end   = int(parsed[2])
                    val   = parsed[self.target]
                    if chrom not in self.chroms:
                        self.chroms.append(chrom)
                    if chrom in self.records:
                        chromlist = self.records[chrom]
                    else:
                        self.records[chrom] = chromlist = {}
                    if start in chromlist:
                        rec = chromlist[start]
                    else:
                        chromlist[start] = rec = ColMergerRecord(chrom, start, end, self.nsources, missing=self.missing)
                    rec.addValue(idx, val)
        sys.stderr.write("done.\n")

    def parseAllSources(self):
        for idx in range(self.nsources):
            self.parseSource(idx)

    def writeMerged(self, out):
        out.write("#Chrom\tStart\tEnd\t" + "\t".join(self.colnames) + "\n")
        for chrom in self.chroms:
            chromlist = self.records[chrom]
            for pos in sorted(chromlist.keys()):
                rec = chromlist[pos]
                out.write("{}\t{}\t{}".format(rec.chrom, rec.start, rec.end))
                for v in rec.data:
                    out.write("\t{}".format(v))
                out.write("\n")

### Averager
### This command reads two -mat files and computes the average methylation of each
### sample across both conditions. It then compares the two sets of averages using
### a two-sample t-test.

class Averager(Script.Command):
    """report and compare genome-wide methylation levels in two sets of replicates"""
    _cmd = "avgmeth"
    matfile1 = None
    nreps1 = 0
    matfile2 = None
    nreps2 = 0
    outfile = None
    pval = 0.01

    def parseArgs(self, args):
        next = ""
        for a in args:
            if a in ['-p']:
                next = a
            elif next == '-p':
                self.pval = P.toFloat(a)
                next = ""
            if self.matfile1 == None:
                self.matfile1 = P.isFile(a)
            elif self.matfile2 == None:
                self.matfile2 = P.isFile(a)
            elif self.outfile == None:
                self.outfile = a
        return (self.matfile1 and self.matfile2)

    def usage(self, parent, out=sys.stdout):
        out.write("""dmaptools.py - Operate on methylation data.

Usage: dmaptools.py avgmeth matfile1 matfile2 [outfile]

Compute the average global methylation rates for all replicates in `matfile1' and `matfile2'
and report them, then test the significance of the difference between the two groups using
a two-sample t-test. The output (written to standard output or to `outfile' if specified)
has the following format:

File1: (name of matfile1)
Replicates1: (number of replicates in matfile1)
MethRates1: (comma-separated average methylation value for each replicate in matfile1)
File2: (name of matfile2)
Replicates2: (number of replicates in matfile2)
MethRates2: (comma-separated average methylation value for each replicate in matfile2)
Tstatistics: (T statistics from t-test)
P-value: (P-value from t-test)
Significant: (Y or N indicating if P-value < 0.01)

""")
                
    def methAvg(self, filename):
        nrows = 0
        sys.stderr.write("Reading {}... ".format(filename))
        with open(filename, "r") as f:
            hdr = readDelim(f)
            nreps = len(hdr) - 4
            counts = [0]*nreps
            sums = [0]*nreps

            for line in f:
                nrows += 1
                fields = line.split("\t")
                for i in range(nreps):
                    d = fields[i+4]
                    if d != "NA":
                        v = float(d)
                        if v >= 0:
                            counts[i] += 1
                            sums[i] += v

        goodreps = 0            # replicates for which we have data
        avgs = []
        for i in range(nreps):
            if counts[i] > 0:
                avgs.append(sums[i] / counts[i])
                goodreps += 1
        sys.stderr.write("{} rows.\n".format(nrows))
        return (goodreps, avgs)

    def report(self, out, avgs1, avgs2):
        out.write("File1: " + self.matfile1 + "\n")
        out.write("Replicates1: " + str(self.nreps1) + "\n")
        out.write("MethRates1: " + ",".join([str(x) for x in avgs1]) + "\n")
        out.write("File2: " + self.matfile2 + "\n")
        out.write("Replicates2: " + str(self.nreps2) + "\n")
        out.write("MethRates2: " + ",".join([str(x) for x in avgs2]) + "\n")
        (tstat, pval) = scipy.stats.ttest_ind(avgs1, avgs2)
        out.write("Tstatistics: " + str(tstat) + "\n")
        out.write("P-value: " + str(pval) + "\n")
        out.write("Significant: " + ("Y" if pval <= self.pval else "N") + "\n")

    def run(self):
        (nreps, avgs1) = self.methAvg(self.matfile1)
        self.nreps1 = nreps
        (nreps, avgs2) = self.methAvg(self.matfile2)
        self.nreps2 = nreps
        if self.outfile:
            with open(self.outfile, "w") as out:
                self.report(out, avgs1, avgs2)
        else:
            self.report(sys.stdout, avgs1, avgs2)

### Histcomparer
### This command reads two -mat files and reports the fraction of sites in each bin of % methylation
### for each replicate in the two conditions. It then determines whether the fractions are significantly
### different in each bin.

class Histcomparer(Averager):
    """compare % methylation rates in two sets of replicates"""
    _cmd = "histmeth"

    def usage(self, parent, out=sys.stdout):
        out.write("""dmaptools.py - Operate on methylation data.

Usage: dmaptools.py histmeth matfile1 matfile2 [outfile]

(more to come)

""")
        
    def methHist(self, filename):
        nrows = 0
        sys.stderr.write("Reading {}... ".format(filename))
        with open(filename, "r") as f:
            hdr = f.readline().split("\t")
            nreps = len(hdr) - 4
            counts = [0]*nreps
            hist = [ [0]*nreps for i in range(10) ]

            for line in f:
                nrows += 1
                fields = line.split("\t")
                for i in range(nreps):
                    d = fields[i+4]
                    if d != "NA":
                        v = float(d)
                        if v >= 0:
                            bin = int(v*10)
                            if bin == 10:
                                bin = 9
                            counts[i] += 1
                            hist[bin][i] += 1
        # print(hist)
        # print(counts)
        fracs = [ [ 1.0*hist[b][i] / counts[i] if counts[i] else 0 for i in range(nreps) ] for b in range(10) ]
        sys.stderr.write("{} rows.\n".format(nrows))
        return (nreps, fracs)

    def report(self, out, fracs1, fracs2):
        out.write("#Bin\tPval\tSignificant\tAvg1\tAvg2\n")
        for b in range(10):
            label = "{}-{}%".format(b*10, (b+1)*10)
            row1 = fracs1[b]
            row2 = fracs2[b]
            (tstat, pval) = scipy.stats.ttest_ind(row1, row2)
            out.write("\t".join([label, str(pval), ("Y" if pval <= self.pval else "N"), str(sum(row1) / self.nreps1), str(sum(row2) / self.nreps2)]) + "\n")

    def run(self):
        (nreps, fracs1) = self.methHist(self.matfile1)
        self.nreps1 = nreps
        (nreps, fracs2) = self.methHist(self.matfile2)
        self.nreps2 = nreps
        if self.outfile:
            with open(self.outfile, "w") as out:
                self.report(out, fracs1, fracs2)
        else:
            self.report(sys.stdout, fracs1, fracs2)

## DMRs

class DMRwriter():
    out = None
    maxdist = 0                 # Maximum distance between DMRs for joining
    samedir = True              # If true, only DMRs in the same direction will be joined
    chrom = ""
    nd = 0
    growing = []
    
    def __init__(self, out, maxdist, samedir=True):
        self.out = out
        self.maxdist = maxdist
        self.samedir = samedir
        self.growing = []
        out.write("#Chrom\tStart\tEnd\tDiff\tPval\n")

    def writeDMR(self):
        """Write out the DMRs in `growing'."""
        if self.nd == 0:
            return
        elif self.nd == 1:
            self.out.write("{}\t{}\t{}\t{}\t{}\n".format(*self.growing[0]))
        else:
            totdiff = 0.0
            maxpval = 0.0
            for d in self.growing:
                end = d[2]
                if self.samedir:
                    totdiff += d[3]
                else:
                    totdiff += abs(d[3])
                maxpval = max(maxpval, d[4])
            self.out.write("{}\t{}\t{}\t{}\t{}\n".format(self.growing[0][0], self.growing[0][1], end, totdiff/self.nd, maxpval))
        self.growing = []
        self.nd = 0

    def finish(self):
        self.writeDMR()

    def addDMR(self, data):
        if data[0] != self.chrom: # If chrom is different...
            self.writeDMR()       # write current DMRs and start over
            self.chrom = data[0]
        elif self.nd == 0:
            pass
        elif (data[1] - self.growing[-1][2]) <= self.maxdist: # these can be joined...
            if self.samedir and (data[3] * self.growing[-1][3]) < 0: # but not if they go in different directions and we want samedir
                self.writeDMR()
        else:
            self.writeDMR()       # can't join, write current DMRs and start over
        self.growing.append(data)
        self.nd += 1

class DMR(Script.Command):
    """detect differentially methylated regions"""
    _cmd = "dmr"
    winsize = 100
    minsites1 = 0   # Minimum number of sites in test condition
    minsites2 = 4   # Minimum number of sites in control condition
    mincov = 4
    methdiff = 0.2
    pval = 0.01
    gap = 0         # Maximum number of non-DMR regions that can be joined
    samedir = True  # Only join DMRs in same direction?
    bedfile1 = None # BED file for test condition
    bedfile2 = None # BED file for control condition
    outfile = None
    growing = None              # buffer for DMRs that can potentially be joined
    jump = False                # Skip to this chrom?
    one = False                 # Do a single chromosome?

    def parseArgs(self, args):
        next = ""
        for a in args:
            if next == '-w':
                self.winsize = P.toInt(a)
                next = ""
            elif next == '-s':
                self.minsites2 = P.toInt(a)
                next = ""
            elif next == '-t':
                self.minsites1 = P.toInt(a)
                next = ""
            elif next == '-c':
                self.mincov = P.toInt(a)
                next = ""
            elif next == '-d':
                self.methdiff = P.toFloat(a)
                next = ""
            elif next == '-p':
                self.pval = P.toFloat(a)
                next = ""
            elif next == '-g':
                self.gap = P.toInt(a)
                next = ""
            elif next == '-o':
                self.outfile = a
                next = ""
            elif next == '-j':
                self.jump = a
                next = ""
            elif next == "-J":
                self.jump = a
                self.one = True
                next = ""
            elif a in ['-w', '-s', '-t', '-c', '-d', '-p', '-o', '-g', '-j', '-J']:
                next = a
            elif a == '-a':
                self.samedir = False
            elif self.bedfile1 == None:
                self.bedfile1 = P.isFile(a)
            else:
                self.bedfile2 = P.isFile(a)
        if self.bedfile1 == None or self.bedfile2 == None:
            P.errmsg(P.NOFILE)
        return True

    def usage(self, parent, out=sys.stdout):
        out.write("""dmaptools.py - Operate on methylation data.

Usage: dmaptools.py dmr [options] testbed ctrlbed

This command compares two BED files containing methylation rates for two different conditions, 
(test and control respectively) and detects regions of differential methylation (DMRs). The 
genome is divided into consecutive windows of `winsize' nucleotides. A window is identified as 
a DMR if: it contains at least `minsites1' sites in the test dataset and `minsites2' sites in
the control dataset having coverage higher than `mincov'; if the difference between the methylation
rates in test and control is higher than `methdiff'; and if the P-value of this difference, 
computed using Fisher's exact test, is smaller than `pval'.

Consecutive DMRs will be joined if the are separated by not more than `gap' windows. Only DMRs 
with a differential methylation in the same direction (ie, both positive or both negative) will
be joined, unless the -a option is specified. When regions are joined, the differential methylation
value is the average of those of the joined regions (or the average of their absolute values if
-a was specified) and the P-value is the highest of those in the joined regions.

Options:

 -o outfile   | Write output to `outfile' instead of standard output.
 -w winsize   | Set window size (default: {}).
 -t minsites1 | Minimum number of sites from test in DMR (default: {}).
 -s minsites2 | Minimum number of sites from control in DMR (default: {}).
 -c mincov    | Minimum coverage of sites for -t and -s (default: {}).
 -d methdiff  | Minimum difference of methylation rates (default: {}).
 -p pval      | P-value threshold (default: {}).
 -g gap       | Maximum gap for DMR joining (default: {}).
 -a           | Allow joining of DMRs in different directions.

""".format(self.winsize, self.minsites1, self.minsites2, self.mincov, self.methdiff, self.pval, self.gap))

    def isDMR(self, data1, data2):
        """data1 = test, data2 = control."""
        ngood1 = 0
        totC1 = 0
        totT1 = 0
        ngood2 = 0
        totC2 = 0
        totT2 = 0

        pos1 = [ d[2] for d in data1 ]
        pos2 = [ d[2] for d in data2 ]

        for d in data1:
            if d[0] >= self.mincov and d[2] in pos2:
                ngood1 += 1
                totC1 += d[1]
                totT1 += (d[0]-d[1])
        for d in data2:
            if d[0] >= self.mincov and d[2] in pos1:
                ngood2 += 1
                totC2 += d[1]
                totT2 += (d[0]-d[1])

        # Do we have a sufficient number of sites?
        if ngood1 < self.minsites1 or ngood2 < self.minsites2:
            return None

        # Compute methdiff ratio. Is it above our threshold?
        ratio1 = 1.0 * totC1 / (totC1 + totT1)
        ratio2 = 1.0 * totC2 / (totC2 + totT2)
        diff = ratio1 - ratio2
        if abs(diff) < self.methdiff:
            return None
        #print(diff, [[totC1, totT1], [totC2, totT2]])

        # Compute p-value
        (odds, pval) = scipy.stats.fisher_exact([[totC1, totT1], [totC2, totT2]])
        #print(odds, pval, diff, [[totC1, totT1], [totC2, totT2]])
        if pval <= self.pval:
            return (pval, diff)
        else:
            return None
        
    def findDMRs(self, out, avgout=None):
        DW = DMRwriter(out, self.gap*self.winsize, samedir=self.samedir)
        BR1 = BEDreader(self.bedfile1, jump=self.jump)
        BR2 = BEDreader(self.bedfile2, jump=self.jump)
        chrom = BR1.chrom       # assume that both BED files start with the same chrom - should check this!
        
        start = 0
        end = self.winsize
        nfound = 0
        totfound = 0

        while BR1.stream != None and BR2.stream != None:
            data1 = BR1.readUntil(chrom, end)
            data2 = BR2.readUntil(chrom, end)

            if stringp(data1): # BR1 is at new chrom?
                sys.stderr.write("{}: {} DMRs\n".format(chrom, nfound))
                if self.one:
                    break
                BR2.skipToChrom(data1)        # BR2 should join it
                totfound += nfound
                nfound = 0
                chrom = data1
                start = 0
                end = self.winsize

            elif stringp(data2): # BR2 is at new chrom?
                sys.stderr.write("{}: {} DMRs\n".format(chrom, nfound))
                if self.one:
                    break
                BR1.skipToChrom(data2)        # BR1 should join it
                totfound += nfound
                nfound = 0
                chrom = data2
                start = 0
                end = self.winsize
            
            else:               # Both are lists, OK
                if len(data1) > 0 and len(data2) > 0:
                    result = self.isDMR(data1, data2)
                    if result:
                        nfound += 1
                        (pval, diff) = result
                        DW.addDMR([chrom, start, end, diff, pval])
                    if avgout:
                        ratios1 = [ (1.0*v[1])/v[0] for v in data1 ]
                        ratios2 = [ (1.0*v[1])/v[0] for v in data2 ]
            # Move forward
            start += self.winsize                                                    
            end += self.winsize
        sys.stderr.write("{}: {} DMRs\n".format(chrom, nfound))
        sys.stderr.write("Total: {} DMRs\n".format(totfound))
        DW.finish()
        # BR1.close()
        # BR2.close()

    def run(self):
        if self.outfile:
            with open(self.outfile, "w") as out:
                self.findDMRs(out)
        else:
            self.findDMRs(sys.stdout)

class DMR2writer():
    out = None
    minsites = 3                # Minimum number of sites in a window (-t)
    maxsitedist = 1000          # Maximum distance between sites in a window (-s)
    mindmrsize = 1000           # Minimum size of a DMR (-w)
    chrom = ""                  # Chrom of current DMR
    pos = 0                     # Position of last site added
    direction = ""              # Either "+" or "-"
    positions = []              # Position of sites in this DMR
    data = []                   # Diffmeth values for sites in this DMR
    insig_max = None            # Maximum number of insignificant sites (--insig)
    insig_count = 0             # Number of insignificant sites in this DMR
   
    
    def start(self, chrom, pos, diff):
        self.maybeWriteDMR()
        self.chrom = chrom
        self.pos = pos
        self.direction = "+" if diff > 0 else "-"
        self.positions = [pos]
        self.data = [diff]

    def maybeWriteDMR(self):
        ns = len(self.positions)
        if ns < self.minsites:
            #sys.stderr.write("skipping due to sites\n")
            return
        start = self.positions[0]
        end = self.positions[-1]
        if end - start < self.mindmrsize:
            #sys.stderr.write("skipping due to length\n")
            return
        #sys.stderr.write("{}\n".format(self.data))
        #sys.stderr.write("Writing DMR: {}\n".format((self.chrom, start, end, end-start, sum(self.data) / ns, ns)))
        self.out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(self.chrom, start, end, end-start, sum(self.data) / ns, ns))

    def add(self, chrom, pos, diff):
        # If we are tracking insignificant sites...
        if self.insig_max is not None:
            # ...and the # of insig. sites is greater than max
            if self.insig_count > self.insig_max:
                # ... start new DMR
                self.start(chrom, pos, diff)
                self.insig_count = 0
                return
        #sys.stderr.write("Adding: {}\n".format((chrom, pos, diff)))

        # If we changed chromosomes...
        if chrom != self.chrom:
            # ... start new DMR
            self.start(chrom, pos, diff)
            return
        thisDir = "+" if diff > 0 else "-"
        # If we changed chromosomes...
        if thisDir != self.direction:
            # ... start new DMR
            self.start(chrom, pos, diff)
            return

        # If we're too far away from previous site...
        if (pos - self.pos) > self.maxsitedist:
            # ... start new DMR
            self.start(chrom, pos, diff)
            return

        
        # Otherwise, extend current DMR
        self.pos = pos
        self.positions.append(pos)
        self.data.append(diff)
        
class DMR2(Script.Command):
    """detect differentially methylated regions (method #2)"""
    _cmd = "dmr2"
    DW       = None          # DMR2writer
    mincov   = 4             # Minimum coverage of sites considered (-c)
    methdiff = 0.2           # Minimum diff meth (-d)
    bedfile1 = None          # BED file for test condition
    bedfile2 = None          # BED file for control condition
    outfile  = None
    track_insig = False      # True if insignificant sites should be tracked

    def parseArgs(self, args):
        self.DW = DMR2writer()
        next = ""
        for a in args:
            if next == '-w':
                self.DW.mindmrsize = P.toInt(a)
                next = ""
            elif next == '-t':
                self.DW.minsites = P.toInt(a)
                next = ""
            elif next == '-s':
                self.DW.maxsitedist = P.toInt(a)
                next = ""
            elif next == '-c':
                self.mincov = P.toInt(a)
                next = ""
            elif next == '-d':
                self.methdiff = P.toFloat(a)
                next = ""
            elif next == '-p':
                self.pval = P.toFloat(a)
                next = ""
            elif next == '-o':
                self.outfile = a
                next = ""
            elif next == '--insig':
                self.DW.insig_max = int(a)
                self.track_insig = True
                next = ""
            elif a in ['-w', '-s', '-t', '-c', '-d', '-p', '-o', '--insig']:
                next = a
            elif self.bedfile1 == None:
                self.bedfile1 = P.isFile(a)
            else:
                self.bedfile2 = P.isFile(a)
        if self.bedfile1 == None or self.bedfile2 == None:
            P.errmsg(P.NOFILE)
        return True

    def usage(self, parent, out=sys.stdout):
        out.write("""dmaptools.py - Operate on methylation data.

Usage: dmaptools.py dmr2 [options] testbed ctrlbed

This command compares two BED files containing methylation rates for two different conditions, 
(test and control respectively) and detects regions of differential methylation (DMRs). A DMR
is defined as a sequence of at least `minsites' consecutive sites having coverage higher than
`mincov', and all showing differential methylation in the same direction, over `methdiff'.

The output file contains chromosome, start, and end position of each DMR, the average diffmeth
of all sites in the DMR, and the number of sites.

Options:

 -o outfile   | Write output to `outfile' instead of standard output.
 -t minsites  | Minimum number of sites from test in DMR (default: {}).
 -c mincov    | Minimum coverage of sites for -t and -s (default: {}).
 -d methdiff  | Minimum difference of methylation rates (default: {}).
 -s maxdist   | Maximum distance between sites in a window (default: {}).
 -w minsize   | Minimum size of a DMR (default: {}).
 --insig num  | Maximum number of sites below methdiff allowed in DMR (default: {}).

""".format(DMR2writer.minsites, self.mincov, self.methdiff, DMR2writer.maxsitedist, DMR2writer.mindmrsize, DMR2writer.insig_max))
            
    def findDMRs(self, out):
        self.DW.out = out
        out.write("#Chrom\tStart\tEnd\tLen\tDiffmeth\tNsites\n")
        BR = DualBEDreader(self.bedfile1, self.bedfile2)
        while True:
            if BR.readNext():
                c1 = BR.current1
                c2 = BR.current2
                if c1[0] < self.mincov or c2[0] < self.mincov:
                    continue
                m1 = 1.0 * c1[1] / c1[0]
                m2 = 1.0 * c2[1] / c2[0]
                d = m1 - m2
                if abs(d) > self.methdiff:
                    self.DW.add(BR.chrom, BR.pos, d)
                # if we are tracking
                elif self.track_insig:
                    self.DW.insig_count += 1
            else:
                return
        self.DW.maybeWriteDMR()
        
    def run(self):
        try:
            if self.outfile:
                with open(self.outfile, "w") as out:
                    self.findDMRs(out)
            else:
                self.findDMRs(sys.stdout)
        except KeyboardInterrupt:
            return
        except IOError:
            return

### WINAVG

class WINAVG(Script.Command):
    """generate a BED file containing average methylation in consecutive windows"""
    _cmd = "winavg"
    winsize = 100
    minsites = 0   # Minimum number of sites in window
    mincov = 4     # Minimum coverage of sites counted
    bedfile = None # Input file
    outfile = None # Output file

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == '-w':
                self.winsize = P.toInt(a)
                prev = ""
            elif prev == '-t':
                self.minsites = P.toInt(a)
                prev = ""
            elif prev == '-c':
                self.mincov = P.toInt(a)
                prev = ""
            elif prev == '-o':
                self.outfile = a
                prev = ""
            elif a in ['-w', '-t', '-c', '-o']:
                prev = a
            elif self.bedfile is None:
                self.bedfile = P.isFile(a)
        if self.bedfile is None:
            P.errmsg(P.NOFILE)
        else:
            return True

    def usage(self, parent, out=sys.stdout):
        out.write("""dmaptools.py - Operate on methylation data.

Usage: dmaptools.py winavg [options] bedfile

This command generates a BED file containing average metylation in consecutive windows from a BED
file containing methylation rates. The input file `bedfile' should be in the format produced by
mcall or cscall; in particular the first four columns should contain chromosome, start of site,
end of site, site % methylation respectively. The output file will also have four columns: 
chromosome, start of window, end of window, average % methylation in window.

Options:
 
 -o outfile   | Write output to `outfile' instead of standard output.
 -w winsize   | Set window size (default: {}).
 -t minsites  | Minimum number of sites in window (default: {}).
 -c mincov    | Minimum coverage of sites for -t (default: {}).

""".format(self.winsize, self.minsites, self.mincov))

    def getAvg(self, data):
        ngood = 0
        totCov = 0
        totC = 0

        for d in data:
            if d[0] >= self.mincov:
                ngood += 1
                totCov += d[0]
                totC += d[1]
        if ngood >= self.minsites:
            return 1.0*totC/totCov
        else:
            return False

    def winAvg(self, out):
        BR = BEDreader(self.bedfile)
        chrom = BR.chrom
        start = 0
        end = self.winsize
        nwins = 0
        totwins = 0

        while BR.stream != None:
            data = BR.readUntil(chrom, end)
            if stringp(data): # New chrom?
                sys.stderr.write("{}: {} windows\n".format(chrom, nwins))
                totwins += nwins
                nwins = 0
                chrom = data
                start = 0
                end = self.winsize
            else:
                if len(data) > 0:
                    avg = self.getAvg(data)
                    if avg:
                        out.write("{}\t{}\t{}\t{}\n".format(chrom, start, end, avg))
                        nwins += 1
                start += self.winsize
                end += self.winsize
        sys.stderr.write("{}: {} windows\n".format(chrom, nwins))
        sys.stderr.write("Total: {} windows\n".format(totwins))

    def run(self):
        if self.outfile:
            with open(self.outfile, "w") as out:
                self.winAvg(out)
        else:
            self.winAvg(sys.stdout)

### WINMAT

class WINMAT(Script.Command):
    """like winavg, but using a -mat file as input"""
    _cmd = "winmat"
    winsize = 100
    minsites = 0   # Minimum number of sites in window
    matfile = None # Input file
    outfile = None # Output file

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == '-w':
                self.winsize = P.toInt(a)
                prev = ""
            elif prev == '-t':
                self.minsites = P.toInt(a)
                prev = ""
            elif prev == '-o':
                self.outfile = a
                prev = ""
            elif a in ['-w', '-t', '-o']:
                prev = a
            elif self.matfile == None:
                self.matfile = P.isFile(a)
        if self.matfile == None:
            P.errmsg(P.NOFILE)

    def usage(self, parent, out=sys.stdout):
        out.write("""dmaptools.py - Operate on methylation data.

Usage: dmaptools.py winmat [options] matfile

This command generates a BED file containing average metylation in consecutive windows from a -mat
file.

Options:
 
 -o outfile   | Write output to `outfile' instead of standard output.
 -w winsize   | Set window size (default: {}).
 -t minsites  | Minimum number of sites in window (default: {}).

""".format(self.winsize, self.minsites))
    
    def winMat(self, out):
        BR = MATreader(self.matfile)
        chrom = BR.chrom
        start = 0
        end = self.winsize
        nwins = 0
        totwins = 0

        out.write("#Chrom\tStart\tEnd\t" + "\t".join(BR.hdr[4:]) + "\n")
        while BR.stream != None:
            data = BR.readUntil(chrom, end)
            if stringp(data): # New chrom?
                sys.stderr.write("{}: {} windows\n".format(chrom, nwins))
                totwins += nwins
                nwins = 0
                chrom = data
                start = 0
                end = self.winsize
            else:
                nd = len(data)
                if nd > 0:
                    sums = [0]*BR.nreps
                    cnts = [0]*BR.nreps
                    for d in data:
                        for i in range(BR.nreps):
                            v = d[i]
                            if v != 'NA':
                                v = float(v)
                                if v > 0:
                                    sums[i] += float(d[i])
                                    cnts[i] += 1
                    avgs = [ str(sums[i]/cnts[i] if cnts[i] > 0 else 0) for i in range(BR.nreps) ]
                    out.write("{}\t{}\t{}\t{}\n".format(chrom, start, end, "\t".join(avgs)))
                    nwins += 1
                start += self.winsize
                end += self.winsize
        sys.stderr.write("{}: {} windows\n".format(chrom, nwins))
        sys.stderr.write("Total: {} windows\n".format(totwins))

    def run(self):
        if self.outfile:
            with open(self.outfile, "w") as out:
                self.winMat(out)
        else:
            self.winMat(sys.stdout)

### CMERGE

class CMERGE(Script.Command):
    """merge columns from multiple files into a single output file"""
    _cmd = "cmerge"
    colmerger = None
    filenames = []
    targetcol = 3
    outfile = None
    missing = "NA"
    skiphdr = False

    def parseArgs(self, args):
        self.filenames = []
        prev = ""
        for a in args:
            if prev == '-c':
                self.targetcol = P.toInt(a) - 1
                prev = ""
            elif prev == '-o':
                self.outfile = a
                prev = ""
            elif prev == '-x':
                self.missing = a
                prev = ""
            elif a in ['-c', '-o', '-x']:
                prev = a
            elif a == '-s':
                self.skiphdr = True
            else:
                self.filenames.append(P.isFile(a))
        if len(self.filenames) == 0:
            P.errmsg(P.NOFILE)
        return True

    def usage(self, parent, out=sys.stdout):
        out.write("""dmaptools.py - Operate on methylation data.

Usage: dmaptools.py cmerge [options] files...

This command merges columns from multiple input files into a single output file. The input files
should have chromosome, start, and end position in the first three columns, and a value associated
with each region in an additional column (the fourth one by default). For each region in the input
files, the output file will contain chromosome, start, end, and the target columns from all input
files in order.

Options:
 
 -o outfile   | Write output to `outfile' instead of standard output.
 -c targetcol | Specify column containing values to be merged (default: {}).
 -x missing   | Value to use for missing elements (default: {}).
 -s           | Skip first row in each input file (default: False).

""".format(self.targetcol, self.missing))

    def run(self):
        self.colmerger = ColMerger(self.filenames, self.targetcol, self.skiphdr)
        self.colmerger.parseAllSources()
        if self.outfile:
            sys.stderr.write("Writing merged file to {}\n".format(self.outfile))
            with open(self.outfile, "w") as out:
                self.colmerger.writeMerged(out)
        else:
            self.colmerger.writeMerged(sys.stdout)

### REGAVG

class REGAVG(Script.Command):
    """compute average methylation over a set of regions"""
    _cmd = "regavg"
    bedfile   = None
    regfile   = None
    outfile   = None
    genesfile = None            # File containing names of genes represented in output
    label     = None

    # Sizes
    vectsize = 1000             # Size of vector representing transcript
    margsize = 100              # Size of vectors representing up/downstream
    upsizebp = 2000             # Size (in bp) of up/downstream region in non-scaled mode
    winsize  = 40               # Smooting window size
    scaled   = True             # If false, up/down regions are not scaled

    # Computed
    margfact = 0.0
    totsize = 0                 # Length of averages vector
    vector = None
    regreader = None
    bedreader = None

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == '-o':
                self.outfile = a
                prev = ""
            elif prev == "-w":
                self.vectsize = P.toInt(a)
                prev = ""
            elif prev == "-m":
                self.margsize = P.toInt(a)
                prev = ""
            elif prev == "-u":
                self.upsizebp = P.toInt(a)
                prev = ""
            elif prev == "-s":
                self.winsize = P.toInt(a)
                prev = ""
            elif prev == "-l":
                self.label = a
                prev = ""
            elif prev == "-g":
                self.genesfile = a
                prev = ""
            elif a in ["-o", "-w", "-m", "-u", "-s", "-l", "-g"]:
                prev = a
            elif a == "-f":
                self.scaled = False
            elif self.bedfile is None:
                self.bedfile = P.isFile(a)
            else:
                self.regfile = P.isFile(a)
        if self.scaled:
            self.margfact = 1.0 * self.margsize / self.vectsize
        else:
            self.margfact = 1.0 * self.margsize / self.upsizebp
        self.totsize = self.margsize + self.vectsize + self.margsize
            
        self.vector = np.zeros((2, self.totsize))

    def usage(self, parent, out=sys.stdout):
        out.write("""dmaptools.py - Operate on methylation data.

Usage: dmaptools.py regavg [options] bedfile regionsfile

Compute average methylation over a set of regions (e.g. gene transcripts) mapping
site positions to a fixed-size vector. File `bedfile' should have at least four
columns: chromosome, site start, site end, methylation. File `regionsfile' should 
have at least four columns: chromsome, region start, region end, strand (+ or -). 

Output is tab-delimited with three columns: position, average methylation at that 
position, moving average at that position (computed over a window extending S 
positions in both directions, where S is specified with the -s option).

Options:

  -o O | Write results to file O (default: stdout).
  -w W | Map regions to a vector of size W (default: {}).
  -m M | Map up/downstream of region to a vector of size M (default: {}).
  -s S | Use smooting window of S positions (default: {}).
  -u U | Size of up/downstream regions in bp in non-scaled mode (default: {}).
  -f   | Do not scale up/downstream regions (in this case, -u is used to specify
         size of up/downstream regions in bp).

""".format(self.vectsize, self.margsize, self.winsize, self.upsizebp))

    def run(self):
        self.bedreader = METHreader(self.bedfile, skipHdr=False)
        self.bedreader.readNext()
        self.regreader = REGreader(self.regfile, skipHdr=False)
        self.regreader.readNext()
        totst = 0               # Total number of sites
        tottx = 0               # Total number of transcripts
        totns = 0               # Total number of sites detected
        totnt = 0               # Total number of transcripts with at least one site
        updn = self.upsizebp    # Assume we're in fixed mode

        sys.stderr.write("#Chromosome\tTranscripts\tSites\tFound Transcripts\tFound Sites\n")
        if self.genesfile:
            genesout = open(self.genesfile, "w")
        else:
            genesout = None
        try:
            while True:
                regions = self.regreader.readChromosome()
                if regions is None:
                    break
                sites = self.bedreader.readChromosome()
                if not sites:
                    continue
                chrom = regions[0][0]
                nst = len(sites)
                ntx = len(regions)
                totst += nst
                tottx += ntx
                ns = 0
                nt = 0
                
                for tx in regions:
                    found = False
                    txlen = tx[2] - tx[1]
                    if self.scaled:
                        updn = int(txlen * self.margfact)

                    p1 = tx[1] - updn
                    p2 = tx[2] + updn
                    for site in sites:
                        if site[1] > p2:
                            break
                        if site[1] >= p1:
                            if self.scaled:
                                self.storeSite(site, p1, p2, tx[3])
                            else:
                                self.storeSiteFixed(site, tx)
                            ns += 1
                            if not found:
                                found = True
                                if genesout:
                                    genesout.write("{}\n".format(tx[4]))
                    if found:
                        nt += 1
                sys.stderr.write("{}\t{}\t{}\t{}\t{}\n".format(chrom, ntx, nst, nt, ns))
                totns += ns
                totnt += nt
        finally:
            if genesout:
                genesout.close()
        sys.stderr.write("Total\t{}\t{}\t{}\t{}\n".format(tottx, totst, totnt, totns))

        if self.outfile:
            with open(self.outfile, "w") as out:
                self.writeResults(out)
        else:
            self.writeResults(sys.stdout)

    def writeResults(self, stream):
        if self.label:
            stream.write("#Position\tTX\t{} - average\t{} - moving average\n".format(self.label, self.label))
        else:
            stream.write("#Position\tTX\taverage\tmoving average\n")
        avgs = self.vector[0] / self.vector[1]
        mavgs = np.zeros(self.totsize)
        b = self.totsize - self.margsize
        for idx in range(self.winsize, self.totsize - self.winsize):
            mavgs[idx] = sum(self.vector[0][idx-self.winsize:idx+self.winsize]) / sum(self.vector[1][idx-self.winsize:idx+self.winsize]) 
        for idx in range(self.totsize):
            if self.margsize < idx < b:
                tx = "0"
            else:
                tx = "."
            if mavgs[idx] == 0.0:
                stream.write("{}\t{}\t{}\t\n".format(idx - self.margsize, tx, avgs[idx]))
            else:
                stream.write("{}\t{}\t{}\t{}\n".format(idx - self.margsize, tx, avgs[idx], mavgs[idx]))

    def storeSite(self, site, p1, p2, strand):
        size = p2 - p1
        frac = 1.0 * (site[1] - p1) / size
        if strand == '-':
            frac = 1.0 - frac
        idx = int(round(frac * (self.totsize - 1)))
        self.vector[0][idx] += site[2]
        self.vector[1][idx] += 1

    def storeSiteFixed(self, site, tx):
        p = site[1]
        if p < tx[1]:
            d = (tx[1] - p)
            # sys.stderr.write("d = {}\n".format(d))
            if tx[3] == '+':
                idx = self.margsize - int(round(d * self.margfact))
            else:
                idx = self.margsize + self.vectsize + int(round(d * self.margfact)) - 1
        elif p > tx[2]:
            d = (p - tx[2])
            # sys.stderr.write("d = {}\n".format(d))
            if tx[3] == '+':
                idx = self.margsize + self.vectsize + int(round(d * self.margfact)) - 1
            else:
                idx = self.margsize - int(round(d * self.margfact))
        else:
            size = tx[2] - tx[1]
            frac = 1.0 * (site[1] - tx[1]) / size
            if tx[3] == '-':
                frac = 1.0 - frac
            idx = self.margsize + int(round(frac * (self.vectsize - 1)))
        if idx == self.totsize:
            idx = self.totsize - 1 # hack
        # sys.stderr.write("Idx = {}\n".format(idx))
        self.vector[0][idx] += site[2]
        self.vector[1][idx] += 1
        
### Correlation of methylation values

class Colpair():                # Used by CORR
    col1 = None
    col2 = None
    name1 = ""
    name2 = ""
    v1 = []
    v2 = []

    def __init__(self, c1, c2):
        self.col1 = c1 + 3
        self.col2 = c2 + 3
        self.v1 = []
        self.v2 = []

class CORR(Script.Command):
    """compute correlation between methylation levels of replicates of a condition"""
    _cmd = "corr"
    matfile = None
    outfile = None
    colpairs = []

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-o":
                self.outfile = a
                prev = ""
            elif a == "-o":
                prev = a
            elif self.matfile is None:
                self.matfile = P.isFile(a)
            else:
                pair = self.splitCols(a)
                if pair:
                    cp = Colpair(pair[0], pair[1])
                    self.colpairs.append(cp)

    def usage(self, parent, out=sys.stdout):
        out.write("TBW\n")
        
    def splitCols(self, s):
        pieces = s.split(",")
        if len(pieces) == 2:
            try:
                a = int(pieces[0])
                b = int(pieces[1])
                return (a, b)
            except:
                sys.stderr.write("Warning: argument `{}' should be of the form P,Q where P and Q are column numbers.\n")
        return None

    def run(self):
        name1 = ""
        name2 = ""
        v1 = []
        v2 = []

        with open(self.matfile, "r") as f:
            c = csv.reader(f, delimiter='\t')
            hdr = next(c)
            for cp in self.colpairs:
                cp.name1 = hdr[cp.col1]
                cp.name2 = hdr[cp.col2]
            for line in c:
                for cp in self.colpairs:
                    x1 = float(line[cp.col1])
                    x2 = float(line[cp.col2])
                    if x1 >= 0 and x2 >= 0:
                        cp.v1.append(x1)
                        cp.v2.append(x2)

        for cp in self.colpairs:
            a1 = np.array(cp.v1, dtype=float)
            a2 = np.array(cp.v2, dtype=float)
            cc = np.corrcoef(a1, a2)
            sys.stdout.write("Sample1:\t{}\n".format(cp.name1))
            sys.stdout.write("Sample2:\t{}\n".format(cp.name2))
            sys.stdout.write("Num sites:\t{}\n".format(len(a1)))
            sys.stdout.write("Mean1:\t{}\n".format(np.mean(a1)))
            sys.stdout.write("Mean2:\t{}\n".format(np.mean(a2)))
            sys.stdout.write("Correlation:\t{}\n\n".format(cc[0,1]))

### Difference of differential methylation rates.

class DIFF(Script.Command):
    """Compute the difference between differential methylation rates in different contrasts"""
    _cmd = "dodmeth"
    bedfile1 = None
    bedfile2 = None
    outfile = None
    column = 4
    threshold = 0.0

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-o":
                self.outfile = a
                prev = ""
            elif prev == "-c":
                self.column = P.toInt(a) - 1
                prev = ""
            elif prev == "-t":
                self.threshold = P.toFloat(a)
                prev = ""
            elif a in ["-o", "-c", "-t"]:
                prev = a
            elif self.bedfile1 is None:
                self.bedfile1 = P.isFile(a)
            else:
                self.bedfile2 = P.isFile(a)
        if self.bedfile1 == None or self.bedfile2 == None:
            P.errmsg(P.NOFILE)

    def usage(self, parent, out=sys.stdout):
        out.write("TBW.\n")

    def run(self):
        if self.outfile:
            with open(self.outfile, "w") as out:
                self.do_difference(out)
        else:
            self.do_difference(sys.stdout)

    def do_difference(self, out):
        with open(self.bedfile1, "r") as f1:
            with open(self.bedfile2, "r") as f2:
                f1.readline()
                f2.readline()
                r1 = csv.reader(f1, delimiter='\t')
                r2 = csv.reader(f2, delimiter='\t')
                (nout, nbad) = self.do_difference_aux(out, r1, r2)
                sys.stderr.write("{} differences written.\n".format(nout))
                sys.stderr.write("{} outliers skipped.\n".format(nbad))

    def do_difference_aux(self, out, r1, r2):
        nout = 0                # Site differences written
        nbad = 0                # Outliers

        line1 = r1.next()
        line2 = r2.next()
        chrom = line1[0]        # Assume both start with the same chrom
        pos1 = int(line1[1])
        pos2 = int(line2[1])
        while True:
            try:
                if pos1 == pos2:
                    v1 = float(line1[self.column])
                    v2 = float(line2[self.column])
                    if v1 >= self.threshold and v2 >= self.threshold:
                        diff = v1 - v2
                        out.write("{}\t{}\t{}\t{}\n".format(chrom, pos1, pos1 + 1, diff))
                        nout += 1
                    else:
                        nbad += 1
                    line1 = r1.next()
                    line2 = r2.next()
                    pos1 = int(line1[1])
                    pos2 = int(line2[1])
                elif pos1 < pos2:
                    line1 = r1.next()
                    pos1 = int(line1[1])
                else:
                    line2 = r2.next()
                    pos2 = int(line2[1])

                if line1[0] != chrom:
                    chrom = line1[0]
                    line2 = self.readUntil(r2, chrom)
                    pos2 = int(line2[1])
                elif line2[0] != chrom:
                    chrom = line2[0]
                    line1 = self.readUntil(r1, chrom)
                    pos1 = int(line1[1])
            except StopIteration:
                return (nout, nbad)

    def readUntil(self, r, chrom):
        while True:
            line = r.next()
            if line[0] == chrom:
                return line

### For Metilene output

def makePayload(ncols):

    def payload_maker():
        return [ [0.0, 0] for _ in range(ncols) ]
    return payload_maker

class DMRAVG(Script.Command):
    """Generate per-sample averages in DMR regions from Metilene output."""
    _cmd = "dmravg"
    dmrfile = None
    matfile = None
    outfile = "/dev/stdout"
    data = {}
    ncols = 0
    header = None

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-o":
                self.outfile = a
                prev = ""
            elif a in ["-o"]:
                prev = a
            elif self.dmrfile is None:
                self.dmrfile = a
            elif self.matfile is None:
                self.matfile = a
        if self.dmrfile and self.matfile:
            return True
        else:
            P.errmsg(P.nofile)

    def run(self):
        self.readDMRs()
        self.readMatrix()
        self.writeOutput()

    def readDMRs(self):
        self.data = BEDdict(self.dmrfile)
        sys.stderr.write("{} regions read from DMR file {}.\n".format(self.data.nregs, self.dmrfile))

    def readMatrix(self):
        nsites = 0
        nfound = 0
        with open(self.matfile, "r") as f:
            c = csv.reader(f, delimiter='\t')
            self.header = next(c)
            self.ncols = len(self.header) - 2 
            self.data.setPayloads(makePayload(self.ncols))
            for line in c:
                nsites += 1
                chrom = line[0]
                pos = int(line[1])
                reg = self.data.findPos(chrom, pos)
                if not reg:
                    continue

                # We found a region containing this position
                nfound += 1
                col = 1
                for idx in range(self.ncols):
                    col += 1
                    if line[col] == "NA":
                        continue
                    v = float(line[col])
                    reg.payload[idx][0] += v
                    reg.payload[idx][1] += 1
        sys.stderr.write("{} sites read from file {}.\n{} sites in DMRs.\n".format(nsites, self.matfile, nfound))

    def writeOutput(self):
        with open(self.outfile, "w") as out:
            out.write("#Chrom\tStart\tEnd\t" + "\t".join(self.header[2:]) + "\n")
            for chrom in self.data.allChroms():
                for reg in self.data.chromRegions(chrom):
                    out.write("{}\t{}\t{}".format(reg.chrom, reg.start, reg.end))
                    for p in reg.payload:
                        out.write("\t{}".format(p[0] / p[1] if p[1] > 0 else 0))
                    out.write("\n")
                
### Report average methylation rates by chromosome

class BYCHROM(Script.Command):
    """generate a table of average methylation rate by chromosome in each sample"""
    _cmd = "bychr"
    bedfiles = []
    labels = []
    nsamples = 0
    chroms = []
    data = {}
    outfile = "/dev/stdout"

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-o":
                self.outfile = a
                prev = ""
            elif prev == "-l":
                self.labels = a.split(",")
                prev = ""
            elif a in ["-o", "-l"]:
                prev = a
            else:
                self.bedfiles.append(a)
        self.nsamples = len(self.bedfiles)
        if not self.labels:
            self.labels = [ filenameNoExt(f) for f in self.bedfiles ]
        self.chroms = []
        self.data = {}

    def usage(self, parent, out=sys.stdout):
        out.write("TBW.\n")
        
    def run(self):
        self.readFirst()
        i = 1
        for bed in self.bedfiles[1:]:
            self.readOther(bed, i)
            i += 1

        with open(self.outfile, "w") as out:
            out.write("#Chrom\t" + "\t".join(self.labels) + "\n")
            for ch in self.chroms:
                cdata = self.data[ch]
                out.write(ch + "\t" + "\t".join([str(x) for x in cdata]) + "\n")

    def readFirst(self):
        chrom = ""
        nr = 0
        nc = 0
        sys.stderr.write("Parsing {}... ".format(self.bedfiles[0]))
        with open(self.bedfiles[0], "r") as f:
            c = csv.reader(f, delimiter='\t')
            for line in c:
                this = line[0]
                if "_" in this: # Skip "fake" chromosomes
                    continue
                if this != chrom: # new chromosome?
                    if chrom:
                        # sys.stderr.write("Creating cdata for {}\n".format(chrom))
                        cdata = [0.0]*self.nsamples
                        cdata[0] = float(nc)/nr
                        self.data[chrom] = cdata
                        nr = 0
                        nc = 0
                    self.chroms.append(this)
                    chrom = this
                nr += int(line[4])
                nc += int(line[5])
            cdata = [0.0]*self.nsamples
            cdata[0] = float(nc)/nr
            self.data[chrom] = cdata
        sys.stderr.write("done.\n")

    def readOther(self, filename, i):
        chrom = ""
        nr = 0
        nc = 0
        sys.stderr.write("Parsing {}... ".format(filename))
        with open(filename, "r") as f:
            c = csv.reader(f, delimiter='\t')
            for line in c:
                this = line[0]
                if "_" in this:
                    continue
                if this != chrom: # new chromosome?
                    if chrom:
                        # sys.stderr.write("Creating cdata for {}\n".format(chrom))
                        cdata = self.data[chrom]
                        cdata[i] = float(nc)/nr
                        nr = 0
                        nc = 0
                    chrom = this
                nr += int(line[4])
                nc += int(line[5])
            cdata = self.data[chrom]
            cdata[i] = float(nc)/nr
        sys.stderr.write("done.\n")

# CpG Islands stats

class CpGisland():
    start  = 0                  # island
    end    = 0
    start2 = 0                  # shore
    end2   = 0
    start3 = 0                  # shelf
    end3   = 0
    shore  = 2000
    shelf  = 4000
    
    def __init__(self, start, end):
        self.start  = start
        self.end    = end
        self.start2 = start - self.shore
        self.end2   = end + self.shore
        self.start3 = start - self.shelf
        self.end3   = end + self.shelf

    def classify(self, position):
        if self.start <= position <= self.end:
            return "I"
        elif self.start2 <= position <= self.end2:
            return "S"
        elif self.start3 <= position <= self.end3:
            return "F"
        else:
            return ""

class CPGCOUNT(Script.Command):
    """Count number of sites in CpG islands, shores, shelves in one or more BED file."""
    _cmd = "cpg"
    cpgfile = None
    bedfiles = []
    cpgdb = {}
    reportfile = "/dev/stdout"
    outfile = True

    def __init__(self):
        self.bedfiles = []
        self.cpgdb = {}

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-r":
                self.reportfile = a
                prev = ""
            elif prev == "-f":
                CpGisland.shelf = int(a)
                prev = ""
            elif prev == "-e":
                CpGisland.shore = int(a)
                prev = ""
            elif a in ["-r", "-f", "-e"]:
                prev = a
            elif a == "-x":
                self.outfile = False
            elif self.cpgfile is None:
                self.cpgfile = a
            else:
                self.bedfiles.append(a)
        return (self.cpgfile and self.bedfiles)

    def usage(self, parent, out=sys.stdout):
        out.write("""dmaptools.py - Operate on methylation data.

Usage: dmaptools.py cpg [options] cpgfile bedfiles...

This command reports the number of methylation sites in the specified BED files 
that fall within or close to the CpG islands defined in `cpgfile'. Sites are classified
as being internal to the island, in the island shore, or in the island shelf.

Output is written to a file with the same name as each bedfile, with extension .cpg.csv,
unless -x is specified, in which case no output is written.

Options:

  -x   | Do not write output file(s).
  -r R | Write overall report to file R (default: no report)
  -e E | Specify size of CpG island shore (default: {})
  -f F | Specify size of CpG island shelf (default: {})

""".format(CpGisland.shore, CpGisland.shelf))

    def run(self):
        self.readCpGfile()
        with open(self.reportfile, "w") as out:
            out.write("#Filename\tSites\tNisland\tNshore\tNshelf\tSitelist\n")
            for bed in self.bedfiles:
                (n0, n1, n2, n3, of) = self.countSites(bed)
                out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(bed, n0, n1, n2, n3, of))

    def readCpGfile(self):
        sys.stderr.write("Reading CpG locations from file {}... ".format(self.cpgfile))
        ncpg = 1
        with open(self.cpgfile, "r") as f:
            c = csv.reader(f, delimiter='\t')
            for line in c:
                chrom = line[0]
                if chrom not in self.cpgdb:
                    self.cpgdb[chrom] = []
                s = CpGisland(int(line[1]), int(line[2]))
                self.cpgdb[chrom].append(s)
                ncpg += 1
        sys.stderr.write("done, {} sites read.\n".format(ncpg))

    def findIsland(self, chrom, pos):
        if chrom not in self.cpgdb:
            return None
        ilist = self.cpgdb[chrom]
        for i in ilist:
            if i.start3 <= pos <= i.end3:
                return i
            if i.start3 > pos:
                return None

    def countSites(self, bedfile):
        n0 = n1 = n2 = n3 = 0
        out = None
        if self.outfile:
            outfilename = filenameNoExt(bedfile) + ".cpg.csv"
            out = open(outfilename, "w")
        try:
            sys.stderr.write("Reading file {}... ".format(bedfile))
            with open(bedfile, "r") as f:
                c = csv.reader(f, delimiter='\t')
                for line in c:
                    chrom = line[0]
                    pos = safeInt(line[1])
                    if pos is None:
                        continue
                    n0 += 1
                    island = self.findIsland(chrom, pos)
                    if island:
                        x = island.classify(pos)
                        if x == "I":
                            n1 += 1
                        elif x == "S":
                            n2 += 1
                        elif x == "F":
                            n3 += 1
                        if x and out:
                            out.write("{}\t{}\t{}\t{}\n".format(bedfile, chrom, pos, x))
            sys.stderr.write("done, {} sites read.\n".format(n0))
        finally:
            if out:
                out.close()
        return (n0, n1, n2, n3, outfilename)

class MatMerger(Script.Command):
    """merge methylation matrix files into a single file"""
    _cmd = "matmerge"
    matfiles = []
    headers = []
    wanted = None
    nona = False
    outfile = "/dev/stdout"

    def __init__(self):
        self.matfiles = []
        self.headers = []
        self.wanted = defaultdict(int)

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-o":
                self.outfile = a
                prev = ""
            elif a in ["-o"]:
                prev = a
            elif a == "-n":
                self.nona = True
            else:
                self.matfiles.append(a)
        return self.matfiles

    def usage(self, parent, out=sys.stdout):
        out.write("""dmaptools.py - Operate on methylation data.

Usage: dmaptools.py matmerge [options] matfiles...

Merge two or more matfiles (SAMPLE-mat.csv) into a single matrix, suitable for MDS plots.
Only sites appearing in all input files will be written to the output.

Options:

  -o O | Write output to file O (default: standard output).
  -n   | Remove lines containing NA values.
""")

    def run(self):
        nfiles = len(self.matfiles)
        for mf in self.matfiles:
            self.updateWanted(mf)
        self.findWanted(nfiles)

        for i in range(nfiles):
            self.readMatfile(self.matfiles[i], i)

        self.writeMatrix()

    def updateWanted(self, filename):
        sys.stderr.write("Reading sites from {}... ".format(filename))
        ns = 0
        with open(filename, "r") as f:
            f.readline()
            c = csv.reader(f, delimiter='\t')
            for row in c:
                ns += 1
                key = row[0] + ":" + row[1]
                self.wanted[key] += 1
        sys.stderr.write("{} sites\n".format(ns))

    def findWanted(self, nfiles):
        newdic = {}
        for k in self.wanted.keys():
            if self.wanted[k] == nfiles:
                newdic[k] = [None]*nfiles
        self.wanted = newdic
        sys.stderr.write("Common sites: {}\n".format(len(self.wanted)))

    def readMatfile(self, filename, i):
        sys.stderr.write("Reading methylation data from {}... ".format(filename))
        with open(filename, "r") as f:
            c = csv.reader(f, delimiter='\t')
            hdr = next(c)
            names = [ n.rstrip(".C") for n in hdr[4:] ]
            self.headers.append("\t".join(names))
            for row in c:
                key = row[0] + ":" + row[1]
                if key in self.wanted:
                    data = row[4:]
                    data = [ "NA" if x == "-1.0" else str(int(1000*float(x))) for x in data ]
                    self.wanted[key][i] = "\t".join(data)
        sys.stderr.write("done.\n")

    def writeMatrix(self):
        sys.stderr.write("Writing merged matrix to {}... ".format(self.outfile))
        with open(self.outfile, "w") as out:
            out.write("\t".join(self.headers) + "\n")
            for k in self.wanted:
                row = self.wanted[k]
                good = True
                if self.nona:
                    for r in row:
                        if "NA" in r:
                            good = False
                            break
                if good:
                    out.write(k + "\t" + "\t".join(row) + "\n")
        sys.stderr.write("done.\n")

## window-based DMR analysis:
### Params:
### window size
### min number of sites in window
### min coverage of each site in each sample
### methylation rate difference
### P-value
### number of non-DMRs for merging
### how to combine diffmeth in joined DMRs
###
### Output:
### chrom  start  end  diffmeth
###
### Additional: 
### average methylation level for each sample in each DMR

### T-test to compare methylation rates in replicates (from mat files)
### include methavg.py - DONE
### in histogram: compare percent in each bin for replicates of each condition - DONE

### in bedplotter: allow for larger window, using line graph instead of bar
### multiple samples overlayed on single graph

### check that mbed plotter works, generate figure 1A

### Re-run pipeline for genediffmeth in more regions - DONE
### Add links to mat files - DONE

### Main

class Prog(Script.Script):

    def main(self, args):
        self.standardOpts(args)
        cmd = args[0]
        cl = self.findCommand(cmd)
        if cl:
            M = cl()
            if M.parseArgs(args[1:]):
                M.run()
            else:
                M.usage(None)
        else:
            self.usage(None)
        
P = Prog("dmaptools.py", version="1.0",
         errors=[('NOCMD', 'Missing command', 'The first argument should be a command name'),
                 ('NOFILE', 'Missing input file(s)', 'One or more input file(s) is missing')])
P.addCommands([Merger, Averager, Histcomparer, DMR, DMR2, WINAVG, WINMAT, CMERGE, REGAVG, CORR, DIFF, BYCHROM, CPGCOUNT, DMRAVG, MatMerger])

cmdlist = ""
maxl = max([len(c) for c in P._commandNames])
fstr = "{:" + str(maxl) + "} - {}.\n"
for cmd in P._commandNames:
    cl = P.findCommand(cmd)
    cmdlist += fstr.format(cmd, cl.__doc__)

P.setDocstrings({'main': """dmaptools.py - Operate on methylation data.

Usage: dmaptools.py command arguments...

where command is one of: {}

{}
Use `dmaptools.py -h command' to display usage for `command'.

""".format(", ".join(P._commandNames), cmdlist)})

if __name__ == "__main__":
    args = sys.argv[1:]
    nargs = len(args)
    if nargs == 0:
        P.errmsg(P.NOCMD)
    P.main(args)
