#!/usr/bin/env python

import sys
import scipy.stats

import Script

# Main

def usage(what=None):
    if what == None:
        sys.stderr.write("""dmaptools.py - Operate on methylation data.

Usage: dmaptools.py command arguments...

where command is one of: merge, avgmeth, dmr

merge - report per-replicate methylation values at differentially methylated sites
avgmeth - report and compare genome-wide methylation levels in two sets of replicates
dmr - detect differentially methylated regions

Use `dmaptools.py -h command' to display usage for `command'.

""")
    elif what == 'merge':
        sys.stderr.write("""dmaptools.py - Operate on methylation data.

Usage: dmaptools.py merge mcompfile matfile1 matfile2 [outfile]

Merge methylation data from two "mat" files `matfile1' and `matfile2' at the sites
listed in `mcompfile', containing differentially-methylated C positions. Write
the results to standard output or to `outfile' if specified.

""")
    elif what == 'avgmeth':
        sys.stderr.write("""dmaptools.py - Operate on methylation data.

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

""")
    else:
        P.usage()

P = Script.Script("dmaptools.py", version="1.0", usage=usage,
                  errors=[('NOCMD', 'Missing command', 'The first argument should be one of: merge, avgmeth, dmr.')])

# Merger

class Merger():
    mcompfile = None
    matfile1 = None
    matfile2 = None
    outfile = None

    def __init__(self, args):
        nfiles = 0
        for a in args:
            if nfiles == 0:
                mcompfile = P.isFile(a)
                nfiles += 1
            elif nfiles == 1:
                matfile1 = P.isFile(a)
                nfiles += 1
            elif nfiles == 2:
                matfile2 = P.isFile(a)
                nfiles += 1
            elif nfiles == 3:
                outfile = a
        if nfiles < 3:
            P.errmsg(P.NOFILE)

    def makeMatMap(self):
        hdr = []
        mmap = {}
        with open(self.matfile, "r") as f:
            line = f.readline().rstrip("\r\n").split("\t")
            hdr = line[4:]
            p = f.tell()
            while True:
                line = f.readline()
                if line == '':
                    break
                line = line.rstrip("\r\n").split("\t")
                key = line[0] + ":" + line[1]
                # print "{} -> {}".format(key, p)
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
        (hdr1, map1) = makeMatMap(self.matfile1)
        sys.stderr.write("done, {} sites found.\n".format(len(map1)))
        sys.stderr.write("Indexing sites in `{}'... ".format(self.matfile2))
        (hdr2, map2) = makeMatMap(self.matfile2)
        sys.stderr.write("done, {} sites found.\n".format(len(map2)))

        nhdr1 = len(hdr1)
        nhdr2 = len(hdr2)

        # print "Header lengths: {}, {}".format(nhdr1, nhdr2)
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
                            dl2 = m2.readline().rstrip("\r\n").split("\t")
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

### Averager
### This command reads two -mat files and computes the average methylation of each
### sample across both conditions. It then compares the two sets of averages using
### a two-sample t-test.

class Averager():
    matfile1 = None
    nreps1 = 0
    matfile2 = None
    nreps2 = 0
    outfile = None

    def __init__(self, args):
        for a in args:
            if self.matfile1 == None:
                self.matfile1 = P.isFile(a)
            elif self.matfile2 == None:
                self.matfile2 = P.isFile(a)
            elif self.outfile == None:
                self.outfile = a

    def methAvg(self, filename):
        nrows = 0
        sys.stderr.write("Reading {}... ".format(filename))
        with open(filename, "r") as f:
            hdr = f.readline().split("\t")
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
        avgs = [ sums[i] / counts[i] for i in range(nreps) ]
        sys.stderr.write("{} rows.\n".format(nrows))
        return (nreps, avgs)

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

### window-based DMR analysis:
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
### include methavg.py
### in histogram: compare percent in each bin for replicates of each condition

### in bedplotter: allow for larger window, using line graph instead of bar
### multiple samples overlayed on single graph

### check that mbed plotter works, generate figure 1A

### Re-run pipeline for genediffmeth in more regions
### Add links to mat files

if __name__ == "__main__":
    args = sys.argv[1:]
    nargs = len(args)
    if nargs == 0:
        P.errmsg(P.NOCMD)

    P.standardOpts(args)

    cmd = args[0]
    
    if cmd == 'merge':
        M = Merger(args[1:])
        M.run()
    elif cmd == 'avgmeth':
        M = Averager(args[1:])
        M.run()
    else:
        P.usage()

