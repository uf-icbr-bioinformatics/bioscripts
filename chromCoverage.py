#!/usr/bin/env python

import sys
import pysam

import Script

def usage(what=None):
    if what == "bedmode":
        sys.stderr.write("""chromCoverage.py - Bedmode: report coverage in a set of regions

Usage: chromCoverage.py [-o F] -b B [-a] bamfile

Read regions from BED file B, and write to standard output (or to F if -o is specified)
a tab-delimited file with one row for each region and the following columns:

  chrom start end nreads nbases avgcov

where nreads is the number of reads that (partially) overlap this region; nbases is the 
sum of all bases from the reads overlapping this region; avgcov is the nbases divided by
the length of the region. If a read partially overlaps the region, it will be counted 
fractionally: for example, if a 100bp read has an overlap of 40 with the region, it will
be counted as 0.4 reads.

If the -a option is specified, each row in the output file will consist of the entire 
row from the original BED file followed by the three columns nreads, nbases, avgcov.

""")
    else:
        sys.stderr.write("""chromCoverage.py - Report per-chromosome coverage.

Usage: chromCoverage.py [-o F] [-m MIN] [-c CHROM] bamfile
       chromCoverage.py [-o F] [-m MIN] [coveragefile]
       chromCoverage.py [-o F] -b B [-a] bamfile
       chromCoverage.py [-o F] -x coveragefiles...

Read coverage data from a BAM file or from a coveragefile (or standard input
if not filename is provided) and write by-chromosome coverage data to standard 
output. The coveragefile should be in the format produced by the bamtools 'coverage' 
command:

  chrom   position   depth

The output file contains six columns:

  chrom  total  length  coverage  efflen  effcov

total    - sum of depth at all positions
length   - observed chromosome length (highest position)
coverage - ratio between total and length
efflen   - number of bases having depth >= MIN
effperc  - percent of bases having depth >= MIN
effcov   - ratio between total and efflen

Options:

 -h, --help | Print this usage message.
 -v         | Print version number.
 -o F       | Write output to file F (default: stdout).
 -m MIN     | only positions with depth over MIN are considered when
              computing 'total' and 'efflen' (default: {}).
 -c CHROM   | only output data for chromosome CHROM (default: all chromosomes).
 -x         | combine multiple output files into a single coverage file.
 -b B       | Annotate BED file B (see -h bedmode)

""".format(Cov.mincov))
    
P = Script.Script("chromCoverage.py", version="1.0", usage=usage)

class Cov():
    chrom = ""
    wanted = None
    outfile = None
    out = sys.stdout
    mincov = 0
    total = 0
    maxpos = 0
    effbases = 0
    bedfile = None
    allcols = False             # If true, output all columns from BED file

    def dump(self):
        print "total={}, maxpos={}, effbases={}".format(self.total, self.maxpos, self.effbases)

    def report(self):
        if self.chrom != "":
            if self.maxpos > 0 and self.effbases > 0:
                self.out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.chrom, self.total, self.maxpos, 1.0*self.total/self.maxpos, 
                                                                     self.effbases, round(100.0*self.effbases/self.maxpos, 1), 
                                                                     1.0*self.total/self.effbases))
            else:
                self.out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.chrom, self.total, self.maxpos, 0, self.effbases, 0, 0))
                
        self.total = 0
        self.maxpos = 0
        self.effbases = 0

    def add(self, chrom, pos, cov, tot):
        if chrom != self.chrom:
            self.update(tot)
            self.report()
            # print "Started chrom {}".format(chrom)
            self.chrom = chrom
        self.maxpos = pos
        self.effbases += 1
        self.total += cov

    def update(self, other):
        other.maxpos += self.maxpos
        other.effbases += self.effbases
        other.total += self.total

    def doCovFile(self, stream, Tot):
        try:
            self.out.write("#Chrom\tTotal\tLength\tCoverage\tEfflen\tEffperc\tEffcov\n")
            while True:
                line = stream.readline()
                if not line:
                    break
                parsed = line.split("\t")
                if len(parsed) > 2:
                    chrom = parsed[0]
                    pos = int(parsed[1])
                    cov = int(parsed[2])
                    if cov > self.mincov:
                        self.add(chrom, pos, cov, Tot)
            self.update(Tot)
            self.report()
            Tot.report()
        finally:
            stream.close()

    def doBAMfile(self, filename, Tot):
        bf = pysam.AlignmentFile(filename, "rb" )
        self.out.write("#Chrom\tTotal\tLength\tCoverage\tEfflen\tEffperc\tEffcov\n")
        try:
            try:
                for pileupcolumn in bf.pileup(self.wanted):
                    if pileupcolumn.n > self.mincov:
                        self.add(pileupcolumn.reference_name, pileupcolumn.pos, pileupcolumn.n, Tot)
            except KeyboardInterrupt:
                pass
            except ValueError:
                pass
            self.update(Tot)
            self.report()
            Tot.report()
        finally:
            bf.close()

    def doBEDfile(self, bamfile, bedfile):
        bf = pysam.AlignmentFile(bamfile, "rb")
        self.out.write("#Chrom\tStart\tEnd\tReads\tBases\tAvgCov\n")
        try:
            with open(bedfile, "r") as bed:
                for bedline in bed:
                    if len(bedline) == 0:
                        continue
                    if bedline[0] == 0:
                        continue
                    parsed = bedline.rstrip("\r\n").split("\t")
                    if len(parsed) < 3:
                        continue
                    nreads = 0
                    nbases = 0
                    start = int(parsed[1])
                    end = int(parsed[2])
                    for read in bf.fetch(parsed[0], start, end):
                        if read.is_unmapped:
                            continue
                        ov = read.get_overlap(start, end)
                        readlen = read.query_length
                        nreads += (1.0 * ov) / readlen
                        nbases += ov
                        # print "s={}, e={}, l={}, ov={}, nr={}, nb={}".format(read.reference_start, read.reference_end, readlen, ov, nreads, nbases)
                        # raw_input()
                    if self.allcols:
                        self.out.write("\t".join(parsed) + "\t{:.2f}\t{}\t{:.2f}\n".format(nreads, nbases, (1.0 * nbases) / (end - start)))
                    else:
                        self.out.write("{}\t{}\t{}\t{:.2f}\t{}\t{:.2f}\n".format(parsed[0], start, end, nreads, nbases, (1.0 * nbases) / (end - start)))
        finally:
            bf.close()

def parseArgs(C, args):
    next = ""
    filenames = []

    for a in args:
        if next == '-m':
            C.mincov = P.toInt(a)
            next = ""
        elif next == "-o":
            C.outfile = a
            next = ""
        elif next == "-c":
            C.wanted = a
            next = ""
        elif next == "-b":
            C.bedfile = P.isFile(a)
            next = ""
        elif a in ['-m', "-o", "-c", "-b"]:
            next = a
        elif a == "-a":
            C.allcols = True
        else:
            filenames.append(P.isFile(a))
    return filenames

if __name__=="__main__":

    C = Cov()
    Tot = Cov()
    Tot.chrom = "Total"
    args = sys.argv[1:]
    P.standardOpts(args)

    filenames = parseArgs(C, args)
    try:
        if C.outfile:
            C.out = open(C.outfile, "w")
        if filenames == []:
            C.doCovFile(sys.stdin, Tot)
        elif filenames[0].endswith(".bam"):
            if C.bedfile:
                C.doBEDfile(filenames[0], C.bedfile)
            else:
                C.doBAMfile(filenames[0], Tot)
        else:
            C.doCovFile(open(filenames[0], "r"), Tot)
    finally:
        C.out.close()

