#!/usr/bin/env python

import sys
import pysam

import Script

def usage():
    sys.stderr.write("""regionsCount.py - Compute coverage in specified regions.

Usage: regionsCount.py [-z] [-q qual] [-o outfile] bamfile bedfile

This program examines a BAM file `bamfile' and computes coverage in all intervals
contained in BED file `bedfile'. For each line in the BED file the output (sent to 
standard output, or to `outfile' if specified) consists of the original contents of
the line followed by two additional columns: the total coverage in the interval and
its RPKM (number of reads in the interval in millions divided by the size of the 
interval in kB).

Options:
 -h         | print this usage message.
 -o outfile | write output to `outfile' (default: standard output)
 -q qual    | discard reads with quality score below `qual' (default: {})
 -z         | do not discard intervals with coverage of 0

""".format(BAMreader.qual))

P = Script.Script("regionsCount.py", version="1.0", usage=usage)

B = None
bedfile = None
outfile = None

class BAMreader(object):
    bamfile = None
    aln = None                  # Alignment object
    zeros = True
    qual = 10
    nalignments = 0

    def setBamfile(self, bamfile):
        self.bamfile = bamfile
        self.aln = pysam.AlignmentFile(self.bamfile, "rb")

    def countAlignments(self):
        sys.stderr.write("Counting alignments in {}... ".format(self.bamfile))
        self.nalignments = self.aln.count(read_callback='all')
        sys.stderr.write("done, {} alignments.\n".format(self.nalignments))
        return self.nalignments

    def countAlignmentsInRegion(self, chrom, start, end):
        aiter = self.aln.fetch(chrom, start, end)
        c = 0
        for a in aiter:
            if a.mapping_quality >= self.qual:
                c += 1
        return c

    def addCountsToBED(self, bedfile, out, zeros=True):
        self.countAlignments()
        rpkm_factor = 1000000000.0/self.nalignments
        with open(bedfile, "r") as f:
            for line in f:
                line = line.rstrip("\r\n")
                if len(line) > 0 and line[0] != '#':
                    parsed = line.split("\t")
                    start = int(parsed[1])
                    end = int(parsed[2])
                    c = self.countAlignmentsInRegion(parsed[0], start, end)
                    rl = end-start
                    if rl > 0:
                        rpkm = c * rpkm_factor / rl
                        if c > 0 or self.zeros:
                            out.write(line + "\t" + str(c) + "\t" + str(rpkm) + "\n")
                else:
                    out.write(line)

def parseArgs(args):
    global B
    global bedfile
    global outfile
    nra = 0
    next = ''

    P.standardOpts(args)
    B = BAMreader()
    for a in args:
        if a == '-h':
            usage()
        elif next == '-o':
            outfile = a
            next = ''
        elif next == '-q':
            B.qual = P.toInt(a)
            next = ''
        elif a == '-z':
            B.zeros = False
        elif a in ['-o', '-q']:
            next = a
        elif nra == 0:
            B.setBamfile(P.isfile(a))
            nra += 1
        elif nra == 1:
            bedfile = P.isFile(a)
            nra += 1
        else:
            usage()
    if nra != 2:
        usage()
    return B

if __name__ == "__main__":
    parseArgs(sys.argv[1:])
    if outfile:
        with open(outfile, "w") as out:
            B.addCountsToBED(bedfile, out)
    else:
        B.addCountsToBED(bedfile, sys.stdout)
