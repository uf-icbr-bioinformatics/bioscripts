#!/usr/bin/env python

import sys
import pysam

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

    B = BAMreader()
    for a in args:
        if a == '-h':
            usage()
            return B
        elif next == '-o':
            outfile = a
            next = ''
        elif next == '-q':
            B.qual = int(a)
            next = ''
        elif a == '-z':
            B.zeros = False
        elif a in ['-o', '-q']:
            next = a
        elif nra == 0:
            B.setBamfile(a)
            nra += 1
        elif nra == 1:
            bedfile = a
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
