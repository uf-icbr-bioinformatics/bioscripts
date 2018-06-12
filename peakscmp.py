#!/usr/bin/env python

import sys
import csv
import math

import Script

def usage(what=None):
    global P
    sys.stderr.write("""{} - compare peak files from MACS2

Usage: {} [options] peaks1 peaks2

Options:

  -o O | Write output to file O (default: standard output)

  -v V | Use V as the overlap fraction to determine if two
         peaks match. The overlap fraction is defined as the
         ration between the intersection of the two peak regions
         divided by their union. Should be a number between 0 
         and 1. Default: {}.

  -f F | Use F as the normalization factor between the two peak
         files. If not specified, this is computed automatically
         based on the numbered of filtered tags.

  -a L | Only return rows in which the absolute value of the fold 
         change is greater than L.

  -A L | Only return rows in which the fold change is greater than
         L (if positive) or smaller than L (if negative). 
         Use "-h limits" for examples.

  -b L | Only return rows in which the absolute value of the fold 
         change is smaller than L.

  -B L | Only return rows in which the fold change is between 0 and
         L (if positive) or between L and 0 (if negative). 
         Use "-h limits" for examples.

""".format(P.name, P.name, P.minOverlap))

class macsDiff(Script.Script):
    filename1 = ""
    filename2 = ""
    outfile = None
    ntags1 = 0
    ntags2 = 0
    factor = None
    regions = {}
    minOverlap = 0.5

    def init(self):
        self.regions   = {}

    def parseArgs(self, args):
        next = ""

        self.standardOpts(args)
        for a in args:
            if next == "-f":
                self.factor = self.toFloat(a)
                next = ""
            elif next == "-v":
                self.minOverlap = self.toFloat(a)
                next = ""
            elif next == "-o":
                self.outfile = a
                next = ""
            elif a in ["-f", "-v", "-o"]:
                next = a
            elif not self.filename1:
                self.filename1 = self.isFile(a)
            elif not self.filename2:
                self.filename2 = self.isFile(a)
        # Validate
        if not (self.filename1 and self.filename2):
            self.errmsg(self.NOFILE)
        if self.minOverlap <= 0 or self.minOverlap > 1:
            self.errmsg(self.BADOVERLAP)

    def readHeader(self, stream):
        ntags = 0

        for line in stream:
            if not line:
                return ntags
            if line.startswith("chr"):
                return ntags
            if line.startswith("# tags after filtering"):
                ntags = int(line.rstrip("\r\n")[37:])
    
    def readRegions(self, stream):
        nregs = 0
        for fields in csv.reader(stream, delimiter='\t'):
            chrom = fields[0]
            if chrom not in self.regions:
                self.regions[chrom] = []
            self.regions[chrom].append((int(fields[1]), int(fields[2]), float(fields[5])))
            nregs += 1
        return nregs

    def readFirst(self):
        with open(self.filename1, "r") as f:
            self.ntags1 = self.readHeader(f)
            nregs = self.readRegions(f)
        sys.stderr.write("Tags 1: {}\n".format(self.ntags1))
        return nregs

    def findRegion(self, regions, start, end):
        for r in regions:
            if r[0] > end:
                return None
            if (start <= r[0] and end >= r[1]) or (r[0] <= start <= r[1]) or (r[0] <= end <= r[1]):
                return r
        return None

    def overlap(self, reg, start, end):
        p1 = min(reg[0], start)
        q1 = max(reg[0], start)
        p2 = max(reg[1], end)
        q2 = min(reg[1], end)
        return 1.0 * (q2-q1) / (p2-p1)

    def readSecond(self, out):
        curr = ""
        currRegs = []

        out.write("#Chrom\tStart1\tEnd1\tPileup1\tStart2\tEnd2\tPileup2\tlog2(FC)\n")
        with open(self.filename2, "r") as f:
            nin   = 0
            nout  = 0
            ninc  = 0
            noutc = 0
            nup   = 0
            ndn   = 0
            self.ntags2 = self.readHeader(f)
            if not self.factor:
                self.factor = 1.0 * self.ntags1 / self.ntags2
            sys.stderr.write("Tags 2: {}\n".format(self.ntags2))
            sys.stderr.write("Factor: {}\n".format(self.factor))
            for fields in csv.reader(f, delimiter='\t'):
                nin += 1
                ninc += 1
                chrom = fields[0]
                if chrom != curr:
                    curr  = chrom
                    currRegs = self.regions[chrom]
                    if noutc > 0:
                        sys.stderr.write("{}: {}/{} common peaks.\n".format(chrom, noutc, ninc))
                    ninc = 0
                    noutc = 0
                start = int(fields[1])
                end   = int(fields[2])
                reg = self.findRegion(currRegs, start, end)
                if reg:
                    ov = self.overlap(reg, start, end)
                    if ov > self.minOverlap:
                        x = float(fields[5]) * self.factor
                        if x == 0:
                            fc = 0
                        else:
                            fc = math.log(reg[2] / x, 2)
                            if fc > 0:
                                nup += 1
                            else:
                                ndn += 1
                        out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, reg[0], reg[1], reg[2], start, end, x, fc))
                        nout += 1
                        noutc += 1
            sys.stderr.write("Peaks 2: {}\n".format(nin))
            sys.stderr.write("Common : {}\n".format(nout))
        return (nup, ndn)

    def doOutput(self):
        if self.outfile:
            with open(self.outfile, "w") as out:
                return self.readSecond(out)
        else:
            return self.readSecond(sys.stdout)

if __name__ == "__main__":
    global P
    P = macsDiff("peakscmp", version="1.0", usage=usage,
                  errors=[("BADOVERLAP", "Bad overlap value", "The value for -v ({}) should be between 0 and 1.")])
    P.parseArgs(sys.argv[1:])
    n1 = P.readFirst()
    sys.stderr.write("Peaks 1: {}\n".format(n1))
    (nup, ndn) = P.doOutput()
    sys.stderr.write("Up  : {}\nDown: {}\n".format(nup, ndn))
