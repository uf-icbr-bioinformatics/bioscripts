#!/usr/bin/env python

import sys

class PTB():
    # Streams
    infile = None
    outfile = None
    bedfile = None
    repfile = None
    instream = None
    outstream = None
    bedstream = None
    repstream = None

    # Output control
    trackName = "coverage"
    maxGap = 10
    minCov = 0
    minScore = 0.0

    # Runtime vars
    currentChrom = ""
    startPos = 0                # Start of current region
    lastPos = 0                 # Last pos seen in current region
    prevEnd = 0                 # End of previous region
    coverage = 0                # Incremental coverage in current region
    covered = 0                 # Size of covered regions in this chromosome
    uncovered = 0               # Size of uncovered regions in this chromosome
    totCovered = 0
    totUncovered = 0

    def init(self, args):
        next = ""
        for a in args:
            if a == '-h':
                return self.usage()
            elif next == '-i':
                self.infile = a
                next = ""
            elif next == '-o':
                self.outfile = a
                next = ""
            elif next == '-r':
                self.repfile = a
                next = ""
            elif next == '-bf':
                self.bedfile = a
                next = ""
            elif next == '-bn':
                self.trackName = a
                next = ""
            elif next == '-g':
                self.maxGap = int(a)
                next = ""
            elif next == '-s':
                self.minScore = float(a)
                next = ""
            elif next == '-c':
                self.minCov = int(a)
                next = ""
            elif a in ['-i', '-o', '-r', '-bf', '-bn', '-g', '-s', '-c']:
                next = a
            else:
                sys.stderr.write("Unknown option '{}'.\n".format(a))
                return self.usage()
        if self.infile:
            self.instream = open(self.infile, "r")
        else:
            self.instream = sys.stdin
        if self.outfile:
            self.outstream = open(self.outfile, "w")
        else:
            self.outstream = sys.stdout
        if self.bedfile:
            self.bedstream = open(self.bedfile, "w")
        if self.repfile:
            self.repstream = open(self.repfile, "w")
        else:
            self.repstream = sys.stderr

    def usage(self):
        sys.stderr.write("""pileupToBED.py - Convert a samtools pileup to a BED file.

Usage: pileupToBED.py [options]

Options:

-i infile     | samtools pileup file (default: stdin)
-o outfile    | Filename for output in tab-delimited format (default: stdout)
-bf bedfile   | Filename for output in BED format
-bn name      | Track name for BED file
-r reportfile | Filename for by-chromosome report
-g gap        | Set gap between regions (default: {})
-s score      | Specify minimum score (default: {})
-c cov        | Specify minimum coverage (default: {})

(c) 2016, A. Riva, DiBiG, ICBR Bioinformatics, University of Florida
""".format(self.maxGap, self.minScore, self.minCov))
        sys.exit(-1)

    def close(self):
        if self.infile:
            self.instream.close()
        if self.outfile:
            self.outstream.close()
        if self.bedfile:
            self.bedstream.close()
        if self.repfile:
            self.repstream.close()

    def writeRegion(self):
        if self.currentChrom != "":
            l = self.lastPos - self.startPos + 1
            s = 1.0*self.coverage/l
            self.covered += l
            self.uncovered += (self.startPos - self.prevEnd - 1)
            if s > self.minScore:
                if self.bedfile:
                    self.bedstream.write("{}\t{}\t{}\t{}\n".format(self.currentChrom, self.startPos, self.lastPos, s))
                if self.outfile:
                    self.outstream.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(self.currentChrom, self.startPos, self.lastPos, l, self.coverage, s))

    def newChrom(self, chrom, pos, cov):
        # print "Found chrom {}".format(chrom)
        if self.repfile and self.currentChrom!= "":
            covpct = (100.0 * self.covered) / (self.covered + self.uncovered)
            uncpct = (100.0 * self.uncovered) / (self.covered + self.uncovered)
            self.repstream.write("{}\t{}\t{:.2f}\t{}\t{:.2f}\n".format(self.currentChrom, self.covered, covpct, self.uncovered, uncpct))
        if chrom:
            self.currentChrom = chrom
            self.startPos = pos
            self.lastPos = pos
            self.prevEnd = 0
            self.coverage = cov
            self.totCovered += self.covered
            self.totUncovered += self.uncovered
            self.covered = 0
            self.uncovered = 0

    def newRegion(self, pos, cov):
        self.prevEnd = self.lastPos + 1
        self.startPos = pos
        self.lastPos = pos
        self.coverage = cov

    def main(self):
        instream = self.instream
        outstream = self.outstream
        if self.bedfile:
            self.bedstream.write("track type=bedGraph name={}\n".format(self.trackName))
        if self.outfile:
            outstream.write("Chrom\tStart\tEnd\tLength\tCoverage\tAvgCov\n")

        while True:
            line = instream.readline()
            if line == '':
                break
            parsed = line.rstrip("\r\n").split("\t")
            chrom = parsed[0]
            pos = int(parsed[1])
            cov = int(parsed[3])
            if cov > self.minCov:
                if chrom != self.currentChrom:             # New chromosome?
                    self.writeRegion()
                    self.newChrom(chrom, pos, cov)
                elif pos > self.lastPos + self.maxGap:     # Same chromosome but too far away?
                    self.writeRegion()
                    self.newRegion(pos, cov)
                else:                                      # Still in region
                    self.lastPos = pos
                    self.coverage += cov
        self.writeRegion()
        self.newChrom(False, 0, 0) # To write final row in report

if __name__ == "__main__":
    P = PTB()
    P.init(sys.argv[1:])
    try:
        P.main()
    finally:
        P.close()