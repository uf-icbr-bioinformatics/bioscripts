#!/usr/bin/env python

import sys

import Script

def usage():
    sys.stderr.write("TBW\n")

P = Script.Script("inserts.py", version="1.0", usage=usage)

class GTFregions():
    entries = {}
    currentChrom = ""
    currentRegions = []
    span = 1000

    def __init__(self, span=1000):
        self.entries = {}
        self.currentChrom = ""
        self.currentRegions = []
        self.span = span

    def addRegion(self, chrom, position):
        if chrom != self.currentChrom:
            if self.currentRegions != []:
                self.entries[self.currentChrom] = self.currentRegions
            self.currentChrom = chrom
            if chrom in self.entries:
                self.currentRegions = self.entries[self.currentChrom]
            else:
                self.currentRegions = []
        self.currentRegions.append((position - self.span, position + self.span))

    def doneAdding(self):
        if self.currentRegions != []:
            self.entries[self.currentChrom] = self.currentRegions
        for (chrom, regions) in self.entries.iteritems():
            regions.sort(key=lambda s: s[0])
            print "{}: {} regions.".format(chrom, len(regions))

    def parseGTF(self, filename):
        with open(filename, "r") as f:
            for line in f:
                if not line.startswith('#'):
                    fields = line.split('\t')
                    if fields[2] == 'transcript':
                        chrom = fields[0]
                        if fields[6] == '+':
                            pos = int(fields[3])
                        else:
                            pos = int(fields[4])
                        self.addRegion(chrom, pos)
        self.doneAdding()

    def posInRegion(self, chrom, pos):
        regions = self.entries[chrom]
        for r in regions:
            if r[0] > pos:
                return False
            elif r[0] <= pos <= r[1]:
                return True

class Params():
    maxsize = 1000
    outfile = None
    gtffile = None
    regsize = 1000
    gtfregions = None

    def __init__(self, args):
        P.standardOpts(args)
        next = ""
        for a in args:
            if next == '-o':
                self.outfile = a
                next = ""
            elif next == '-gtf':
                self.gtffile = P.isFile(a)
                next = ""
            elif next == '-s':
                self.maxsize = P.toInt(a)
                next = ""
            elif next == '-r':
                self.regsize = P.toInt(a)
                next = ""
            elif a in ['-o', '-gtf', '-s', '-r']:
                next = a
        if self.gtffile:
            G = GTFregions()
            G.parseGTF(self.gtffile)
            self.gtfregions = G

def main(PA):
    totlines = 0
    data = [0]*PA.maxsize
    G = PA.gtfregions

    while True:
        line = sys.stdin.readline()
        if line == '':
            break
        totlines += 1
        fields = line.rstrip("\r\n").split("\t")
        if fields[6] == '=' and fields[8][0] != '-':
            c = int(fields[8])
            if c < PA.maxsize:
                good = True
                if G:
                    chrom = fields[2]
                    pos = int(fields[3])
                    good = G.posInRegion(chrom, pos)
                if good:
                    data[c] += 1

    if PA.outfile:
        out = open(PA.outfile, "w")
    else:
        out = sys.stdout
    try:
        for i in range(PA.maxsize):
            out.write("{}\t{}\t{}\n".format(i, data[i], 1.0*data[i]/totlines))
    finally:
        if PA.outfile:
            out.close()

if __name__ == "__main__":
    PA = Params(sys.argv[1:])
    main(PA)

