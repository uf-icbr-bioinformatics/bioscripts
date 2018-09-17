#!/usr/bin/env python

import sys
import random
import numpy as np
from Script import Script

class CovStats():
    maxCov = 0
    numCovBases = 0
    counts = []

    def __init__(self):
        self.counts = [ [0, 0],   # 25th percentile
                        [0, 0],   # 50th percentile
                        [0, 0],    # 75th percentile
                        [1, 0],
                        [5, 0],
                        [10, 0],
                        [20, 0],
                        [30, 0],
                        [50, 0],
                        [100, 0] ]

    def pass1(self, vector, size):
        self.vectorSize = size
        for i in range(size):
            x = vector[i]
            if x > 0:
                self.numCovBases += 1
                if x > self.maxCov:
                    self.maxCov = x
        self.counts[0][0] = int(np.around(self.maxCov * 0.25))
        self.counts[1][0] = int(np.around(self.maxCov * 0.50))
        self.counts[2][0] = int(np.around(self.maxCov * 0.75))

    def pass2(self, vector, size):
        for i in range(size):
            x = vector[i]
            for c in self.counts:
                if x >= c[0]:
                    c[1] += 1

    def report(self, scale, genomesize):
        for i in range(3, len(self.counts)):
            bc = self.counts[i][1] / scale
            pct = bc / genomesize
            sys.stdout.write("Bases covered at {:3}X: {:,}bp ({:.1f}%)\n".format(self.counts[i][0], int(bc), 100.0 * pct))
            
        bc = self.counts[0][1] / scale
        pct = bc / genomesize
        sys.stdout.write("25th percentile - bases covered at {}X: {:,}bp ({:.1f}%)\n".format(self.counts[0][0], int(bc), 100.0 * pct))
        bc = self.counts[1][1] / scale
        pct = bc / genomesize
        sys.stdout.write("50th percentile - bases covered at {}X: {:,}bp ({:.1f}%)\n".format(self.counts[1][0], int(bc), 100.0 * pct))
        bc = self.counts[2][1] / scale
        pct = bc / genomesize
        sys.stdout.write("75th percentile - bases covered at {}X: {:,}bp ({:.1f}%)\n".format(self.counts[2][0], int(bc), 100.0 * pct))

class SNPsite():
    all1freq = 0.0
    coverage = 0
    all1count = 0
    all2count = 0

    def __init__(self):
        self.all1freq = random.random()

    def seen(self):
        self.coverage += 1
        if random.random() <= self.all1freq:
            self.all1count += 1
        else:
            self.all2count += 1

    def sqerror(self):
        x = self.all1freq - 1.0 * self.all1count / self.coverage
        return x * x
        
class SNPstats():
    positions = []
    snps = {}
    counts = []
    ndetected = 0
    
    def __init__(self, size, sitefreq):
        nsites = int(sitefreq * size)
        sys.stderr.write("Simulating {} SNP positions.\n".format(nsites))
        self.positions = np.random.choice(size, nsites, replace=False)
        self.positions.sort()
        for p in self.positions:
            self.snps[p] = SNPsite()
        self.counts = [ [1, 0],
                        [5, 0],
                        [10, 0],
                        [20, 0],
                        [30, 0],
                        [50, 0],
                        [100, 0] ]

    def addread(self, start, end):
        for i in range(start, end):
            if i in self.snps:
                self.snps[i].seen()

    def snpCoverage(self):
        for p in self.positions:
            snp = self.snps[p]
            x = snp.coverage
            if x:
                self.ndetected += 1
                for c in self.counts:
                    if x >= c[0]:
                        c[1] += 1

    def avgError(self):
        errs = []
        for p in self.positions:
            snp = self.snps[p]
            if snp.coverage > 0:
                errs.append(snp.sqerror())
        if errs:
            return np.mean(errs)
        else:
            return None
                        
    def report(self):
        self.snpCoverage()
        sys.stdout.write("\n=== SNPs ===\n")
        nsnps = len(self.positions)
        sys.stdout.write("Simulated: {}\n".format(nsnps))
        sys.stdout.write("Detected: {} ({:.1f}%)\n".format(self.ndetected, 100.0 * self.ndetected / nsnps))
        for c in self.counts:
            sys.stdout.write("  at {}X: {} ({:.1f}%)\n".format(c[0], c[1], 100.0 * c[1] / nsnps))
        sys.stdout.write("Average squared error: {}\n".format(self.avgError()))
        sys.stdout.write("\n")

def usage():
    sys.stderr.write("""simcov.py - Simulate short read coverage.

Usage: simcov.py [options]

Options:

  -l L | Set read length to L (defalt: {})
  -n N | Set number of reads to N (default: {})
  -g G | Set target genome size to G (default: {})
  -p   | Enable paired-end mode (default: single-end)
  -t T | Set insert size to T in paired-end mode (default: {})
  -f F | Set site frequency to F (default: no site analysis)
  -s S | Set simulation vector size to S (default: {})

Values for -l, -n, -g and -t can be followed by G (for billion) or M (for million).
The value for -f can be expressed as a float or a fraction (e.g. 1/8)

""".format(Simcov.readlen, Simcov.nreads, Simcov.genomeSize, Simcov.insertSize, Simcov.vectorSize))
    
class Simcov(Script):
    readlen = 150
    paired = False
    insertSize = 400
    nreads = 10000000
    
    genomeSize = 3100000000
    vectorSize = 1000000
    vector = None
    scale = 1.0

    snpfreq = None
    snpstats = None
    
    def parseArgs(self, args):
        self.standardOpts(args)
        prev = ""
        for a in args:
            if prev == "-l":
                self.readlen = self.toInt(a)
                prev = ""
            elif prev == "-t":
                self.insertSize = self.toInt(a)
                prev = ""
            elif prev == "-n":
                self.nreads = self.toInt(a, units=True)
                prev = ""
            elif prev == "-g":
                self.genomeSize = self.toInt(a, units=True)
                prev = ""
            elif prev == "-s":
                self.vectorSize = self.toInt(a, units=True)
                prev = ""
            elif prev == "-f":
                self.snpfreq = self.toFloat(a)
                prev = ""
            elif a in ["-l", "-t", "-n", "-g", "-s", "-f"]:
                prev = a
            elif a == "-p":
                self.paired = True
                
    def run(self):
        self.vector = np.zeros(self.vectorSize, dtype=int)
        self.scale = 1.0 * self.vectorSize / self.genomeSize
        effReads = int(np.around(self.nreads * self.scale))
        self.writeConfiguration()
        if self.snpfreq:
            self.snpstats = SNPstats(self.vectorSize, self.snpfreq)
        if self.paired:
            self.simulatePaired(effReads)
        else:
            self.simulateUnpaired(effReads)
        CS = CovStats()
        CS.pass1(self.vector, self.vectorSize)
        CS.pass2(self.vector, self.vectorSize)
        self.report(CS)

    def simulateUnpaired(self, nreads):
        low = 1 - self.readlen
        high = self.vectorSize
        for i in range(nreads):
            pos = np.random.randint(low, high)
            start = max(pos, 0)
            end = min(pos + self.readlen, self.vectorSize)
            for p in range(start, end):
                self.vector[p] += 1
            if self.snpstats:
                self.snpstats.addread(start, end)

    def simulatePaired(self, nreads):
        low = 1 - self.insertSize
        high = self.vectorSize
        for i in range(nreads):
            pos = np.random.randint(low, high) # position of insert

            # left read
            start = max(pos - self.readlen, 0)
            end   = max(pos, 0)
            for p in range(start, end):
                self.vector[p] += 1
            if self.snpstats:
                self.snpstats.addread(start, end)

            # right read
            start = min(pos + self.insertSize, self.vectorSize)
            end   = min(start + self.readlen, self.vectorSize)
            for p in range(start, end):
                self.vector[p] += 1
            if self.snpstats:
                self.snpstats.addread(start, end)
            
    def writeConfiguration(self):
        sys.stdout.write("=== Configuration ===\n")
        sys.stdout.write("Genome size: {:,}bp\n".format(self.genomeSize))
        sys.stdout.write("Number of reads: {:,}\n".format(self.nreads))
        sys.stdout.write("Read length: {:,}bp\n".format(self.readlen))
        sys.stdout.write("Mode: {}\n".format("paired" if self.paired else "unpaired"))
        if self.paired:
            nbp = self.nreads * 2 * self.readlen
            sys.stdout.write("Insert size: {:,}bp\n".format(self.insertSize))
        else:
            nbp = self.nreads * self.readlen
        ecov = 1.0 * nbp / self.genomeSize
        sys.stdout.write("Expected average coverage: {:.1f}X\n".format(ecov))
        
    def report(self, CS):
        rcov = 1.0 * np.sum(self.vector) / self.vectorSize
        sys.stdout.write("\n=== Results ===\n")
        sys.stdout.write("Effective average coverage: {:.1f}X\n".format(rcov))
        sys.stdout.write("Max coverage: {:,}\n".format(CS.maxCov))
        pctcov = 1.0 * CS.numCovBases / self.vectorSize
        #sys.stdout.write("Genome covered at 1X: {:,}bp ({:.1f}%)\n".format(int(self.genomeSize * pctcov), pctcov * 100.0))
        CS.report(self.scale, self.genomeSize)

        if self.snpfreq:
            self.snpstats.report()
            #SS = SNPstats(self.vectorSize, self.snpfreq)
            #SS.pass1(self.vector)
            #print SS.counts
        
if __name__ == "__main__":
    S = Simcov("simcov.py", version="1.0", usage=usage)
    S.parseArgs(sys.argv[1:])
    S.run()
    
