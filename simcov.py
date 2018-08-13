#!/usr/bin/env python

import sys
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
        self.counts[0][0] = np.around(self.maxCov * 0.25)
        self.counts[1][0] = np.around(self.maxCov * 0.50)
        self.counts[2][0] = np.around(self.maxCov * 0.75)

    def pass2(self, vector, size):
        for i in range(size):
            x = vector[i]
            for c in self.counts:
                if x >= c[0]:
                    c[1] += 1

    def report(self, scale, genomesize):
        for i in range(3, 10):
            bc = self.counts[i][1] / scale
            pct = bc / genomesize
            sys.stdout.write("Bases covered at {}X: {:,}bp ({:.1f}%)\n".format(self.counts[i][0], int(bc), 100.0 * pct))

class Simcov(Script):
    readlen = 150
    paired = False
    insertSize = 400
    nreads = 50000000
    
    genomeSize = 3100000000
    vectorSize = 1000000
    vector = None
    scale = 1.0

    def run(self):
        self.vector = np.zeros(self.vectorSize, dtype=int)
        self.scale = 1.0 * self.vectorSize / self.genomeSize
        effReads = int(np.around(self.nreads * self.scale))
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

    def report(self, CS):
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

        rcov = 1.0 * np.sum(self.vector) / self.vectorSize
        sys.stdout.write("\n=== Results ===\n")
        sys.stdout.write("Effective average coverage: {:.1f}X\n".format(rcov))
        sys.stdout.write("Max coverage: {:,}\n".format(CS.maxCov))
        pctcov = 1.0 * CS.numCovBases / self.vectorSize
        sys.stdout.write("Genome covered at 1X: {:,}bp ({:.1f}%)\n".format(int(self.genomeSize * pctcov), pctcov * 100.0))
        CS.report(self.scale, self.genomeSize)

if __name__ == "__main__":
    S = Simcov("simcov.py")
    S.run()
    
