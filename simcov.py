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
        for i in range(size):
            x = vector[i]
            if x > 0:
                self.numCovBases += 1
                if x > self.maxCov:
                    self.maxCov = x
        self.counts[0][0] = np.around(self.maxCov * 0.25)
        self.counts[1][0] = np.around(self.maxCov * 0.50)
        self.counts[2][0] = np.around(self.maxCov * 0.75)
        print self.numCovBases

    def pass2(self, vector, size):
        for i in range(size):
            x = vector[i]
            for c in self.counts:
                if x >= c[0]:
                    c[1] += 1
        print self.counts
        
class Simcov(Script):
    readlen = 150
    paired = False
    insertSize = 400
    nreads = 10000000
    
    genomeSize = 3100000000
    vectorSize = 1000000
    vector = None

    def run(self):
        self.vector = np.zeros(self.vectorSize, dtype=int)
        scale = 1.0 * self.vectorSize / self.genomeSize
        effReads = int(np.around(self.nreads * scale))
        print effReads
        if self.paired:
            self.simulatePaired(effReads)
        else:
            self.simulateUnpaired(effReads)
        print np.max(self.vector)
        CS = CovStats()
        CS.pass1(self.vector, self.vectorSize)
        CS.pass2(self.vector, self.vectorSize)
        #print self.vector[:100]

    def simulateUnpaired(self, nreads):
        low = 1 - self.readlen
        high = self.vectorSize
        for i in range(nreads):
            pos = np.random.randint(low, high)
            start = max(pos, 0)
            end = min(pos + self.readlen, self.vectorSize)
            for p in range(start, end):
                self.vector[p] += 1

if __name__ == "__main__":
    S = Simcov("simcov.py")
    S.run()
    
