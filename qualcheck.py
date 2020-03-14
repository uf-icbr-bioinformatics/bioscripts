#!/usr/bin/env python

import sys, gzip
from math import pow

def decode(q):
    return ord(q) - 33

class Stats(object):
    nreads = 0
    nbases = 0                  # Total bases seen
    sumqual = 0                 # Sum of qualities (for average)
    minqual = 100               # Minumum quality seen
    maxqual = 0                 # Maximum quality seen
    minreadqual = 100           # Minimum per-read quality seen
    maxreadqual = 0             # Maximum per-read quality seen

    def add(self, quals):
        nb = 0                  # number of bases in this read
        sq = 0                  # sum of qualities in this read

        self.nreads += 1
        for q in quals:
            qv = decode(q)
            self.nbases += 1
            self.sumqual += qv
            self.minqual = min(self.minqual, qv)
            self.maxqual = max(self.maxqual, qv)
            nb += 1
            sq += qv
        readqual = 1.0 * sq / nb
        self.minreadqual = min(self.minreadqual, readqual)
        self.maxreadqual = max(self.maxreadqual, readqual)

    def report(self):
        avgqual = 1.0 * self.sumqual / self.nbases
        err_rate = pow(10, -avgqual/10.0)
        out = sys.stdout
        out.write("Number of reads: {}\n".format(self.nreads))
        out.write("Number of bases: {}\n".format(self.nbases))
        out.write("Average read length: {:.3f}\n".format(1.0 * self.nbases / self.nreads))
        out.write("Overall average quality: {:.3f}\n".format(avgqual))
        out.write("Quality range: {} - {}\n".format(self.minqual, self.maxqual))
        out.write("Per-read quality range: {} - {}\n".format(self.minreadqual, self.maxreadqual))
        out.write("Average error rate: {}\n".format(err_rate))
        out.write("Expected number of errors: {:.2f}\n".format(self.nbases * err_rate))

def collect(fastqfile, maxreads=0):
    nr = 0
    S = Stats()
    with gzip.open(fastqfile, "r") as f:
        while True:
            try:
                read = f.readline()
                if not read:
                    break
                read = f.readline()
                read = f.readline()
                quals = f.readline().strip("\r\n")
                S.add(quals)
                nr += 1
                if nr == maxreads:
                    break
            except KeyboardInterrupt:
                break
    S.report()

def main(args):
    fastqfilea = []
    maxreads = 0
    prev = ""
    for a in args:
        if prev == "-n":
            maxreads = int(a)
            prev = ""
        elif a in ["-n"]:
            prev = a
        else:
            fastqfiles.append(a)

    if fastqfiles:
        for fq in fastqfiles:
            collect(fq, maxreads=maxreads)
    else:
        usage()

def usage():
    sys.stdout.write("""qualcheck.py - Compute quality stats on fastq files.

Usage: qualcheck.py [options] fastqfiles...

Options:

  -n N | Examine the first N reads

This program reads one or more fastq files, printing out the following information for each:

Number of reads
Number of bases
Average read length
Overall average quality
Range of quality scores
Per-read quality range
Average error rate
Expected number of errors

"Per-read quality" is the average quality score over the whole read.
"Average error rate" is the probability of error corresponding to the
overall average quality Q (= 10^(-Q/10)).

Processing can be interrupted with Ctrl-C, and the program will print
the statistics computed up to that point.

"""

if __name__ == "__main__":
    args = sys.argv[1:]
    main(args)
