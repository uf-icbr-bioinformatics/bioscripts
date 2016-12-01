#!/usr/bin/env python

import sys

class Cov():
    chrom = ""
    mincov = 0
    total = 0
    maxpos = 0
    effbases = 0

    def dump(self):
        print "total={}, maxpos={}, effbases={}".format(self.total, self.maxpos, self.effbases)

    def report(self):
        if self.chrom != "":
            sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.chrom, self.total, self.maxpos, 1.0*self.total/self.maxpos, self.effbases, round(100.0*self.effbases/self.maxpos, 1), 1.0*self.total/self.effbases))
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
        
def parseArgs(C, args):
    next = ""
    filename = None

    for a in args:
        if next == '-m':
            C.mincov = int(a)
        elif a == '-h':
            sys.stderr.write("""Usage: chromCoverage.py [-h] [-m min] [coveragefile]
       chromCoverage.py -c combinedfile coveragefiles...

Read coverage data from standard input (or from coveragefile if provided)
and write by-chromosome coverage data to standard output. The input file 
should be in the format produced by the bamtools 'coverage' command:

  chrom   position   depth

The output file contains six columns:

  chrom  total  length  coverage  efflen  effcov

total    - sum of depth at all positions
length   - chromosome length (highest position)
coverage - ratio between total and length
efflen   - number of bases having non-zero depth
effperc  - percent of bases having non-zero depth
effcov   - ratio between total and efflen

If -m is specified, only positions with depth over 'min' are considered
for the computation of 'total' and 'efflen'.

With -c, combine multiple 'coveragefiles' into a single 'combinedfile'.

The -h option prints this usage message.
""")
            exit(-1)
        elif a == '-m':
            next = a
        else:
            filename = a
    return filename

if __name__=="__main__":

    C = Cov()
    Tot = Cov()
    Tot.chrom = "Total"
    filename = parseArgs(C, sys.argv[1:])

    if filename == None:
        stream = sys.stdin
    else:
        stream = open(filename, "r")

    sys.stdout.write("Chrom\tTotal\tLength\tCoverage\tEfflen\tEffperc\tEffcov\n")
    while True:
        line = stream.readline()
        if not line:
            break
        parsed = line.split("\t")
        if len(parsed) > 2:
            chrom = parsed[0]
            pos = int(parsed[1])
            cov = int(parsed[2])
            if cov > C.mincov:
                C.add(chrom, pos, cov, Tot)
    C.update(Tot)
    C.report()
    Tot.report()
    stream.close()

