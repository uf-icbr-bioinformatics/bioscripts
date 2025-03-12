#!/usr/bin/env python

import sys
import csv

class Peak(object):
    chrom = ""
    start = 0
    end = 0
    weight = 0
    summit = 0

    def write(self):
        if self.chrom:
            sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}:{}-{}\n".format(self.chrom, self.start, self.end, self.weight, self.summit, self.chrom, self.start, self.end))

def main():
    P = Peak()
    c = csv.reader(sys.stdin, delimiter='\t')
    for row in c:
        chrom = row[0]

        astart = int(row[1])
        aend   = int(row[2])
        bstart = int(row[7])
        bend   = int(row[8])
        start  = max(astart, bstart)
        end    = min(aend, bend)
        aweight = float(row[3]) * (end-start) / (aend-astart)
        bweight = float(row[3]) * (end-start) / (bend-bstart)
        weight  = (aweight + bweight) / 2
        summit  = (float(row[4]) + float(row[10])) / 2

        if chrom == P.chrom and start == P.start and end == P.end:
            P.weight += weight
            P.summit = max(P.summit, summit)
        else:
            try:
                P.write()
            except BrokenPipeError:
                break
            P.chrom = chrom
            P.start = start
            P.end   = end
            P.weight = weight
            P.summit = summit

main()
