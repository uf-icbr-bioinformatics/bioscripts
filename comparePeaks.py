#!/usr/bin/env python

import sys

import Utils

def readPeaks(filename):
    d = {}
    nread = 0
    with open(filename, "r") as f:
        for line in f:
            if line[0] == '#':
                continue
            parsed = line.split("\t")
            if len(line) < 4:
                continue
            chrom = parsed[0]
            if chrom not in d:
                d[chrom] = []
            start = Utils.safeInt(parsed[1])
            if start:
                d[chrom].append((start, Utils.safeInt(parsed[2]), Utils.safeInt(parsed[3])))
            nread += 1
    sys.stderr.write("{}: {} peaks in {} chromosomes.\n".format(filename, nread, len(d)))
    return d

def findMatchingPeak(peaks, start, end, size, over=0.5):
    """Find a peak in `peaks' that overlaps (`start', `end') by at least `over' (as a fraction of `size')."""
    for p in peaks:
        if p[0] > end:
            return None
        if start <= p[0] < p[1] <= end:
            return ("I", p)     # Increase
        elif p[0] <= start < end <= p[1]:
            return ("D", p)     # Decrease
        elif start <= p[0] < end <= p[1]:
            ov = 1.0 * (end - p[0]) / size
            if ov >= over:
                return ("O", p, ov)
        elif p[0] <= start < p[1] <= end:
            ov = 1.0 * (p[1] - start) / size
            if ov >= over:
                return ("O", p, ov)
    return None

def comparePeaks(filename1, filename2, out):
    dict1 = readPeaks(filename1)
    dict2 = readPeaks(filename2)
    chroms = dict1.keys()
    chroms.sort()

    nI = 0
    nD = 0
    nO = 0

    for chrom in chroms:
        l1 = dict1[chrom]
        l2 = dict2[chrom]
        for p2 in l2:
            res = findMatchingPeak(l1, p2[0], p2[1], p2[2])
            if res:
                key = res[0]
                p1 = res[1]
                if key == 'I':
                    nI += 1
                    out.write("I\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, p1[0], p1[1], p2[0], p2[1], p2[2]-p1[2]))
                elif key == 'D':
                    nD += 1
                    out.write("D\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, p1[0], p1[1], p2[0], p2[1], p2[2]-p1[2]))
                elif key == 'O':
                    nO += 1
                    out.write("O\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, p1[0], p1[1], p2[0], p2[1], p2[2]-p1[2]))

    sys.stderr.write("Increased: {}\nDecreased: {}\nShifted: {}\n".format(nI, nD, nO))
            

if __name__ == "__main__":
    args = sys.argv
    comparePeaks(args[1], args[2], sys.stdout)
