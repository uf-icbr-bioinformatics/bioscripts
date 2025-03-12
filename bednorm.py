#!/usr/bin/env python

import sys
import csv

def main(args):
    factors = []
    scale = 1
    col = 3

    prev = ""
    for a in args:
        if prev == "-s":
            scale = int(a)
            prev = ""
        elif prev == "-c":
            col = int(a) - 1
            prev = ""
        elif a in ["-s", "-c"]:
            prev = a
        else:
            factors.append(int(a))

    nf = len(factors)
    c = csv.reader(sys.stdin, delimiter='\t')
    for row in c:
        for i in range(nf):
            row[i+col] = int(scale * int(row[i+col]) / factors[i])
        sys.stdout.write("\t".join([str(r) for r in row]) + "\n")

if __name__ == "__main__":
  args = sys.argv[1:]
  main(args)

