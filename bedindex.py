#!/usr/bin/env python

import sys
import csv
import os.path

def bidx_filename(filename):
    return os.path.splitext(filename)[0] + ".bidx"

def loadindex(filename):
    bidx = bidx_filename(filename)
    if os.path.isfile(bidx):
        with open(bidx, "r") as f:
            data = f.read().split("\n")
        idx = {}
        for d in data:
            if "\t" in d:
                [chrom, fp] = d.split("\t")
                idx[chrom] = int(fp)
        return idx
    else:
        sys.stderr.write("Error: index file {} not found.\n".format(bidx))
        return None

def bedindex(filename):
    bidx = bidx_filename(filename)
    tag = ""
    sys.stderr.write("{} => {}\n".format(filename, bidx))
    with open(bidx, "w") as out:
        with open(filename, "r") as f:
            while True:
                pos = f.tell()
                line = f.readline()
                if not line:
                    break
                if line[0] == '#':
                    continue
                tp = line.find("\t")
                if tp < 0:
                    continue
                v = line[:tp]
                if v != tag:
                    out.write("{}\t{}\n".format(v, pos))
                    tag = v

def usage():
    sys.stderr.write("""Usage: bedindex.py file.bed

Write an index of the specified BED file to standard output.

""")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        for f in sys.argv[1:]:
            bedindex(f)
    else:
        usage()

