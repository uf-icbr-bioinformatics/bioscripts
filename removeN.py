#!/usr/bin/env python

import sys

MAXL = 70
isfirst = True

def main(instream, out):
    nout = 0
    for line in instream:
        if len(line) > 0 and line[0] == '>':
            if 0 < nout < MAXL:
                out.write("\n")
            out.write(line)
            nout = 0
        else:
            for c in line:
                if c not in ['N', 'n', '\r', '\n']:
                    out.write(c)
                    nout += 1
                    if nout == MAXL:
                        out.write('\n')
                        nout = 0

if __name__ == "__main__":
    nargs = len(sys.argv) - 1
    if nargs == 0:
        sys.stderr.write("Usage: removeN.py infile [outfile]\n")
        exit(1)
    infile = sys.argv[1]
    with open(infile, "r") as f:
        if nargs == 1:
            main(f, sys.stdout)
        elif nargs == 2:
            with open(sys.argv[2], "w") as out:
                main(f, out)
