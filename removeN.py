#!/usr/bin/env python

import sys

### Error messages

class Err():
    HELP = 1
    NOFILE = 2

    MSGS = {1: "Help requested.",
            2: "Input file not specified."}

def errmsg(code, exit=False):
    sys.stderr.write("Error: " + Err.MSGS[code] + "\n")
    if exit:
        sys.exit(code)

MAXL = 70

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

def usage():
    sys.stderr.write("""removeN.py - remove Ns from sequences in FASTA file

Usage: removeN.py [-l len] infile [outfile]

Remove all occurrences of N or n from the sequences in input multi-FASTA 
file `infile'. Output is written to standard output or to `outfile' if 
specified, in FASTA format.

Options:

 -l len | set line length in output file to `len' (default: {})

(c) 2016, A. Riva, DiBiG, ICBR Bioinformatics, University of Florida
""".format(MAXL))
    sys.exit(Err.HELP)

def parseArgs(args):
    infile = None
    outfile = None
    next = ""
    global MAXL

    if '-h' in args:
        usage()
    for a in args:
        if next == "-E":
            errmsg(int(a), True)
        elif next == "-l":
            MAX = int(a)
            next = ""
        elif a in ["-E", "-l"]:
            next = a
        elif infile == None:
            infile = a
        else:
            outfile = a
    if infile == None:
        errmsg(Err.NOFILE, True)
    return (infile, outfile)

if __name__ == "__main__":
    (infile, outfile) = parseArgs(sys.argv[1:])
    with open(infile, "r") as f:
        if outfile:
            with open(sys.argv[2], "w") as out:
                main(f, out)
        else:
            main(f, sys.stdout)
