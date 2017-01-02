#!/usr/bin/env python

import sys
import Script

MAXL = 70

def usage():
    sys.stderr.write("""removeN.py - Remove Ns from sequences in FASTA file.

Usage: removeN.py [-l len] infile [outfile]

Remove all occurrences of N or n from the sequences in input multi-FASTA 
file `infile'. Output is written to standard output or to `outfile' if 
specified, in FASTA format.

Options:

 -l len | set line length in output file to `len' (default: {})

""".format(MAXL))

### Program object

P = Script.Script("removeN", version="1.0", usage=usage)

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

def parseArgs(args):
    infile = None
    outfile = None
    next = ""
    global MAXL

    P.standardOpts(args)
    for a in args:
        if next == "-l":
            MAX = P.toInt(a)
            next = ""
        elif a in ["-l"]:
            next = a
        elif infile == None:
            infile = a
        else:
            outfile = a
    if infile == None:
        P.errmsg(P.NOFILE)
    return (infile, outfile)

if __name__ == "__main__":
    (infile, outfile) = parseArgs(sys.argv[1:])
    try:
        with open(infile, "r") as f:
            if outfile:
                with open(sys.argv[2], "w") as out:
                    main(f, out)
            else:
                main(f, sys.stdout)
    except Exception as e:
        P.errmsg(P.ERR, e)
