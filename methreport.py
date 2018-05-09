#!/usr/bin/env python

## (c) 2014, Alberto Riva (ariva@ufl.edu)
## DiBiG, ICBR Bioinformatics, University of Florida

import sys
from Bio import SeqIO
import Script

# Script object

def usage():
    sys.stderr.write("""methreport.py - Report methylation rate at CG and GC positions.

Usage: methreport.py [-gcg] infile [outfile]

`Infile' should be a multi-FASTA file in which the first sequence is assumed to
be the reference. All other sequences should have the same length as the reference
and be aligned to it. This program will identify all CG and GC positions in the 
reference sequence (including GCG positions if the -gcg option is specified) and 
report to the output the nummber and fraction of unconverted Cs at each position.

Output is written to `outfile', is specified, or to standard output. Entries for
CG positions are written first, followed by those for GC positions. The output is
in tab-delimited format with three columns: position, number of unconverted Cs, 
fraction of unconverted Cs (over total number of sequences examined).

Options:
 -gcg | Do not exclude GCG positions from analysis.

""")

P = Script.Script("methreport.py", version="1.0", usage=usage)

EXCLGCG=False

# Utility classes

class refDesc():
    """A class containing the reference sequence, its length, and a list of CG and GC positions."""
    sequence = None
    length = 0
    CGpositions = []
    GCpositions = []
    numCGs = 0
    numGCs = 0
    global EXCLGCG

    def __init__(self, ref):
        length = len(ref)
        self.sequence = ref
        self.length = length
        self.CGpositions = detectCG(ref, length, EXCLGCG)
        self.GCpositions = detectGC(ref, length, EXCLGCG)
        self.numCGs = len(self.CGpositions)
        self.numGCs = len(self.GCpositions)

### General

def loadSequences(filename):
    return SeqIO.parse(filename, "fasta")

def detectCG(seq, length, excludeGCG=False):
    """Returns the list of C positions in CG dinucleotides in sequence `seq'.
If `excludeGCG' is True, ignores GCG positions."""
    result = []
    candidate = False

    for i in range(length):
        if (seq[i] == 'C'):
            candidate = i
        elif (seq[i] == 'G') and candidate:
            if excludeGCG:
                if (i < 2) or seq[i-2] != 'G':
                    result.append(candidate)
                candidate = False
            else:
                result.append(candidate)
                candidate = False
        else:
            candidate = False
    return result

def detectGC(seq, length, excludeGCG=False):
    """Returns the list of C positions in GC dinucleotides in sequence `seq'.
If `excludeGCG' is True, ignores GCG positions."""
    result = []
    candidate = False

    for i in range(1, length):
        if (seq[i] == 'C') and (seq[i-1] == 'G'):
            # This is a GC position. Now check position i+1 if excluding GCGs
            if excludeGCG:
                if (i == length-1) or seq[i+1] != 'G':
                    result.append(i)
            # Otherwise, simply add the position
            else:
                result.append(i)
    return result

def formatTabDelim(stream, l):
    stream.write("\t".join(l) + "\n")


def main():
    global EXCLGCG
    infile = None
    outfile = ""

    # Parse arguments
    args = sys.argv[1:]
    P.standardOpts(args)
    for arg in args:
        if arg == "-gcg":
            EXCLGCG = True
        elif infile == "":
            infile = P.isFile(arg)
        else:
            outfile = arg
    
    if not infile:
        P.errmsg(P.NOFILE)

    nreads = 0

    seqs = loadSequences(infile)
    rd = refDesc(seqs.next())   # reference sequence

    print("Reference sequence loaded from file `{}'.".format(infile))
    print("{}bp, {} CG positions, {} GC positions.".format(rd.length, rd.numCGs, rd.numGCs))

    CGarr = [ [p, 0] for p in rd.CGpositions ]
    GCarr = [ [p, 0] for p in rd.GCpositions ]

    print("Reading sequences...")
    for s in seqs:
        nreads += 1
        for p in CGarr:
            if s[p[0]] == 'C':
                p[1] += 1
        for p in GCarr:
            if s[p[0]] == 'C':
                p[1] += 1
                
    if outfile:
        out = open(outfile, "w")
    else:
        out = sys.stdout
    try:
        formatTabDelim(out, ["#CG", "Sites:", str(len(CGarr))])
        formatTabDelim(out, ["CG pos", "Num C", "Perc C"])
        for pair in CGarr:
            formatTabDelim(out, [str(pair[0]), str(pair[1]), str(1.0 * pair[1] / nreads)])
        out.write("\n")
        formatTabDelim(out, ["#GC", "Sites:", str(len(GCarr))])
        formatTabDelim(out, ["GC pos", "Num C", "Perc C"])
        for pair in GCarr:
            formatTabDelim(out, [str(pair[0]), str(pair[1]), str(1.0 * pair[1] / nreads)])
    finally:
        if outfile:
            out.close()

if __name__ == "__main__":
    main()

