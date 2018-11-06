#!/usr/bin/env python

import sys
import gzip
import os.path

import Script
import Utils

__doc__ = """This is doc."""

### Program definition

def usage():
    sys.stderr.write("""countseqs.py - Count sequences in fasta/fastq files.

Usage: countseqs.py [-h] [-m] [-t] [-o outfile] files...

Prints the number of sequences contained in the specified files. Files can be 
in fasta or fastq format, optionally compressed with gzip. Output is in four 
columns (tab-delimited): filename, number of sequences, total number of 
bases, average sequence length.

Options:

-h, --help | Print this usage message
-o outfile | Write output to outfile (instead of standard output)
-t         | Print total of all files at the end
-m         | Print number of reads in millions
-p         | Verify-paired mode: prints statistics on (in)correct pairing.
-c C       | Trim all reads to lenght C, or to a subsequence if C is in 
           | `slice' notation: A:B = from position A to B (1-based, inclusive),
           | A: = from position A to the end, etc.

""")

P = Script.Script("countseqs.py", "1.0", usage=usage)

OUTPUT = sys.stdout
TOTAL = False
MILLIONS = False
CUT = False

def printReads(r):
    if MILLIONS:
        return "{:.1f}M".format(r/1000000.0)
    else:
        return r

def getReadName(h):
    p = h.find(" ")
    if p == -1:
        return h
    else:
        return h[:p]

def genOpen(filename, mode):
    """Generalized open() function - works on both regular files and .gz files."""
    (name, ext) = os.path.splitext(filename)
    if ext == ".gz":
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def countSeqs(filename):
    nseqs = 0
    nbases = 0
    with genOpen(filename, "r") as f:
        line = f.readline()
        if len(line) > 0:
            if line[0] == '>':
                (nseqs, nbases) = countSeqsFasta(f)
                OUTPUT.write("{}\t{}\t{}\t{:.1f}\n".format(filename, printReads(nseqs), nbases, 1.0*nbases/nseqs))
            elif line[0] == '@':
                (nseqs, nbases) = countSeqsFastq(f)
                OUTPUT.write("{}\t{}\t{}\t{:.1f}\n".format(filename, printReads(nseqs), nbases, 1.0*nbases/nseqs))
            else:
                sys.stderr.write("Error: file `{}' is not in Fasta or FastQ format.\n".format(filename))
    return (nseqs, nbases)

def countSeqsFasta(f):
    nseqs = 1
    nbases = 0
    for line in f:
        if line[0] == '>':
            nseqs += 1
        else:
            nbases += len(line) -1
    return (nseqs, nbases)

def countSeqsFastq(f):
    nseqs = 1
    nbases = 0
    while True:
        # Skip rest of first record
        nbases += len(f.readline())-1
        f.readline()
        f.readline()
        line = f.readline()
        if line == '':
            break
        if line[0] == '@':
            nseqs += 1
    return (nseqs, nbases)

def verifyPaired(filename1, filename2):
    total = 0
    zero1 = 0
    zero2 = 0
    qdiff1 = 0
    qdiff2 = 0
    lendiff = 0
    qualdiff = 0
    namemismatch = 0
    with genOpen(filename1, "r") as f1:
        with genOpen(filename2, "r") as f2:
            while True:
                h1 = f1.readline()
                r1 = f1.readline()
                d1 = f1.readline()
                q1 = f1.readline()
                h2 = f2.readline()
                r2 = f2.readline()
                d2 = f2.readline()
                q2 = f2.readline()
                if h1 == '' and h2 == '':
                    break
                total += 1
                if getReadName(h1) != getReadName(h2):
                    namemismatch += 1
                l1 = len(r1)
                ql1 = len(q1)
                l2 = len(r2)
                ql2 = len(q2)
                if l1 == 0:
                    zero1 += 1
                if l2 == 0:
                    zero2 += 1
                if l1 != ql1:
                    qdiff1 += 1
                if l2 != ql2:
                    qdiff2 += 1
                if l1 != l2:
                    lendiff += 1
                if ql1 != ql2:
                    qualdiff += 1
    sys.stdout.write("""Reads: {}
Zero length in 1: {}
Zero length in 2: {}
Read/qual length mismatch in 1: {}
Read/qual length mismatch in 2: {}
Read length mismatch: {}
Qual length mismatch: {}
Name mismatch: {}
    """.format(total, zero1, zero2, qdiff1, qdiff2, lendiff, qualdiff, namemismatch))

def cutReads(filename1, filename2, outfile1, outfile2):
    nin = 0
    nout = 0
    f1 = genOpen(filename1, "r")
    f2 = genOpen(filename2, "r")
    o1 = genOpen(outfile1, "w")
    o2 = genOpen(outfile2, "w")
    try:
        while True:
            h1 = f1.readline().rstrip("\r\n")
            r1 = f1.readline().rstrip("\r\n")
            d1 = f1.readline().rstrip("\r\n")
            q1 = f1.readline().rstrip("\r\n")
            h2 = f2.readline().rstrip("\r\n")
            r2 = f2.readline().rstrip("\r\n")
            d2 = f2.readline().rstrip("\r\n")
            q2 = f2.readline().rstrip("\r\n")
            if h1 == '' and h2 == '':
                break
#            if len(r1) < CUT or len(r2) < CUT:
#                continue
            o1.write("{}\n{}\n{}\n{}\n".format(h1, r1[CUT], d1, q1[CUT]))
            o2.write("{}\n{}\n{}\n{}\n".format(h2, r2[CUT], d2, q2[CUT]))
    finally:
        f1.close()
        f2.close()
        o1.close()
        o2.close()
    sys.stderr.write("""{} reads in,
{} reads written.
""".format(nin, nout))

def checkQuality(filename):
    minq = 1000
    maxq = 0
    with genOpen(filename, "r") as f:
        while True:
            h = f.readline().rstrip("\r\n")
            r = f.readline().rstrip("\r\n")
            d = f.readline().rstrip("\r\n")
            q = f.readline().rstrip("\r\n")
            if h == '':
                break
            for v in q:
                x = ord(v)
                if x > maxq:
                    maxq = x
                if x < minq:
                    minq = x
    sys.stderr.write("{}\t{}\t{}\n".format(filename, minq, maxq))

if __name__ == "__main__":
    mode = "n"
    files = []
    next = ""
    args = sys.argv[1:]
    P.standardOpts(args)
    for a in args:
        if next == '-o':
            OUTPUT = open(a, "w")
            next = ""
        elif next == '-c':
            mode = 'c'
            CUT = Utils.parseSlice(a)
            next = ""
        elif a in ['-o', '-c']:
            next = a
        elif a == '-m':
            MILLIONS = True
        elif a == '-t':
            TOTAL = True
        elif a == '-p':
            mode = 'p'
        elif a == '-q':
            mode = 'q'
        else:
            #files.append(P.isFile(a))
            files.append(a)

    if len(files) == 0:
        P.errmsg(P.NOFILE)

    try:
        if mode == 'p':
            verifyPaired(files[0], files[1])
        elif mode == 'c':
            cutReads(files[0], files[1], files[2], files[3])
        elif mode == 'q':
            for f in files:
                checkQuality(f)
        else:
            total = 0
            totbases = 0
            for filename in files:
                (nseqs, nbases) = countSeqs(filename)
                total += nseqs
                totbases += nbases
            if TOTAL:
                OUTPUT.write("Total\t{}\t{}\t{:.1f}\n".format(printReads(total), totbases, 1.0*totbases/total))
    finally:
        OUTPUT.close()
