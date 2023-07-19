#!/usr/bin/env python

import sys
import gzip
import os.path

import Script
import Utils

__doc__ = """Count the number of sequences in the specified fasta or fastq files (and more)."""

### Program definition

def usage(what=None):
    sys.stderr.write("""countseqs.py - Count sequences in fasta/fastq files.

    Usage: countseqs.py [options] files...

Prints the number of sequences contained in the specified files. Files can be 
in fasta or fastq format, optionally compressed with gzip. Output is in four 
columns (tab-delimited): filename, number of sequences, total number of 
bases, average sequence length.

If -p is specified, checks for proper read pairing in a pair of fastq files.

If -c is specified, trims all reads to a specified interval.

Options:

  -h, --help | Print this usage message.
  -o outfile | Write output to outfile (instead of standard output).
  -t         | Print total of all files at the end.
  -m         | Print number of reads in millions.
  -q         | Report minimum and maximum quality score in each fastq file.
  -p         | Verify-paired mode: prints statistics on (in)correct pairing.
  -c C       | Trim all reads in a pair of fastq files to length C, or to a 
             | subsequence if C is in `slice' notation: 
             |   A:B = from position A to B (1-based, inclusive),
             |   A: = from position A to the end, etc.
  -C C       | Like -c, but for a single fastq file.

""")

OUTPUT = "/dev/stdout"
TOTAL = False
MILLIONS = False
CUT = False

class CountSeqs(Script.Script):
    output = "/dev/stdout"
    total = False
    millions = False
    cut = False

    def run(self, args):
        mode = "n"
        self.standardOpts(args)
        self.parseArgs(args, "+o,t,m,p,q,+c,+C")
        self.output = self.getOpt("o") or OUTPUT
        self.millions = self.getOpt("m")
        self.total = self.getOpt("t")
        files = self.getArgs()

        with open(self.output, "w") as out:
            if self.getOpt("p"):
                verifyPaired(files[0], files[1], out)
            elif self.getOpt("q"):
                for f in files:
                    checkQuality(f, out)
            elif self.getOpt("c"):
                self.cut = Utils.parseSlice(self.getOpt("c"))
                cutReads(self.cut, files[0], files[1], files[2], files[3])
            elif self.getOpt("C"):
                self.cut = Utils.parseSlice(self.getOpt("C"))
                cutReads1(self.cut, files[0], self.output)
            else:
                total = 0
                totbases = 0
                for filename in files:
                    (nseqs, nbases) = self.countSeqs(filename, out)
                    total += nseqs
                    totbases += nbases
                if self.total:
                    out.write("Total\t{}\t{}\t{:.1f}\n".format(self.printReads(total), totbases, 1.0*totbases/total))

    def printReads(self, r):
        if self.millions:
            return "{:.1f}M".format(r/1000000.0)
        else:
            return r

    def countSeqs(self, filename, out):
        nseqs = 0
        nbases = 0
        with genOpen(filename, "rt") as f:
            line = f.readline()
            if len(line) > 0:
                if line[0] == '>':
                    (nseqs, nbases, minlen, maxlen) = self.countSeqsFasta(f)
                    out.write("{}\t{}\t{}\t{:.1f}\t{}\t{}\n".format(filename, self.printReads(nseqs), nbases, 1.0*nbases/nseqs, minlen, maxlen))
                elif line[0] == '@':
                    (nseqs, nbases, minlen, maxlen) = self.countSeqsFastq(f)
                    out.write("{}\t{}\t{}\t{:.1f}\t{}\t{}\n".format(filename, self.printReads(nseqs), nbases, 1.0*nbases/nseqs, minlen, maxlen))
                else:
                    sys.stderr.write("Error: file `{}' is not in Fasta or FastQ format.\n".format(filename))
        return (nseqs, nbases)

    def countSeqsFasta(self, f):
        nseqs = 1
        nbases = 0
        minlen = 1000000
        maxlen = 0
        for line in f:
            if line[0] == '>':
                nseqs += 1
            else:
                rl = len(line) -1
                nbases += rl
                if rl < minlen:
                    minlen = rl
                if rl > maxlen:
                    maxlen = rl
        return (nseqs, nbases, minlen, maxlen)

    def countSeqsFastq(self, f):
        nseqs = 1
        nbases = 0
        minlen = 1000000
        maxlen = 0
        while True:
            # Skip rest of first record
            rl = len(f.readline())-1
            nbases += rl
            if rl < minlen:
                minlen = rl
            if rl > maxlen:
                maxlen = rl
            f.readline()
            f.readline()
            line = f.readline()
            if line == '':
                break
            if line[0] == '@':
                nseqs += 1
        return (nseqs, nbases, minlen, maxlen)

P = CountSeqs("countseqs.py", "1.0", usage=usage)

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

def verifyPaired(filename1, filename2, out):
    total = 0
    zero1 = 0
    zero2 = 0
    qdiff1 = 0
    qdiff2 = 0
    lendiff = 0
    qualdiff = 0
    namemismatch = 0
    with genOpen(filename1, "rt") as f1:
        with genOpen(filename2, "rt") as f2:
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
    out.write("""R1 fastq: {}
R2 fastq: {}
Reads: {}
Zero length in 1: {}
Zero length in 2: {}
Read/qual length mismatch in 1: {}
Read/qual length mismatch in 2: {}
Read length mismatch: {}
Qual length mismatch: {}
Name mismatch: {}
""".format(filename1, filename2, total, zero1, zero2, qdiff1, qdiff2, lendiff, qualdiff, namemismatch))

def cutReads(CUT, filename1, filename2, outfile1, outfile2):
    nin = 0
    nout = 0
    f1 = genOpen(filename1, "rt")
    f2 = genOpen(filename2, "rt")
    o1 = genOpen(outfile1, "wt")
    o2 = genOpen(outfile2, "wt")
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
            nin += len(r1) + len(r2)
            rc1 = r1[CUT]
            rc2 = r2[CUT]
            nout += len(rc1) + len(rc2)
            o1.write("{}\n{}\n{}\n{}\n".format(h1, rc1, d1, q1[CUT]))
            o2.write("{}\n{}\n{}\n{}\n".format(h2, rc2, d2, q2[CUT]))
    finally:
        f1.close()
        f2.close()
        o1.close()
        o2.close()
    sys.stderr.write("""{} bases in,
{} bases written.
""".format(nin, nout))

def cutReads1(CUT, filename1, outfile1):
    nin = 0
    nout = 0
    f1 = genOpen(filename1, "rt")
    o1 = genOpen(outfile1, "wt")
    try:
        while True:
            h1 = f1.readline().rstrip("\r\n")
            r1 = f1.readline().rstrip("\r\n")
            d1 = f1.readline().rstrip("\r\n")
            q1 = f1.readline().rstrip("\r\n")
            if h1 == '':
                break
            nin += len(r1)
            rc = r1[CUT]
            nout += len(rc)
            o1.write("{}\n{}\n{}\n{}\n".format(h1, rc, d1, q1[CUT]))
    finally:
        f1.close()
        o1.close()
    sys.stderr.write("""{} bases in,
{} bases written.
""".format(nin, nout))

def checkQuality(filename, out):
    minq = 1000
    maxq = 0
    with genOpen(filename, "r") as f:
        while True:
            try:
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
            except KeyboardInterrupt:
                break
    out.write("{}\t{}\t{}\n".format(filename, minq, maxq))

if __name__ == "__main__":
    args = sys.argv[1:]
    P.run(args)
    sys.exit(0)

    mode = "n"
    files = []
    next = ""
    P.standardOpts(args)
    for a in args:
        if next == '-o':
            OUTPUT = open(a, "w")
            next = ""
        elif next == '-c':
            mode = 'c'
            CUT = Utils.parseSlice(a)
            next = ""
        elif next == '-C':
            mode = 'C'
            CUT = Utils.parseSlice(a)
            next = ""
        elif a in ['-o', '-c', '-C']:
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
        elif mode == 'C':
            cutReads1(files[0], files[1])
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
