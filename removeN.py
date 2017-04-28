#!/usr/bin/env python

import sys
import Script

def usage():
    sys.stderr.write("""removeN.py - Remove Ns and/or gaps from sequences in FASTA file.

Remove all occurrences of N or n from the sequences in input multi-FASTA 
file `infile'. Output is written to standard output or to `outfile' if 
specified, in FASTA format. The -l option specifies the line length of the 
resulting output file.

In `clean' mode (enabled with -c), removes from output all sequences with
an excessive number of gaps (-). The decision whether to remove a sequence 
is controlled by the two parameters -s (length of longest allowed stretch
of '-') and -f (highest allowed fraction of - vs sequence length).

Usage: removeN.py [options] infile

Options:
  -o outfile | Write output to outfile (default: stdout).
  -l L       | Set line length in output FASTA file to L (default: {})
  -c         | Enable 'clean' mode.
  -s S       | Do not output sequence if it contains a stretch of '-' 
               longer than S (default: {}).
  -f F       | Do not output sequence if the fraction of '-' it contains
               is larger than F (default: {}).

""".format(Remn.MAXL, Remn.STRETCH, Remn.FRAC))

### Program object

class Remn(Script.Script):
    infile = None
    outfile = None
    MAXL = 70
    STRETCH = 10
    FRAC = 0.1
    CLEAN = False
    nread = 0                   # For clean mode
    nwritten = 0                # For clean mode

    def runMain(self, out):
        with open(self.infile, "r") as instream:
            nout = 0
            for line in instream:
                if len(line) > 0 and line[0] == '>':
                    if 0 < nout < self.MAXL:
                        out.write("\n")
                    out.write(line)
                    nout = 0
                else:
                    for c in line:
                        if c not in ['N', 'n', '\r', '\n']:
                            out.write(c)
                            nout += 1
                            if nout == self.MAXL:
                                out.write('\n')
                                nout = 0

    def runClean(self, out):
        # sys.stderr.write("S={}, F={}\n".format(self.STRETCH, self.FRAC))
        seqname = ""
        seq = ""
        with open(self.infile, "r") as instream:
            for line in instream:
                line = line.rstrip("\r\n")
                if len(line) > 0 and line[0] == '>':
                    self.maybeWriteSeq(out, seqname, seq)
                    seq = ""
                    seqname = line
                else:
                    seq += line
        self.maybeWriteSeq(out, seqname, seq)
        sys.stderr.write("{} sequences read, {} written.\n".format(self.nread, self.nwritten))

    def maybeWriteSeq(self, out, seqname, seq):
        if len(seq) == 0:
            return              # Don't write empty sequences
        self.nread += 1
        nbases = 0
        ndashes = 0
        strstart = 0
        strend = 0
        strmax = 0
        state = 'seq'
        pos = 0

        for b in seq:
            nbases += 1
            if state == 'seq':
                if b == '-':
                    ndashes += 1
                    strstart = pos
                    state = 'dash'
            else:
                if b == '-':
                    ndashes += 1
                else:
                    strend = pos
                    l = strend-strstart
                    if l > strmax:
                        strmax = l
                    state = 'seq'
            pos += 1

        if state == 'dash':
            l = pos - strstart
            if l > strmax:
                strmax = l
        # sys.stderr.write("strmax={}, frac={}  ".format(strmax, 1.0*ndashes/nbases))
        if (strmax <= self.STRETCH) and (1.0*ndashes/nbases <= P.FRAC):
            self.nwritten += 1
            self.writeSeq(out, seqname, seq)
            # sys.stderr.write("Y\n")
        else:
            pass
            # sys.stderr.write("N\n")

    def writeSeq(self, out, seqname, seq):
        out.write(seqname + "\n")
        nout = 0
        for b in seq:
            out.write(b)
            nout += 1
            if nout == self.MAXL:
                out.write("\n")
                nout = 0
        if nout != 0:
            out.write("\n")

P = Remn("removeN", version="1.0", usage=usage)

def parseArgs(args):
    next = ""

    P.standardOpts(args)
    for a in args:
        if next == "-l":
            P.MAXL = P.toInt(a)
            next = ""
        elif next == "-o":
            P.outfile = a
            next = ""
        elif next == "-s":
            P.STRETCH = P.toInt(a)
            next = ""
        elif next == "-f":
            P.FRAC = P.toFloat(a)
            next = ""
        elif a in ["-l", "-o", "-s", "-f"]:
            next = a
        elif a == "-c":
            P.CLEAN = True
        elif P.infile == None:
            P.infile = P.isFile(a)
    if P.infile == None:
        P.errmsg(P.NOFILE)

if __name__ == "__main__":
    parseArgs(sys.argv[1:])
    try:
        if P.outfile:
            with open(P.outfile, "w") as out:
                if P.CLEAN:
                    P.runClean(out)
                else:
                    P.runMain(out)
        else:
            if P.CLEAN:
                P.runClean(sys.stdout)
            else:
                P.runMain(sys.stdout)

    except Exception as e:
        P.errmsg(P.ERR, e)
