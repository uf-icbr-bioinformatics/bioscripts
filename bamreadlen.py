#!/usr/bin/env python

import sys
import math
import os.path
import pysam

class BAMflagAnalyzer():
    infiles = []
    nreads = 0
    bits = [("PAIRED", 1),
            ("PROPER_PAIR", 2),
            ("UNMAP", 4),
            ("MUNMAP", 8),
            ("REVERSE", 16),
            ("MREVERSE", 32),
            ("READ1", 64),
            ("READ2", 128),
            ("SECONDARY", 256),
            ("QCFAIL", 512),
            ("DUP", 1024),
            ("SUPPLEMENTARY", 2048)]
    counts = {}

    def __init__(self):
        for b in self.bits:
            self.counts[b[0]] = 0

    def parseArgs(self, args):
        self.infiles = args

    def run(self):
        for f in self.infiles:
            self.parseBAM(f)
        self.report()

    def parseBAM(self, bamfile):
        bf = pysam.AlignmentFile(bamfile, "rb")
        try:
            for rec in bf.fetch():
                self.nreads += 1
                if (self.nreads % 1000000) == 0:
                    sys.stderr.write(chr(13) + "{:,} reads processed...".format(self.nreads) + "\033[K")
                v = rec.flag
                for b in self.bits:
                    if b[1] > v:
                        break
                    if v & b[1] != 0:
                        self.counts[b[0]] += 1
        except:
            bf.close()
            sys.stdout.write("\n")

    def report(self):
        sys.stdout.write("{:14} {:10d}\n".format("TOTAL:", self.nreads))
        for b in self.bits:
            c = self.counts[b[0]]
            sys.stdout.write("{:14} {:10d} ({:.2f}%)\n".format(b[0] + ":", c, 100.0 * c / self.nreads))

class BAManalyzer():
    infiles = []
    nreadsin = 0
    sumlen = 0
    maxlen = 0
    minlen = 1000000
    mode = 'avg' # or max, min
    howmany = 1000
    skip = 10

    def usage(self):
        sys.stdout.write("""bamreadlen.py - Determine read length from BAM file.

Usage: bamreadlen.py [options] bamfiles...

Compute read length from one or more BAM files. By default, the program will
examine one read every 10, until 1000 reads are reached in each file. Options:

  -h   | Print usage message.
  -n N | Number of reads to test in each BAM file (default: {}). Set this and
         -s to 0 to check all reads.
  -s S | Number of reads to skip between tested reads (default: {}).
  -m M | Use mode M, one of 'avg' (show average read length, rounded up), 
         'min' or 'max' (shortest or longest read respectively) or 'all'
         (all three values on separate lines).

""")

    def parseArgs(self, args):
        if "-h" in args:
            return self.usage()
        prev = ""
        for a in args:
            if prev == "-n":
                self.howmany = int(a)
                prev = ""
            elif prev == "-s":
                self.skip = int(a)
                prev = ""
            elif prev == "-m":
                self.mode = a
                prev = ""
            elif a in ["-n", "-s", "-m"]:
                prev = a
            elif os.path.isfile(a):
                self.infiles.append(a)
            else:
                sys.stderr.write("Unrecognized option: {}\n".format(a))
        if not self.infiles:
            return self.usage()

    def getReadlen(self):
        for f in self.infiles:
            self.getOneReadlen(f)
        if self.mode == 'avg':
            sys.stdout.write("{}\n".format(int(math.ceil(1.0*self.sumlen/self.nreadsin))))
        elif self.mode == 'max':
            sys.stdout.write("{}\n".format(self.maxlen))
        elif self.mode == 'min':
            sys.stdout.write("{}\n".format(self.minlen))
        elif self.mode == 'all':
            sys.stdout.write("{}\n".format(int(math.ceil(1.0*self.sumlen/self.nreadsin))))
            sys.stdout.write("{}\n".format(self.minlen))
            sys.stdout.write("{}\n".format(self.maxlen))

    def getOneReadlen(self, filename):
        seen = 0
        nreads = 0
        bf = pysam.AlignmentFile(filename, "rb")
        for rec in bf.fetch():
            seen += 1
            if seen >= self.skip:
                l = rec.rlen
                nreads += 1
                self.nreadsin += 1
                self.sumlen += l
                if l > self.maxlen:
                    self.maxlen = l
                if l < self.minlen:
                    self.minlen = l
                seen = 0
                if nreads == self.howmany:
                    break
        bf.close()

if __name__ == "__main__":
    cmd = os.path.split(sys.argv[0])[1]
    args = sys.argv[1:]
    if cmd == "bamreadlen.py":
        B = BAManalyzer()
        B.parseArgs(args)
        B.getReadlen()
    elif cmd == "bamstats.py":
        B = BAMflagAnalyzer()
        B.parseArgs(args)
        B.run()


