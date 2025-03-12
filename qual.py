#!/usr/bin/env python

import sys
import gzip
import os.path
import Utils

def convert(qstring):
    quals = []
    for c in qstring:
        q = ord(c) - 33
        quals.append(q)
    return quals

def average(quals):
    return 1.0 * sum(quals) / len(quals)

def doTrimming(qstring, quals, score):
    l = len(qstring)
    s = sum(quals)
    i = l - 1
    while True:
        a = 1.0 * s / l
        if a >= score:
            return (qstring[:l], quals[:l])
        l -= 1
        if l == 0:
            return (None, None)
        s -= quals[l]

class Main(object):
    score = None
    clip = None
    filenames = []
    readstdin = True
    maxreads = 0
    q30 = False

    def __init__(self):
        self.filenames = []

    def usage(self):
        sys.stdout.write("""qual.py - Analyze quality scores

Usage: qual.py [options] [args...]

Where options are:

  -t T   | Trim quality scores from end until average quality reaches T.
  -c A:B | Clip quality score string to specified interval (one-based).
  -m M   | When reading from files or stdin, examine an most M reads.
  -q Q   | Disable regular output, print fraction of reads with average quality >= Q.

If no args... are specified, reads quality strings from standard input.

If args are existing filenames, they are assumed to be in compressed fastq
format, and the program reads quality scores from there.

Otherwise, args are assumed to be quality score strings.

""")

    def parseArgs(self, args):
        if "-h" in args or "--help" in args:
            return False
        prev = ""
        for a in args:
            if prev == "-t":
                self.score = int(a)
                prev = ""
            elif prev == "-c":
                self.clip = Utils.parseSlice(a)
                prev = ""
            elif prev == "-q":
                self.q30 = float(a)
                prev = ""
            elif prev == "-m":
                self.maxreads = int(a)
                prev = ""
            elif a in ["-t", "-c", "-q", "-m"]:
                prev = a
            elif os.path.isfile(a):
                self.filenames.append(a)
            else:
                self.showQuality(a)
                self.readstdin = False
        return True

    def run(self):
        if self.filenames:
            for fn in self.filenames:
                with gzip.open(fn, "rt") as f:
                    (nr, gq) = self.showQualityStream(f)
                    if self.q30:
                        sys.stdout.write("{}\t{}\t{:.1f}\n".format(fn, nr, 100.0*gq/nr))
        elif self.readstdin:
            (nr, gq) = self.showQualityStream(sys.stdin)
            if self.q30:
                sys.stdout.write("stdin\t{}\t{:.1f}\n".format(nr, 1.0*gq/nr))

    def showQualityStream(self, f):
        nr = 0                  # number of reads seen
        gq = 0                  # number of reads with good quality scores
        n = 0
        while True:
            line = f.readline()
            if not line:
                break
            n += 1
            if n == 4:
                aq = self.showQuality(line.rstrip("\r\n"))
                if aq and self.q30 and aq >= self.q30:
                    gq += 1
                n = 0
                nr += 1
                if nr == self.maxreads:
                    break
        return (nr, gq)

    def showQuality(self, qstring):
        quals = convert(qstring)
        if self.score:
            (qstring, quals) = doTrimming(qstring, quals, self.score)
        if self.clip:
            qstring = qstring[self.clip]
            quals = quals[self.clip]
        if quals:
            aq = average(quals)
            if not self.q30:
                sys.stdout.write("{}\t{}\t{:.2f}\t{}\t{}\n".format(qstring, len(qstring), aq, min(quals), max(quals)))
            return aq
        return None

if __name__ == "__main__":
    M = Main()
    if M.parseArgs(sys.argv[1:]):
        M.run()
    else:
        M.usage()
