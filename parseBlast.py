#!/usr/bin/env python

import sys
import Script

doc = """A simple script to parse BLAST output filtering hits for:
- % Identity
- Alignment length
- Mismatches
- E-value
- Bit score
"""

# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
# M02337:42:000000000-ANMUK:1:1101:2619:10376     NC_002551.1     88.00   25      3       0       2       76      7854    7928    3e-08   55.6

class BLASTlimits():
    identity = None
    alnlength = None
    mismatches = None
    evalue = None
    bitscore = None

    def checkRecord(self, record):
        """Return True if this blast record (a list) matches the limits in this object."""
        if self.identity:
            if float(record[2]) < self.identity:
                return False
        if self.alnlength:
            if int(record[3]) < self.alnlength:
                return False
        if self.mismatches:
            if int(record[4]) > self.mismatches:
                return False
        if self.evalue:
            if float(record[10]) > self.evalue:
                return False
        if self.bitscore:
            if float(record[11]) < self.bitscore:
                return False
        return True

def usage():
    sys.stderr.write("Coming soon\n")

class ParseBlast(Script.Script):
    limits = None
    infiles = []
    outfile = None

    def parseArgs(self, args):
        if not args or "-h" in args or "--help" in args:
            return usage()
        prev = ""
        self.limits = BLASTlimits()
        self.infiles = []
        for a in args:
            if prev == "-i":
                self.limits.identity = self.toFloat(a)
                prev = ""
            elif prev == "-l":
                self.limits.alnlength = self.toInt(a)
                prev = ""
            elif prev == "-m":
                self.limits.mismatches = self.toInt(a)
                prev = ""
            elif prev == "-e":
                self.limits.evalue = self.toFloat(a)
                prev = ""
            elif prev == "-b":
                self.limits.bitscore = self.toFloat(a)
                prev = ""
            elif prev == "-o":
                self.outfile = a
                prev = ""
            elif a in ["-i", "-l", "-m", "-e", "-b", "-o"]:
                prev = a
            else:
                self.infiles.append(self.isFile(a))

    def parseOne(self, f, out):
        for line in f:
            if line[0] == '#':
                continue
            parsed = line.rstrip("\r\n").split("\t")
            if self.limits.checkRecord(parsed):
                out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(parsed[1], parsed[8], parsed[9], parsed[2], parsed[3], parsed[10], parsed[11], parsed[0]))

    def run(self):
        if self.outfile:
            out = open(self.outfile, "w")
        else:
            out = sys.stdout
        try:
            if self.infiles:
                for filename in self.infiles:
                    with open(filename, "r") as f:
                        self.parseOne(f, out)
            else:
                self.parseOne(sys.stdin, out)
        finally:
            if self.outfile:
                out.close()

P = ParseBlast("parseBlast", version="1.0", usage=usage)

if __name__ == "__main__":
    P.parseArgs(sys.argv[1:])
    P.run()
