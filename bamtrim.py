#!/usr/bin/env python

import sys
import gzip
import pysam

class Main(object):
    bamfile = ""
    fq1file = ""
    fq2file = ""
    aln = None
    fq1 = None
    fq2 = None

    read1_name = ""
    read1_len = 0
    read2_name = ""
    read2_len = 0

    def __init__(self, bamfile, fq1file, fq2file):
        self.bamfile = bamfile
        self.fq1file = fq1file
        self.fq2file = fq2file

    def init(self):
        self.fq1 = gzip.open(self.fq1file, "rt")
        self.fq2 = gzip.open(self.fq2file, "rt")

    def get_read(self):
        r = self.fq1.readline()
        name = r.split(" ")[0][1:]
        self.read1_name = name
        seq = self.fq1.readline().rstrip()
        self.read1_len = len(seq)
        self.fq1.readline()
        self.fq1.readline()

        r = self.fq2.readline()
        name = r.split(" ")[0][1:]
        self.read2_name = name
        seq = self.fq2.readline().rstrip()
        self.read2_len = len(seq)
        self.fq2.readline()
        self.fq2.readline()

        if self.read1_name != self.read2_name:
            sys.stderr.write("Error: fastq files are not properly paired! ({} / {})\n".format(self.read1_name, self.read2_name))
            sys.exit(1)

    def run(self, outfile):
        self.aln = pysam.AlignmentFile(self.bamfile, "rb", check_sq=False)
        out = pysam.AlignmentFile(outfile, "wb", template=self.aln)
        self.get_read()
        #print((self.read1_name, self.read1_len, self.read2_name, self.read2_len))
        try:
            while True:
                try:
                    r1 = next(self.aln)
                    r2 = next(self.aln)
                except StopIteration:
                    break
                if r1.query_name != r2.query_name:
                    sys.stderr.write("Error: reads in BAM file are not properly paired!\n")
                    sys.exit(2)
                if r1.query_name == self.read1_name:
                    r1.set_tag("Yt", self.read1_len, value_type="i")
                    r2.set_tag("Yt", self.read2_len, value_type="i")
                    out.write(r1)
                    out.write(r2)
                    self.get_read()
        finally:
            out.close()

class Converter(object):
    bamfile = ""

    def __init__(self, bamfile):
        self.bamfile = bamfile

    def run(self, fq1file, fq2file):
        with gzip.open(fq1file, "wt") as o1, gzip.open(fq2file, "wt") as o2:
            aln = pysam.AlignmentFile(self.bamfile, "rb", check_sq=False)
            while True:
                try:
                    r1 = next(aln)
                    r2 = next(aln)
                except StopIteration:
                    break
                    
                o1.write("@" + r1.query_name + "\n")
                o2.write("@" + r2.query_name + "\n")

                if r1.has_tag("Yt") and r2.has_tag("Yt"):
                    n1 = r1.get_tag("Yt")
                    n2 = r2.get_tag("Yt")
                    o1.write(r1.query_sequence[:n1])
                    o1.write("\n+\n")
                    o1.write(r1.query_qualities[:n1])
                    o1.write("\n")
                    o2.write(r2.query_sequence[:n2])
                    o2.write("\n+\n")
                    o2.write(r2.query_qualities[:n2])
                    o2.write("\n")

            aln.close()

def main(args):
    if args[0] == "-f":
        C = Converter(args[1])
        C.run(args[2], args[3])
    else:
        M = Main(args[0], args[1], args[2])
        M.init()
        M.run(args[3])

if __name__ == "__main__":
    main(sys.argv[1:])
