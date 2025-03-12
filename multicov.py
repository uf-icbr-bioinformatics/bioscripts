#!/usr/bin/env python

import sys
import csv
import os.path
import subprocess as sp

class Multicov(object):
    bed = None
    bams = []
    ncols = 0
    scale = 1000000.0
    outfile = "/dev/stdout"
    _mode = "a"
    nreads = []

    def usage(self):
        sys.stdout.write("""multicov.py - compute normalized coverage of one or more BAM files on a set of regions.

Usage: multicov.py [options] file.bed files.bam...

Where options are:

  -o O | Write output to file O (default: stdout).
  -s S | Multiply all normalized values by S (default: {}).

This program is similar to `samtools bedcov', but coverage values are automatically normalized to the total
number of mapped reads in each BAM file.

Requires BAM files to be indexed, and `samtools' to be in PATH.

""".format(self.scale))

    def parseArgs(self, args):
        if "-h" in args or "--help" in args:
            return self.usage()
        self.bams = []
        prev = ""
        for a in args:
            if prev == "-o":
                self.outfile = a
                self._mode = "w"
                prev = ""
            elif prev == "-s":
                self.scale = float(a)
                prev = ""
            elif a in ["-o", "-s"]:
                prev = a
            elif self.bed is None:
                self.bed = a
            else:
                self.bams.append(a)
        if not os.path.isfile(self.bed):
            sys.stderr.write("Error: BED file `{}' does not exist.\n".format(self.bed))
            sys.exit(1)
        for b in self.bams:
            if not os.path.isfile(b):
                sys.stderr.write("Error: BAM file `{}' does not exist.\n".format(self.bed))
                sys.exit(1)
        self.ncols = len(self.bams)
        self.names = [ os.path.splitext(os.path.split(b)[1])[0] for b in self.bams ]
        return self.bed and self.bams

    def getNreads(self):
        sys.stderr.write("Computing number of reads:\n")
        for b in self.bams:
            cmdline = "samtools idxstats " + b
            out = sp.run(cmdline, shell=True, capture_output=True)
            nr = 0
            for row in out.stdout.decode().split("\n"):
                parts = row.split("\t")
                if len(parts) < 3:
                    continue
                if parts[0] != "*":
                    nr += int(parts[2])
            sys.stderr.write("{}\t{}\n".format(b, nr))
            self.nreads.append(nr)

    def getCoverage(self):
        sys.stderr.write("Generating coverage table...\n")
        with open(self.outfile, self._mode) as out:
            out.write("Region\t" + "\t".join(self.names) + "\n")
            cmdline = "samtools bedcov {} {}".format(self.bed, " ".join(self.bams))
            proc = sp.Popen(cmdline, shell=True, stdout=sp.PIPE)
            for line in proc.stdout:
                row = line.decode().strip().split("\t")
                out.write("{}:{}-{}".format(row[0], row[1], row[2]))
                for i in range(self.ncols):
                    x = int(row[i+3]) * self.scale / self.nreads[i]
                    out.write("\t" + str(x))
                out.write("\n")

def main(args):
    M = Multicov()
    if M.parseArgs(args):
        M.getNreads()
        M.getCoverage()
    else:
        M.usage()

if __name__ == "__main__":
    main(sys.argv[1:])
