#!/usr/bin/env python

import sys
import csv
import math

def parseSlice(s):
    if "-" in s:
        parts = s.split("-")
        return slice(int(parts[0]) - 1, int(parts[1]))
    else:
        p = int(s)
        return slice(p-1, p)

class SimpleDiff():
    filename = None
    outfile = "/dev/stdout"
    labels = None
    colname1 = "avg1"
    colname2 = "avg2"
    alpha = 1.0
    slice1 = None
    slice2 = None

    def process(self, f, out, header=True):
        nin = 0
        nout = 0
        na = self.slice1.stop - self.slice1.start
        nb = self.slice2.stop - self.slice2.start
        if header:
            f.readline()
        c = csv.reader(f, delimiter='\t')
        for line in c:
            nin += 1
            data1 = line[self.slice1]
            data2 = line[self.slice2]
            data1 = [ float(v) for v in data1 ]
            data2 = [ float(v) for v in data2 ]
            amin = min(data1)
            amax = max(data1)
            bmin = min(data2)
            bmax = max(data2)
            if amin > bmax:
                # A over B
                r1 = amax - amin
                r2 = bmax - bmin
                d = self.alpha * max(r1, r2)
                if (amin - bmax) > d:
                    avg1 = sum(data1) / na
                    avg2 = sum(data2) / nb
                    if avg1 > 0 and avg2 > 0:
                        out.write("{}\t{}\t{}\t{}\n".format(line[0], avg1, avg2, math.log(avg1/avg2, 2.0)))
                        nout += 1
            elif bmin > amax:
                # B over A
                r1 = amax - amin
                r2 = bmax - bmin
                d = self.alpha * max(r1, r2)
                if (bmin - amax) > d:
                    avg1 = sum(data1) / na
                    avg2 = sum(data2) / nb
                    if avg1 > 0 and avg2 > 0:
                        out.write("{}\t{}\t{}\t{}\n".format(line[0], avg1, avg2, math.log(avg1/avg2, 2.0)))
                        nout += 1
        return (nin, nout)

    def parseArgs(self, args):
        prev = ""
        if "-h" in args or "--help" in args:
            return self.usage()
        for a in args:
            if prev == "-a":
                self.alpha = float(a)
                prev = ""
            elif prev == "-o":
                self.outfile = a
                prev = ""
            elif prev == "-l":
                self.labels = parseSlice(a)
                prev = ""
            elif prev == "-c1":
                self.colname1 = a
                prev = ""
            elif prev == "-c2":
                self.colname2 = a
                prev = ""
            elif a in ["-a", "-o", "-l", "-c1", "-c2"]:
                prev = a
            elif self.filename is None:
                self.filename = a
            elif self.slice1 is None:
                self.slice1 = parseSlice(a)
            elif self.slice2 is None:
                self.slice2 = parseSlice(a)

        if (self.filename and self.slice1 and self.slice2):
            return True
        else:
            return self.usage()

    def usage(self):
        sys.stdout.write("""Usage: simplediff.py [options] exprfile slice1 slice2

This program performs "simple" differential analysis on gene expression data. `exprfile'
should be a file containing gene expression values with genes on the rows and samples
in the columns. `slice1' and `slice2' should be expressions of the form P-Q indicating
which columns contain the data for the two conditions being compared (e.g., if the first
condition is represented by three columns starting at column 5, use 5-7).

Options:

  -a A  | Set the alpha parameter to A (see below). Default: {}.
  -o O  | Write output to file O.
  -c1 C | Set label for average of condition 1 values to C. Default: {}.
  -c1 C | Set label for average of condition 2 values to C. Default: {}.

A gene is considered to be differentially expressed between two groups of samples (A and B)
if the two following conditions hold:

  * The two sets of expression values are totally separated, ie:
  
      the minimum expression values for the samples in A is larger than the maximum in B
      -OR-
      the minimum expression values for the samples in B is larger than the maximum in A

  * The distance between the two sets of values (the difference between the maximum of 
    the "lower" one and the minimum of the "upper" one) is larger than the largest of the
    two ranges of values in A and B, multiplied by the alpha parameter.

Example: A = {{10, 12, 16}}
         B = {{20, 21, 22}}

The two sets are separated, because min(B) > max(A). The distance between the two sets is
4 (20-16), range(A) = 6, range(B) = 2. If alpha is set to 1.0 (the default) then this
gene would NOT be considered significantly different, because the largest range is 6, 
and 6 * alpha > 4. If alpha was set to 0.5, the gene would be called as different.

""".format(self.alpha, self.colname1, self.colname2))


    def run(self):
        with open(self.outfile, "w") as out:
            with open(self.filename, "r") as f:
                (nin, nout) = self.process(f, out)
                sys.stderr.write("{} in, {} out\n".format(nin, nout))

if __name__ == "__main__":
    SD = SimpleDiff()
    if SD.parseArgs(sys.argv[1:]):
        SD.run()
