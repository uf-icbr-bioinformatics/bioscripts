#!/usr/bin/env python

import sys
import numpy

import Script

def parseFixed(f):
    data = f.split(" ")
    d = {}
    for w in data:
        if '=' in w:
            pieces = w.split("=")
            d[pieces[0]] = pieces[1]
    return d

def getChromName(filename):
    with open(filename, "r") as f:
        while True:
            line = f.readline().rstrip("\r\n")
            if line.startswith("fixed"):
                fd = parseFixed(line)
                return fd['chrom']
    
def fillVector(filename, length):
    start = 0
    step = 0
    v = numpy.zeros(length)
    with open(filename, "r") as f:
        f.readline()
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith("fixed"):
                fd = parseFixed(line)
                start = int(fd['start'])
                step = int(fd['step'])
            else:
                x = float(line)
                for i in range(start, start+step):
                    v[i] = x
                start += step
    return v

def writeVector(out, vd, name):
    out.write("track type=wiggle_0\n")
    out.write("fixedStep chrom={} start=0 step=1 span=1\n".format(name))
    for i in range(len(vd)):
        out.write("{}\n".format(vd[i]))

class SubWig(Script.Script):
    wig1 = None
    wig2 = None
    outfile = None
    length = 0
    clip = False
    reverse = False

    def parseArgs(self, args):
        prev = ""
        self.standardOpts(args)
        for a in args:
            if prev == "-o":
                self.outfile = a
                prev = ""
            elif prev == "-l":
                self.length = self.toInt(a)
                prev = ""
            elif a in ["-o", "-l"]:
                prev = a
            elif a == '-c':
                self.clip = True
            elif a == '-r':
                self.reverse = True
            elif self.wig1 is None:
                self.wig1 = self.isFile(a)
            else:
                self.wig2 = self.isFile(a)

    def run(self):
        name = getChromName(self.wig1)
        vd = self.wigDifference()
        if self.outfile:
            with open(self.outfile, "w") as out:
                writeVector(out, vd, name)
            sys.stderr.write("{} written.\n".format(self.outfile))
        else:
            writeVector(sys.stdout, vd, name)

    def wigDifference(self):
        v1 = fillVector(self.wig1, self.length)
        v2 = fillVector(self.wig2, self.length)
        d = v1 - v2
        if self.clip:
            if self.reverse:
                d = d.clip(max=0.0)
            else:
                d = d.clip(min=0.0)
        return d

def usage():
    sys.stdout.write("""subwig.py - Subtract two WIG files.

Usage: subwig.py [options] file1.wig file2.wig

This command will subtract the values in file2.wig from those
in file1.wig, and create a new WIG file of the difference. The 
result will be written to the file specified with the -o option, 
or to standard output.

Limitations:
  - Currently only works on single-chromosome WIG files;
  - The length of the sequence needs to be specified;
  - The output WIG file is not optimized (ie, zero regions
    are not removed).
These will be fixed in a future version.

Options:

  -l L | Length of sequence (required)
  -o O | Write output to file O (default: stdout).
  -c   | Clip mode: don't output values below 0.
  -r   | Reverse-strand mode (all values should be negative).
""")



if __name__ == "__main__":
    args = sys.argv[1:]
    S = SubWig("subwig.py", version="1.0", usage=usage)
    S.parseArgs(args)
    if S.wig1 and S.wig2 and S.length:
        S.run()
    else:
        usage()
