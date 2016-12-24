#!/usr/bin/env python

import sys
import os.path
import subprocess

import utils

def usage():
    sys.stderr.write("""bamToWig.py - Convert BAM file to WIG track for the UCSC genome browser.

Usage: bamToWig.py [-dntwfp] [-o wigfile] bamfile

Options: 
 -o FILE  | Name of output file (default: stdout)
 -w W     | Set window size to W (default: {})
 -n N     | Normalize values by number of reads N
 -s SCALE | Scale values by SCALE (y = r * SCALE / N)
 -t S     | Set track title to S
 -d D     | Set track description to D
 -f       | Converts diff file to bedGraph (value column: 8)
 -m       | Convert diff meth file to bedGraph (value column: 5)
 -p       | Convert a Homer peaks.txt file to bedGraph

""".format(trackdata.window))

### Program object

P = utils.Prog("bamToWig", version="1.0", usage=usage)

class trackdata():
    outfile = None
    window = 100
    normalize = 1
    scale = 1
    name = None
    description = None
    diff = False
    meth = False
    homer = False

    def trackHeader(self):
        if self.diff or self.homer:
            l = "track type=bedGraph"
        else:
            l = "track type=wiggle_0"
        if self.name:
            l += " name=" + self.name
        if self.description:
            l += " description=" + self.description
        return l + "\n"

    def trackFirstLine(self, chrom, pos):
        return "fixedStep chrom={} start={} step={} span={}\n".format(chrom, pos, self.window, self.window)

def bamToWig(bamfile, trackdata):
    currChrom = ""
    windowstart = 0
    windowend = 0
    windowsum = 0
    window = trackdata.window
    normalize = trackdata.normalize
    scale = trackdata.scale
    sys.stderr.write("Normalizing on {} reads, scale={}\n".format(normalize, scale))
    f = 1.0 * scale / normalize

    if trackdata.outfile:
        out = open(trackdata.outfile, "w")
    else:
        out = sys.stdout

    try:
        p = subprocess.Popen(['samtools', 'depth', bamfile], stdout=subprocess.PIPE)
        pin = p.stdout

        out.write(trackdata.trackHeader())

        while True:
            row = pin.readline().split("\t")
            if len(row) < 3:
                break
            chrom = row[0]
            pos = int(row[1])
            dp = int(row[2])
            if chrom != currChrom:
                if currChrom != "":
                    out.write("{}\n".format((1.0 * windowsum / window) * f))
                currChrom = chrom
                windowstart = pos
                windowend = windowstart + window
                windowsum = dp
                out.write(trackdata.trackFirstLine(chrom, pos))
            if pos > windowend:
                out.write("{}\n".format((1.0 * windowsum / window) * f))
                if pos > windowend + window: # gap - need to start a new track
                    out.write(trackdata.trackFirstLine(chrom, pos))
                    windowstart = pos
                    windowend = windowstart + window
                    windowsum = dp
                else:           #  still in current track
                    windowstart = windowend
                    windowend = windowstart + window
                    windowsum = dp
            else:
                windowsum += dp
    finally:
        if trackdata.outfile:
            out.close()

def splitCoords(c):
    p1 = c.find(":")
    p2 = c.find("-")
    return (c[0:p1], c[p1+1:p2], c[p2+1:])

def diffToBedGraph(infile, trackdata):
    if trackdata.outfile:
        out = open(trackdata.outfile, "w")
    else:
        out = sys.stdout
    out.write(trackdata.trackHeader())
    try:
        with open(infile, "r") as f:
            f.readline()
            while True:
                line = f.readline()
                if line == '':
                    break
                parsed = line.rstrip("\r\n").split("\t")
                (chrom, start, end) = splitCoords(parsed[1])
                fc = parsed[7]
                out.write("{}\t{}\t{}\t{}\n".format(chrom, start, end, fc))
    finally:
        if trackdata.outfile:
            out.close()

def methToBedGraph(infile, trackdata):
    if trackdata.outfile:
        out = open(trackdata.outfile, "w")
    else:
        out = sys.stdout
    out.write(trackdata.trackHeader())
    try:
        with open(infile, "r") as f:
            f.readline()
            while True:
                line = f.readline()
                if line == '':
                    break
                parsed = line.rstrip("\r\n").split("\t")
                chrom = parsed[0]
                pos = int(parsed[1])
                fc = float(parsed[4])
                out.write("{}\t{}\t{}\t{}\n".format(chrom, pos, pos+1, fc))
    finally:
        if trackdata.outfile:
            out.close()

def homerToBedGraph(infile, trackdata):
    if trackdata.outfile:
        out = open(trackdata.outfile, "w")
    else:
        out = sys.stdout
    out.write(trackdata.trackHeader())
    try:
        with open(infile, "r") as f:
            for line in f:
                if line != '' and line[0] != '#':
                    parsed = line.rstrip("\r\n").split("\t")
                    chrom = parsed[1]
                    start = int(parsed[2])
                    end = int(parsed[3])
                    fc = float(parsed[5])
                    out.write("{}\t{}\t{}\t{}\n".format(chrom, start, end, fc))
    finally:
        if trackdata.outfile:
            out.close()
    
def parseArgs(args, td):
    infile = None
    next = ""

    P.standardOpts(args)
    for a in args:
        if a == "-f":
            td.diff = True
        elif a == "-m":
            td.meth = True
        elif a == "-p":
            td.homer = True
        elif next == "-n":
            td.normalize = P.toInt(a)
            next = ""
        elif next == "-s":
            td.scale = P.toFloat(a)
            next = ""
        elif next == "-o":
            td.outfile = a
            next = ""
        elif next == "-w":
            td.window = P.toInt(a)
            next = ""
        elif next == "-t":
            td.name = a
            next = ""
        elif next == "-d":
            td.description = a
            next = ""
        elif a in ["-n", "-o", "-w", "-t", "-d", "-p", "-s"]:
            next = a
        else:
            infile = P.isFile(a)
    return infile

def main(args):
    td = trackdata()
    infile = parseArgs(args, td)
    if not infile:
        P.errmsg(P.NOFILE)
    if td.diff:
        diffToBedGraph(infile, td)
    elif td.meth:
        methToBedGraph(infile, td)
    elif td.homer:
        homerToBedGraph(infile, td)
    else:
        bamToWig(infile, td)

if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) > 0:
        main(args)
    else:
        usage();

