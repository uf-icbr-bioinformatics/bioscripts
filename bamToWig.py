#!/usr/bin/env python

import sys
import os.path
import subprocess
import pysam

import Script
from Utils import LinkedList, LinkedListPool, GenomicWindower

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
 -e       | Do not remove ERCC controls
 -f       | Converts diff file to bedGraph (value column: 8)
 -m       | Convert diff meth file to bedGraph (value column: 5)
 -p       | Convert a Homer peaks.txt file to bedGraph
 -a       | ATAC mode (build pileup on read starts only).

""".format(trackdata.window))

### Program object

P = Script.Script("bamToWig", version="1.0", usage=usage)

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
    ercc = False                # If True, preserve ERCC reads
    atac = False

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
    wanted = True
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
                wanted = (not chrom.startswith("ERCC") or trackdata.ercc)
                windowstart = pos
                windowend = windowstart + window
                windowsum = dp
                if wanted:
                    out.write(trackdata.trackFirstLine(chrom, pos))
            if pos > windowend:
                if wanted: 
                    out.write("{}\n".format((1.0 * windowsum / window) * f))
                if pos > windowend + window: # gap - need to start a new track
                    if wanted:
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
    data = {}

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
                    if chrom not in data:
                        data[chrom] = []
                    start = int(parsed[2])
                    end = int(parsed[3])
                    fc = float(parsed[5])
                    data[chrom].append((start, end, fc))
            
        # Write regions in order, taking care of handling overlaps
        for chrom in sorted(data):
            data[chrom].sort(key=lambda r: r[0])
            prev = data[chrom][0]
            for row in data[chrom][1:]:
                if row[0] <= prev[1]:
                    prev[1] = max(prev[1], row[1])
                    prev[2] += row[2]
                else:
                    out.write("{}\t{}\t{}\t{}\n".format(chrom, prev[0], prev[1], prev[2]))
                    prev = row
            out.write("{}\t{}\t{}\t{}\n".format(chrom, prev[0], prev[1], prev[2]))
    finally:
        if trackdata.outfile:
            out.close()

# ATAC mode
    
class PosLinkedList(LinkedListPool):
    """Specialize the smaller() method so it works on objects of the form (chr, pos)."""

    def smaller(self, a, b):
        return a[1] <= b[1]

def bamToWigA_old(bamfile, trackdata):
    wanted = True
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

    L = PosLinkedList()
    L.preallocate(200)
    G = GenomicWindower(window, out)
    try:
        #a = pysam.AlignmentFile(bamfile, "rb")
        #for r in a.fetch():
        for line in sys.stdin:
            line = line.rstrip("\n").split("\t")
            #chrom = r.reference_name
            chrom = line[0]

            # First, let's check if we've changed chromosomes.
            # If so, all elements in the linked list need to be printed out.
            if chrom != currChrom:
                curr = L.current()
                if curr and (chrom != curr[0]):
                    #out.write("# New chromosome, emptying list of length {}.\n".format(L.length()))
                    #L.dump(out)
                    while True:
                        #out.write("{}\t{}\t-\n".format(curr[0], curr[1]))
                        G.add(curr)
                        x = L.pop()
                        L.deallocate(x)
                        curr = L.current()
                        if curr is None or curr[0] == chrom:
                            break

            #if r.is_read1
            if line[2] == "1":

                # Check if the linked list contains positions that are before 
                # the one we're seeing. If so, output them.
                #p = r.pos
                p = int(line[1])
                curr = L.current()
                if curr and (curr[1] <= p):
                    while True:
                        #out.write("{}\t{}\t-\n".format(curr[0], curr[1]))
                        #print("Output: 2 @ {}".format(curr[1]))
                        G.add(curr)
                        x = L.pop()
                        L.deallocate(x)
                        curr = L.current()
                        if curr is None or curr[1] > p:
                            break
                #print("Output: 1 @ {}".format(p))
                # Finally output entry for this read
                G.add((chrom, p))
                #out.write("{}\t{}\t+\n".format(chrom, p))
                currChrom = chrom

            else:
                #print("Storing {}".format(r.pos + r.qlen - 1))
                #L.insert((r.reference_name, r.pos + r.qlen - 1))
                #print("Storing: 2 @ {}".format(line[1]))
                L.insert((chrom, int(line[1])))
    finally:
        G.close()
        sys.stderr.write("Total conses: {}\n".format(L.nconsed))
        if trackdata.outfile:
            out.close()
            
def bamToWigA(trackdata):
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

    G = GenomicWindower(window, out)
    try:
        for line in sys.stdin:
            line = line.rstrip("\n").split("\t")
            G.add((line[0], int(line[1])))
    finally:
        G.close()
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
        elif a == "-e":
            td.ercc = True
        elif a == "-a":
            td.atac = True
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
    if not infile and not td.atac:
        P.errmsg(P.NOFILE)
    if td.diff:
        diffToBedGraph(infile, td)
    elif td.meth:
        methToBedGraph(infile, td)
    elif td.homer:
        homerToBedGraph(infile, td)
    elif td.atac:
        bamToWigA(td)
    else:
        bamToWig(infile, td)

if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) > 0:
        main(args)
    else:
        P.errmsg(P.NOFILE)


