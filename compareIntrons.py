#!/usr/bin/env python

import sys
import math
import os.path

import Script

### Program object

def usage():
    progname = os.path.split(sys.argv[0])[1]
    sys.stderr.write("""{} - analyze intron retention.

Usage: {} [options] intronsdb

Options:
 -i1 FILE | Name of BED file for introns in sample 1 (required)
 -i2 FILE | Name of BED file for introns in sample 2 (required)
 -j1 FILE | Name of BED file for junctions in sample 1
 -j2 FILE | Name of BED file for junctions in sample 2
 -o  FILE | Set output file to FILE (default: stdout)
 -fc F    | Set fold change threshold to F (default: {})
 -t  T    | Set coverage threshold to T (default: {})

""".format(progname, progname, Params.fc, Params.thr))

P = Script.Script("compareIntrons", version="1.0", usage=usage)

def readBEDfile(bedfile):
    dict = {}
    sys.stderr.write("Reading `{}'... ".format(bedfile))
    with open(bedfile, "r") as f:
        for line in f:
            parsed = line.rstrip("\r\n").split("\t")
            dict[parsed[3]] = float(parsed[7])
    sys.stderr.write("done, {} entries.\n".format(len(dict)))
    return dict

def dget(key, dict):
    if key in dict:
        return dict[key]
    else:
        return 0.0

class Params():
    bedfile = None
    introns1file = None
    juncs1file = None
    introns2file = None
    juncs2file = None
    outfile = None
    fc = 1
    thr = 0.00001
    intr1 = {}
    junc1 = {}
    intr2 = {}
    junc2 = {}

    def __init__(self):
        self.intr1 = {}
        self.junc1 = {}
        self.intr2 = {}
        self.junc2 = {}

    def parseArgs(self, args):
        P.standardOpts(args)
        next = ""
        for a in args:
            if next == "-i1":
                self.introns1file = P.isFile(a)
                next = ""
            elif next == "-j1":
                self.juncs1file = P.isFile(a)
                next = ""
            elif next == "-i2":
                self.introns2file = P.isFile(a)
                next = ""
            elif next == "-j2":
                self.juncs2file = P.isFile(a)
                next = ""
            elif next == "-o":
                self.outfile = a
                next = ""
            elif next == "-fc":
                self.fc = P.toFloat(a)
                next = ""
            elif next == "-t":
                self.thr = P.toFloat(a)
                next = ""
            elif a in ["-i1", "-j1", "-i2", "-j2", "-fc", "-o", "-t"]:
                next = a
            else:
                self.bedfile = a
        if self.bedfile == None or self.introns1file == None or self.introns2file == None:
            P.errmsg(P.NOFILE)
    
    def readFiles(self):
        self.intr1 = readBEDfile(self.introns1file)
        if self.juncs1file != None:
            self.juncs1 = readBEDfile(self.juncs1file)
        self.intr2 = readBEDfile(self.introns2file)
        if self.juncs2file != None:
            self.juncs2 = readBEDfile(self.juncs2file)

    def compare(self, out):
        nin = 0
        nup = 0
        ndown = 0
        maxup = 0
        maxdn = 0
        with open(self.bedfile, "r") as f:
            for line in f:
                nin += 1
                parsed = line.rstrip("\r\n").split("\t")
                intron = parsed[3]
                gene = parsed[4]
                sp = intron.split("_")
                tx = sp[0]
                intid = sp[1]
                iv1 = dget(intron, self.intr1)
                iv2 = dget(intron, self.intr2)
                if self.juncs1file:
                    iv1 = (iv1 + dget(intron + "_a", self.juncs1) + dget(intron + "_b", self.juncs1)) / 3.0
                if self.juncs2file:
                    iv2 = (iv2 + dget(intron + "_a", self.juncs2) + dget(intron + "_b", self.juncs2)) / 3.0

                if iv1 == 0 and iv2 == 0:
                    pass
                elif iv1 == 0:
                    if iv2 >= self.thr:
                        out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(gene, tx, intid, iv1, iv2, "+inf"))
                elif iv2 == 0:
                    if iv1 >= self.thr:
                        out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(gene, tx, intid, iv1, iv2, "-inf"))
                else:
                    l2fc = math.log(iv2/iv1, 2)
                    if abs(l2fc) > self.fc:
                        out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(gene, tx, intid, iv1, iv2, l2fc))
                        if l2fc > 0:
                            nup += 1
                            if l2fc > maxup:
                                maxup = l2fc
                        else:
                            ndown += 1
                            if l2fc < maxdn:
                                maxdn = l2fc
        return (nin, nup, ndown, maxup, maxdn)

if __name__ == "__main__":
    PA = Params()
    PA.parseArgs(sys.argv[1:])
    PA.readFiles()
    if PA.outfile:
        with open(PA.outfile, "w") as out:
            (nin, nup, ndown, maxup, maxdn) = PA.compare(out)
    else:
        (nin, nup, ndown, maxup, maxdn) = PA.compare(sys.stdout)

    sys.stderr.write("{}\t{}\t{}\t{}\t{}\n".format(nin, nup, ndown, maxup, maxdn))
