#!/usr/bin/env python

import sys
import csv
import os.path

def parseSpec(s):
    name = None
    i0 = s.find("%")
    if i0 > 0:
        name = s[:i0]
        s = s[i0+1:]
    i1 = s.find(":")
    i2 = s.find("-")
    if i1 > 0:
        try:
            chrom = s[:i1]
            if i2 > i1:
                start = int(s[i1+1:i2])
                end   = int(s[i2+1:])
                return (chrom, start, end, name)
            else:
                start = int(s[i1+1:])
                return (chrom, start, None, name)
        except ValueError:
            return None

class Region():
    name = ""
    chrom = ""
    start = 0
    end = 0

    def __init__(self, chrom, start, end=None, name=None):
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end

    def __repr__(self):
        w = "{}:{}".format(self.chrom, self.start)
        if self.name:
            w = self.name + "%" + w
        if self.end:
            w = w + "-" + str(self.end)
        return w

    def matchRegion(self, start, end):
        #print ((start, end, self.start, self.end))
        if start <= self.start <= end:
            return True
        if self.end:
            if start <= self.end <= end:
                return True
            if self.start <= start <= end <= self.end:
                return True
        return False

def makeRegion(chrom, start, end, name=None):
    try:
        start = int(start)
        if end:
            end = int(end)
        return Region(chrom, start, end, name=name)
    except ValueError:
        return None

class BedGrep():
    regions = []
    filenames = []

    showRegions = False
    showLineNumbers = False

    def __init__(self):
        self.regions = []
        self.filenames = []

    def readRegionsFromFile(self, filename):
        with open(filename, "r") as f:
            c = csv.reader(f, delimiter='\t')
            for line in c:
                reg = None
                if len(line) > 3:
                    reg = makeRegion(line[0], line[1], line[2], line[3])
                if len(line) > 2:
                    reg = makeRegion(line[0], line[1], line[2])
                elif len(line) == 2:
                    reg = makeRegion(line[0], line[1])
                if reg:
                    self.regions.append(reg)

    def parseArgs(self, args):
        for a in args:
            if a[0] == '@':
                atfile = a[1:]
                if os.path.isfile(atfile):
                    self.readRegionsFromFile(atfile)
            elif os.path.isfile(a):
                self.filenames.append(a)
            elif a == '-v':
                self.showRegions = True
            elif a == '-n':
                self.showLineNumbers = True
            else:
                spec = parseSpec(a)
                if spec:
                    reg = makeRegion(*spec)
                    self.regions.append(reg)
        if self.showRegions and self.regions:
            sys.stderr.write("Regions:\n")
            for r in self.regions:
                sys.stderr.write(str(r) + "\n")
    
    def grepOne(self, filename):
        with open(filename, "r") as f:
            c = csv.reader(f, delimiter='\t')
            ln = 0
            for line in c:
                ln += 1
                if len(line) == 0:
                    continue
                if line[0][0] == '#':
                    continue
                try:
                    chrom = line[0]
                    start = int(line[1])
                    end   = int(line[2])
                    for r in self.regions:
                        if r.chrom == chrom and r.matchRegion(start, end):
                            #print "match!"
                            if r.name:
                                postfix = "\t" + r.name
                            else:
                                postfix = ""
                            if self.showLineNumbers:
                                sys.stdout.write("{}:".format(ln))
                            sys.stdout.write("\t".join(line) + postfix + "\n")
                except ValueError:
                    pass
                except IndexError:
                    pass
                        
    def grepAll(self):
        for f in self.filenames:
            self.grepOne(f)

if __name__ == "__main__":
    args = sys.argv[1:]
    BG = BedGrep()
    BG.parseArgs(args)
    if BG.regions and BG.filenames:
        BG.grepAll()
