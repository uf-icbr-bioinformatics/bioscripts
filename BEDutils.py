#!/usr/bin/env python

### Utilities and commands for working with BED files

import sys
import csv
import random
import os.path

### Error classes

class BEDNotIndexed(Exception):
    filename = ""

    def __init__(self, filename):
        self.filename = filename
    
### BED indexer

class BEDindexer():
    filename = ""

    def __init__(self, filename):
        self.filename = filename
        
    def bidx_filename(self):
        return os.path.splitext(self.filename)[0] + ".bidx"

    def loadindex(self):
        bidx = self.bidx_filename()
        if os.path.isfile(bidx):
            with open(bidx, "r") as f:
                data = f.read().split("\n")
            idx = {}
            for d in data:
                if "\t" in d:
                    [chrom, fp] = d.split("\t")
                    idx[chrom] = int(fp)
            return idx
        else:
            raise BEDNotIndexed(self.filename)

    def bedindex(self):
        bidx = self.bidx_filename()
        tag = ""
        sys.stderr.write("{} => {}\n".format(self.filename, bidx))
        with open(bidx, "w") as out:
            with open(self.filename, "r") as f:
                while True:
                    pos = f.tell()
                    line = f.readline()
                    if not line:
                        break
                    if line[0] == '#':
                        continue
                    tp = line.find("\t")
                    if tp < 0:
                        continue
                    v = line[:tp]
                    if v != tag:
                        out.write("{}\t{}\n".format(v, pos))
                        tag = v

### BED reader

class BEDreader():
    filename = None
    stream   = None
    reader   = None
    current  = None             # Payload
    chrom    = ""
    pos      = 0
    idx      = None             # For index, if it exists
    
    def __init__(self, filename, header=False, jump=False):
        self.filename = filename
        self.stream = open(self.filename, "r")
        self.reader = csv.reader(self.stream, delimiter='\t')
        if jump:
            self.jumpTo(jump)
        elif header:
            self.next()
        self.init()

    def __iter__(self):
        return self

    def init(self):
        pass                    # Can be specialized
    
    def close(self):
        self.stream.close()

    def storeCurrent(self, data):
        self.current = data

    def jumpTo(self, chrom):
        if not self.idx:
            BI = BEDindexer(self.filename)
            self.idx = BI.loadindex()
        if chrom in self.idx:
            fp = self.idx[chrom]
            self.stream.seek(fp)
            self.next()
        else:
            sys.stderr.write("Warning: `{}' not found in BED index.\n".format(jump))
    
    def next(self):
        """Read one line from stream and store it in the `current' attribute. Also sets `chrom' 
and `pos' to its first and second elements."""
        if self.stream == None:
            raise StopIteration
        while True:
            data = self.reader.next()
            if data == None:
                self.stream.close()
                self.stream = None
                raise StopIteration
            if data[0][0] == '#':
                continue
            else:
                self.chrom = data[0]
                self.pos   = int(data[1])
                self.storeCurrent(data)
                return self.current

    def skipToChrom(self, chrom):
        """Read lines until finding one that starts with `chrom'."""
        # print("Skipping to chrom {} for {}".format(chrom, self.filename))
        while self.chrom != chrom:
            self.next()
            if self.stream == None:
                break

    def readUntil(self, chrom, limit):
        """Read lines until reaching one that is after `pos' or is on a different chromosome. Returns:
- None if the BED file is finished,
- The new chromosome, if different from chrom,
- The list of records read otherwise.
"""
        result = []
        if self.stream is None:
            return None
        if chrom != self.chrom:
            return self.chrom
        while True:
            if not limit or self.pos < limit:
                result.append(self.current)
                self.next()
                if self.stream == None:
                    break
                if chrom != self.chrom:
                    break
            else:
                break
        return result

    def readChromosome(self):
        """Call next() before the first call to this method."""
        if self.stream is None:
            return None
        self.next()
        result = [self.current]
        thisChrom = self.chrom

        try:
            while True:
                if self.chrom != thisChrom:
                    return result
                result.append(self.current)
                self.next()
            return result
        except StopIteration:
            return result

### Basic BEDreader - a BEDreader that only stores chrom and coords

class BasicBEDreader(BEDreader):

    def init(self):
        self.current = [0, 0]
        
    def storeCurrent(self, data):
        self.current[0] = int(data[1])
        self.current[1] = int(data[2])

### Dual BED reader - for two bed files at once

class DualBEDreader():
    bed1 = None
    bed2 = None
    chrom = ""
    pos = 0
    current1 = None
    current2 = None
    stop = False
    
    def __init__(self, filename1, filename2, readerclass=BEDreader):
        self.bed1 = readerclass(filename1)
        self.bed2 = readerclass(filename2)
        self.bed1.next()
        self.bed2.next()
        if self.bed1.chrom != self.bed2.chrom:
            sys.stderr.write("Error: BED files start on different chromosomes ({}, {}).\n".format(self.bed1.chrom, self.bed2.chrom))
        else:
            self.chrom = self.bed1.chrom

    def __iter__(self):
        return self
    
    def close(self):
        self.bed1.close()
        self.bed2.close()
        
    def next(self):
        """Read the two BED files looking for the next site present in both. Returns True if a site is found,
in which case the chrom, pos, and current attributes are set. Returns False if one of the two files ends."""
        if self.stop:
            raise StopIteration
        while True:
            # Check if one of the two BEDs is now at a different chrom
            if self.bed1.chrom != self.bed2.chrom:
                # If so, advance the other one to the same chrom
                if self.bed1.chrom == self.chrom:
                    self.bed1.skipToChrom(self.bed2.chrom)
                else:
                    self.bed2.skipToChrom(self.bed1.chrom)
                self.chrom = self.bed1.chrom

            # Now look at positions. If they are the same, add to region
            if self.bed1.pos == self.bed2.pos:
                self.chrom = self.bed1.chrom
                self.pos = self.bed1.pos
                self.current1 = self.bed1.current
                self.current2 = self.bed2.current
                try:
                    self.bed1.next()
                    self.bed2.next()
                except StopIteration:
                    self.stop = True # delay StopIteration to next cycle
                return True
            elif self.bed1.pos < self.bed2.pos:
                self.bed1.next()
            else:
                self.bed2.next()
                
### BED grep

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
        if start <= self.start <= end:
            return True
        if self.end:
            if start <= self.end <= end:
                return True
            if self.start <= start <= end <= self.end:
                return True
        return False

def parseSpec(s):
    """Parse a specification of the form: chrom:start-end%name, returning the
four components in a tuple. End and name are optional."""
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

def makeRegion(chrom, start, end, name=None):
    """Create a Region object given the supplied components (start and end can be strings)."""
    try:
        start = int(start)
        if end:
            end = int(end)
        return Region(chrom, start, end, name=name)
    except ValueError:
        return None

class BEDgrep():
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
        B = BEDreader(filename)
        ln = 0
        for line in B:
            ln += 1
            chrom = B.chrom
            try:
                start = int(B.current[1])
                end   = int(B.current[2])
            except ValueError:
                continue
            except IndexError:
                continue
            for r in self.regions:
                if r.chrom == B.chrom and r.matchRegion(start, end):
                    if r.name:
                        postfix = "\t" + r.name
                    else:
                        postfix = ""
                    if self.showLineNumbers:
                        sys.stdout.write("{}:".format(ln))
                    sys.stdout.write("\t".join(line) + postfix + "\n")
                        
    def grepAll(self):
        for f in self.filenames:
            self.grepOne(f)

### BEDreduce - output lines from input BED file with specified probability

class BEDreduce():
    filename = ""
    prob = 1.0

    def __init__(self, filename, prob):
        self.filename = filename
        self.prob = prob

    def run(self):
        nin = 0
        nout = 0
        B = BEDreader(self.filename)
        for line in B:
            nin += 1
            if random.random() <= self.prob:
                sys.stdout.write("\t".join(line) + "\n")
                nout += 1
        sys.stderr.write("{} records read\n{} records written\n".format(nin, nout))

### Main

def run_bedindex():
    BI = BEDindexer(sys.argv[1])
    BI.bedindex()

def run_bedgrep():
    BG = BEDgrep()
    BG.parseArgs(sys.argv[1:])
    if BG.regions and BG.filenames:
        BG.grepAll()

def run_bedreduce():
    BR = BEDreduce(sys.argv[1], float(sys.argv[2]))
    BR.run()
                   
if __name__ == "__main__":
    cmd = os.path.split(sys.argv[0])[1]
    if cmd == "bedindex.py":
        run_bedindex()
    elif cmd == "bedgrep.py":
        run_bedgrep()
    elif cmd == "bedreduce.py":
        run_bedreduce()
    else:
        sys.stderr.write("Please call this script as bedindex.py or bedgrep.py.\n")
        
        
