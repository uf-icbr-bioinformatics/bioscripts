#!/usr/bin/env python

import csv
import Utils

class Region(object):
    chrom = ""
    start = 0
    end = 0
    payload = None

    def __init__(self, chrom, start=0, end=0):
        if start:
            try:
                self.start = int(start)
                self.end   = int(end)
                self.chrom = chrom
            except ValueError:
                return None
        else:
            self.parseSpec(chrom)

    def parseSpec(self, chrom):
        parts = chrom.split(":")
        self.chrom = parts[0]
        parts = parts[1].split("-")
        self.start = int(parts[0])
        self.end = int(parts[1])

    def overlap(self, other):
        """Returns True if this region (partially) overlaps region `other'."""
        return (self.start <= other.start <= self.end) or (self.start <= other.end <= self.end) or (other.start <= self.start and other.end >= self.end)
        
class RegionSource(object):
    sources = []
    chromcol = 0
    endcol = 2
    _nsources = 0
    _idx = 0
    _reader = None

    def __init__(self, sources):
        self.sources = sources
        self._nsources = len(sources)

    def next(self):
        if self._idx >= self._nsources:
            return None
        if self._reader:
                row = self._reader.readline()
                if row:
                    reg = Region(row[self._reader.col], row[self._reader.col+1], row[self._reader.col+self.endcol])
                    if reg:
                        reg.payload = row
                        return reg
                    else:
                        return self.next()
                else:
                    self._reader.close()
                    self._reader = None
                    self._idx += 1
                    return self.next()
        else:
            src = self.sources[self._idx]
            if src[0] == '@':
                self._reader = Utils.AtFileReader(src)
                row = self._reader.readline()
                if row:
                    reg = Region(row[self._reader.col], row[self._reader.col+1], row[self._reader.col+2])
                    reg.payload = row
                    return reg
                else:
                    self._reader.close()
                    self._reader = None
                    self._idx += 1
                    return self.next()
            else:
                self._idx += 1
                return Region(src)

class BEDdict(object):
    bedfile = ""
    d = {}
    nregs = 0

    def __init__(self, bedfile):
        self.bedfile = bedfile
        self.d = {}
        if bedfile:
            with open(bedfile, "r") as f:
                c = csv.reader(f, delimiter='\t')
                for line in c:
                    if line[0][0] == '#':
                        continue
                    chrom = line[0]
                    reg = Region(chrom, line[1], line[2])
                    if chrom in self.d:
                        self.d[chrom].append(reg)
                    else:
                        self.d[chrom] = [reg]
                    self.nregs += 1
            for k in self.d.keys():
                self.d[k].sort(key=lambda r: r.start)

    def allChroms(self):
        return sorted(self.d.keys())

    def chromRegions(self, chrom):
        if chrom in self.d:
            return self.d[chrom]
        else:
            return []
        
    def find(self, reg):
        """Returns the first region in this BEDdict overlapping `reg', or None if not found."""
        if reg.chrom in self.d:
            for r in self.d[reg.chrom]:
                if r.start > reg.end:
                    return None
                elif reg.overlap(r):
                    return r
        return None

    def findPos(self, chrom, pos):
        """Returns the first region in this BEDdict containing position `pos' on chromosome `chrom', or None if not found."""
        if chrom in self.d:
            for r in self.d[chrom]:
                if r.start > pos:
                    return None
                if r.start <= pos <= r.end:
                    return r
        return None

    def setPayloads(self, func):
        """Set all payloads to the value returned by function `func'."""
        for reglist in self.d.values():
            for reg in reglist:
                reg.payload = func()
