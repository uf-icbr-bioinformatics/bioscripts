#!/usr/bin/env python

import sys
import csv
import Script

def usage():
    sys.stderr.write("""to be written""")

# Utils

def f2dd(x):
    return "{:.2f}".format(x)

def parseCoords(c):
    """Parse a pair of coordinates in the form X..Y and return them as ints."""
    dp = c.find(".")
    if dp > 0:
        return (int(c[0:dp]), int(c[dp+2:]))
    else:
        return None

def parseStartEnd(s):
    """Parse a range specification in Genbank format. It can contain complement and join operations."""
    cl = len('complement')
    jl = len('join')
    strand = 1
    if s[0:cl] == 'complement':
        strand = -1
        s = s[cl+1:-1]
    if s[0:jl] == 'join':
        s = s[jl+1:-1]
    cp = s.find(",")
    if cp > 0:
        pairs = [ parseCoords(z) for z in s.split(",") ]
        introns = []
        for i in range(len(pairs)-1):
            introns.append((pairs[i][1], pairs[i+1][0]))
        return (pairs[0][0], pairs[-1][1], strand, introns)
    else:
        (st, en) = parseCoords(s)
        return (st, en, strand, None)

# Classes

class CSVreader():
    _reader = None
    ignorechar = '#'

    def __init__(self, stream, delimiter='\t'):
        self._reader = csv.reader(stream, delimiter=delimiter)

    def __iter__(self):
        return self

    def next(self):
        row = self._reader.next()
        if len(row) == 0 or row[0][0] == self.ignorechar:
            return self.next()
        else:
            return row

class Genelist():
    chroms = []
    genes = {}
    ngenes = 0
    indexes = {}
    btFlags = {}
    currentChrom = ""
    currentGenes = ""

    def __init__(self):
        self.chroms = []
        self.genes = {}
        self.btFlags = {"": True, "*": True}

    def setWanted(self, wanted):
        """Set all biotypes in `wanted' to True, all others to False."""
        for w in wanted:
            self.btFlags[w] = True
        self.btFlags["*"] = False

    def setNotWanted(self, notwanted):
        """Set all biotypes in `notwanted' to False, all others to True."""
        for w in notwanted:
            self.btFlags[w] = False
        self.btFlags["*"] = True

    def add(self, gene, chrom):
        bt = gene.biotype
        if bt in self.btFlags:
            w = self.btFlags[bt]
        else:
            w = self.btFlags["*"]
        if w:
            if chrom not in self.chroms:
                self.chroms.append(chrom)
                self.genes[chrom] = []
            self.genes[chrom].append(gene)
            self.ngenes += 1

    def selectChrom(self, chrom):
        if chrom <> self.currentChrom and chrom in self.chroms:
            self.currentChrom = chrom
            self.currentGenes = self.genes[chrom]
        return self.currentGenes

    def genesOnChrom(self, chrom):
        if chrom in self.genes:
            return self.selectChrom(chrom)
        else:
            return []

    def allGeneNames(self):
        result = []
        for (chrom, cgenes) in self.genes.iteritems():
            for cg in cgenes:
                result.append(cg.name)
        return result

    def findGene(self, name, chrom=None):
        if chrom:
            cgenes = self.genes[chrom]
            for g in cgenes:
                if g.mrna == name or g.name == name:
                    return g
            return None
        else:
            for ch in self.chroms:
                g = self.findGene(name, chrom=ch)
                if g:
                    return g
            return None

    def sortGenes(self):
        for chrom in self.chroms:
            self.genes[chrom].sort(key=lambda g:g.start)

    def buildIndexes(self):
        step = 100
        idxs = {}
        for chrom, genes in self.genes.iteritems():
            ng = len(genes)
            d = []
            # print "chr={}, genes={}".format(chrom, len(genes))
            for i in range(0, len(genes), step):
                i2 = min(i+step, ng) - 1
                # print "i1={}, i2={}".format(i, i+step-1)
                gfirst = genes[i]
                glast  = genes[i2]
                d.append([gfirst.start, glast.end, i, i2])
            idxs[chrom] = d
            # print d
            # raw_input()
        self.indexes = idxs

    def positionsToRange(self, chrom, start, end):
        """Returns the first and last index for genes in the range `start' to `end'."""
        first = last = 0
        if chrom in self.indexes:
            idxs = self.indexes[chrom]
            for i in range(0, len(idxs)):
                iblock = idxs[i]
                if iblock[0] <= start <= iblock[1]:
                    sb = iblock
                    eb = iblock
                    if i > 0:
                        sb = idxs[i-1]
                    if i < len(idxs) - 1:
                        eb = idxs[i+1]
                    first = sb[2]
                    last  = eb[3]
                    break
        return (first, last)

    def classifyIntersection(self, astart, aend, g):
        """A is mine, B is other."""
        bstart = g.txstart
        bend   = g.txend
        if astart < bstart:
            if bstart <= aend <= bend:
                how = 'left'
                alen = aend-astart+1
                blen = bend-bstart+1
                isiz = aend-bstart+1
                return (g, how, 1.0*isiz/alen, 1.0*isiz/blen)
            elif aend > bend:
                how = 'contains'
                alen = aend-astart+1
                isiz = bend-bstart+1
                return (g, how, 1.0*isiz/alen, 1.0)
        elif bstart <= astart <= bend:
            if aend <= bend:
                how = 'contained'
                blen = bend-bstart+1
                isiz = aend-astart+1
                return (g, how, 1.0, 1.0*isiz/blen)
            else:
                how = 'right'
                alen = aend-astart+1
                blen = bend-bstart+1
                isiz = bend-astart+1
                return (g, how, 1.0*isiz/alen, 1.0*isiz/blen)
        return None

    def allIntersecting(self, chrom, start, end):
        """Returns all genes in `chrom' that intersect the `start-end' region."""
        result = []
        self.selectChrom(chrom)
        genes = self.currentGenes
        (first, last) = self.positionsToRange(chrom, start, end)
        # print "Looking at {}-{}".format(first, last)
        for i in range(first, last+1):
            ix = self.classifyIntersection(start, end, genes[i])
            if ix:
                result.append(ix)
        return result
        
    def intersectFromBED(self, bedfile, outfile):
        with open(outfile, "w") as out:
            with open(bedfile, "r") as f:
                f.readline()
                for line in f:
                    parsed = line.rstrip("\n\r").split("\t")
                    allint = self.allIntersecting(parsed[0], int(parsed[1]), int(parsed[2]))
                    for a in allint:
                        g = a[0]
                        out.write("\t".join(parsed + [g.name, g.biotype, a[1], f2dd(a[2]*100), f2dd(a[3]*100)]) + "\n")

    def notIntersectFromBED(self, bedfile, outfile):
        with open(outfile, "w") as out:
            with open(bedfile, "r") as f:
                f.readline()
                for line in f:
                    parsed = line.rstrip("\n\r").split("\t")
                    s = int(parsed[1])
                    e = int(parsed[2])
                    allint = self.allIntersecting(parsed[0], s, e)
                    if allint == []:
                        out.write("\t".join(parsed) + "\t{}\n".format(1+e-s))

    def intronsToBED(self, bedfile):
        with open(bedfile, "w") as out:
            for chrom in self.chroms:
                for gene in self.genes[chrom]:
                    intid = 1
                    prevexon = gene.exons[0]
                    for ex in gene.exons[1:]:
                        start = prevexon[1] + 1
                        end = ex[0] - 1
                        if gene.strand == 1:
                            out.write("{}\t{}\t{}\t{}_{}\t{}\t{}\n".format(chrom, start, end, gene.mrna, intid, gene.name, gene.strand))
                        else:
                            out.write("{}\t{}\t{}\t{}_{}\t{}\t{}\n".format(chrom, end, start, gene.mrna, intid, gene.name, gene.strand))
                            # print [chrom, end, start, gene.mrna, intid, gene.name, gene.strand]
                            # raw_input()
                        intid += 1
                        prevexon = ex

    def junctionsToBED(self, bedfile, size=5):
        with open(bedfile, "w") as out:
            for chrom in self.chroms:
                for gene in self.genes[chrom]:
                    intid = 1
                    for ex in gene.exons:
                        out.write("{}\t{}\t{}\t{}_{}_a\t{}\t{}\n".format(chrom, ex[0]-10, ex[0]+10, gene.mrna, intid-1, gene.name, gene.strand))
                        out.write("{}\t{}\t{}\t{}_{}_b\t{}\t{}\n".format(chrom, ex[1]-10, ex[1]+10, gene.mrna, intid, gene.name, gene.strand))
                        intid += 1

# Gene class

class Transcript():
    ID = ""
    name = ""
    chrom = ""
    strand = ""
    accession = ""
    enst = ""
    txstart = 0
    txend = 0
    cdsstart = None
    cdsend = None
    strand = None
    exons = []
    smallrects = []
    largerects = []

    def __init__(self, ID, chrom, strand, txstart, txend):
        self.ID = ID
        self.chrom = chrom
        self.strand = strand
        self.txstart = txstart
        self.txend = txend
        self.exons = [(txstart, txend)] # We initially assume transcript has no introns
        self.smallrects = []
        self.largerects = []

    def dump(self, prefix="", short=False):
        if short:
            print "{}{} {}:{}-{} {}".format(prefix, self.ID, self.chrom, self.txstart, self.txend, self.exons)
        else:
            print prefix + "ID: " + self.ID
            print prefix + "Chrom: " + self.chrom
            print "{}Strand: {}".format(prefix, self.strand)
            print "{}Transcript: {}-{}".format(prefix, self.txstart, self.txend)
            print "{}CDS: {}-{}".format(prefix, self.cdsstart, self.cdsend)
            print "{}Exons: {}".format(prefix, self.exons)

    def addExon(self, start, end):
        self.exons.append((start, end))

    def setIntrons(self, introns):
        self.exons = [(self.txstart, introns[0][0])]
        for i in range(len(introns)-1):
            self.exons.append((introns[i][1], introns[i+1][0]))
        self.exons.append((introns[-1][1], self.txend))

    def setCDS(self, cdsstart, cdsend):
        """Set the CDS of this transcript to `cdsstart' and `cdsend'. This also sets the
smallrects and largerects lists."""
        self.cdsstart = cdsstart
        self.cdsend   = cdsend
        self.smallrects = []    # we're recomputing them
        self.largerects = []    # from scratch
        small = True
        for e in self.exons:
            if (e[0] <= self.cdsstart < self.cdsend <= e[1]): # CDS entirely contained in exon?
                self.smallrects.append((e[0], self.cdsstart))
                self.largerects.append((self.cdsstart, self.cdsend))
                self.smallrects.append((self.cdsend, e[1]))
            elif (e[0] <= self.cdsstart <= e[1]):             # Exon contains start of CDS? 
                self.smallrects.append((e[0], self.cdsstart))
                self.largerects.append((self.cdsstart, e[1]))
                small = False
            elif (e[0] <= self.cdsend <= e[1]):               # Exon contains end of CDS?
                self.largerects.append((e[0], self.cdsend))
                self.smallrects.append((self.cdsend, e[1]))
                small = True
            elif small:
                self.smallrects.append(e)
            else:
                self.largerects.append(e)

    def posInExon(self, pos):
        """Return True if position `pos' is in one of the exons of this transcript."""
        for e in self.exons:
            if e[0] <= pos <= e[1]:
                return True
        return False

    def positionMatch(self, pos, mode, distance):
        """Return True if position `pos' matches transcript according to `mode'.
`mode' can be one of b, P, d, e, or i. `distance' is used when `mode' includes p or d."""
        for m in mode:
            if m == 'b':
                if self.txstart <= pos <= self.txend:
                    return True
            elif m == 'p':
                if self.strand == 1:
                    if self.txstart - distance <= pos < self.txstart:
                        return True
                else:
                    if self.txend < pos <= self.txend + distance:
                        return True
            elif m == 'd':
                if self.strand == -1:
                    if self.txstart - distance <= pos < self.txstart:
                        return True
                else:
                    if self.txend < pos <= self.txend + distance:
                        return True
            else:
                ex = self.posInExon(pos)
                if m == 'e' and ex:
                    return True
                if m == 'i' and not ex:
                    return True
        return False

    def classifyPosition(self, pos, distance):
        """Returns a single character that classifies position `pos' within this transcript.
Possible return values are 'p' (up to `distance' bp upstream of the transcript), 'd'
(up to `distance' bp downstream of the transcript), 'E' (in a coding exon), 'e' (in a
non-coding exon), 'i' (in an intron)."""
        if pos < self.txstart - distance or pos > self.txend + distance:
            return False

        if pos < self.txstart:
            if self.strand == 1:
                return 'p'
            else:
                return 'd'
        if pos > self.txend:
            if self.strand == 1:
                return 'd'
            else:
                return 'p'
        if self.posInExon(pos):
            if self.cdsstart <= pos <= self.cdsend:
                return 'E'
            else:
                return 'e'
        return 'i'

    def geneLimits(self, upstream, downstream, ref='B'):
        """Returns the start and end coordinates of a region extending from `upstream' bases
upstream of the TSS to `downstream' bases downstream of the TSE. If `ref' is equal to S, both
coordinates are computed according to the TSS, and if it is 'E', both coordinates are computed
according to the TSE. Takes strand of transcript into account."""
        if ref == 'B':
            if self.strand == 1:
                return (self.txstart - upstream, self.txend + downstream)
            else:
                return (self.txstart - downstream, self.txend + upstream)
        elif ref == 'S':
            if self.strand == 1:
                return (self.txstart - upstream, self.txstart + downstream)
            else:
                return (self.txend - downstream, self.txend + upstream)
        elif ref == 'E':
            if self.strand == 1:
                return (self.txend - upstream, self.txend + downstream)
            else:
                return (self.txend - downstream, self.txend + upstream)

    def distanceFrom(self, position, ref='S'):
        """Returns the distance from `position' to the TSS (if ref=S) or the TSE (if ref=E).
The distance is negative if `position' is upstream of the reference point, positive if
downstream."""
        if ref == 'S':
            if self.strand == 1:
                return position - self.txstart
            else:
                return self.txend - position
        else:
            if self.strand == 1:
                return position - self.txend
            else:
                return self.txend - position

class Gene():
    ID = ""
    name = ""
    geneid = ""
    ensg = ""
    biotype = ""
    chrom = ""
    strand = ""
    start = None                # leftmost txstart
    end = None                  # rightmost txend
    transcripts = []
    data = []

    def __init__(self, ID, chrom, strand):
        self.ID          = ID
        self.chrom       = chrom
        self.strand      = strand
        self.transcripts = []
        self.data        = []

    def dump(self):
        print "ID: " + self.ID
        print "Gene: " + self.name
        print "GeneID: " + self.geneid
        print "Chrom: " + self.chrom
        print "Strand: {}".format(self.strand)
        print "Transcripts:"
        for t in self.transcripts:
            t.dump(prefix="  ", short=True)

    def addTranscript(self, transcript):
        self.transcripts.append(transcript)
        if self.start:
            self.start = min(self.start, transcript.txstart)
        else:
            self.start = transcript.txstart
        if self.end:
            self.end = max(self.end, transcript.txend)
        else:
            self.end = transcript.txend

class GeneLoader():
    gl = None
    filename = ""
    currGene = None
    currTranscript = None

    def __init__(self, filename):
        self.filename = filename
        self.gl = Genelist()

    def validateChrom(self, chrom):
        if chrom.find("_") > 0:
            return False
        if chrom.find("random") > 0:
            return False
        if not chrom.startswith('chr') or chrom.startswith('Chr'):
            chrom = "chr" + chrom
        return chrom

    def load(self, sort=True, index=True):
        self._load()
        if sort:
            self.gl.sortGenes()
        if index:
            self.gl.buildIndexes()
        return self.gl

class refGeneLoader(GeneLoader):
    genes = {}                  # Dictionary of genes by name

    def _load(self):
        self.genes = {}
        with open(self.filename, "r") as f:
            reader = CSVreader(f)
            for line in reader:
                chrom = self.validateChrom(line[2])
                if chrom:
                    name = line[12]
                    if name in self.genes:
                        self.currGene = self.genes[name]
                    else:
                        self.currGene = self.genes[name] = Gene(name, chrom, line[3])
                        self.currGene.name = name
                        self.gl.add(self.currGene, chrom)
                    transcript = Transcript(line[1], chrom, line[3], int(line[4]), int(line[5]))
                    transcript.accession = line[1]
                    transcript.exons = zip(line[9].rstrip(",").split(","), line[10].rstrip(",").split(","))
                    transcript.setCDS(int(line[6]), int(line[7]))
                    self.currGene.addTranscript(transcript)

class GenbankLoader(GeneLoader):

    def _load(self):
        chrom = ""
        thisGene = None
        start = None
        end = None
        strand = None
        cdsstart = None
        cdsend = None
        section = ""

        infeatures = False
        with open(self.filename, "r") as f:
            for line in f:
                line = line.rstrip("\r\n")
                if infeatures:
                    key = line[0:20].strip()
                    if key == '':   # still in previous section?
                        if section == 'gene':
                            data = line[21:]
                            # print "In gene section: {}".format(data)
                            if data[0:5] == '/gene':
                                # print "Setting name to {}".format(data[7:-1])
                                thisGene.name = thisGene.ID = data[7:-1]
                            elif data[0:10] == '/locus_tag':
                                if thisGene.name == '':
                                    thisGene.name = data[12:-1]
                                if thisGene.ID == '':
                                    thisGene.ID = data[12:-1]
                    elif key == 'gene':
                        # print "Found new gene"
                        # raw_input()
                        section = key
                        start, end, strand, ignore = parseStartEnd(line[21:])
                        # print "New gene at ({}, {})".format(start, end)
                        thisGene = Gene("", chrom, strand)
                        thisGene.chrom = chrom
                        self.gl.add(thisGene, chrom)
                        thisGene.addTranscript(Transcript("", chrom, strand, start, end))
                    elif key == 'CDS':
                        cdsstart, cdsend, ignore, introns = parseStartEnd(line[21:])
                        # print "Setting CDS to ({}, {})".format(cdsstart, cdsend)
                        if introns:
                            # print "Setting introns: {}".format(introns)
                            thisGene.transcripts[0].setIntrons(introns)

                        thisGene.transcripts[0].setCDS(cdsstart, cdsend)
                        # print thisGene.smallrects
                        # print thisGene.largerects
                    else:
                        section = ''
                elif line[0:9] == 'ACCESSION': # 'LOCUS':
                    chrom = line[12:line.find(" ", 12)]
                elif line[0:8] == 'FEATURES':
                    infeatures = True

class GTFloader(GeneLoader):

    def parseAnnotationsGTF(self, ann):
        """Parse GTF annotations `ann' and return them as a dictionary."""
        anndict = {}
        pieces = [ s.strip(" ") for s in ann.split(";") ]
        for p in pieces:
            f = p.find(" ")
            if f > 0:
                key = p[0:f]
                val = p[f+1:].strip('"')
                anndict[key] = val
        return anndict

    def _load(wanted=[], notwanted=[]):
        """Read genes from a DTF file `filename'. If `wanted' is specified, only include genes whose biotype
is in this list (or missing). If `notwanted' is specified, only include genes whose biotype is not in this list (or missing)."""
        g = None                # Current gene
        gt = None               # Current transcript

        if wanted <> []:
            self.gl.setWanted(wanted)
        elif notwanted <> []:
            self.gl.setNotWanted(notwanted)

        with open(filename, "r") as f:
            reader = CSVreader(f)
            for line in f:
                if line[6] == '+':
                    strand = 1
                else:
                    strand = -1
                chrom = self.validateChrom(line[0])
                btype = line[2]
                if btype == 'gene':
                    if False: # gt:
                        gt.setCDS(g.cdsstart, g.cdsend)
                        g.dump()
                        raw_input()
                        genes.add(gt, gt.chrom)
                    g = Gene([], chrom=chrom, strand=strand) 
                    gt = None
                    ann = parseAnnotationsGTF(line[8])
                    if 'gene_name' in ann:
                        g.name = ann['gene_name']
                    else:
                        g.name = ann['gene_id']
                    g.geneid = ann['gene_id']
                    g.biotype = ann['gene_biotype']
                elif btype == 'transcript':
                    if gt:
                        gt.setCDS(gt.cdsstart, gt.cdsend)
                        # gt.dump()
                        # raw_input()
                    gt = Gene([], chrom=g.chrom, strand=g.strand) # clone gene into transcript
                    gt.name = g.name
                    gt.geneid = g.geneid
                    gt.biotype = g.biotype
                    gt.txstart = int(line[3])
                    gt.txend   = int(line[4])
                    ann = parseAnnotationsGTF(line[8])
                    gt.mrna = ann['transcript_id']
                    if 'transcript_name' in ann:
                        gt.txname = ann['transcript_name']
                    genes.add(gt, gt.chrom)
                elif btype == 'CDS':
                    start = int(line[3])
                    end   = int(line[4])
                    if gt.cdsstart == None:
                        gt.cdsstart = start
                    gt.cdsend = end
                elif btype == 'exon':
                    start = int(line[3])
                    end   = int(line[4])
                    # print (start, end)
                    gt.exons.append((start, end))
        genes.add(gt, gt.chrom)

## Read genes from a Genbank file

def readGenesFromGenbank(filename):
    """Parse a genbank file `filename' and return a list of all the genes it contains."""
    genes = Genelist()
    chrom = ""
    thisGene = None
    start = None
    end = None
    strand = None
    cdsstart = None
    cdsend = None
    section = ""

    infeatures = False
    with open(filename, "r") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if infeatures:
                key = line[0:20].strip()
                if key == '':   # still in previous section?
                    if section == 'gene':
                        data = line[21:]
                        # print "In gene section: {}".format(data)
                        if data[0:5] == '/gene':
                            # print "Setting name to {}".format(data[7:-1])
                            thisGene.name = data[7:-1]
                        elif data[0:10] == '/locus_tag' and thisGene.name == '':
                            thisGene.name = data[12:-1]
                elif key == 'gene':
                    # print "Found new gene"
                    # raw_input()
                    section = key
                    start, end, strand, ignore = parseStartEnd(line[21:])
                    # print "New gene at ({}, {})".format(start, end)
                    thisGene = Gene([(start, end)], strand=strand)
                    thisGene.chrom = chrom
                    genes.add(thisGene, chrom)
                elif key == 'CDS':
                    cdsstart, cdsend, ignore, introns = parseStartEnd(line[21:])
                    # print "Setting CDS to ({}, {})".format(cdsstart, cdsend)
                    if introns:
                        # print "Setting introns: {}".format(introns)
                        thisGene.setIntrons(introns)
                    
                    thisGene.setCDS(cdsstart, cdsend)
                    # print thisGene.smallrects
                    # print thisGene.largerects
                else:
                    section = ''
            elif line[0:9] == 'ACCESSION': # 'LOCUS':
                chrom = line[12:line.find(" ", 12)]
            elif line[0:8] == 'FEATURES':
                infeatures = True
    genes.sortGenes()
    genes.buildIndexes()
    return genes

## Read genes from a refFlat file (UCSC genome browser)

def parseRefExons(starts, ends):
    starts = starts.rstrip(",").split(",")
    ends   = ends.rstrip(",").split(",")
    return [ (int(s), int(e)) for (s,e) in zip(starts, ends) ]

def readGenesFromRefFlat(filename):
    genes = Genelist()

    with open(filename, "r") as f:
        for line in f:
            line = line.rstrip("\r\n").split("\t")
            if line[3] == '+':
                strand = 1
            else:
                strand = -1
            g = Gene(parseRefExons(line[9], line[10]), strand=strand)
            g.name = line[0]
            g.mrna = line[1]
            g.chrom = line[2]
            g.setCDS(int(line[6]), int(line[7]))
            genes.add(g, g.chrom)
    genes.sortGenes()
    genes.buildIndexes()
    return genes

## Read genes from a GFF3 file

def parseAnnotations(ann):
    anndict = {}
    pieces = [ s.strip(" ") for s in ann.split(";") ]
    for p in pieces:
        pair = p.split("=")
        if len(pair) == 2:
            anndict[pair[0]] = pair[1].strip('"')
    return anndict

def readGenesFromGFF3(filename):
    genes = Genelist()

    with open(filename, "r") as f:
        g = None
        for line in f:
            line = line.rstrip("\r\n").split("\t")
            if line[6] == '+':
                strand = 1
            else:
                strand = -1
            chrom = line[0]
            btype = line[2]
            if btype == 'gene':
                if g:
                    g.setCDS(g.cdsstart, g.cdsend)
                    # g.dump()
                    # raw_input()
                    if g.txstart > 0:
                        genes.add(g, g.chrom)
                g = Gene([], chrom=chrom, strand=strand) 
                ann = parseAnnotations(line[8])
                g.name = ann['ID'][5:]
            elif btype == 'mRNA':
                start = int(line[3])
                end   = int(line[4])
                g.txstart = start
                g.txend = end
                ann = parseAnnotations(line[8])
                g.mrna = ann['ID'][11:]
            elif btype == 'CDS':
                start = int(line[3])
                end   = int(line[4])
                if g.cdsstart == None:
                    g.cdsstart = start
                g.cdsend = end
            elif btype == 'exon':
                start = int(line[3])
                end   = int(line[4])
                g.exons.append((start, end))
        if g.txstart > 0:
            genes.add(g, g.chrom)
    genes.buildIndexes()
    return genes

def parseAnnotationsGTF(ann):
    anndict = {}
    pieces = [ s.strip(" ") for s in ann.split(";") ]
    for p in pieces:
        f = p.find(" ")
        if f > 0:
            key = p[0:f]
            val = p[f+1:].strip('"')
            anndict[key] = val
    return anndict

def readGenesFromGTF(filename, wanted=[], notwanted=[]):
    """Read genes from the GTF file `filename' and return a Genelist object
containing them. If `wanted' is specified, only include genes whose biotype
is in this list (or missing). If `notwanted' is specified, only include genes
whose biotype is not in this list (or missing)."""
    genes = Genelist()
    if wanted <> []:
        genes.setWanted(wanted)
    elif notwanted <> []:
        genes.setNotWanted(notwanted)

    with open(filename, "r") as f:
        g = None                # Current gene
        gt = None               # Current transcript
        for line in f:
            if len(line) > 0 and line[0] <> '#':
                line = line.rstrip("\r\n").split("\t")
                if line[6] == '+':
                    strand = 1
                else:
                    strand = -1
                chrom = line[0]
                if not chrom.startswith('chr') or chrom.startswith('Chr'):
                    chrom = "chr" + chrom
                btype = line[2]
                if btype == 'gene':
                    if False: # gt:
                        gt.setCDS(g.cdsstart, g.cdsend)
                        g.dump()
                        raw_input()
                        genes.add(gt, gt.chrom)
                    g = Gene([], chrom=chrom, strand=strand) 
                    gt = None
                    ann = parseAnnotationsGTF(line[8])
                    if 'gene_name' in ann:
                        g.name = ann['gene_name']
                    else:
                        g.name = ann['gene_id']
                    g.geneid = ann['gene_id']
                    g.biotype = ann['gene_biotype']
                elif btype == 'transcript':
                    if gt:
                        gt.setCDS(gt.cdsstart, gt.cdsend)
                        # gt.dump()
                        # raw_input()
                    gt = Gene([], chrom=g.chrom, strand=g.strand) # clone gene into transcript
                    gt.name = g.name
                    gt.geneid = g.geneid
                    gt.biotype = g.biotype
                    gt.txstart = int(line[3])
                    gt.txend   = int(line[4])
                    ann = parseAnnotationsGTF(line[8])
                    gt.mrna = ann['transcript_id']
                    if 'transcript_name' in ann:
                        gt.txname = ann['transcript_name']
                    genes.add(gt, gt.chrom)
                elif btype == 'CDS':
                    start = int(line[3])
                    end   = int(line[4])
                    if gt.cdsstart == None:
                        gt.cdsstart = start
                    gt.cdsend = end
                elif btype == 'exon':
                    start = int(line[3])
                    end   = int(line[4])
                    # print (start, end)
                    gt.exons.append((start, end))
        genes.add(gt, gt.chrom)
    genes.buildIndexes()
    return genes

