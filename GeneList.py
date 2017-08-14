#!/usr/bin/env python

import os
import sys
import os.path
import sqlite3 as sql

import Utils

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

def parseCoords(c):
    """Parse a pair of coordinates in the form X..Y and return them as ints."""
    dp = c.find(".")
    if dp > 0:
        return (int(c[0:dp]), int(c[dp+2:]))
    else:
        return None

# Classes

class Genelist():
    chroms = []
    genes = {}
    ngenes = 0
    indexes = {}
    btFlags = {}
    currentChrom = ""
    currentGenes = ""
    source = ""

    def __init__(self):
        self.chroms = []
        self.genes = {}
        self.btFlags = {"": True, "*": True}

    def saveAllToDB(self, filename):
        """Save all genes to the database represented by connection `conn'."""
        tot  = 0                # Total number of genes written
        conn = sql.connect(filename)
        with conn:
            for chrom in self.chroms:
                sys.stderr.write("  {}... ".format(chrom))
                n = 0
                for g in self.genes[chrom]:
                    g.saveToDB(conn)
                    n += 1
                sys.stderr.write("{} genes written.\n".format(n))
                tot += n
            conn.execute("INSERT INTO Source (filename) values (?);", (self.source,))
        return tot

    def getGenesTable(self):
        table = {}
        for chrom in self.chroms:
            for gene in self.genes[chrom]:
                table[gene.ID] = {'gene_id': gene.ID,
                                  'gene_biotype': gene.biotype or "???",
                                  'gene_name': gene.name}
        return table

    def getTranscriptsTable(self):
        table = {}
        for chrom in self.chroms:
            for gene in self.genes[chrom]:
                for tx in gene.transcripts:
                    table[tx.ID] = {'transcript_id': tx.ID,
                                    'gene_biotype': gene.biotype or "???",
                                    'transcript_name': tx.name,
                                    'gene_name': gene.name}
        return table

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
                if g.ID == name or g.name == name:
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
                        out.write("\t".join(parsed + [g.name, g.biotype, a[1], Utils.f2dd(a[2]*100), Utils.f2dd(a[3]*100)]) + "\n")

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

class GenelistDB(Genelist):
    conn = None                 # DB connection
    preloaded = False           # Did we preload all genes into the Genelist?

    def allGeneNames(self):
        curr = self.conn.cursor()
        curr.execute("SELECT name FROM Genes ORDER BY name")
        names = [ r[0] for r in curr.fetchall() ]
        curr.close()
        return names

    def findGene(self, name, chrom=None):
        gcur = self.conn.cursor()
        tcur = self.conn.cursor()
        ecur = self.conn.cursor()
        gcur.execute("SELECT ID, name, geneid, ensg, biotype, chrom, strand, start, end FROM Genes WHERE ID=? OR name=? OR geneid=? OR ensg=?", 
                     (name, name, name, name))
        row = gcur.fetchone()
        if row:
            gid = row[0]
            g = Gene(gid, row[5], row[6])
            for pair in zip(['ID', 'name', 'geneid', 'ensg', 'biotype', 'chrom', 'strand', 'start', 'end'], row):
                setattr(g, pair[0], pair[1])
            for trow in tcur.execute("SELECT ID, name, accession, enst, chrom, strand, txstart, txend, cdsstart, cdsend FROM Transcripts WHERE parentID=?", (gid,)):
                tid = trow[0]
                tr = Transcript(tid, trow[4], trow[5], trow[6], trow[7])
                tr.exons = []
                for pair in zip(['ID', 'name', 'accession', 'enst', 'chrom', 'strand', 'txstart', 'txend', 'cdsstart', 'cdsend'], trow):
                    setattr(tr, pair[0], pair[1])
                for erow in ecur.execute("SELECT start, end FROM Exons WHERE ID=? ORDER BY idx", (tid,)):
                    tr.addExon(erow[0], erow[1])
                g.addTranscript(tr)
            return g
        else:
            return None

    def allTranscriptNames(self):
        curr = self.conn.cursor()
        curr.execute("SELECT ID FROM Transcripts ORDER BY ID")
        names = [ r[0] for r in curr.fetchall() ]
        curr.close()
        return names

    def findTranscript(self, name, chrom=None):
        tcur = self.conn.cursor()
        ecur = self.conn.cursor()

        tcur.execute("SELECT ID, name, accession, enst, chrom, strand, txstart, txend, cdsstart, cdsend FROM Transcripts WHERE ID=? OR name=? OR accession=? OR enst=?",
                     (name, name, name, name))
        trow = tcur.fetchone()
        if trow:
            tid = trow[0]
            tr = Transcript(tid, trow[4], trow[5], trow[6], trow[7])
            tr.exons = []
            for pair in zip(['ID', 'name', 'accession', 'enst', 'chrom', 'strand', 'txstart', 'txend', 'cdsstart', 'cdsend'], trow):
                setattr(tr, pair[0], pair[1])
            for erow in ecur.execute("SELECT start, end FROM Exons WHERE ID=? ORDER BY idx", (tid,)):
                tr.addExon(erow[0], erow[1])
            return tr
        else:
            return None

    def getAllTranscripts(self):
        tcur = self.conn.cursor()
        ecur = self.conn.cursor()

        for trow in tcur.execute("SELECT t.ID, t.name, accession, enst, t.chrom, t.strand, txstart, txend, cdsstart, cdsend, g.name FROM Transcripts t, Genes g WHERE t.parentID = g.ID ORDER BY t.chrom, txstart"):
            if trow:
                tid = trow[0]
                tr = Transcript(tid, trow[4], trow[5], trow[6], trow[7])
                tr.gene = trow[10]
                tr.exons = []
                for pair in zip(['ID', 'name', 'accession', 'enst', 'chrom', 'strand', 'txstart', 'txend', 'cdsstart', 'cdsend'], trow):
                    setattr(tr, pair[0], pair[1])
                for erow in ecur.execute("SELECT start, end FROM Exons WHERE ID=? ORDER BY idx", (tid,)):
                    tr.addExon(erow[0], erow[1])
                yield tr

    def findGenes(self, query, args):
        result = []
        qcur = self.conn.cursor()
        gcur = self.conn.cursor()
        tcur = self.conn.cursor()
        ecur = self.conn.cursor()
        for geneIDrow in qcur.execute(query, args):
            geneID = geneIDrow[0]
            row = gcur.execute("SELECT ID, name, geneid, ensg, biotype, chrom, strand, start, end FROM Genes WHERE ID=?", (geneID,)).fetchone()
            if row:
                gid = row[0]
                g = Gene(gid, row[5], row[6])
                for pair in zip(['ID', 'name', 'geneid', 'ensg', 'biotype', 'chrom', 'strand', 'start', 'end'], row):
                    setattr(g, pair[0], pair[1])
                for trow in tcur.execute("SELECT ID, name, accession, enst, chrom, strand, txstart, txend, cdsstart, cdsend FROM Transcripts WHERE parentID=?", (gid,)):
                    tid = trow[0]
                    tr = Transcript(tid, trow[4], trow[5], trow[6], trow[7])
                    tr.exons = []
                    for pair in zip(['ID', 'name', 'accession', 'enst', 'chrom', 'strand', 'txstart', 'txend', 'cdsstart', 'cdsend'], trow):
                        setattr(tr, pair[0], pair[1])
                    for erow in ecur.execute("SELECT start, end FROM Exons WHERE ID=? ORDER BY idx", (tid,)):
                        tr.addExon(erow[0], erow[1])
                    g.addTranscript(tr)
                result.append(g)
        return result

    def allIntersecting(self, chrom, start, end):
        """Returns all genes in `chrom' that intersect the `start-end' region."""
        return self.findGenes("SELECT ID from Genes where chrom=? and ((? <= start) and (start <= ?) or ((? <= end) and (end <= ?)) or ((start <= ?) and (end >= ?)))",
                              (chrom, start, end, start, end, start, end))

# Transcript class

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
    gene = ""                   # Not normally used...

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

    def saveToDB(self, conn, parentID):
        idx = 0
        for ex in self.exons:
            conn.execute("INSERT INTO Exons(ID, idx, chrom, start, end) VALUES (?, ?, ?, ?, ?);",
                         (self.ID, idx, self.chrom, ex[0], ex[1]))
            idx += 1

        try:
            conn.execute("INSERT INTO Transcripts (ID, parentID, name, accession, enst, chrom, strand, txstart, txend, cdsstart, cdsend) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);",
                         (self.ID, parentID, self.name, self.accession, self.enst, self.chrom, self.strand, self.txstart, self.txend, self.cdsstart, self.cdsend))
        except sql.IntegrityError:
            sys.stderr.write("Error: transcript ID {} is not unique.\n".format(self.ID))

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

    def classifyPosition(self, pos, updistance, dndistance):
        """Returns a single character that classifies position `pos' within this transcript.
Possible return values are 'u' (up to `distance' bp upstream of the transcript), 'd'
(up to `distance' bp downstream of the transcript), 'E' (in a coding exon), 'e' (in a
non-coding exon), 'i' (in an intron)."""
        if self.txstart <= pos <= self.txend:
            if self.posInExon(pos):
                if self.cdsstart <= pos <= self.cdsend:
                    return 'E'
                else:
                    return 'e'
            return 'i'
        elif self.strand == 1:
            if self.txstart - updistance <= pos < self.txstart:
                return 'u'
            elif self.txend < pos <= self.txend + dndistance:
                return 'd'
            else:
                return False
        else:
            if self.txstart - dndistance <= pos < self.txstart:
                return 'd'
            elif self.txend < pos <= self.txend + updistance:
                return 'u'
            else:
                return False
            
    def getRegion(self, params):
        """Returns the region for this transcript as (start, end) according to the 'regwanted', 
'updistance, and 'dndistance' attributes in the `params' object."""
        if self.txstart and self.txend:
            if params.regwanted == 'b':
                return (self.txstart, self.txend)
            elif params.regwanted == 'u':
                if self.strand == 1:
                    return (self.txstart - params.updistance, self.txstart + params.dndistance)
                else:
                    return (self.txend - params.dndistance, self.txend + params.updistance)
            elif params.regwanted == 'd':
                if self.strand == 1:
                    return (self.txend - params.updistance, self.txend + params.dndistance)
                else:
                    return (self.txstart - params.dndistance, self.txstart + params.updistance)
        else:
            return None


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

### Gene class

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
        print "ENSG: " + self.ensg
        print "Biotype: " + self.biotype
        print "Region: {}:{}-{}".format(self.chrom, self.start, self.end)
        print "Strand: {}".format(self.strand)
        print "Transcripts:"
        for t in self.transcripts:
            t.dump(prefix="  ", short=True)

    def saveToDB(self, conn):
        for tr in self.transcripts:
            tr.saveToDB(conn, self.ID)
        conn.execute("INSERT INTO Genes (ID, name, geneid, ensg, biotype, chrom, strand, start, end) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);",
                     (self.ID, self.name, self.geneid, self.ensg, self.biotype, self.chrom, self.strand, self.start, self.end))

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

    def classifyPosition(self, position, updistance, dndistance):
        """Returns a string containing all possible classifications of `position' for the transcript of this gene."""
        result = []
        for tr in self.transcripts:
            c = tr.classifyPosition(position, updistance, dndistance)
            if c and c not in result:
                result.append(c)
        if len(result) == 0:
            return "o"
        else:
            return "".join(sorted(result))

    def getRegion(self, params):
        """Returns the region for this gene as (start, end) according to the 'regwanted', 
'updistance, and 'dndistance' attributes in the `params' object."""
        #print self.name, self.start, self.end
        if self.start and self.end:
            if params.regwanted == 'b':
                return (self.start, self.end)
            elif params.regwanted == 'u':
                if self.strand == 1:
                    return (self.start - params.updistance, self.start + params.dndistance)
                else:
                    return (self.end - params.dndistance, self.end + params.updistance)
            elif params.regwanted == 'd':
                if self.strand == 1:
                    return (self.end - params.updistance, self.end + params.dndistance)
                else:
                    return (self.start - params.dndistance, self.start + params.updistance)
        else:
            return None
    
def getLoader(filename):
    """Returns the name of the loader class to be used to load the gene database in `filename'.
Currently recognizes the following extensions: .gb (genbank), .gtf (GTF), .gff/.gff3 (GFF), .db
(sqlite3 database)."""
    loadermap = {'.gb':   GenbankLoader,
                 '.gtf':  GTFloader,
                 '.gff':  GFFloader,
                 '.gff3': GFFloader,
                 '.GTF':  GTFloader,
                 '.GFF':  GFFloader,
                 '.GFF3': GFFloader,
                 '.db':   DBloader,
                 '.csv':  refGeneLoader}

    ext = os.path.splitext(filename)[1]
    if ext in loadermap:
        return loadermap[ext]
    else:
        return None

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

    def load(self, sort=True, index=True, preload=True):
        self._load(preload=preload)
        if sort:
            self.gl.sortGenes()
        if index:
            self.gl.buildIndexes()
        self.gl.source = self.filename
        return self.gl

class refGeneLoader(GeneLoader):
    genes = {}                  # Dictionary of genes by name

    def _load(self, preload=True):
        self.genes = {}
        with open(self.filename, "r") as f:
            reader = Utils.CSVreader(f)
            for line in reader:
                chrom = self.validateChrom(line[2])
                if chrom:
                    name = line[0]
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

    def _load(self, preload=True):
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

    def parseAnnotations(self, ann):
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

    def _load(self, wanted=[], notwanted=[], preload=True):
        """Read genes from a GTF file `filename'. If `wanted' is specified, only include genes whose biotype
is in this list (or missing). If `notwanted' is specified, only include genes whose biotype is not in this list (or missing)."""

        if wanted <> []:
            self.gl.setWanted(wanted)
        elif notwanted <> []:
            self.gl.setNotWanted(notwanted)

        with open(self.filename, "r") as f:
            reader = Utils.CSVreader(f)
            for line in reader:
                strand = 1 if line[6] == '+' else -1
                chrom = self.validateChrom(line[0])
                if not chrom:
                    continue    # skip this entry
                btype = line[2]
                ann = self.parseAnnotations(line[8])

                if btype == 'gene':
                    self.currGene = Gene(ann['gene_id'], chrom, strand) 
                    self.currTranscript = None
                    if 'gene_name' in ann:
                        self.currGene.name = ann['gene_name']
                    else:
                        self.currGene.name = ann['gene_id']
                        self.currGene.biotype = ann['gene_biotype']
                    # print "Adding {} to {}".format(self.currGene.ID, chrom)
                    # raw_input()
                    self.gl.add(self.currGene, chrom)

                elif btype == 'transcript':
                    txid = ann['transcript_id']
                    self.currTranscript = Transcript(txid, chrom, strand, int(line[3]), int(line[4])) # clone gene into transcript
                    #self.currTranscript.name    = self.currGene.name
                    self.currTranscript.geneid  = self.currGene.geneid
                    self.currTranscript.biotype = Utils.dget('transcript_biotype', ann)
                    self.currTranscript.name  = Utils.dget('transcript_name', ann)
                    if txid.startswith("ENST"):
                        self.currTranscript.enst = txid
                    self.currTranscript.exons   = []
                    self.currGene.addTranscript(self.currTranscript)
                elif btype == 'CDS':
                    start = int(line[3])
                    end   = int(line[4])
                    self.currTranscript.setCDS(start, end)
                elif btype == 'exon':
                    start = int(line[3])
                    end   = int(line[4])
                    self.currTranscript.addExon(start, end)

class GFFloader(GeneLoader):

    def parseAnnotations(self, ann):
        anndict = {}
        pieces = [ s.strip(" ") for s in ann.split(";") ]
        for p in pieces:
            pair = p.split("=")
            if len(pair) == 2:
                anndict[pair[0]] = pair[1].strip('"')
        return anndict

    def cleanID(self, idstring):
        if not idstring:
            return ""
        if idstring.startswith("gene:"):
            return idstring[5:]
        if idstring.startswith("transcript:"):
            return idstring[11:]
        return idstring

    def _load(self, wanted=[], notwanted=[], preload=True):
        chrom = ""
        strand = 0
        orphans = 0             # Entries not following their parents

        if wanted <> []:
            self.gl.setWanted(wanted)
        elif notwanted <> []:
            self.gl.setNotWanted(notwanted)

        with open(self.filename, "r") as f:
            reader = Utils.CSVreader(f)
            for line in reader:
                if len(line) < 8:
                    print "|"+line+"|"
                if line[6] == '+':
                    strand = 1
                else:
                    strand = -1
                chrom = self.validateChrom(line[0])
                if not chrom:
                    continue
                tag = line[2]

                if tag in ['gene', 'miRNA_gene', 'lincRNA_gene']:
                    ann   = self.parseAnnotations(line[8])
                    gid   = self.cleanID(Utils.dget('ID', ann))
                    self.currGene = Gene(gid, chrom, strand)
                    self.currGene.name = Utils.dget('Name', ann, "")
                    self.currGene.biotype = Utils.dget('biotype', ann)
                    self.gl.add(self.currGene, chrom)
                    
                elif tag in ['mRNA', 'transcript', 'processed_transcript', 'pseudogenic_transcript', 'pseudogene', 'processed_pseudogene', 'miRNA', 'lincRNA']:
                    ann = self.parseAnnotations(line[8])
                    tid = self.cleanID(Utils.dget('ID', ann))
                    pid = self.cleanID(Utils.dget('Parent', ann))
                    self.currTranscript = Transcript(tid, chrom, strand, int(line[3]), int(line[4]))
                    self.currTranscript.name = Utils.dget('Name', ann, "")
                    self.currTranscript.biotype = Utils.dget('biotype', ann)
                    self.currTranscript.exons = [] # Exons come later in the file
                    if pid == self.currGene.ID:
                        self.currGene.addTranscript(self.currTranscript)
                    else:
                        orphans += 1

                elif tag == 'exon':
                    ann = self.parseAnnotations(line[8])
                    pid = self.cleanID(Utils.dget('Parent', ann))
                    if pid == self.currTranscript.ID:
                        start = int(line[3])
                        end   = int(line[4])
                        self.currTranscript.addExon(start, end)
                    else:
                        orphans += 1

                elif tag == 'CDS':
                    ann = self.parseAnnotations(line[8])
                    pid = self.cleanID(Utils.dget('Parent', ann))
                    if pid == self.currTranscript.ID:
                        start = int(line[3])
                        end   = int(line[4])
                        self.currTranscript.setCDS(start, end)
                    else:
                        orphans += 1

        sys.stderr.write("Orphans: {}\n".format(orphans))

class DBloader(GeneLoader):
    conn = None

    def _load(self, preload=True, wanted=[], notwanted=[]):
        self.gl = GenelistDB()
        self.gl.preloaded = preload
        self.conn = sql.connect(self.filename)
        self.gl.conn = self.conn
        if preload:
            gcur = self.conn.cursor()
            tcur = self.conn.cursor()
            ecur = self.conn.cursor()
            for row in gcur.execute("SELECT ID, name, geneid, ensg, biotype, chrom, strand, start, end FROM Genes"):
                gid = row[0]
                g = Gene(gid, row[5], row[6])
                for pair in zip(['ID', 'name', 'geneid', 'ensg', 'biotype', 'chrom', 'strand', 'start', 'end'], row):
                    setattr(g, pair[0], pair[1])
                self.gl.add(g, g.chrom)
                for trow in tcur.execute("SELECT ID, name, accession, enst, chrom, strand, txstart, txend, cdsstart, cdsend FROM Transcripts WHERE parentID=?", (gid,)):
                    tid = trow[0]
                    tr = Transcript(tid, trow[4], trow[5], trow[6], trow[7])
                    for pair in zip(['ID', 'name', 'accession', 'enst', 'chrom', 'strand', 'txstart', 'txend', 'cdsstart', 'cdsend'], trow):
                        setattr(tr, pair[0], pair[1])
                    for erow in ecur.execute("SELECT start, end FROM Exons WHERE ID=? ORDER BY idx", (tid,)):
                        tr.addExon(erow[0], erow[1])
                    g.addTranscript(tr)
        else:
            row = self.conn.execute("SELECT count(*) FROM Genes").fetchone()
            self.gl.ngenes = row[0]

### Database stuff

def initializeDB(filename):
    """Create a new database in 'filename' and write the Genes, Transcripts, and Exons tables to it."""
    conn = sql.connect(filename)
    try:
        conn.execute("DROP TABLE IF EXISTS Genes;")
        conn.execute("CREATE TABLE Genes (ID varchar primary key, name varchar, geneid varchar, ensg varchar, biotype varchar, chrom varchar, strand int, start int, end int);")
        for field in ['name', 'geneid', 'ensg', 'chrom', 'start', 'end']:
            conn.execute("CREATE INDEX Genes_{} on Genes({});".format(field, field))
        conn.execute("DROP TABLE IF EXISTS Transcripts;")
        conn.execute("CREATE TABLE Transcripts (ID varchar primary key, parentID varchar, name varchar, accession varchar, enst varchar, chrom varchar, strand int, txstart int, txend int, cdsstart int, cdsend int);")
        for field in ['parentID', 'name', 'accession', 'enst', 'chrom', 'txstart', 'txend']:
            conn.execute("CREATE INDEX Trans_{} on Transcripts({});".format(field, field))
        conn.execute("DROP TABLE IF EXISTS Exons;")
        conn.execute("CREATE TABLE Exons (ID varchar, idx int, chrom varchar, start int, end int);")
        for field in ['ID', 'chrom', 'start', 'end']:
            conn.execute("CREATE INDEX Exon_{} on Exons({});".format(field, field))
        conn.execute("DROP TABLE IF EXISTS Source;")
        conn.execute("CREATE TABLE Source (filename VARCHAR, ts TIMESTAMP DEFAULT CURRENT_TIMESTAMP NOT NULL);")
    finally:
        conn.close()

