#!/usr/bin/env python

import os
import sys
import os.path
import sqlite3 as sql

import Utils

## DB Utilities

def selectValue(conn, query, *args):
    r = conn.execute(query.format(*args))
    v = r.fetchone()
    if v:
        return v[0]
    else:
        return None

def selectColumn(conn, query, *args):
    result = []
    r = conn.execute(query.format(*args))
    while True:
        v = r.fetchone()
        if v:
            result.append(v[0])
        else:
            break
    return result

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

# Classification of positions

class Classif(object):
    idx = 0
    name = ""
    subject = None
    extra = ""

    def __init__(self, subject=None, extra=""):
        self.subject = subject
        self.extra = extra

class CL_NONE(Classif):
    idx = 0
    name = 'o'

class CL_CODING(Classif):
    idx = 90
    name = 'E'

class CL_EXON(Classif):
    idx = 80
    name = 'e'

class CL_UPSTREAM(Classif):
    idx = 70
    name = 'u'

class CL_INTRON(Classif):
    idx = 60
    name = 'i'

class CL_DNSTREAM(Classif):
    idx = 50
    name = 'd'

class CL_ENHANCER(Classif):
    idx = 65
    name = 'h'

def addClassification(cls, clist):
    """Add classification `cls' to list `clist' only if
it does not already appear."""
    for x in clist:
        if x.name == cls.name:
            if cls.name in "Ee":
                if cls.extra == x.extra:
                    return clist
            else:
                return clist
    clist.append(cls)
    return clist

def sortClassifications(cls):
    cls.sort(key=lambda c: c.idx, reverse=True)
    return cls

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

    # Next two methods allow a genelist to be used in a with statement
    # Useful for GenelistDB.
    def __enter__(self):
        return self

    def __exit__(self, a, b, c):
        pass

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
            self.setTSS(conn)
            self.setCanonical(conn)
            self.setCounts(conn)
            conn.execute("INSERT INTO Source (filename) values (?);", (self.source,))
        return tot

    def setTSS(self, conn):
        sys.stderr.write("Setting TSS...\n")
        conn.execute("UPDATE Genes SET tss=start WHERE strand='1';")
        conn.execute("UPDATE Genes SET tss=end WHERE strand='-1';")
        conn.execute("UPDATE Transcripts SET tss=txstart WHERE strand='1';")
        conn.execute("UPDATE Transcripts SET tss=txend WHERE strand='-1';")

    def setCanonical(self, conn):
        sys.stderr.write("Setting Canonical...\n")
        c = conn.cursor()
        for r in conn.execute("SELECT ID FROM Genes;"):
            gid = r[0]
            maxtr = 0
            best  = ""
            for tr in c.execute("SELECT ID, txstart, txend FROM Transcripts WHERE parentID=?;", (gid,)):
                m = int(tr[2]) - int(tr[1])
                if m > maxtr:
                    maxtr = m
                    best = tr[0]
            c.execute("UPDATE Genes SET canonical=? WHERE ID=?;", (best, gid))
            c.execute("UPDATE Transcripts SET canonical='Y' WHERE ID=?;", (best,))

    def setCounts(self, conn):
        sys.stderr.write("Updating Counts...\n")
        ngenes = selectValue(conn, "SELECT count(*) from Genes;")
        ntrans = selectValue(conn, "SELECT count(*) from Transcripts;")
        nexons = selectValue(conn, "SELECT count(*) from Exons;")
        conn.execute("INSERT INTO Counts (ngenes, ntranscripts, nexons) VALUES (?, ?, ?);", (ngenes, ntrans, nexons))

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
        if chrom != self.currentChrom and chrom in self.chroms:
            self.currentChrom = chrom
            self.currentGenes = self.genes[chrom]
        return self.currentGenes

    def genesOnChrom(self, chrom):
        if chrom in self.genes:
            l = self.selectChrom(chrom)
            return l
        else:
            return []

    def allGeneNames(self):
        result = []
        for (chrom, cgenes) in Utils.get_iterator(self.genes):
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
        for chrom, genes in Utils.get_iterator(self.genes):
            ng = len(genes)
            d = []
            # print("chr={}, genes={}".format(chrom, len(genes)))
            for i in range(0, len(genes), step):
                i2 = min(i+step, ng) - 1
                # print("i1={}, i2={}".format(i, i+step-1))
                gfirst = genes[i]
                glast  = genes[i2]
                d.append([gfirst.start, glast.end, i, i2])
            idxs[chrom] = d
            # print(d)
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
        # print("Looking at {}-{}".format(first, last))
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
                            # print([chrom, end, start, gene.mrna, intid, gene.name, gene.strand])
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
    dbname = None
    dbconn = None
    _level = 0
    dbcurs = None
    preloaded = False           # Did we preload all genes into the Genelist?
    genesTable = {}
    transcriptsTable = {}

    def __enter__(self):
        self._level += 1
        if not self.dbconn:
            # sys.stderr.write("Opening connection to db {}\n".format(self.dbname))
            self.dbconn = sql.connect(self.dbname)
        return self.dbconn

    def __exit__(self, exc_type, exc_value, exc_tb):
        if self.dbconn:
            self._level -= 1
            if self._level == 0:
                # sys.stderr.write("Closing db connection.\n")
                self.dbconn.close()
                self.dbconn = None

    def getGenesTable(self):
        if not self.genesTable:
            with self:
                for row in self.dbconn.execute("SELECT ID, biotype, name FROM genes;"):
                    self.genesTable[row[0]] = {'gene_id': row[0],
                                               'gene_biotype': row[1] or "???",
                                               'gene_name': row[2]}
        return self.genesTable

    def getTranscriptsTable(self):
        if not self.transcriptsTable:
            with self:
                for row in self.dbconn.execute("SELECT transcripts.ID, genes.biotype, transcripts.name, genes.name FROM transcripts, genes WHERE genes.ID=transcripts.parentID;"):
                    self.transcriptsTable[row[0]] = {'transcript_id': row[0],
                                                     'gene_biotype': row[1] or "???",
                                                     'transcript_name': row[2],
                                                     'gene_name': row[3]}
        return self.transcriptsTable

    def allGeneNames(self):
        names = []
        with self:
            #conn = sql.connect(self.dbname)
#        try:
            curr = self.dbconn.execute("SELECT name FROM Genes ORDER BY name")
            names = [ r[0] for r in curr.fetchall() ]
#        finally:
#            conn.close()
        return names

    def findGene(self, name, chrom=None):
        """Returns the gene called `name'. The name is matched against the ID, name, geneid, and ensg fields."""
        with self:
            row = self.dbconn.execute("SELECT ID, name, geneid, ensg, biotype, chrom, strand, start, end FROM Genes WHERE ID=? OR name=? OR geneid=? OR ensg=?", 
                                      (name, name, name, name)).fetchone()
            if row:
                gid = row[0]
                g = Gene(gid, row[5], row[6])
                for pair in zip(['ID', 'name', 'geneid', 'ensg', 'biotype', 'chrom', 'strand', 'start', 'end'], row):
                    setattr(g, pair[0], pair[1])
                for trow in self.dbconn.execute("SELECT ID, name, accession, enst, chrom, strand, txstart, txend, cdsstart, cdsend FROM Transcripts WHERE parentID=?", (gid,)):
                    tid = trow[0]
                    tr = Transcript(tid, trow[4], trow[5], trow[6], trow[7])
                    tr.exons = []
                    for pair in zip(['ID', 'name', 'accession', 'enst', 'chrom', 'strand', 'txstart', 'txend', 'cdsstart', 'cdsend'], trow):
                        setattr(tr, pair[0], pair[1])
                    for erow in self.dbconn.execute("SELECT start, end FROM Exons WHERE ID=? ORDER BY idx", (tid,)):
                        tr.addExon(erow[0], erow[1])
                    g.addTranscript(tr)
                return g
            else:
                return None

    def allTranscriptNames(self):
        names = []
        with self:
            curr = self.dbconn.execute("SELECT ID FROM Transcripts ORDER BY ID")
            names = [ r[0] for r in curr.fetchall() ]
        return names

    def findTranscript(self, name, chrom=None):
        """Returns the transcript called `name'. The name is matched against the ID, name, accession, and enst fields."""
        with self:
            trow = self.dbconn.execute("SELECT ID, name, accession, enst, chrom, strand, txstart, txend, cdsstart, cdsend FROM Transcripts WHERE ID=? OR name=? OR accession=? OR enst=?",
                                       (name, name, name, name)).fetchone()
            if trow:
                tid = trow[0]
                tr = Transcript(tid, trow[4], trow[5], trow[6], trow[7])
                tr.exons = []
                for pair in zip(['ID', 'name', 'accession', 'enst', 'chrom', 'strand', 'txstart', 'txend', 'cdsstart', 'cdsend'], trow):
                    setattr(tr, pair[0], pair[1])
                for erow in self.dbconn.execute("SELECT start, end FROM Exons WHERE ID=? ORDER BY idx", (tid,)):
                    tr.addExon(erow[0], erow[1])
                return tr
            else:
                return None

    def getAllTranscripts(self):
        """Returns an iterator that lopps over all transcripts."""
        with self:
            tcur = self.dbconn.cursor()
            ecur = self.dbconn.cursor()

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

    def findGenes(self, query, args=[]):
        """Returns the list of all genes that satisfy the `query'. `query' should be a SQL statement that returns
a single column of values from the ID field."""
        result = []
        with self:
            gcur = self.dbconn.cursor()
            for geneIDrow in self.dbconn.execute(query, args):
                geneID = geneIDrow[0]
                row = gcur.execute("SELECT ID, name, geneid, ensg, biotype, chrom, strand, start, end FROM Genes WHERE ID=?", (geneID,)).fetchone()
                if row:
                    g = Gene(geneID, row[5], row[6])
                    for pair in zip(['ID', 'name', 'geneid', 'ensg', 'biotype', 'chrom', 'strand', 'start', 'end'], row):
                        setattr(g, pair[0], pair[1])
                    for trow in self.dbconn.execute("SELECT ID, name, accession, enst, chrom, strand, txstart, txend, cdsstart, cdsend, canonical FROM Transcripts WHERE parentID=?", (geneID,)):
                        tid = trow[0]
                        tr = Transcript(tid, trow[4], trow[5], trow[6], trow[7])
                        tr.exons = []
                        for pair in zip(['ID', 'name', 'accession', 'enst', 'chrom', 'strand', 'txstart', 'txend', 'cdsstart', 'cdsend'], trow):
                            setattr(tr, pair[0], pair[1])
                        tr.canonical = (trow[10] == 'Y')
                        for erow in self.dbconn.execute("SELECT start, end FROM Exons WHERE ID=? ORDER BY idx", (tid,)):
                            tr.addExon(erow[0], erow[1])
                        g.addTranscript(tr)
                    result.append(g)
            return result

    def allIntersecting(self, chrom, start, end):
        """Returns all genes in `chrom' that intersect the `start-end' region."""
        return self.findGenes("SELECT ID from Genes where chrom=? and ((? <= start) and (start <= ?) or ((? <= end) and (end <= ?)) or ((start <= ?) and (end >= ?)))",
                              (chrom, start, end, start, end, start, end))

    def findClosestGene(self, chrom, start, end, transcripts=True, biotype=None, tss=False, canonical=False):
        """Find the closest gene to the region chrom:start-end, in either direction. Returns a tuple: (gene, distance).
Distance can be positive (downstream of `end') or negative (upstream of `start'). If `transcripts' is True, look at
transcript instead of genes; in this case the return value is (transcript ID, distance). If `biotype' is specified,
restrict the search to genes with that biotype (e.g., protein_coding). If `tss' is True, distances are computed relative
to the gene's (or transcript's) TSS.
"""
        g1 = None
        g2 = None
        p1 = 0
        p2 = 0
        d1 = 0
        d2 = 0
        args = {'table': "Transcripts" if transcripts else "Genes",
                'chrom': chrom,
                'start': start,
                'end': end,
                'fstart': "start",
                'fend': "end",
                'biotype': "AND biotype='{}' ".format(biotype) if biotype else "",
                'canon': "AND canonical='Y' " if (transcripts and canonical) else ""}
        if transcripts:
            args['fstart'] = 'txstart'
            args['fend'] = 'txend'
        elif tss:
            args['fstart'] = 'tss'
            args['fend'] = 'tss'

        with self:
            query0 = "SELECT ID, {fstart}, {fend} FROM {table} WHERE chrom='{chrom}' AND {fstart} <= {start} AND {fend} >= {end} {biotype} {canon} ORDER BY {fstart} DESC LIMIT 1;".format(**args) # containing
            query1 = "SELECT ID, {fend}   FROM {table} WHERE chrom='{chrom}' AND {fend} <= {start} {biotype} {canon} ORDER BY {fend} DESC LIMIT 1;".format(**args) # gene is upstream of region
            query2 = "SELECT ID, {fstart} FROM {table} WHERE chrom='{chrom}' AND {fstart} >= {end} {biotype} {canon} ORDER BY {fstart} LIMIT 1;".format(**args) # gene is downstream of region
            r0 = self.dbconn.execute(query0).fetchone()
            if r0:
                g0 = r0[0]
                p1 = r0[1]
                p2 = r0[2]
                d1 = start - p1
                d2 = p2 - end
                if d1 < d2:
                    return (g0, d1)
                else:
                    return (g0, -d2)
            r1 = self.dbconn.execute(query1).fetchone()
            if r1:
                g1 = r1[0]
                p1 = r1[1]
                d1 = start - p1
                # print (r1, d1)
            r2 = self.dbconn.execute(query2).fetchone()
            if r2:
                g2 = r2[0]
                p2 = r2[1]
                d2 = p2 - end
                # print (r2, d2)
            if p1 == 0 and p2 == 0: # No results?? This should not happen...
                return (None, 0)
            elif p1 == 0:           # No upstream, we only have downstream
                return (g2, -d2)
            elif p2 == 0:           # No downstream, we only have upstream
                return (g1, d1)
            else:
                #print (d1, d2)
                if d1 < d2:         # We have both, but upstream is closer
                    return (g1, d1)
                else:
                    return (g2, -d2)

    def getGeneInfo(self, geneid, query):
        with self:
            gcur = self.dbconn.cursor()
            row = gcur.execute(query, (geneid,)).fetchone()
            return row

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
    canonical = False
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
            print("{}{} {}:{}-{} {}".format(prefix, self.ID, self.chrom, self.txstart, self.txend, self.exons))
        else:
            print(prefix + "ID: " + self.ID)
            print(prefix + "Chrom: " + self.chrom)
            print("{}Strand: {}".format(prefix, self.strand))
            print("{}Transcript: {}-{}".format(prefix, self.txstart, self.txend))
            print("{}CDS: {}-{}".format(prefix, self.cdsstart, self.cdsend))
            print("{}Exons: {}".format(prefix, self.exons))

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

    def updateCDS(self, cdsstart, cdsend):
        if self.cdsstart:
            self.cdsstart = min(self.cdsstart, cdsstart)
        else:
            self.cdsstart = cdsstart
        if self.cdsend:
            self.cdsend = max(self.cdsend, cdsend)
        else:
            self.cdsend = cdsend

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
        """If position `pos' is in one of the exons of this transcript, return the exon number, otherwise return False."""
        nexons = len(self.exons) + 1
        ne = 1
        for e in self.exons:
            if e[0] <= pos <= e[1]:
                if self.strand == 1:
                    return ne
                else:
                    return nexons - ne
            ne += 1
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
        """Returns a Classif instance that classifies position `pos' within this transcript,
including `updistance' bases upstream of the TSS and `dndistance' bases downstream of the TES.
If the class is CL_CODING or CL_EXON, the `extra' attribute contains the exon number.
"""
        if self.txstart <= pos <= self.txend:
            ne = self.posInExon(pos)
            if ne:
                if self.cdsstart <= pos <= self.cdsend:
                    return CL_CODING(subject=self, extra=ne)
                else:
                    return CL_EXON(subject=self, extra=ne)
            return CL_INTRON(subject=self)
        elif self.strand == 1:
            if self.txstart - updistance <= pos < self.txstart:
                return CL_UPSTREAM(subject=self)
            elif self.txend < pos <= self.txend + dndistance:
                return CL_DNSTREAM(subject=self)
            else:
                return None
        else:
            if self.txstart - dndistance <= pos < self.txstart:
                return CL_DNSTREAM(subject=self)
            elif self.txend < pos <= self.txend + updistance:
                return CL_UPSTREAM(subject=self)
            else:
                return None
            
    def getRegion(self, params):
        """Returns the region for this transcript as (start, end) according to the 'regwanted', 
'updistance, and 'dndistance' attributes in the `params' object."""
        if self.txstart and self.txend:
            if params.regwanted == 'b':
                return (self.txstart, self.txend)
            elif params.regwanted == 'u':
                if self.strand == 1:
                    return (max(0, self.txstart - params.updistance), self.txstart + params.dndistance)
                else:
                    return (max(0, self.txend - params.dndistance), self.txend + params.updistance)
            elif params.regwanted == 'd':
                if self.strand == 1:
                    return (max(0, self.txend - params.updistance), self.txend + params.dndistance)
                else:
                    return (max(0, self.txstart - params.dndistance), self.txstart + params.updistance)
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
    description = ""
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
        print("ID: " + self.ID)
        print("Gene: " + self.name)
        print("GeneID: " + self.geneid)
        print("ENSG: " + self.ensg)
        print("Biotype: " + self.biotype)
        print("Region: {}:{}-{}".format(self.chrom, self.start, self.end))
        print("Strand: {}".format(self.strand))
        print("Transcripts:")
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

    def classifyPosition(self, position, updistance, dndistance, best=False):
        """Returns all possible classifications of `position' for the transcript of this gene.
If `best' is True, only returns a single classification, ie the first one to appear in this list: """
        result = []
        for tr in self.transcripts:
            c = tr.classifyPosition(position, updistance, dndistance)
            if c:
                c.subject = self
                addClassification(c, result)
        if len(result) == 0:
            return [CL_NONE()]
        else:
            sortClassifications(result)
        if best:
            return [result[0]]
        else:
            return result

    def getRegion(self, params):
        """Returns the region for this gene as (start, end) according to the 'regwanted', 
'updistance, and 'dndistance' attributes in the `params' object."""
        #print(self.name, self.start, self.end)
        if self.start and self.end:
            if params.regwanted == 'b':
                return (self.start, self.end)
            elif params.regwanted == 'B':
                if self.strand == 1:
                    return (max(0, self.start - params.updistance), self.end + params.dndistance)
                else:
                    return (max(0, self.start - params.dndistance), self.end + params.updistance)
            elif params.regwanted == 'p':
                if self.strand == 1:
                    return (max(0, self.start - params.updistance), self.start)
                else:
                    return (self.end, self.end + params.dndistance)
            elif params.regwanted == 'u':
                if self.strand == 1:
                    return (max(0, self.start - params.updistance), self.start + params.dndistance)
                else:
                    return (max(0, self.end - params.dndistance), self.end + params.updistance)
            elif params.regwanted == 'd':
                if self.strand == 1:
                    return (max(0, self.end - params.updistance), self.end + params.dndistance)
                else:
                    return (max(0, self.start - params.dndistance), self.start + params.updistance)
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
        if chrom.startswith("scaff"):
            return chrom
        if chrom.find("_") > 0:
            return False
        if chrom.find("random") > 0:
            return False
        if chrom.startswith('chr') or chrom.startswith('Chr'):
            return chrom
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
                            # print("In gene section: {}".format(data))
                            if data[0:5] == '/gene':
                                # print("Setting name to {}".format(data[7:-1]))
                                thisGene.name = thisGene.ID = data[7:-1]
                            elif data[0:10] == '/locus_tag':
                                if thisGene.name == '':
                                    thisGene.name = data[12:-1]
                                if thisGene.ID == '':
                                    thisGene.ID = data[12:-1]
                    elif key == 'gene':
                        # print("Found new gene")
                        # raw_input()
                        section = key
                        start, end, strand, ignore = parseStartEnd(line[21:])
                        # print("New gene at ({}, {})".format(start, end))
                        thisGene = Gene("", chrom, strand)
                        thisGene.chrom = chrom
                        self.gl.add(thisGene, chrom)
                        thisGene.addTranscript(Transcript("", chrom, strand, start, end))
                    elif key == 'CDS':
                        cdsstart, cdsend, ignore, introns = parseStartEnd(line[21:])
                        # print("Setting CDS to ({}, {})".format(cdsstart, cdsend))
                        if introns:
                            # print("Setting introns: {}".format(introns))
                            thisGene.transcripts[0].setIntrons(introns)

                        thisGene.transcripts[0].setCDS(cdsstart, cdsend)
                        # print(thisGene.smallrects)
                        # print(thisGene.largerects)
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

        if wanted != []:
            self.gl.setWanted(wanted)
        elif notwanted != []:
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
                    # print("Adding {} to {}".format(self.currGene.ID, chrom))
                    # raw_input()
                    self.gl.add(self.currGene, chrom)

                elif btype == 'transcript':
                    txid = ann['transcript_id']
                    if self.currTranscript:
                        self.currTranscript.cdsstart=int(line[3])
                        self.currTranscript.cdsend=int(line[4])
                        self.currTranscript.setCDS(self.currTranscript.cdsstart, self.currTranscript.cdsend) # Seems redundant, but we can only do this after all exons have been collected
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
                    #self.currTranscript.setCDS(start, end)
                    self.currTranscript.updateCDS(start, end)
                elif btype == 'exon':
                    start = int(line[3])
                    end   = int(line[4])
                    self.currTranscript.addExon(start, end)
        self.currTranscript.setCDS(self.currTranscript.cdsstart, self.currTranscript.cdsend) # Seems redundant, but we can only do this after all exons have been collected

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

        if wanted != []:
            self.gl.setWanted(wanted)
        elif notwanted != []:
            self.gl.setNotWanted(notwanted)

        with open(self.filename, "r") as f:
            reader = Utils.CSVreader(f)
            for line in reader:
#                if len(line) < 8:
#                    print("|"+line+"|")
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
                    self.currGene.start = int(line[3])
                    self.currGene.end   = int(line[4])
                    self.currGene.name = Utils.dget('Name', ann, "")
                    self.currGene.biotype = Utils.dget('biotype', ann)
                    self.gl.add(self.currGene, chrom)
                    
                elif tag in ['mRNA', 'transcript', 'processed_transcript', 'pseudogenic_transcript', 'pseudogene', 'processed_pseudogene', 'miRNA', 'lincRNA']:
                    ann = self.parseAnnotations(line[8])
                    tid = self.cleanID(Utils.dget('ID', ann))
                    pid = self.cleanID(Utils.dget('Parent', ann))
                    if self.currTranscript:
                        self.currTranscript.setCDS(self.currTranscript.cdsstart, self.currTranscript.cdsend) # Seems redundant, but we can only do this after all exons have been collected
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
                        #self.currTranscript.setCDS(start, end)
                        self.currTranscript.updateCDS(start, end)
                    else:
                        orphans += 1
        self.currTranscript.setCDS(self.currTranscript.cdsstart, self.currTranscript.cdsend) # Seems redundant, but we can only do this after all exons have been collected
        sys.stderr.write("Orphans: {}\n".format(orphans))

class DBloader(GeneLoader):
    conn = None

    def _load(self, preload=True, wanted=[], notwanted=[]):
        self.gl = GenelistDB()
        self.gl.dbname = self.filename
        self.gl.preloaded = preload
        self.conn = sql.connect(self.filename)
        try:
            if preload:
                gcur = self.conn.cursor()
                tcur = self.conn.cursor()
                ecur = self.conn.cursor()
                for row in gcur.execute("SELECT ID, name, geneid, ensg, biotype, chrom, strand, start, end, description FROM Genes"):
                    gid = row[0]
                    g = Gene(gid, row[5], row[6])
                    for pair in zip(['ID', 'name', 'geneid', 'ensg', 'biotype', 'chrom', 'strand', 'start', 'end', 'description'], row):
                        setattr(g, pair[0], pair[1])
                    self.gl.add(g, g.chrom)
                    for trow in tcur.execute("SELECT ID, name, accession, enst, chrom, strand, txstart, txend, cdsstart, cdsend, canonical FROM Transcripts WHERE parentID=?", (gid,)):
                        tid = trow[0]
                        tr = Transcript(tid, trow[4], trow[5], trow[6], trow[7])
                        for pair in zip(['ID', 'name', 'accession', 'enst', 'chrom', 'strand', 'txstart', 'txend', 'cdsstart', 'cdsend'], trow):
                            setattr(tr, pair[0], pair[1])
                        tr.canonical = (trow[10] == 'Y')
                        for erow in ecur.execute("SELECT start, end FROM Exons WHERE ID=? ORDER BY idx", (tid,)):
                            tr.addExon(erow[0], erow[1])
                        g.addTranscript(tr)
            else:
                try:
                    ncur = self.conn.execute("SELECT ngenes FROM Counts")
                except:
                    ncur = self.conn.execute("SELECT count(*) FROM Genes") # Fallback method for databases that don't have the Counts table yet...
                self.gl.ngenes = ncur.fetchone()[0]
        finally:
            self.conn.close()

### Database stuff

def initializeDB(filename):
    """Create a new database in 'filename' and write the Genes, Transcripts, and Exons tables to it."""
    conn = sql.connect(filename)
    try:
        conn.execute("DROP TABLE IF EXISTS Genes;")
        conn.execute("CREATE TABLE Genes (ID varchar primary key, name varchar, geneid varchar, ensg varchar, biotype varchar, description text, chrom varchar, strand int, start int, end int, canonical varchar, tss int default null);")
        for field in ['name', 'geneid', 'ensg', 'chrom', 'start', 'end']:
            conn.execute("CREATE INDEX Genes_{} on Genes({});".format(field, field))
        conn.execute("DROP TABLE IF EXISTS Transcripts;")
        conn.execute("CREATE TABLE Transcripts (ID varchar primary key, parentID varchar, name varchar, accession varchar, enst varchar, chrom varchar, strand int, txstart int, txend int, cdsstart int, cdsend int, canonical char(1) default 'N', tss int default null);")
        for field in ['parentID', 'name', 'accession', 'enst', 'chrom', 'txstart', 'txend']:
            conn.execute("CREATE INDEX Trans_{} on Transcripts({});".format(field, field))
        conn.execute("DROP TABLE IF EXISTS Exons;")
        conn.execute("CREATE TABLE Exons (ID varchar, idx int, chrom varchar, start int, end int);")
        for field in ['ID', 'chrom', 'start', 'end']:
            conn.execute("CREATE INDEX Exon_{} on Exons({});".format(field, field))
        conn.execute("DROP TABLE IF EXISTS Source;")
        conn.execute("CREATE TABLE Source (filename VARCHAR, ts TIMESTAMP DEFAULT CURRENT_TIMESTAMP NOT NULL);")
        conn.execute("DROP TABLE IF EXISTS Counts;")
        conn.execute("CREATE TABLE Counts (ngenes INT DEFAULT 0, ntranscripts INT DEFAULT 0, nexons INT DEFAULT 0);")
    finally:
        conn.close()

### Top level

def loadGenes(filename, preload=True):
    loaderClass = getLoader(filename)
    if not loaderClass:
        return None
    loader = loaderClass(filename)
    sys.stderr.write("Loading genes database from {}... ".format(filename))
    gl = loader.load(preload=preload)
    sys.stderr.write("{} genes loaded.\n".format(gl.ngenes))
    return gl

