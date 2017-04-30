#!/usr/bin/env python

import sys

import Utils
import Script

def usage(what=None):
    if what == 'rank':
        sys.stderr.write("""Usage: sorter.py rank [options...] infile [outfile]

Options:

  -o F | Write output to file F (default: stdout).
  -i I | Row identifiers are in column I (default: 1).
  -c C | Columns to use for score computation are specified by colspec C.
  -r   | Reverse output order: low-to-high instead of high-to-low.
  -s S | Use method S to compute row score. Possible values are: 'A' (average),
         'a' (average, missing values ignored), 'S' (sum), 's' (sum, missing
         values ignored), 'm' (min), 'M' (max). Default: {}.

A column specification C has the form C1,C2,...,Cn where each C can be either a 
number (1-based column number), X-Y (from column X to column Y inclusive), X+K 
(K columns starting at X).

""".format(Sorter.score))

    elif what == 'order':
        sys.stderr.write("""Usage: sorter.py order [options...] infile [outfile]

Options:

  -o F | Write output to file F (default: stdout).
  -i I | Row identifiers are in column I (default: 1).
  -k K | Read row order from file K (default: stdin). This file should be in the
         format produced by the `rank' command.

""")
    else:
        sys.stderr.write("""sorter.py - Sort tables based on values in one or more columns.

Usage sorter.py command [options... ] infile [outfile]

where `command' is one of: rank, order.

If command is `rank', computes a score for each row, and will output a tab-delimited file 
containing row identifiers in the first column and scores in the second one. The file 
will be sorted by decreasing values of the score.

If command is `order', sorts the input file according to the ordering in a rank file (produced
by the `rank' command).

Use the -h option followed by the name of the command for command-specific options.

""")

### Functions to compute scores

def score_ASas(cols, values, avg, default):
    n = 0
    s = 0
    for c in cols:
        v = Utils.safeFloat(values[c], default=default)
        if v != None:
            s += v
            n += 1

    if avg:
        return s/n
    else:
        return s

def score_A(cols, values):
    return score_ASas(cols, values, True, 0)

def score_a(cols, values):
    return score_ASas(cols, values, True, None)

def score_S(cols, values):
    return score_ASas(cols, values, False, 0)

def score_s(cols, values):
    return score_ASas(cols, values, False, None)

def score_m(cols, values):
    m = sys.float_info.max
    for c in cols:
        v = Utils.safeFloat(values[c], default=None)
        if v and v < m:
            m = v
    return m

def score_M(cols, values):
    m = -sys.float_info.max
    for c in cols:
        v = Utils.safeFloat(values[c], default=None)
        if v and v > m:
            m = v
    return m

SCOREFUNCS = {'A': score_A,
              'a': score_a,
              'S': score_S,
              's': score_s,
              'm': score_m,
              'M': score_M}

### Main class

class Sorter(Script.Script):
    mode = ""
    infile = None
    outfile = None
    rankfile = None
    orderfile = ""
    idcol = 0
    score = 'a'
    scorecols = None
    reverse = True

    def parseArgs(self, args):
        self.standardOpts(args)
        next = ''

        if len(args) == 0:
            self.errmsg(self.NOCMD)

        self.mode = args[0]
        for a in args[1:]:
            if next == '-o':
                self.outfile = a
                next = ""
            elif next == '-i':
                self.idcol = self.toInt(a) - 1
                next = ""
            elif next == '-c':
                self.scorecols = Utils.parseColspec(a)
                next = ""
            elif next == '-s':
                self.score = a
                next = ""
            elif next == '-k':
                self.rankfile = self.isFile(a)
                next = ""
            elif a in ['-o', '-i', '-c', '-s', '-k']:
                next = a
            elif a == '-r':
                self.reverse = False
            else:
                self.infile = self.isFile(a)

        if self.infile == None:
            self.errmsg(self.NOFILE)
        if self.mode == 'rank' and self.scorecols == None:
            self.errmsg(self.NOCOLS)
        return True

    def run(self):
        if self.mode == 'rank':
            self.doRank()
        elif self.mode == 'order':
            self.doOrder()

    def doRank(self):
        scorefunc = SCOREFUNCS[self.score]
        scores = []

        with open(self.infile, "r") as f:
            r = Utils.CSVreader(f)
            for line in r:
                gene = line[self.idcol]
                score = scorefunc(self.scorecols, line)
                scores.append([gene, score])

        scores.sort(key=lambda r: r[1], reverse=self.reverse)

        if self.outfile:
            out = open(self.outfile, "w")
        else:
            out = sys.stdout
        try:
            for sc in scores:
                out.write("{}\t{}\n".format(*sc))
        finally:
            out.close()

    def readRanking(self):
        ranking = []
        if self.rankfile:
            f = open(self.rankfile, "r")
        else:
            f = sys.stdin
        try:
            r = Utils.CSVreader(f)
            for line in r:
                ranking.append(line[0])
        finally:
            f.close()
        return ranking

    def readRows(self):
        rows = {}
        with open(self.infile, "r") as f:
            r = Utils.CSVreader(f)
            for line in r:
                rows[line[self.idcol]] = line
        return rows

    def doOrder(self):
        ranking = self.readRanking()
        rows = self.readRows()
        if self.outfile:
            out = open(self.outfile, "w")
        else:
            out = sys.stdout
        try:
            print ranking
            for g in ranking:
                out.write("\t".join(rows[g]) + "\n")
        finally:
            out.close()

P = Sorter("sorter.py", version="1.0", usage=usage,
           errors=[('NOCMD', 'Missing command', 'The first argument should be one of: rank, order.'),
                   ('NOCOLS', 'Missing column specification', 'Please specify one or more target columns with the -c option.')])

if __name__ == "__main__":
    if P.parseArgs(sys.argv[1:]):
        P.run()
