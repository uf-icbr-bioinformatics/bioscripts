#!/usr/bin/env python

import sys

import Utils
import Script

def usage(what=None):
    if what == 'rank':
        sys.stderr.write("""Usage: sorter.py rank [options...] infile [outfile]

Options:

  -o F | Write output to file F (default: stdout).
  -i I | Row identifiers are in column I (default: 0).
  -c C | Columns to use for score computation are specified by colspec C.
  -s S | Use method S to compute row score. Possible values are: 'a' (average),
         'm' (min), 'M' (max). Default: {}.

A column specification C has the form C1,C2,...,Cn where each C can be either a 
number (1-based column number), X-Y (from column X to column Y inclusive), X+K 
(K columns starting at X).

""".format(Sorter.score))


    elif what == 'order':
        pass
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

class Sorter(Script.Script):
    mode = ""
    infile = None
    outfile = None
    orderfile = ""
    idcol = 0
    score = 'a'
    scorecols = None

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
                self.idcol = self.toInt(a)
                next = ""
            elif next == '-c':
                self.scorecols = Utils.parseColspec(a)
                next = ""
            elif next == '-s':
                self.score = a
            elif a in ['-o', '-i', '-c', '-s']:
                next = a
            else:
                self.infile = self.isFile(a)
        if self.infile == None:
            self.errmsg(self.NOFILE)
        if self.scorecols == None:
            self.errmsg(self.NOCOLS)

P = Sorter("sorter.py", version="1.0", usage=usage,
           errors=[('NOCMD', 'Missing command', 'The first argument should be one of: rank, order.'),
                   ('NOCOLS', 'Missing column specification', 'Please specify one or more target columns with the -c option.')])

if __name__ == "__main__":
    if P.parseArgs(sys.argv[1:]):
        P.run()
