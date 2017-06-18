#!/usr/bin/env python

import sys
import csv
import ast
import math
import gzip
from Utils import genOpen, convertValue

### filter C1=I set f=C1+C2 return avg(f)
### set f=C1+C2 print C1 C2 f

### Exceptions

class SkipEntry(Exception):
    pass

### Functions

ACTIVE = None

def FUNC(value):
    ACTIVE.perform(value)

SUM = FUNC
ADD = FUNC
AVG = FUNC
MEAN = FUNC
STDEV = FUNC
MIN = FUNC
MAX = FUNC

### Terms

class Term():
    termtype = "term"
    parent = None
    source = ""
    code = None

    def __init__(self, source):
        self.init(source)

    def init(self, source):
        self.source = source
        a = ast.parse(source, mode='eval')
        self.code = compile(a, '<ast>', mode='eval')

    def execute(self):
        return eval(self.code, globals(), self.parent.bindings)

    def dump(self):
        return self.source

class FilterTerm(Term):
    termtype = "filter"

    def dump(self):
        return "IF " + self.source

    def execute(self):
        f = eval(self.code, globals(), self.parent.bindings)
        if not f:
            raise SkipEntry

class SetTerm(Term):
    termtype = "set"
    varname = ""

    def __init__(self, source):
        q = source.find("=")
        if q > 0:
            self.varname = source[:q]
            self.init(source[q+1:])
        else:
            sys.stderr.write("An assignment term should have the form V=expr.")
            return None

    def execute(self):
        value = eval(self.code, globals(), self.parent.bindings)
        self.parent.bindings[self.varname] = value
        return value

    def dump(self):
        return "{} = {}".format(self.varname, self.source)

class PrintTerm(Term):
    termtype = "print"
    variables = []

    def init(self, source):
        if source:
            self.addVariables(source)
    
    def addVariables(self, vars):
        for v in [ s.strip(" \t") for s in vars.split(",") ]:
            self.variables.append(v)

    def execute(self):
        if self.variables == []:
            self.variables = self.parent.colnames
        out = self.parent.out
        bdg = self.parent.bindings
        bdg['CM'] += 1
        out.write(str(bdg[self.variables[0]]))
        for v in self.variables[1:]:
            out.write("\t")
            out.write(str(bdg[v]))
        out.write("\n")

    def dump(self):
        return "PRINT " + ",".join(self.variables)

class ReturnTerm(Term):
    termtype = "return"

    def execute(self):
        global ACTIVE
        ACTIVE = self
        eval(self.code, globals(), self.parent.bindings)

    def dump(self):
        return "DO " + self.source

class SumTerm(ReturnTerm):
    sum = 0

    def reset(self):
        self.sum = 0

    def execute(self):
        global ACTIVE
        ACTIVE = self
        eval(self.code, globals(), self.parent.bindings)

    def perform(self, value):
        #sys.stderr.write("Adding {}\n".format(value))
        if not type(value).__name__ == 'str':
            self.sum += value

    def result(self):
        return self.sum

class AvgTerm(ReturnTerm):
    sum = 0
    nvals = 0

    def reset(self):
        self.sum = 0
        self.nvals = 0

    def execute(self):
        global ACTIVE
        ACTIVE = self
        eval(self.code, globals(), self.parent.bindings)

    def perform(self, value):
        if not type(value).__name__ == 'str':
            self.sum += value
            self.nvals += 1

    def result(self):
        return 1.0*self.sum/self.nvals
        
class StdevTerm(ReturnTerm):
    sum = 0
    sumsq = 0
    nvals = 0

    def reset(self):
        self.sum = 0
        self.sumsq = 0
        self.nvals = 0

    def execute(self):
        global ACTIVE
        ACTIVE = self
        eval(self.code, globals(), self.parent.bindings)

    def perform(self, value):
        if not type(value).__name__ == 'str':
            self.sum += value
            self.sumsq += value * value
            self.nvals += 1

    def result(self):
        e = 1.0 * self.sum / self.nvals
        return math.sqrt((1.0 * self.sumsq / self.nvals) - (e * e))
    
class MinTerm(ReturnTerm):
    minvalue = sys.float_info.max

    def reset(self):
        self.minvalue = sys.float_info.max

    def perform(self, value):
        if not type(value).__name__ == 'str':
            self.minvalue = min(self.minvalue, value)
    
    def result(self):
        return self.minvalue

class MaxTerm(ReturnTerm):
    maxvalue = sys.float_info.min

    def reset(self):
        self.minvalue = sys.float_info.min

    def perform(self, value):
        if not type(value).__name__ == 'str':
            self.maxvalue = max(self.maxvalue, value)
    
    def result(self):
        return self.maxvalue

TERMMAP = [ ('SUM', SumTerm),
            ('ADD', SumTerm),
            ('AVG', AvgTerm),
            ('MEAN', AvgTerm),
            ('STDEV', StdevTerm),
            ('MIN', MinTerm),
            ('MAX', MaxTerm) ]

def funcToTerm(source):
    global TERMMAP
    p = source.find("(")
    if p > 0:
        tag = source[:p].upper()
        for tm in TERMMAP:
            if tag == tm[0]:
                return (tm[1], tag + source[p:])
    return (None, source)

### Driver object

class Driver():
    src = sys.stdin
    out = sys.stdout
    infile = None
    outfile = None
    ncols = None
    colnames = None
    terms = []
    bindings = {}
    variables = {}

    # Default terms
    printTerm = None
    returnTerms = []

    # Debugging
    dumpRecipe = False

    def __init__(self):
        self.terms = []
        self.bindings = {'CN': 0, 'CM':0}
        self.returnTerms = []

    def parseArgs(self, args):
        na = 0
        next = ""
        for a in args:
            if next == "-n":
                self.setNcols(int(a))
                next = ""
            elif a in ["-n"]:
                next = a
            elif a == '-d':
                self.dumpRecipe = True
            elif na == 0:
                if a[0] == '@':
                    with open(a[1:], 'r') as f:
                        rec = f.read()
                else:
                    rec = a
                self.parseRecipe(rec)
                na += 1
            elif na == 1:
                self.infile = a
                na += 1
            elif na == 2:
                self.outfile = a
                na += 1

    def initialize(self):
        if self.infile:
            self.src = genOpen(self.infile, "r")
        if self.outfile:
            self.out = open(self.outfile, "w")

    def cleanup(self):
        if self.infile:
            self.src.close()
        if self.outfile:
            self.out.close()

    def addTerm(self, term):
        term.parent = self
        self.terms.append(term)
        if term.termtype == "set":
            self.variables[term.varname] = term

    def parseRecipe(self, recipe):
        words = recipe.split()
        mode = ""
        for w in words:
            if w in ["filter", "if", "set", "print", "do", "return"]:
                mode = w
            elif mode == "filter" or mode == "if":
                # sys.stderr.write("Adding FILTER term `{}'\n".format(w))
                self.addTerm(FilterTerm(w))
            elif mode == "set":
                # sys.stderr.write("Adding SET term `{}'\n".format(w))
                self.addTerm(SetTerm(w))
            elif mode == "print":
                if self.printTerm == None:
                    sys.stderr.write("Adding PRINT term `{}'\n".format(w))
                    pt = PrintTerm(w)
                    self.addTerm(pt)
                    self.printTerm = pt
                else:
                    # sys.stderr.write("Adding variable to PRINT term `{}'\n".format(w))
                    pt.addVariables(w)
            elif mode == "do" or mode == "return":
                (cls, neww) = funcToTerm(w)
                # print (cls, neww)
                if cls:
                    rt = cls(neww)
                    self.addTerm(rt)
                    self.returnTerms.append(rt)
                else:
                    sys.stderr.write("Warning: no function found in term `{}'.\n".format(w))
        if not (self.printTerm or self.returnTerms):
            pt = PrintTerm(None)
            self.addTerm(pt)
            self.printTerm = pt

    def execute(self):
        for term in self.terms:
            term.execute()

    def setNcols(self, ncols):
        self.ncols = ncols
        self.colnames = [ "C" + str(i) for i in range(1, ncols+1) ]

    def bindColumnValues(self, row):
        """Assign elements of `row' to the C1... Cn variables after converting them
to their numeric representation if possible."""
        for colname, val in zip(self.colnames, row):
            self.bindings[colname] = convertValue(val)

    def processAllRows(self):
        f = csv.reader(self.src, delimiter='\t')
        try:
            for row in f:
                if len(row) == 0:
                    continue
                if len(row[0]) > 0 and row[0][0] == '#':
                    continue
                if not self.ncols:
                    self.setNcols(len(row))
                self.bindings['CN'] += 1
                try:
                    self.processRow(row)
                except SkipEntry:
                    pass
        except IOError:
            return

    def processRow(self, row):
        self.bindColumnValues(row)
        self.execute()

    def printReturns(self):
        if len(self.returnTerms) > 0:
            self.out.write("\t".join([rt.source for rt in self.returnTerms]) + "\n")
            results = [str(rt.result()) for rt in self.returnTerms]
            self.out.write("\t".join(results) + "\n")

    def dump(self):
        sys.stderr.write("*** Recipe dump:\n")
        for term in self.terms:
            sys.stderr.write("  " + term.dump() + "\n")
        sys.stderr.write("***\n")
            
## Test

def test():
    d = Driver()
    d.parseRecipe("set x=C1+C2 print C1,C2,x")
    d.dump()
    d.setNcols(3)
    d.processAllRows()
    print d.colnames
    print d.bindings
    print d.variables

if __name__ == "__main__":
    d = Driver()
    d.parseArgs(sys.argv[1:])
    if d.dumpRecipe:
        d.dump()
    d.initialize()
    try:
        d.processAllRows()
        d.printReturns()
    finally:
        d.cleanup()

