#!/usr/bin/env python

import sys
import csv
import ast
import gzip
import os.path

### Utils

def convertValue(v):
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except ValueError:
            return v

### filter C1=I set f=C1+C2 return avg(f)
### set f=C1+C2 print C1 C2 f


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
        self.variables = [ s.strip(" \t") for s in source.split(",") ]
        print self.variables

    def execute(self):
        out = self.parent.out
        bdg = self.parent.bindings
        out.write(str(bdg[self.variables[0]]))
        for v in self.variables[1:]:
            out.write("\t")
            out.write(str(bdg[v]))
        out.write("\n")

class ReturnTerm(Term):
    termtype = "return"

    def execute(self):
        value = eval(self.code, globals(), self.parent.bindings)
        print value

    def dump(self):
        return "DO " + self.source

### Driver object

class Driver():
    out = sys.stdout
    ncols = None
    colnames = None
    terms = []
    bindings = {}
    variables = {}

    def __init__(self):
        self.terms = []
        self.bindings = {}

    def addTerm(self, term):
        term.parent = self
        self.terms.append(term)
        if term.termtype == "set":
            self.variables[term.varname] = term

    def parseRecipe(self, recipe):
        words = recipe.split()
        mode = ""
        for w in words:
            if w in ["filter", "set", "print"]:
                mode = w
            elif mode == "filter":
                self.addTerm(FilterTerm(w))
            elif mode == "set":
                self.addTerm(SetTerm(w))
            elif mode == "print":
                self.addTerm(PrintTerm(w))
            
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

    def processRow(self, row):
        self.bindColumnValues(row)
        self.execute()

    def dump(self):
        for term in self.terms:
            print term.dump()
            
## Test

def test():
    d = Driver()
    d.parseRecipe("set x=C1+C2 print C1,C2,x")
    d.dump()
    d.setNcols(3)
    d.processRow(["5", "2", "pippo"])
    d.processRow(["6", "12", "pluto"])
    print d.colnames
    print d.bindings
    print d.variables
