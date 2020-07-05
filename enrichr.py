#!/usr/bin/env python

import sys
import json
import requests

import Utils
import Script

class Enrichr(Script.Script):
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/'
    genes = []
    genesfile = None
    description = "Generic gene list"
    setid = None
    libraries = []
    outfile = "/dev/stdout"
    html = None                 # HTML output?
    barchart = None
    pval = 1
    qval = 1
    zscore = 0
    cscore = 0
    maxrank = None
    cssfile = None
    css = """BODY {font-family: arial}"""
    
    def parseArgs(self, args):
        self.standardOpts(args)
        self.command = args[0]
        next = ""

        if self.command == "upload":
            for a in args[1:]:
                if next == "-d":
                    self.description = a
                    next = ""
                elif next == "-u":
                    self.ENRICHR_URL = a
                elif a in ["-d", "-u"]:
                    next = a
                elif a.startswith("@"):
                    self.addGenesFromFile(a)
                else:
                    self.genes.append(a)
            if not self.genes:
                sys.stderr.write("No genes specified, reading from stdin...\n")
                self.addGenesFromStdin()

        elif self.command == "run":
            for a in args[1:]:
                if next == "-o":
                    self.outfile = a
                    next = ""
                elif next == "-n":
                    self.maxrank = self.toInt(a)
                    next = ""
                elif next == "-p":
                    self.pval = self.toFloat(a)
                    next = ""
                elif next == "-q":
                    self.qval = self.toFloat(a)
                    next = ""
                elif next == "-z":
                    self.zscore = self.toFloat(a)
                    next = ""
                elif next == "-c":
                    self.cscore = self.toFloat(a)
                    next = ""
                elif next == "-u":
                    self.ENRICHR_URL = a
                elif a in ["-o", "-n", "-p", "-u", "-q", "-z", "-c"]:
                    next = a
                elif a == "-H":
                    self.html = True
                elif a == "-b":
                    self.barchart = True
                elif self.setid:
                    self.libraries.append(a)
                else:
                    self.setid = a
            if not self.setid:
                self.errmsg(self.NOSETID)
            if not self.libraries:
                self.errmsg(self.NOLIBS)

    def addGenesFromFile(self, filename):
        ar = Utils.AtFileReader(filename)
        for g in ar:
            g = g.strip()
            if g:
                self.genes.append(g)

    def addGenesFromStdin(self):
        for g in sys.stdin:
            g = g.strip()
            if g:
                self.genes.append(g.strip())

    def uploadGeneList(self, genes, description):
        """Submit a list of genes to Enrichr, returns the gene set identifier."""
        genestr = "\n".join(genes)
        payload = { 'list': (None, genestr),
                    'description': (None, description) }
#        print payload
        response = requests.post(self.ENRICHR_URL + "addList", files=payload)
#        print response
        if not response.ok:
            self.errmsg(self.BADCOMM)

        data = json.loads(response.text)
        return data

    def getEnrichment(self, userListId, library):
        query_string = 'enrich?userListId={}&backgroundType={}'
        req = self.ENRICHR_URL + query_string.format(userListId, library)
#        print req
        response = requests.get(req)
#        print response.text
        if not response.ok:
            self.errmsg(self.BADCOMM)

        data = json.loads(response.text)
        return data

    # HTML utils

    def writeHead(self, out):
        if self.cssfile:
            with open(self.cssfile, "r") as f:
                self.css = f.read()
        out.write("""<!DOCTYPE html>
<HTML>
  <HEAD>
    <TITLE>Enrichr Output</TITLE>
    <STYLE>
{}
    </STYLE>
  </HEAD>
  <BODY>""".format(self.css))

    def writeTail(self, out):
        out.write("""  </BODY>
</HTML>
""")
        
    # Top-level commands

    def doUpload(self):
        sys.stderr.write("Uploading set of {} genes...\n".format(len(self.genes)))
        result = self.uploadGeneList(self.genes, self.description)
        sys.stdout.write("{}\n".format(result['userListId']))

    def doRun(self):
        with open(self.outfile, "w") as out:
            if self.barchart:
                self.writeBarChart(out)
            else:
                self.writeEnrichment(out)

    def writeEnrichment(self, out):
        headers = ["Library", "Rank", "Term", "P-value", "Q-value", "Z-score", "Combined score", "Overlapping genes"]
        if self.html:
            self.writeHead(out)
            out.write("<TABLE class='enr-table'>\n  <TR class='enr-header-row'>\n")
            for h in headers:
                out.write("    <TH class='enr-header'>" + h + "</TH>\n")
            out.write("  </TR>\n")
        else:
            out.write("#" + "\t".join(headers) + "\n")

        for lib in self.libraries:
            result = self.getEnrichment(self.setid, lib)
            for row in result[lib]:
                if self.maxrank and row[0] > self.maxrank:
                    break       # More than wanted number of results
                if (row[2] > self.pval or row[6] > self.qval):
                    continue        # P-value too high
                if (abs(row[3]) < self.zscore or row[4] < self.cscore):
                    continue
                if self.html:
                    out.write("  <TR class='enr-row'>")
                    out.write("<TD class='enr-text>{}</TD><TD class='enr-num'>{}</TD><TD class='enr-text'>{}</TD><TD class='enr-num'>{}</TD><TD class='enr-num'>{}</TD><TD class='enr-num'>{}</TD><TD class='enr-num'>{}</TD><TD class='enr-text'>{}</TD>\n".format(lib, row[0], row[1], row[2], row[6], row[3], row[4], ",".join(row[5])))
                    out.write("  </TR>\n")
                else:
                    out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(lib, row[0], row[1], row[2], row[6], row[3], row[4], ",".join(row[5])))
        if self.html:
            out.write("</TABLE>\n")
            self.writeTail(out)

    def writeBarChart(self, out):
        self.writeHead(out)
        for lib in self.libraries:
            result = self.getEnrichment(self.setid, lib)
            self.writeBarTable(out, lib, result[lib])
        self.writeTail(out)

    def writeBarTable(self, out, lib, rows):
        maxscore = 0
        goodrows = []
        for row in rows:
            if self.maxrank and row[0] > self.maxrank:
                break       # More than wanted number of results
            if (row[2] > self.pval or row[6] > self.qval):
                continue        # P-value too high
            if (abs(row[3]) < self.zscore or row[4] < self.cscore):
                continue
            goodrows.append(row)
            maxscore = max(maxscore, row[4])
            
        out.write("""<CENTER><H1>{}</H1>\n""".format(lib))
        out.write("""<TABLE width='80%' style='border: 2px solid black; margin: 10px; padding: 10px'>\n""")
        for row in goodrows:
            perc = int(100.0 * row[4] / maxscore)
            out.write("""  <TR>
    <TD style='border: 5px solid white; padding: 10px; background: linear-gradient(to right, #FF0000 {}%, #FFFFFF 0%);' nowrap>{}</TD>
    <TD align='right'>{:.1f}</TD>
  </TR>
""".format(perc, row[1], row[4]))
        out.write("""</CENTER>\n""")
            
    def usage(self, what=None):
        if what == "upload":
            sys.stderr.write("""enrichr.py - Command-line interface to the Enrichr website

Usage: enrichr.py [options] upload [genes...]

This command uploads a list of genes to the Enrichr server, and prints a set identifier (a
number) to standard output. The set identifier can then used in the `run' command.

Genes can be specifed on the command line as separate arguments. If an argument has the 
form @filename, gene names are read from the first column file `filename', one per line 
(a different column can be specified using @filename:col). If no arguments appear after
the command, gene names are read from standard input. For example, to upload the first 
50 genes in file genes.csv, you can do:

  head -50 genes.csv | enrichr.py upload

Options:

  -d D | Associate description D with the set of genes.

""")
        elif what == "run":
            sys.stderr.write("""enrichr.py - Command-line interface to the Enrichr website

Usage: enrichr.py [options] run setid library [libraries...]

This command performs enrichment analysis on the set of genes identified by `setid' (returned
by the `upload' command) against one or more Enrichr libraries. The full list of Enrichr 
libraries can be found at: http://amp.pharm.mssm.edu/Enrichr/#stats.

Output (written to standard output, or to a file if -o is specified is tab-delimited with the
following columns: library, rank, term name, P-value, Z score, combined score, overlapping genes.
See the Enrichr documentation for a description of these values. If multiple libraries are 
specified, the respective results will be concatenated in the same output; use the library id
in the first column to distinguish the different groups.

Options:

  -o O | Write output to file O (default: stdout)
  -H   | Write output as HTML instead of tab-delimited
  -p P | Only return terms with P-value less than or equal to P (default: {})
  -q Q | Only return terms with q-score (BH-corrected P-value) less than or equal to Q (default: {})
  -z Z | Only return terms with Z-score greater than Z (default: {})
  -c C | Only return terms with combined score greater than C (default: {})
  -n N | Only return the top N terms (default: no limit)

""".format(self.pval, self.qval, self.zscore, self.cscore))

        else:
            sys.stderr.write("""enrichr.py - Command-line interface to the Enrichr website

This program can be used to perform enrichment analysis on sets of genes, using the 
Enrichr engine (http://amp.pharm.mssm.edu/Enrichr/). Its basic usage is:

  enrichr.py [options] command [arguments...]

where `command' is one of: upload, run. Use "-h command" to display help for each command.

General options:

  -h   | Print help message
  -v   | Print version number
  -u U | Change the URL to access the Enrichr API to U
         (this should never be necessary).

""")

if __name__ == "__main__":
    E = Enrichr("enrichr.py", version="1.0", usage=Enrichr.usage,
                errors=[('NOSETID', 'Missing set id', "The first argument to the run command should be a set ID created with the upload command."),
                        ('NOLIBS', 'Missing libraries', "Please specify at least one Enrichr library name after the set ID."),
                        ('BADCOMM', 'Communication error', "Error communicating with the Enrichr server.")])

    E.parseArgs(sys.argv[1:])
    if E.command == "upload":
        E.doUpload()
    elif E.command == "run":
        E.doRun()
