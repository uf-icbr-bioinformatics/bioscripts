#!/usr/bin/env python

import sys
import json
import math
import requests

import Utils
import Script

LIBS = ["Genes_Associated_with_NIH_Grants",
        "Cancer_Cell_Line_Encyclopedia",
        "Achilles_fitness_decrease",
        "Achilles_fitness_increase",
        "Aging_Perturbations_from_GEO_down",
        "Aging_Perturbations_from_GEO_up",
        "Allen_Brain_Atlas_down",
        "Allen_Brain_Atlas_up",
        "ARCHS4_Cell-lines",
        "ARCHS4_IDG_Coexp",
        "ARCHS4_Kinases_Coexp",
        "ARCHS4_TFs_Coexp",
        "ARCHS4_Tissues",
        "BioCarta_2013",
        "BioCarta_2015",
        "BioCarta_2016",
        "BioPlanet_2019",
        "BioPlex_2017",
        "CCLE_Proteomics_2020",
        "ChEA_2013",
        "ChEA_2015",
        "ChEA_2016",
        "Chromosome_Location",
        "Chromosome_Location_hg19",
        "ClinVar_2019",
        "CORUM",
        "Data_Acquisition_Method_Most_Popular_Genes",
        "dbGaP",
        "DepMap_WG_CRISPR_Screens_Broad_CellLines_2019",
        "DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019",
        "Disease_Perturbations_from_GEO_down",
        "Disease_Perturbations_from_GEO_up",
        "Disease_Signatures_from_GEO_down_2014",
        "Disease_Signatures_from_GEO_up_2014",
        "DisGeNET",
        "Drug_Perturbations_from_GEO_2014",
        "Drug_Perturbations_from_GEO_down",
        "Drug_Perturbations_from_GEO_up",
        "DrugMatrix",
        "DSigDB",
        "Elsevier_Pathway_Collection",
        "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
        "ENCODE_Histone_Modifications_2013",
        "ENCODE_Histone_Modifications_2015",
        "ENCODE_TF_ChIP-seq_2014",
        "ENCODE_TF_ChIP-seq_2015",
        "Enrichr_Libraries_Most_Popular_Genes",
        "Enrichr_Submissions_TF-Gene_Coocurrence",
        "Epigenomics_Roadmap_HM_ChIP-seq",
        "ESCAPE",
        "Gene_Perturbations_from_GEO_down",
        "Gene_Perturbations_from_GEO_up",
        "GeneSigDB",
        "Genome_Browser_PWMs",
        "GO_Biological_Process_2013",
        "GO_Biological_Process_2015",
        "GO_Biological_Process_2017",
        "GO_Biological_Process_2017b",
        "GO_Biological_Process_2018",
        "GO_Cellular_Component_2013",
        "GO_Cellular_Component_2015",
        "GO_Cellular_Component_2017",
        "GO_Cellular_Component_2017b",
        "GO_Cellular_Component_2018",
        "GO_Molecular_Function_2013",
        "GO_Molecular_Function_2015",
        "GO_Molecular_Function_2017",
        "GO_Molecular_Function_2017b",
        "GO_Molecular_Function_2018",
        "GTEx_Tissue_Sample_Gene_Expression_Profiles_down",
        "GTEx_Tissue_Sample_Gene_Expression_Profiles_up",
        "GWAS_Catalog_2019",
        "HMDB_Metabolites",
        "HMS_LINCS_KinomeScan",
        "HomoloGene",
        "Human_Gene_Atlas",
        "Human_Phenotype_Ontology",
        "HumanCyc_2015",
        "HumanCyc_2016",
        "huMAP",
        "InterPro_Domains_2019",
        "Jensen_COMPARTMENTS",
        "Jensen_DISEASES",
        "Jensen_TISSUES",
        "KEA_2013",
        "KEA_2015",
        "KEGG_2013",
        "KEGG_2015",
        "KEGG_2016",
        "KEGG_2019_Human",
        "KEGG_2019_Mouse",
        "Kinase_Perturbations_from_GEO_down",
        "Kinase_Perturbations_from_GEO_up",
        "L1000_Kinase_and_GPCR_Perturbations_down",
        "L1000_Kinase_and_GPCR_Perturbations_up",
        "Ligand_Perturbations_from_GEO_down",
        "Ligand_Perturbations_from_GEO_up",
        "LINCS_L1000_Chem_Pert_down",
        "LINCS_L1000_Chem_Pert_up",
        "LINCS_L1000_Ligand_Perturbations_down",
        "LINCS_L1000_Ligand_Perturbations_up",
        "lncHUB_lncRNA_Co-Expression",
        "MCF7_Perturbations_from_GEO_down",
        "MCF7_Perturbations_from_GEO_up",
        "MGI_Mammalian_Phenotype_2013",
        "MGI_Mammalian_Phenotype_2017",
        "MGI_Mammalian_Phenotype_Level_3",
        "MGI_Mammalian_Phenotype_Level_4",
        "MGI_Mammalian_Phenotype_Level_4_2019",
        "Microbe_Perturbations_from_GEO_down",
        "Microbe_Perturbations_from_GEO_up",
        "miRTarBase_2017",
        "Mouse_Gene_Atlas",
        "MSigDB_Computational",
        "MSigDB_Oncogenic_Signatures",
        "NCI-60_Cancer_Cell_Lines",
        "NCI-Nature_2015",
        "NCI-Nature_2016",
        "NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions",
        "NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions",
        "NIH_Funded_PIs_2017_Human_AutoRIF",
        "NIH_Funded_PIs_2017_Human_GeneRIF",
        "NURSA_Human_Endogenous_Complexome",
        "Old_CMAP_down",
        "Old_CMAP_up",
        "OMIM_Disease",
        "OMIM_Expanded",
        "Panther_2015",
        "Panther_2016",
        "Pfam_Domains_2019",
        "Pfam_InterPro_Domains",
        "PheWeb_2019",
        "Phosphatase_Substrates_from_DEPOD",
        "PPI_Hub_Proteins",
        "ProteomicsDB_2020",
        "Rare_Diseases_AutoRIF_ARCHS4_Predictions",
        "Rare_Diseases_AutoRIF_Gene_Lists",
        "Rare_Diseases_GeneRIF_ARCHS4_Predictions",
        "Rare_Diseases_GeneRIF_Gene_Lists",
        "Reactome_2013",
        "Reactome_2015",
        "Reactome_2016",
        "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO",
        "SILAC_Phosphoproteomics",
        "SubCell_BarCode",
        "SysMyo_Muscle_Gene_Sets",
        "Table_Mining_of_CRISPR_Studies",
        "TargetScan_microRNA",
        "TargetScan_microRNA_2017",
        "TF-LOF_Expression_from_GEO",
        "TF_Perturbations_Followed_by_Expression",
        "Tissue_Protein_Expression_from_Human_Proteome_Map",
        "Tissue_Protein_Expression_from_ProteomicsDB",
        "Transcription_Factor_PPIs",
        "TRANSFAC_and_JASPAR_PWMs",
        "TRRUST_Transcription_Factors_2019",
        "UK_Biobank_GWAS_v1",
        "Virus-Host_PPI_P-HIPSTer_2020",
        "Virus_Perturbations_from_GEO_down",
        "Virus_Perturbations_from_GEO_up",
        "VirusMINT",
        "WikiPathways_2013",
        "WikiPathways_2015",
        "WikiPathways_2016",
        "WikiPathways_2019_Human",
        "WikiPathways_2019_Mouse"]

def sortResults(results, key):
    """Sort the list of results according to `key', one of: p, q, z, c."""
    for res in results.values():
        if key == 'p':
            res.sort(key=lambda e: e[2])
        elif key == 'q':
            res.sort(key=lambda e: e[6])
        elif key == 'z':
            res.sort(key=lambda e: e[3], reverse=True)
        elif key == 'c':
            res.sort(key=lambda e: e[4], reverse=True)

class EnEntry(object):
    rank = 0
    term = ""
    pvalue = None
    qvalue = None
    zscore = None
    cscore = None
    genes = []
    key = None

    def __init__(self, row):
        self.rank = row[0]
        self.term = row[1]
        self.pvalue = row[2]
        self.qvalue = row[6]
        self.zscore = row[3]
        self.cscore = row[4]
        self.genes = row[5]

    def setKey(self, key):
        if key == 'p':
            self.key = -math.log(self.pvalue, 10.0)
        elif key == 'q':
            self.key = -math.log(self.qvalue, 10.0)
        elif key == 'z':
            self.key = self.zscore
        elif key == 'c':
            self.key = self.cscore

class EnSet(object):
    lib = ""
    entries = []

    def __init__(self, lib):
        self.lib = lib
        self.entries = []
        
    def add(self, entry):
        self.entries.append(entry)
        
    def setKey(self, key):
        for e in self.entries:
            e.setKey(key)
        self.entries.sort(key=lambda e: e.key, reverse=True)
        
    def getMaxScore(self):
#        sys.stderr.write("{}\n".format( [e.key for e in self.entries] ))
        return max([e.key for e in self.entries])

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
    barcolor = "FF0000"
    barchart_var = "c"          # one of p, q, z, c
    barchart_value = None       # or one of p, q, z, c
    sortby = 'p'                # or one of p, q, z, c
    pval = 1
    qval = 1
    zscore = 0
    cscore = 0
    maxrank = None
    cssfile = None
    css = """BODY {font-family: arial; font-size: 32pt}
.tab {border: 2px solid black; margin: 10px; padding: 5px}
.bar {border: 5px solid white; padding: 5px}
"""
    
    def parseCmdline(self, args):
        self.standardOpts(args)
        self.command = args[0]
        args = args[1:]
        if self.command == "upload":
            self.parseArgs(args, "+d,+u")
            self.description = self.getOpt("d", self.description)
            self.ENRICHR_URL = self.getOpt("u", self.ENRICHR_URL)

            for a in self.getArgs():
                if a.startswith("@"):
                    self.addGenesFromFile(a)
                else:
                    self.genes.append(a)
            if not self.genes:
                sys.stderr.write("No genes specified, reading from stdin...\n")
                self.addGenesFromStdin()

        elif self.command == "run":
            self.parseArgs(args, "+o,+#n,+.p,+.q,+.z,+.c,+u,+s,H,b,+bc,+bv")
            self.outfile = self.getOpt("o", self.outfile)
            self.maxrank = self.getOpt("n")
            self.sortby  = self.getOpt("s", 'p')
            self.pval    = self.getOpt("p", self.pval)
            self.qval    = self.getOpt("q", self.qval)
            self.zscore  = self.getOpt("z", self.zscore)
            self.cscore  = self.getOpt("c", self.cscore)
            self.ENRICHR_URL = self.getOpt("u", self.ENRICHR_URL)
            self.html    = self.getOpt("H")
            self.barchart = self.getOpt("b")
            self.barcolor = self.getOpt("bc", self.barcolor)
            self.barchart_value = self.getOpt("bv")
            theArgs = self.getArgs()
            if len(theArgs) > 0:
                self.setid = self.getArgs()[0]
                self.libraries = self.getArgs()[1:]
            if not self.setid:
                self.errmsg(self.NOSETID)
            if not self.libraries:
                self.errmsg(self.NOLIBS)
            for lib in self.libraries:
                if lib not in LIBS:
                    self.errmsg(self.BADLIB, lib)

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
#        sys.stderr.write("{}\n".format(response.text))
        if not response.ok:
            self.errmsg(self.BADCOMM)

        data = json.loads(response.text)
        ES = EnSet(data)
        
        nr = 0
        for row in data[library]:
            if self.maxrank and nr > self.maxrank:
                break       # More than wanted number of results
            E = EnEntry(row)
            if (E.pvalue > self.pval or E.qvalue > self.qval):
                continue        # P-value too high
            if (abs(E.zscore) < self.zscore or E.cscore < self.cscore):
                continue
            ES.add(E)
            nr += 1

        ES.setKey(self.sortby)
        return ES

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
            for e in result.entries:
                if self.html:
                    out.write("  <TR class='enr-row'>")
                    out.write("<TD class='enr-text'>{}</TD><TD class='enr-num'>{}</TD><TD class='enr-text'>{}</TD><TD class='enr-num'>{}</TD><TD class='enr-num'>{}</TD><TD class='enr-num'>{}</TD><TD class='enr-num'>{}</TD><TD class='enr-text'>{}</TD>\n".format(lib, e.rank, e.term, e.pvalue, e.qvalue, e.zscore, e.cscore, ",".join(e.genes)))
                    out.write("  </TR>\n")
                else:
                    out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(lib, e.rank, e.term, e.pvalue, e.qvalue, e.zscore, e.cscore, ",".join(e.genes)))
        if self.html:
            out.write("</TABLE>\n")
            self.writeTail(out)

    def writeBarChart(self, out):
        self.writeHead(out)
        for lib in self.libraries:
            result = self.getEnrichment(self.setid, lib)
            self.writeBarTable(out, lib, result)
        self.writeTail(out)

    def writeBarTable(self, out, lib, result):
        maxscore = 0
        goodrows = []
        valcol = False
        if self.barchart_value:
            try:
                valcol = {'p': 2, 'q': 6, 'q': 3, 'z': 4}[self.barchart_value]
            except KeyError:
                pass

        maxscore = result.getMaxScore()
#        sys.stderr.write("Maxscore = {}\n".format(maxscore))
        out.write("""<CENTER><H1>{}</H1>\n""".format(lib))
        out.write("""<TABLE width='80%' class='tab'>\n""")
        for e in result.entries:
            perc = int(100.0 * e.key / maxscore)
#            sys.stderr.write("{} => {}\n".format(e.key, perc))
            out.write("""  <TR>
    <TD class='bar' style='background: linear-gradient(to right, #{} {}%, #FFFFFF 0%);' nowrap>{}</TD>
  </TR>
""".format(self.barcolor, perc, e.term))
        out.write("""</TABLE></CENTER>\n""")

    def doList(self):
        for lib in LIBS:
            sys.stdout.write(lib + "\n")
            
    def usage(self, what=None):
        if what == "upload":
            sys.stdout.write("""enrichr.py - Command-line interface to the Enrichr website

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
            sys.stdout.write("""enrichr.py - Command-line interface to the Enrichr website

Usage: enrichr.py [options] run setid library [libraries...]

This command performs enrichment analysis on the set of genes identified by `setid' (returned
by the `upload' command) against one or more Enrichr libraries. The full list of Enrichr 
libraries can be found at: http://amp.pharm.mssm.edu/Enrichr/#stats.

Output (written to standard output, or to a file if -o is specified is tab-delimited with the
following columns: library, rank, term name, P-value, Z score, combined score, overlapping genes.
See the Enrichr documentation for a description of these values. If multiple libraries are 
specified, the respective results will be concatenated in the same output; use the library id
in the first column to distinguish the different groups. If -H is specified, the table is written
in HTML format instead of tab-delimited.

Options:

  -o O | Write output to file O (default: stdout)
  -H   | Write output as HTML instead of tab-delimited
  -b   | Produce a barchart in HTML format (written to destination specified with -o)
  -s S | Sort results by S, one of p=pvalue, q=qvalue, z=zscore, c=combined score (default: {}).
  -p P | Only return terms with P-value less than or equal to P (default: {})
  -q Q | Only return terms with q-score (BH-corrected P-value) less than or equal to Q (default: {})
  -z Z | Only return terms with Z-score greater than Z (default: {})
  -c C | Only return terms with combined score greater than C (default: {})
  -n N | Only return the top N terms (default: no limit)

""".format(self.sortby, self.pval, self.qval, self.zscore, self.cscore))

        elif what == "list":
            sys.stdout.write("""enrichr.py - Command-line interface to the Enrichr website

Usage: enrichr.py list

This command prints all known libraries to standard output.

""")

        else:
            sys.stdout.write("""enrichr.py - Command-line interface to the Enrichr website

This program can be used to perform enrichment analysis on sets of genes, using the 
Enrichr engine (http://amp.pharm.mssm.edu/Enrichr/). Its basic usage is:

  enrichr.py [options] command [arguments...]

where `command' is one of: upload, list, run. Use "-h command" to display help for each command.

General options:

  -h   | Print help message
  -v   | Print version number
  -u U | Change the URL to access the Enrichr API to U
         (this should never be necessary).

""")

    def main(self):
        try:
            {"upload": self.doUpload,
             "run": self.doRun,
             "list": self.doList }[self.command]()
        except KeyError:
            self.errmsg(self.BADCMD)

if __name__ == "__main__":
    E = Enrichr("enrichr.py", version="1.0", usage=Enrichr.usage,
                errors=[('NOSETID', 'Missing set id', "The first argument to the run command should be a set ID created with the upload command."),
                        ('NOLIBS', 'Missing libraries', "Please specify at least one Enrichr library name after the set ID."),
                        ('BADLIB', 'Bad library name', "{} is not a valid Enrichr library name. Use the list command to display all libraries."),
                        ('BADCOMM', 'Communication error', "Error communicating with the Enrichr server."),
                        ('BADCMD', 'Bad command', "The first argument should be one of: upload, run, list.")
])

    E.parseCmdline(sys.argv[1:])
    E.main()
