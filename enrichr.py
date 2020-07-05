#!/usr/bin/env python

import sys
import json
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
            self.parseArgs(args, "+o,+#n,+.p,+.q,+.z,+.c,+u,H,b")
            self.outfile = self.getOpt("o", self.outfile)
            self.maxrank = self.getOpt("n")
            self.pval    = self.getOpt("p", self.pval)
            self.qval    = self.getOpt("q", self.qval)
            self.zscore  = self.getOpt("z", self.zscore)
            self.cscore  = self.getOpt("c", self.cscore)
            self.ENRICHR_URL = self.getOpt("u", self.ENRICHR_URL)
            self.html    = self.getOpt("H")
            self.barchart = self.getOpt("b")
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

if __name__ == "__main__":
    E = Enrichr("enrichr.py", version="1.0", usage=Enrichr.usage,
                errors=[('NOSETID', 'Missing set id', "The first argument to the run command should be a set ID created with the upload command."),
                        ('NOLIBS', 'Missing libraries', "Please specify at least one Enrichr library name after the set ID."),
                        ('BADLIB', 'Bad library name', "{} is not a valid Enrichr library name. Use the list command to display all libraries."),
                        ('BADCOMM', 'Communication error', "Error communicating with the Enrichr server."),
])

    E.parseCmdline(sys.argv[1:])
    if E.command == "upload":
        E.doUpload()
    elif E.command == "run":
        E.doRun()
    elif E.command == "list":
        E.doList()

