#!/usr/bin/env python

import sys
import subprocess as sp
import Script
import Utils
from BIutils import BIexperiment

### Utils

def parseContrasts(ctrs):
    return [ c.split("^") for c in ctrs ]

def los(l):
    """Format l as a list of strings separated by commas and enclosed in double quotes."""
    return ", ".join([ '"' + s + '"' for s in l])

### Class

class DESeq2(Script.Script):
    args = []
    infiles = []
    mode = "counts"             # -m, --mode
    conditions = []             # -c, --conditions
    contrasts = []              # -d, --contrasts
    pval = 0.01                 # -p, --pval
    mincount = 5                # -t, --min_count
    fmode = "avg"               # -f, --filter
    write_norm = False          # -wn, --write_norm
    write_full = True           # -wf, --write_full
    write_diff = True           # -wd, --write_diff
    mds_file = None             # -g, --mds_plot
    script_file = "/dev/stdout" # -o, --outfile
    execute = False             # -x, --execute
    _experiment = None

    def init(self):
        args = sys.argv[1:]
        prev = ""
        self.standardOpts(args)
        self._experiment = BIexperiment.Experiment()
        for a in args:
            if prev in [ "-n", "--names" ]:
                self.names = Utils.parseList(a)
                prev = ""
            elif prev in [ "-c", "--conditions" ]:
                if self.isFile(a, error=False):
                    self._experiment.initConditionsFromFile(a)
                else:
                    self.parseConditions(a)
                prev = ""
            elif prev in [ "-d", "--contrasts" ]:
                if self.isFile(a, error=False):
                    self._experiment.initContrastsFromFile(a)
                else:
                    contrasts = parseContrasts(Utils.parseList(a))
                    for contr in contrasts:
                        print(contr)
                        self._experiment.addContrast(contr[0], contr[1])
                prev = ""
            elif prev in ["-m", "--mode"]:
                self.mode = a
                prev = ""
            elif prev in [ "-p", "--pval" ]:
                self.pval = self.toFloat(a)
                prev = ""
            elif prev in [ "-t", "--min_count" ]:
                self.mincount = self.toInt(a)
                prev = ""
            elif prev in [ "-f", "--filter" ]:
                self.fmode = a
                prev = ""
            elif prev in [ "-wn", "--write_norm" ]:
                self.write_norm = a
                prev = ""
            elif prev in [ "-wf", "--write_full" ]:
                self.write_full = ( a in "Yy")
                prev = ""
            elif prev in [ "-wd", "--write_diff" ]:
                self.write_diff = ( a in "Yy")
                prev = ""
            elif prev in [ "-g", "--mds_plot" ]:
                self.mds_file = a
                prev = ""
            elif prev in [ "-o", "--outfile" ]:
                self.script_file = a
                prev = ""
            elif a in [ "-n", "--names", "-c", "--conditions", "-d", "--contrasts", "-m", "--mode", "-p", "--pval", "-t", "--min_count", "-wn", "--write_norm", "-wf", "--write_full", "-wd", "--write_diff",
                        "-g", "--mds_plot", "-o", "--outfile", "-f", "--filter"]:
                prev = a
            elif a in [ "-x", "--execute" ]:
                self.execute = True
            else:
                self.infiles.append(a)
        return self.validate()

    def parseConditions(self, s):
        parts = s.split(":")
        for p in parts:
            words = p.split(",")
            self._experiment.addCondition(words[0], words[1:])
#        sys.stderr.write("{}\n".format(self._experiment.conditions))
#        sys.stderr.write("{}\n".format(self._experiment.condsamples))

    def getColumns(self, filename):
        result = {}
        with open(filename, "r") as f:
            line = f.readline().rstrip().split("\t")
        idx = 1
        for name in line[1:]:
            result[name] = idx
            idx += 1
        return result

    def generate(self):
        sys.stderr.write("Writing R script to {}\n".format(self.script_file))
        with open(self.script_file, "w") as out:
            if self.mode == "edger":
                self.generate_edger(out)
            else:
                self.generate_s(out)

    def generate_s(self, out):
        if self.mode == "counts":
            columns = self.getColumns(self.infiles[0])
            # sys.stderr.write("{}\n".format(columns))
            params = { "inputs": los(self.infiles),
                       "smps": los(self._experiment.samples),
                       "conds": los(self._experiment.getSampleLabels(list(columns.keys()))),
                       "colnumbers": columns,
                       "mincov": self.mincount }
        else:
            params = { "inputs": los(self.infiles),
                       "smps": los(self._experiment.samples),
                       "conds": los(self._experiment.sampleLabels()),
                       "colnumbers": columns,
                       "mincov": self.mincount }

        out.write("""# DESeq2 driver script for {} data
library(DESeq2)
""".format(self.mode))
        if self.mode == "rsem":
            out.write("""library(tximport)
a = c({inputs})
names(a) = c({smps})
data = tximport(a, type="rsem", txIn=FALSE, txOut=FALSE)
sampleTable <- data.frame(condition = factor(c({conds})))
rownames(sampleTable) <- colnames(data$counts)
data$length[data$length == 0] <- 1
dds <- DESeqDataSetFromTximport(data, sampleTable, ~condition)
""".format(**params))
        elif self.mode == "counts":
            out.write("""datafile = {inputs}
counts = round(as.matrix(read.csv(datafile, sep='\\t', row.names=1)))
levels = c({smps})
labels = c({conds})
sampleTable = data.frame(condition=factor(labels))
rownames(sampleTable) = colnames(counts)
dds = DESeqDataSetFromMatrix(countData = counts, colData = sampleTable, design = ~ condition)
""".format(**params))

        out.write("""
# Normalization
dds = estimateSizeFactors(dds)
nc = counts(dds, normalize=TRUE)
""".format(**params))
        if self.write_norm:
            out.write("""
# Write normalized counts
write.table(nc, file="{}", sep='\\t')
""".format(self.write_norm))

        self.writeFiltering(out, params)
        
        out.write("""
# Distance
# sampledist = dist(t(data$counts))

dds = DESeq(dds)

""".format(**params))
        for contr in self._experiment.contrasts:
            params = { "test": contr[0],
                       "ctrl": contr[1],
                       "pval": self.pval}
            out.write("""# {test} vs {ctrl}

res = results(dds, contrast=c("condition", "{test}", "{ctrl}"))
resOrdered = res[order(res$pvalue),]

# Write full results
write.table(resOrdered, file="{test}.vs.{ctrl}.diff.csv", sep="\\t")

# Filter by pval, sort by fold change
resSig = subset(res, padj < {pval})
resSigOrder = resSig[ order(resSig$log2FoldChange), ]

write.table(resSigOrder, file="{test}.vs.{ctrl}.sig.csv", sep="\\t")

""".format(**params))

    def writeFiltering(self, out, params):
        mc = params['mincov']
        if self.fmode[0] == 'a':
            nsamples = len(self._experiment.samples)
            out.write("""
# Filtering - average of counts across all samples
keep = (rowSums(counts(dds)) / {}) >= {}
dds = dds[keep,]
""".format(nsamples, mc))
        elif self.fmode[0] == 'c':
            tests = []
            idx = 1
            colnums = params["colnumbers"]
            for c in self._experiment.conditions:
                cs = self._experiment.condsamples[c]
                cidx = [ str(colnums[x]) for x in cs ]
                #sys.stderr.write("{}\n{}\n".format(cs, cidx))
                col_names = "c(" + ",".join(cidx) + ")"
                tests.append("(rowSums(cts[,{}, drop=FALSE]) > {})".format(col_names, mc*len(cidx)))
                idx += len(cs)
            out.write("""
# Filtering - average of counts in each condition
cts = counts(dds)
keep = {}
dds = dds[keep,]
            """.format(" & ".join(tests)))

    def generate_edger(self, out):
        columns = self.getColumns(self.infiles[0])
        nconds = len(self._experiment.conditions)
        params = { "input": self.infiles[0],
                   "columns": los(columns),
                   "smps": los(self._experiment.samples),
                   "conds": los(self._experiment.getSampleLabels(self._experiment.samples)),
                   "allconds": los(self._experiment.conditions),
                   "mincov": self.mincount }
        params["mds"] = """# Draw MDS plot
png(filename="{}", width=1000, height=1000, units="px")
plotMDS(y, method="bcv", col=as.numeric(y$samples$group))
legend("bottomleft", as.character(unique(y$samples$group)), col=1:groupnum, pch=20)
dev.off()
""".format(self.mds_file) if self.mds_file else ""
        params["norm"] = """write.table(y, file="{}", sep="\\t", row.names=TRUE, col.names=TRUE)""".format(self.write_norm) if self.write_norm else ""
        out.write("""# Differential analysis with edgeR

library(edgeR)

do_contrast <- function(filename, fit, contrast) {{
    fnall = paste(filename, "diff.txt", sep=".")
    fnsig = paste(filename, "sig.txt", sep=".")
    lrt <- glmLRT(fit, contrast=contrast)
    message(paste("Writing", fnall, "and", fnsig))
    data = topTags(lrt, n=nrow(lrt), sort.by="logFC")
    write.table(data, file=fnall, sep='\\t', col.names=TRUE, row.names=TRUE)
    data = topTags(lrt, n=nrow(lrt), sort.by="logFC", p.value=0.05)
    write.table(data, file=fnsig, sep='\\t', col.names=TRUE, row.names=TRUE)
}}

datafile = "{input}"
samples = c({smps})
type = c({conds})
conditions = c({allconds})
groupnum = length(conditions)

counts = round(as.matrix(read.csv(datafile, sep='\\t', row.names=1)))
counts = counts[,samples]       # restrict matrix to samples we actually have

message(paste("Number of samples:", length(samples)))
message(paste("Number of conditions:", length(conditions)))

# Create DGEList object, filter it, and normalize
y = DGEList(counts=counts, genes=rownames(counts), group=type)
keep = rowSums(cpm(y)>1) >= length(samples)/2
y = y[keep, , keep.lib.sizes=FALSE]
y =  calcNormFactors(y, method="upperquartile")
message(paste("Number of genes used for differential analysis: ", nrow(y$genes), sep=""))
{norm}
{mds}
# Create design matrix, and rearrange columns to desired order
design <- model.matrix(~ 0 + type)
colnames(design) = sub("type", "", colnames(design))
design <- design[, conditions]

# Fit linear model
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
message("Fitting done")

# Contrasts
""".format(**params))

        for test, ctrl in self._experiment.contrasts:
            idx1 = self._experiment.conditions.index(test)
            idx2 = self._experiment.conditions.index(ctrl)
            vec = [0]*nconds
            vec[idx1] = 1
            vec[idx2] = -1
            out.write("""do_contrast("{}.vs.{}", fit, c({}))\n""".format(test, ctrl, ",".join([str(x) for x in vec])))
        
    def validate(self):
        if self.execute and not self.script_file:
            sys.stderr.write("Error: -x requires an output file (-o option).\n")
            return False
        if not self._experiment.conditions:
            sys.stderr.write("Error: no conditions defined!.\n")
        return True

        # if self.names:
        #     if len(self.names) != len(self.infiles):
        #         sys.stderr.write("Error: there are {} input files but {} sample names.\n".format(len(self.infiles), len(self.names)))
        #         return False
        # else:
        #     for f in self.infiles:
        #         fn = os.path.split(f)[1]
        #         self.names.append(fn.split(".")[0])
        # for contr in self.contrasts:
        #     if len(contr) != 2:
        #         sys.stderr.write("Error: contrasts should have the form S1^S2 where S1 and S2 are condition names.\n")
        #         return False
        #     if contr[0] not in self.conditions:
        #         sys.stderr.write("Error {} is not a condition name.\n".format(contr[0]))
        #         return False
        #     if contr[1] not in self.conditions:
        #         sys.stderr.write("Error {} is not a condition name.\n".format(contr[1]))
        #         return False
        # return True

    def run_script(self):
        sys.stderr.write("Executing R script...\n")
        with open(self.script_file + ".stdout", "w") as out:
            with open(self.script_file + ".stderr", "w") as err:
                sp.call("Rscript " + self.script_file, shell=True, stdout=out, stderr=err)
        sys.stderr.write("Done.\n")

def usage(what=None):
    if what == "options":
        sys.stdout.write("""Usage: rundeseq2.py [options] filenames...

All options have a short and a long form. In the following table,
S indicates a single string value, L indicates a list of strings
separated by commas (e.g. one,two,three), F indicates a filename,
N indicates a number, and Y indicates a flag (Y and y mean true,
everything else means false). Multiple letters indicate that all
specified option types can be used.

 -h [S]          | Display help.
 --help [S]      |

 -v              | Show version number.
 --version       |

 -m S            | Specify mode. Should be one of "counts" (DESeq2 on counts),
 --mode S        | "rsem" (DESeq2 on RSEM files),  or "edger" (edgeR on counts). 
                 | Default: {}.

 -c SF           | Specify condition and sample names, or
 --conditions SF | read them from a file.

 -d SF           | Specify contrasts or read them from a file.
 --contrasts SF  | On the command line, contrasts should have the form:
                 | A^B,C^D (to compare A vs B and C vs D).

 -t N            | Minimum count (genes having less than this
 --min_count N   | number of reads will be ignored) Default: {}.

 -f S            | Filtering strategy. Can be one of: "avg" (average counts
  --filter S     | for a gene in all samples should be higher than min_count),
                 | "cond" (average counts for a gene in each condition should
                 | be higher than min_count. Default: {}.

 -p N            | Specify P-value threshold for
 --pval          | significance. Default: {}.

 -wn F           | Write normalized count matrix to the
 --write_norm F  | specified file. Default: {}.

 -g G            | Create MDS plot, saving it to 
 --mds_plot G    | file G in PNG format.

 -wf Y           | If true, write differential expression values
 --write_full    | for ALL genes to a tab-delimited file, named
                 | A.vs.B.diff.csv. Default: {}.

 -wd Y           | If true, write differential expression values
 --write_diff Y  | for all SIGNIFICANT genes (those with a P-value less than
                 | the limit specified with -p) to a tab-delimited file, named
                 | A.vs.B.sig.csv. Default: {}.

 -o F            | Write script to this file. If not specified,
 --outfile F     | writes to standard output.

 -x              | Execute script after writing it. Requires -o,
 --execute       | and Rscript in PATH. Default: {}.

""".format(DESeq2.mode, DESeq2.mincount, DESeq2.fmode, DESeq2.pval, DESeq2.write_norm,
           DESeq2.write_full, DESeq2.write_diff, DESeq2.execute))
    elif what == "config":
        sys.stdout.write("""Usage: rundeseq2.py [options] filenames...

The experiment structure can be defined by two files, containing information
about samples and contrasts respectively. Both files are tab-delimited with
two columns. The first one (that should be specified with the -c option) has
the following format:

  C1    S1,S2,S3
  C2    S4,S5,S6

where Cn indicates a condition name and Sn are the samples belonging to each 
condition. The second file (specified with the -d option) has the format:

  C1   C2

where C1 and C2 are two conditions to be compared to each other (C1 is test,
C2 is control / baseline condition).

The experiment definition can also be provided on the command line; use
-h examples for more details.

""")
    elif what == "examples":
        sys.stdout.write("""Usage: rundeseq2.py [options] filenames...

As an example, assume we performed an experiment with a control condition (WT)
and two different knockouts (KO1, KO2), and that we have three replicates for 
each condition. The conditions file will be:

  WT   WT_A,WT_B,WT_C
  KO1  KO1_A,KO1_B,KO1_C
  KO2  KO2_A,KO2_B,KO2_C

And the contrasts file will be:

  KO1  WT
  KO2  WT

The same could be accomplished specifying the experiment definition directly on
the command line:

  -c   WT,WT_A,WT_B,WT_C:KO1,KO1_A,KO1_B,KO1_C:KO2,KO2_A,KO2_B,KO2_C
  -d   KO1^WT,KO2^WT

""")

    else:
        sys.stdout.write("""Usage: rundeseq2.py [options] filenames...

This program writes a script to perform differential RNA-Seq analysis
using DESeq2 and optionally runs it. It handles the general case in
which there are multiple experimental conditions, each represented
by multiple samples (replicates), and you want to perform multiple
contrasts between pairs of conditions.

Use:
 "rundeseq2.py -h options" to get help on command options.
 "rundeseq2.py -h config"  to get help on experiment definition files.
 "rundeseq2.py -h examples" for usage examples.

""")

### 

if __name__ == "__main__":
    D = DESeq2("rundeseq2.py", usage=usage, copyright="(c) 2019, A. Riva, ICBR Bioinformatics Core, University of Florida")
    D.generate()
    if D.execute:
        D.run_script()
