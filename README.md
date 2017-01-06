# bioscripts
## A collection of stand-alone scripts for bioinformatics.

This repository contains a collection of scripts for various tasks in bioinformatics. 
They were developed independently for use in various data analysis pipelines.

### Dependencies

All scripts in **bioscripts** are designed to be executed as stand-alone
command-line programs in a Python 2.7 environment. If your Python
executable is not in a standard location, please adjust the first line
(starting with #!) accordingly.

All scripts require the <tt>Script.py</tt> file, that should therefore be
located in the same directory as the script you are executing. In
addition: 

* <tt>bisconv.py</tt> and <tt>regionscount.py</tt> require the pysam library;
* <tt>methreport.py</tt> and <tt>methylfilter.py</tt> require the SeqIO module from BioPython.

### Common arguments

All scripts in *bioscripts* accept the following command-line
arguments:

Argument | Description
-----------|------------
-h, --help | Print usage message.
-v, --version | Print version number.
-E [n]        | Decode error code. If called with an error code *n*, prints the corresponding error message. If called without *n*, displays all error codes for this script with the associated error message.

The -E argument is designed for automated error handling in
scripts. After calling a script (e.g. <tt>prog.py</tt>), check the return value (in the $?
variable), and if it is not zero, print the corresponding error
message with:

```
prog.py -E $?
```

### List of scripts

The following table lists all scripts in this package with a short
description of their purpose.

Name | Description
-----|------------
<tt>bamToWig.py</tt>        | Convert BAM file to WIG track for the UCSC genome browser.
<tt>bisconv.py</tt>         | Extract aligned reads from BAM file based on conversion strand.
<tt>chromCoverage.py</tt>   | Report per-chromosome coverage.
<tt>compareIntrons.py</tt>  | Analyze intron retention.
<tt>countseqs.py</tt>       | Count sequences in fasta/fastq[.gz] files.
<tt>dmaptools.py</tt>       | Merge methylation data.
<tt>mergeCols.py</tt>       | Merge columns from multiple files.
<tt>methreport.py</tt>      | Report methylation rate at CG and GC positions.
<tt>methylfilter.py</tt>    | Separate sequences by average methylation.
<tt>pileupToBED.py</tt>     | Convert a samtools pileup to a BED file.
<tt>regionsCount.py</tt>    | Compute coverage from BAM file in specified regions.
<tt>removeN.py</tt>         | Remove Ns from sequences in FASTA file.
