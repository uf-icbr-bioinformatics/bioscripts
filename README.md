# bioscripts
## A collection of stand-alone scripts for bioinformatics.

This repository contains a collection of scripts for various tasks in bioinformatics. 
They were developed independently for use in various data analysis pipelines, but they
can be used independently of each other.

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
<tt>bamToWig.py</tt>   | 
<tt>bisconv.py</tt>   | 
<tt>chromCoverage.py</tt>   | 
<tt>compareIntrons.py</tt>   | 
<tt>countseqs.py</tt>   | 
<tt>dmaptools.py</tt>   | 
<tt>inserts.py</tt>   | 
<tt>mergeCols.py</tt>   | 
<tt>methreport.py</tt>   | 
<tt>methylfilter.py</tt>   | 
<tt>pileupToBED.py</tt>   | 
<tt>regionsCount.py</tt>   | 
<tt>removeN.py</tt>   | 
