# bioscripts
## A collection of stand-alone scripts for bioinformatics.

This repository contains a collection of scripts for various tasks in bioinformatics. 
They were developed independently for use in various data analysis pipelines, but they
can be used independently of each other.

### Dependencies

All scripts in *bioscripts* are designed to be executed as stand-alone
command-line programs in a Python 2.7 environment. If your Python
executable is not in a standard location, please adjust the first line
(starting with #!) accordingly.

All scripts require the Script.py file, that should therefore be
located in the same directory as the script you are executing. In
addition: 

* bisconv.py and regionscount.py require the pysam library;
* methreport.py and methylfilter.py require the SeqIO module from BioPython.

### Common arguments
