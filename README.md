# BStoolkit
## Purpose
Tools to work with whole-genome bisulfite sequence data 
## Requirements
Written in Python 3 and bash
Requires:
[bedtools](http://bedtools.readthedocs.io/en/latest/)
## BSparse
This script takes an output file from MethylExtract, and writes a bed file of methylation levels at each cytosine. It can also take one or more annotation files (in GFF format), and returns the mean methylation level for each feature. Basic usage is:
```bash
BSparse.py -i MethyExtract.output
``` 
BSparse.py can also take one or more annotation files (in GFF format), and returns the mean methylation level for each feature. This is done by using with -f argument with one or more filenames, for example:
```bash
BSparse.py -i MethyExtract.output -f features1.gff,features2.gff
```

