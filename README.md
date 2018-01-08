# BStoolkit
## Purpose
Tools to work with whole-genome bisulfite sequence data, written in Python 3 and bash.
## BSQC
This takes two paired bisulfite sequencing read files, and runs FastQC on them. Example usage is:
```bash
BSQC.sh -l LeftReads.fastq -r RightReads.fastq
```
It requires [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
## BSparse
This takes an output file from MethylExtract, and writes a bed file of methylation levels at each cytosine. Basic usage is:
```bash
BSparse.py -i MethyExtract.output
``` 
It can also take one or more annotation files (in GFF format), and returns the mean methylation level for each feature. This is done by using with -f argument with one or more filenames, for example:
```bash
BSparse.py -i MethyExtract.output -f features1.gff,features2.gff
```
It requires [bedtools](http://bedtools.readthedocs.io/en/latest/).

