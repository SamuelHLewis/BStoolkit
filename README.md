# BStoolkit
## Purpose
Tools to work with whole-genome bisulfite sequence data, written in Python 3 and bash.
## BSQC
This takes two paired bisulfite sequencing read files, and runs FastQC on them. Example usage is:
```bash
BSQC.sh -l LeftReads.fastq -r RightReads.fastq
```
It requires [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
## BSspike
This maps paired-end bisulfite reads to the genome of a spike-in (e.g. known unmethylated DNA). Basic usage is:
```bash
BSspike.sh -g ~/Path/To/SpikeIn/Genome/ -l LeftReads.fastq -r RightReads.fastq
```
By default it runs on a single core, but it can be run on multiple cores using the -c argument, for example:
```bash
BSspike.sh -c 24 -g ~/Path/To/SpikeIn/Genome/ -l LeftReads.fastq -r RightReads.fastq
```
It requires [Bismark](https://github.com/FelixKrueger/Bismark).
## BStrim
This trims bases from the 3' end of a set of paired-end bisulfite sequencing reads. The default number of bases to trim is 10. Basic usage is:
```bash
BStrim.sh -l LeftReads.fastq -r RightReads.fastq
```
A different number of reads to trim can be specified with the -b argument, for example to specify 20 bases:
```bash
BStrim.sh -b 20 -l LeftReads.fastq -r RightReads.fastq
```
By default it runs on a single core, but it can be run on multiple cores using the -c argument, for example:
```bash
BStrim.sh -c 24 -l LeftReads.fastq -r RightReads.fastq
```
It requires [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
## MEbatch
This uses MethylExtract to estimate the methylation level at all cytosines on a **specific** chromosome/contig. Basic usage is:
```bash
MEbatch.sh -b MappedToAllChromosomes.bam -c Chromosome1 -g Genome.fasta
```
It requires [samtools](http://www.htslib.org/) and [MethylExtract](http://bioinfo2.ugr.es/MethylExtract/).

Analysing chromosomes/contigs separately is particularly useful when analysing highly-fragmented genomes with large numbers of contigs, which can cause MethylExtract to crash if the entire bam file is analysed at once. In this case, MEbatch can be used to run MethylExtract on all contigs in batches with a command along the lines of:
```bash
# write the name of each contig in the bam file to a contig.names file
samtools view -H MappedToAllChromosomes.bam | awk -F"\t" '/@SQ/{print $2}' |  cut -d":" -f2 > contig.names
# run contigs in parallel
cat contig.names | parallel "MEbatch.sh -b MappedToAllChromosomes.bam -c {} -g Genome.fasta"
```
This type of parallel usage requires [GNU parallel](https://www.gnu.org/software/parallel/).
## UpstreamFinder
This finds the coordinates of upstream regions for each gene in a genome, with reference to the first exon of each gene. It takes a GFF file and an integer length of upstream sequence as input, and outputs a GFF file of upstream sequencefeatures. Basic usage is:
```bash
UpstreamFinder.py -i Annotation.gff -l 1000
```
## CountConcat
This takes an annotation file, a methylation levels file (a .output file from MethylExtract), a bam file of mapped RNA-Seq reads, a bam file of mapped siRNAs, and a bam file of mapped piRNAs. It counts the number of RNA-Seq reads, siRNAs and piRNAs mapping to each feature in the annotation file, and concatenates these counts together with the methylation levels into a final "Concatenated.counts" file. Example usage is:
```bash
CountConcat.py -a Annotations.gff -m MethylationLevels.output -r RNASeq.bam -p piRNA.bam -s siRNA.bam
```
It requires [bedtools](http://bedtools.readthedocs.io/en/latest/).


