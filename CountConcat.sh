#!/bin/bash

# defaults for cores (single)
TrimmedBases=10

# command line option parsing
while getopts ":a:m:r:p:s:h" opt; do
	case $opt in
		a)
			# set annotation file
			printf "Feature annotation file (-a): $OPTARG\n"
			AnnotationFile=$OPTARG
			;;
		m)
			# set methylation levels file
			printf "Methylation levels file (-m): $OPTARG\n"
			MethylationFile=$OPTARG
			;;
		r)
			# set RNA-Seq bam file
			printf "RNA-Seq bam file (-r): $OPTARG\n"
			RNASeqFile=$OPTARG
			;;
		p)
			# set piRNA bam file
			printf "piRNA bam file (-p): $OPTARG\n"
			piRNAFile=$OPTARG
			;;
		s)
			# set siRNA bam file
			printf "siRNA bam file (-s): $OPTARG\n"
			siRNAFile=$OPTARG
			;;
		h)
			printf "Valid parameters are:\n-a (annotation file of features in GFF or BED format)\n-m (methylation levels file)\n-r (RNA-Seq bam file)\n-p (piRNA bam file)\n-s (siRNA bam file)\n-h display help message\n"
			exit 1
			;;
	esac
done

# check that input variables are not empty
if [ -z ${AnnotationFile} ]
then
        printf "ERROR: no feature annotation file (-a) specified\n"
        exit 1
fi
if [ -z ${MethylationFile} ]
then
        printf "ERROR: no methylation levels file (-m) specified\n"
        exit 1
fi
if [ -z ${RNASeqFile} ]
then
        printf "ERROR: no RNA-Seq bam file (-r) specified\n"
        exit 1
fi
if [ -z ${piRNAFile} ]
then
        printf "ERROR: no piRNA bam file (-p) specified\n"
        exit 1
fi
if [ -z ${siRNAFile} ]
then
        printf "ERROR: no siRNA bam file (-s) specified\n"
        exit 1
fi

# count RNA-Seq
echo "Counting RNA-Seq reads"
bedtools coverage -s -counts -a "${AnnotationFile}" -b "${RNASeqFile}" > RNA.count
# count siRNA
echo "Counting siRNAs"
bedtools coverage -s -counts -a "${AnnotationFile}" -b "${siRNAFile}" > siRNA.count
# count piRNA
echo "Counting piRNAs"
bedtools coverage -s -counts -a "${AnnotationFile}" -b "${piRNAFile}" > piRNA.count
# cut out RNAseq, siRNA and piRNA count columns
echo "Combining all counts"
awk '{print $NF}' < RNA.count > RNA.count.vector
awk '{print $NF}' < siRNA.count > siRNA.count.vector
awk '{print $NF}' < piRNA.count > piRNA.count.vector
# paste them onto the methylation levels file
cat "${MethylationFile}" | paste - RNA.count.vector | paste - siRNA.count.vector | paste - piRNA.count.vector > temp.all.count
# create header for file
echo -e "Chromosome\tProgram\tFeature\tStart\tEnd\tINTENTIONALLYBLANK\tStrand\tINTENTIONALLYBLANK\tName\tIndividualCytosineMethylation\tMeanMethylation\tRNAseq\tsiRNA\tpiRNA" > Header.count
# add header to concatenated count file
cat Header.count temp.all.count > Concatenated.counts
# remove intermediate files
rm *.count
rm *.vector
# print end of program message
echo "Concatenated count file written to Concatenated.counts"

