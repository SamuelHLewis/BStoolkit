#!/bin/bash

while getopts ":b:c:g:h" opt; do
	case $opt in
			b)
				# set BAM file name
				BAMfile=$OPTARG
				printf "BAM file (-b): ${BAMfile}\n"
				;;
			c)
				# set contig name
				Contig=$OPTARG
				printf "Contig (-c): ${Contig}\n"
				;;
			g)
				# set genome fasta file name
				GenomeFasta=$OPTARG
				printf "Genome fasta file (-g): ${GenomeFasta}\n"
				;;
			h)
				printf "Valid parameters are:\n-b (BAM file name)\n-c (Contig name)\n-g (genome fasta file name)\n-h display help message\n"
				exit 1
                        ;;
	esac
done

# make new directory for this contig
mkdir ./MEbatch/"${Contig}"
# grab the name and sequence for this contig from the genome fasta file (NB: genome file must be unwrapped)
grep -A1 "${Contig}" "${GenomeFasta}" > ./MEbatch/"${Contig}"/"${Contig}".fasta
# extract the entries for this contig from the bam file, and sort the bam file
samtools view -h -F 0x4 -q 10 "${BAMfile}" "${Contig}" | samtools view -hbS - > ./MEbatch/"${Contig}"/"${Contig}".bam
# run MethylExtract on this contig
~/bin/MethylExtract_1.9.1/MethylExtract.pl seq=./MEbatch/"${Contig}"/"${Contig}".fasta inDir=./MEbatch/"${Contig}"/ flagW=99,147 flagC=83,163 context=ALL bedOut=Y p=1 chromDiv=800 memNumReads=500000 minDepthMeth=10 chromSplitted=Y

