#!/bin/bash

while getopts ":b:c:h" opt; do
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
			h)
				printf "Valid parameters are:\n-b (BAM file name)\n-c (Contig name)\n-g (genome fasta file name)\n-h display help message\n"
				exit 1
                        ;;
	esac
done

# extract the entries for this contig from the bam file and print to a new bam file
samtools view -b "${BAMfile}" "${Contig}" > ./MEbatch/"${Contig}"/"${Contig}".bam
# run MethylExtract on this contig
~/bin/MethylExtract_1.9.1/MethylExtract.pl seq=./MEbatch/"${Contig}"/"${Contig}".fasta inDir=./MEbatch/"${Contig}"/ flagW=99,147 flagC=83,163 context=ALL bedOut=Y p=1 chromDiv=800 memNumReads=500000 minDepthMeth=10 chromSplitted=Y

