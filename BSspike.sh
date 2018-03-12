#!/bin/bash

# defaults for cores (single)
Cores=1

# command line option parsing
while getopts ":c:g:l:r:b:h" opt; do
	case $opt in
		c)
			if [ $OPTARG -le 1 ]
			then
				# exit if <1 cores are specified
				printf "ERROR: cores (-c) must be >0\n"
				exit 1
			else
				# set number of cores to user input
				printf "cores (-c): $OPTARG\n"
				Cores=$OPTARG
			fi
			;;
		g)
			# set spike-in genome fasta file directory
			printf "Directory containing genome of spike-in (-g): $OPTARG\n"
			Genome=$OPTARG
			;;
		l)
			# set lefthand read file
			printf "lefthand read file (-l): $OPTARG\n"
			LeftReads=$OPTARG
			;;
		r)
			# set righthand read file
			printf "righthand read file (-r): $OPTARG\n"
			RightReads=$OPTARG
			;;
		h)
			printf "Valid parameters are:\n-c (number of cores, default=1)\n-l (lefthand read file)\n-r (righthand read file)\n-h display help message\n"
			exit 1
			;;
	esac
done

# delete the "SpikeIn" directory if it already exists
if [ -d SpikeIn ]; then
        rm -r SpikeIn;
fi
# make a "SpikeIn" directory and move into it
mkdir SpikeIn
cd SpikeIn
# use Bismark to map reads to spike-in genome in non-directional mode
bismark --multicore $Cores --un --non_directional --genome $Genome -1 ../$LeftReads -2 ../$RightReads
# generate summary report
bismark_methylation_extractor -p --gzip ../*_bismark_bt2_pe.bam

