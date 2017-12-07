#!/bin/bash

# defaults for cores (single)
Cores=1
TrimmedBases=10

# command line option parsing
while getopts ":c:l:r:b:h" opt; do
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

# check if "FastQC" directory exists within the working directory, and make if not
if [ ! -d FastQC ];
	then mkdir ./FastQC;
fi

# QC lefthand and righthand reads, and move fastQC files to FastQC subdir
fastqc -t $Cores $LeftReads
fastqc -t $Cores $RightReads
mv *fastqc* ./FastQC

