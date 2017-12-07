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
		b)
			if [ $OPTARG -le 0 ]
			then
				# exit if <0 trimmed bases are specified
				printf "ERROR: trimmed bases (-b) must be 0 or higher\n"
				exit 1
			else
				# set number of trimmed bases to user input
				printf "trimmed bases (-b): $OPTARG\n"
				TrimmedBases=$OPTARG
			fi
			;;
		h)
			printf "Valid parameters are:\n-c (number of cores, default=1)\n-l (lefthand read file)\n-r (righthand read file)\n-b (trimmed bases)\n-h display help message\n"
			exit 1
			;;
	esac
done

# trim adapters from reads, and user-specified number of bases from 5' and 3' ends of each read
trim_galore --paired --retain_unpaired --illumina --clip_R1 $TrimmedBases --clip_R2 $TrimmedBases --three_prime_clip_R1 $TrimmedBases --three_prime_clip_R2 $TrimmedBases $LeftReads $RightReads

