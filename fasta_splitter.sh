#!/bin/bash

# command line option parsing
while getopts ":i:h" opt; do
	case $opt in
		i)
			# set input file
			printf "input file (-i): $OPTARG\n"
			InputFile=$OPTARG
			;;
		h)
			# print help message
			printf "Valid parameters are:\n-i (input fasta file)\n-h (display help message)\n"
	esac
done

# split genome fasta file
csplit --prefix=seq --digits=6 --silent "${InputFile}" '%>%' '/^>/' '{*}'

# for each individual chromosome/contig fasta file (excluding the input file)...
for file in $( ls | grep -v -F ${InputFile} ); do
	# grab the header, strip off the ">", and set this as the Name
	Name=$( head -n1 $file | cut -d'>' -f2)
	# make a new directory for this contig
	mkdir $Name
	# rename the fasta file, and move it into the directory (with the same name)
	mv $file ./"${Name}"/"${Name}".fasta
done

printf "Input file split into separate fasta files\n"

