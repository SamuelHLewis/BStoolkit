#!/bin/bash

## split genome fasta file
#csplit --prefix=seq --digits=6 --silent Genomic.fas '%>%' '/^>/' '{*}'
## move genome fasta file out of this directory
#mv Genomic.fas ../

# for each individual chromosome/contig fasta file...
for file in $( ls ); do
	# grab the header, strip off the ">", and set this as the Name
	Name=$( head -n1 $file | cut -d'>' -f2)
	# make a new directory for this contig
	mkdir $Name
	# rename the fasta file, and move it into the directory (with the same name)
	mv $file ./"${Name}"/"${Name}".fasta
done
