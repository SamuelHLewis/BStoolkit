#!/usr/bin/env python3

import os
import sys
import shutil
import re
import argparse

###########################
## USER ARGUMENT PARSING ##
###########################
parser = argparse.ArgumentParser(description='Read arguments')
parser.add_argument('-f', '--fasta', type=str, help='fasta file')
args = parser.parse_args()
# input fasta file parsing
FastaFile = args.fasta
if FastaFile is not None:
	print('Fasta file = ' + FastaFile)
else:
	print('ERROR: no fasta file specified')

##########################
## FASTA FILE SPLITTING ##
##########################

Name = ""
Sequence = ""

# trigger to stop file writing after the first name is read in
StartWriting = False

# for each chromosome/contig, write name and sequence to a separate file
for line in open(FastaFile):
	if line.startswith(">"):	
		if StartWriting == True:
			if os.path.exists(Name):
				shutil.rmtree(Name)
			os.mkdir(Name)
			outfile = open("./" + Name + "/" + Name + ".fasta", "wt")
			outfile.write(">" + Name + "\n" + Sequence + "\n")
			outfile.close()
		Name = line.strip(">").strip("\n")
	else:
		Sequence += line.strip("\n")
		StartWriting = True

# write the last name & sequence to a file
if os.path.exists(Name):
	shutil.rmtree(Name)
os.mkdir(Name)
outfile = open("./" + Name + "/" + Name + ".fasta", "wt")
outfile.write(">" + Name + "\n" + Sequence + "\n")
outfile.close()

print("Sequences written to separate files")
