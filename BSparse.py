#!/usr/bin/env python3

import os
import sys
import subprocess
import random as random
import shutil
import re
import argparse

###########################
## USER ARGUMENT PARSING ##
###########################
parser = argparse.ArgumentParser(description='Read arguments')
parser.add_argument('-i', '--infile', type=str, help='input file (should be a .output file from MethylExtract)')
args = parser.parse_args()
# input fasta file parsing
InputFile = args.infile
if InputFile is not None:
	print('Input file is ' + InputFile)
else:
	print('ERROR: no input file specified')

#################
## ACTUAL CODE ##
#################

## function to take a line from a MethyExtract output file, and return a bed-formatted line with % methylation on both strands
def MEconverter(line):
	Line = line.split("\t")
	Chromosome = Line[0]
	Start = str(int(Line[1])-1)
	End = str(int(Line[1])+1)
	Name = "m"+Line[2]
	# note: "Score" here is the mean of the % methylated reads on both strands
	Score = str((((int(Line[3])/int(Line[4]))*100)+((int(Line[6])/int(Line[7]))*100))/2)
	BedLine = Chromosome + "\t" + Start + "\t" + End + "\t" + Name + "\t" + Score
	return(BedLine)

# read input file, keep only lines with methylation called on both strands
BothStrands = []
for line in open(InputFile, "r"):
	if not line.startswith("#"):
		temp = line.split("\t")
		if not "." in temp:
			BothStrands.append(line.strip("\n"))

# convert lines to bed format
BedFormatLines = []
for i in BothStrands:
        BedFormatLines.append(MEconverter(i))

# format output bed file
BedOutput = ""
for i in BedFormatLines:
	BedOutput += i + "\n"
BedOutfile = open(InputFile.replace(".output",".bed"),"wt")
BedOutfile.write(BedOutput)
BedOutfile.close()
print("Bed file written to "+InputFile.replace(".output",".bed"))

# calculate genome-wide methylation level
Ccount = 0
MethylationTotal = 0
for i in BedFormatLines:
	Ccount += 1
	MethylationTotal += float(i.split("\t")[4])
MethylationLevel = MethylationTotal/Ccount
print("Genome-wide methylation level = " + str(MethylationLevel))

# use bedtools to intersect bed with gene models, generating mean methylation level for each gene
# sort bedfile and gff file
# use map to generate methylation levels for each feature (write into function)
