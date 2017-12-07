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

## function to take a line from a MethyExtract output file, and return a bed-formatted line with % methylation on both strands
def MEconverter(line):
	Line = line.split("\t")
	Chromosome = Line[0]
	Start = str(int(Line[1])-1)
	End = str(int(Line[1])+1)
	Name = "m"+Line[2]
	Score = str((int(Line[5])+int(Line[8]))/2)
	BedLine = Chromosome + "\t" + Start + "\t" + End + "\t" + Name + "\t" + Score
	return(BedLine)

# read input file, keep only lines with methylation called on both strands
BothStrands = []
for line in open(InputFile, "r"):
	temp = line.split("\t")
	if not "." in temp:
		BothStrands.append(line.strip("\n"))

for i in BothStrands:
	print(i)

for i in BothStrands:
	if not i.startswith("#"):
		print(MEconverter(i))
