#!/usr/bin/env python3

import os
import sys
import subprocess
import random as random
import shutil
import re
import argparse

#############################
## DEFAULT ARGUMENT VALUES ##
#############################
FeatureFiles = None

###########################
## USER ARGUMENT PARSING ##
###########################
parser = argparse.ArgumentParser(description='Read arguments')
parser.add_argument('-i', '--infile', type=str, help='input file (should be a .output file from MethylExtract)')
parser.add_argument('-f', '--featurefiles', type=str, help='list of feature files in gff format (comma-separated)')
args = parser.parse_args()
# input fasta file parsing
InputFile = args.infile
if InputFile is not None:
	print('Input file = ' + InputFile)
else:
	print('ERROR: no input file specified')
# feature file parsing
if args.featurefiles is not None:
	FeatureFiles = args.featurefiles.split(",")
	for i in FeatureFiles:
		print('Feature file = ' + i)

##########################
## FUNCTION DEFINITIONS ##
##########################

## function to take a line from a MethyExtract output file, and return a bed-formatted line with % methylation on both strands
def MEconverter(line):
	Line = line.split()
	Chromosome = Line[0]
	Start = str(int(Line[1])-1)
	End = str(int(Line[1])+1)
	Name = "m"+Line[2]
	# note: "Score" here is the mean of the % methylated reads on both strands
	Score = str((((int(Line[3])/int(Line[4]))*100)+((int(Line[6])/int(Line[7]))*100))/2)
	BedLine = Chromosome + "\t" + Start + "\t" + End + "\t" + Name + "\t" + Score
	return(BedLine)

## function to take a MethylExtract output file and output a bed file
def MEtoBED(MEfile):
	# read input file, keep only lines with methylation called on both strands
	BothStrands = []
	for line in open(MEfile, "r"):
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
	# output bed file to current working dir, but keep file name
	BedOutfileName = MEfile.split("/")[-1].replace(".output",".bed")
	BedOutfile = open(BedOutfileName,"wt")
	BedOutfile.write(BedOutput)
	BedOutfile.close()
	return(BedOutfileName)

## function to take a bed file of methylation levels for cytosines and a gff feature file, and generate mean methylation levels for each feature
def FeatureMeth(MethBed,GFF):
	# remove path from input filenames (to allow adjusted versions to be written to working directory)
	MethBedInput = MethBed.split("/")[-1]
	GFFInput = GFF.split("/")[-1]
	# sort bedfile and feature files
	cmd = "bedtools sort -i " + MethBed + " > " + MethBedInput.replace(".bed",".sorted.bed")
	subprocess.call(cmd, shell=True)
	cmd = "bedtools sort -i " + GFF + " > " + GFFInput.replace(".gff",".sorted.gff")
	subprocess.call(cmd, shell=True)
	# generate mean methylation levels for each feature, and print the levels for each cytosine used to calculate mean
	cmd = "bedtools map -a " + GFFInput.replace(".gff",".sorted.gff") + " -b " + MethBedInput.replace(".bed",".sorted.bed") + " -c 5 -o collapse,mean > " + GFFInput.replace(".gff",".CG.bed")
	subprocess.call(cmd, shell=True)
	# remove intermediate files
	os.remove(MethBedInput.replace(".bed",".sorted.bed"))
	os.remove(GFFInput.replace(".gff",".sorted.gff"))
	print("Methylation levels for " + GFFInput + " written to " + GFFInput.replace(".gff",".CG.bed"))
	return()

##################
## ACTUAL CALLS ##
##################

# convert MethylExtract file to BED file
InputBed = MEtoBED(MEfile=InputFile)

# if feature files have been specified, calculate mean methylation level for each feature in each feature file, and delete whole-genome bed file
if FeatureFiles is not None:
	for i in FeatureFiles:
		FeatureMeth(MethBed=InputBed,GFF=i)
	os.remove(InputBed)

