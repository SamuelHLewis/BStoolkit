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
parser.add_argument('-f', '--featurefiles', type=str, help='list of feature files in gff format (comma-separated)')
args = parser.parse_args()
# input fasta file parsing
InputFile = args.infile
if InputFile is not None:
	print('Input file is ' + InputFile)
else:
	print('ERROR: no input file specified')
# feature file parsing
FeatureFiles = args.featurefiles.split(",")
if len(FeatureFiles) > 0:
	for i in FeatureFiles:
		print('Feature file = ' + i)

##########################
## FUNCTION DEFINITIONS ##
##########################

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
	BedOutfile = open(MEfile.replace(".output",".bed"),"wt")
	BedOutfile.write(BedOutput)
	BedOutfile.close()
	print("Bed file written to "+MEfile.replace(".output",".bed"))
	# calculate genome-wide methylation level
	BedFile = MEfile.replace(".output",".bed")
	Ccount = 0
	MethylationTotal = 0
	for i in BedFormatLines:
		Ccount += 1
		MethylationTotal += float(i.split("\t")[4])
	MethylationLevel = MethylationTotal/Ccount
	print("Genome-wide methylation level for " + BedFile + "= " + str(MethylationLevel) + "%")
	GenomeMethylationOutput = "Genome-wide methylation for " + BedFile + "\n" + str(MethylationLevel)
	GenomeMethylationOutfile = open(BedFile.replace(".bed",".MethSummary"),"wt")
	GenomeMethylationOutfile.write(GenomeMethylationOutput)
	GenomeMethylationOutfile.close()	
	return(MEfile.replace(".output",".bed"))

## function to take a bed file of methylation levels for cytosines and a gff feature file, and generate mean methylation levels for each feature
def FeatureMeth(MethBed,GFF):
	# sort bedfile and feature files
	cmd = "bedtools sort -i " + MethBed + " > " + MethBed.replace(".bed",".sorted.bed")
	subprocess.call(cmd, shell=True)
	cmd = "bedtools sort -i " + GFF + " > " + GFF.replace(".gff",".sorted.gff")
	subprocess.call(cmd, shell=True)
	# generate methylation levels for each feature
	cmd = "bedtools map -a " + GFF.replace(".gff",".sorted.gff") + " -b " + MethBed.replace(".bed",".sorted.bed") + " -c 5 -o mean > " + GFF.replace(".gff",".CG.bed")
	subprocess.call(cmd, shell=True)
	return()

##################
## ACTUAL CALLS ##
##################

# convert MethylExtract file to BED file
InputBed = MEtoBED(MEfile=InputFile)

# calculate mean methylation level for each feature in each feature file
for i in FeatureFiles:
	FeatureMeth(MethBed=InputBed,GFF=i)
