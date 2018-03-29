#!/usr/bin/env python3

import os
import sys
import subprocess
import shutil
import re
import argparse

###########################
## USER ARGUMENT PARSING ##
###########################
parser = argparse.ArgumentParser(description='Read arguments')
parser.add_argument('-i', '--infile', type=str, help='Input CDS annotation file (gff format)')
parser.add_argument('-l', '--length', type=str, help='Length of upstream sequence to annotate upstream')
args = parser.parse_args()
# input fasta file parsing
InputFile = args.infile
if InputFile is not None:
	print("Input file is " + InputFile)
else:
	print("ERROR: no input file specified")
	sys.exit(0)
# upstream length parsing
if str.isdigit(args.length):
	Length = int(args.length)
	if Length>0:
		print("Upstream length = "+str(Length))
	else:
		print("ERROR: upstream length (-l) must be >0")
else:
	print("ERROR: upstream length (-l) not specified as an integer")
	sys.exit(0)

#################
## ACTUAL CODE ##
#################
# function to find start position of 1st exon for each gene - takes a GFF file, returns a dict of gene:start position pairs
def Exon1Finder(GFF,UpstreamLength):
	Genes = []
	Chromosomes = {}
	ExonStarts = {}
	Strands = {}
	# go through input, recording gene names in original order and the Chromosome, Strand and Exon Starts for each gene
	for line in open(GFF,"r"):
		if not line.startswith("#"):
			temp = line.split("\t")
			if temp[2] == "exon":
				# extract data from GFF fields
				Chromosome = temp[0]
				Strand = temp[6]
				if Strand == "+":
					Start = temp[3]
				elif Strand == "-":
					Start = temp[4]
				NameMatch = re.search("Parent\=[^\;]*",line)
				Name = NameMatch.group().replace("Parent=","")
				# add gene name to list
				if Name not in Genes:
					Genes.append(Name)
				# add Chromosome name to dict if it isn't there already
				if Name not in Chromosomes:
					Chromosomes[Name] = Chromosome
				# add Exon start to dict, either as a new list or appended to the existing list
				if Name not in ExonStarts:
					ExonStarts[Name] = [int(Start)]
				else:
					ExonStarts[Name].append(int(Start))
				# add strand to dict if it isn't there already
				if Name not in Strands:
					Strands[Name] = Strand
	# sort start coordinates for each gene into ascending order
	for i in ExonStarts:
		ExonStarts[i].sort()
	# calculate upstream region coordinates for each gene
	UpstreamStarts = {}
	for i in Genes:
		if Strands[i] == "+":
			UpstreamStarts[i]=ExonStarts[i][0]-UpstreamLength-1
		elif Strands[i] == "-":
			UpstreamStarts[i]=ExonStarts[i][-1]+1
	UpstreamEnds = {}
	for i in Genes:
		if Strands[i] == "+":
			UpstreamEnds[i]=ExonStarts[i][0]-1
		elif Strands[i] == "-":
			UpstreamEnds[i]=ExonStarts[i][-1]+UpstreamLength+1
	# make lists from dicts, preserving original gene input order
	ChromosomesList = []
	StrandsList = []
	UpstreamStartsList = []
	UpstreamEndsList = []
	for i in Genes:
		ChromosomesList.append(Chromosomes[i])
		StrandsList.append(Strands[i])
		UpstreamStartsList.append(UpstreamStarts[i])
		UpstreamEndsList.append(UpstreamEnds[i])
	# clean up all lists, removing entries where the upstream start coordinate is <1
	ToDelete = []
	for i in range(len(Genes)):
		if UpstreamStarts[Genes[i]]<1:
			ToDelete.append(i)
	# NB: deleting elements from lists in reverse order to preserve element numbering until it doesn't matter
	for i in reversed(ToDelete):
		del Genes[i]
		del ChromosomesList[i]
		del UpstreamStartsList[i]
		del UpstreamEndsList[i]
		del StrandsList[i]
	# make list of lists
	Combined = list(zip(ChromosomesList,UpstreamStartsList,UpstreamEndsList,StrandsList,Genes))
	return(Combined)

# final call
result=Exon1Finder(GFF=InputFile,UpstreamLength=Length)
OutputString = ""
for i in result:
	OutputString += i[0] + "\tUpstreamFinder\tUpstreamRegion\t" + str(i[1]) + "\t" + str(i[2]) + "\t.\t" + i[3] + "\t.\tID=" + i[4]  + "\n"
outfile = open("Upstream.gff","wt")
outfile.write(OutputString)
outfile.close()
print("Upstream regions written to Upstream.gff")
