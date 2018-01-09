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
	Chromosomes = []
	Genes = []
	Starts = []
	Strands = []
	for line in open(GFF,"r"):
		if not line.startswith("#"):
			temp = line.split("\t")
			if temp[2] == "exon":
				Chromosome = temp[0]
				Strand = temp[6]
				if Strand == "+":
					Start = temp[3]
				elif Strand == "-":
					Start = temp[4]
				Name = temp[8].split(";")[1].replace("Parent=","")
				if Name not in Genes:
					Chromosomes.append(Chromosome)
					Genes.append(Name)
					Starts.append(Start)
					Strands.append(Strand)
	# calculate upstream region coordinates for each gene
	UpstreamStarts = []
	for i in range(len(Starts)):
		if Strands[i] == "+":
			UpstreamStarts.append(int(Starts[i])-UpstreamLength-1)
		elif Strands[i] == "-":
			UpstreamStarts.append(int(Starts[i])+1)
	UpstreamEnds = []
	for i in range(len(Starts)):
		if Strands[i] == "+":
			UpstreamEnds.append(int(Starts[i])-1)
		elif Strands[i] == "-":
			UpstreamEnds.append(int(Starts[i])+UpstreamLength+1)
	# clean up all lists, removing entries where the upstream start coordinate is <1
	ToDelete = []
	for i in range(len(UpstreamStarts)):
		if UpstreamStarts[i]<1:
			ToDelete.append(i)
	# NB: deleting the elements in reverse order, to ensure that the element numbering isn't disrupted while we're still relying on it
	for i in reversed(ToDelete):
		del Chromosomes[i]
		del UpstreamStarts[i]
		del UpstreamEnds[i]
		del Strands[i]
		del Genes[i]
	# make list of lists
	Combined = list(zip(Chromosomes,UpstreamStarts,UpstreamEnds,Strands,Genes))
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
