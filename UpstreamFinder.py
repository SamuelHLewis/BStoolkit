#!/usr/bin/env python3

import os
import sys
import subprocess
import shutil
import re
import argparse
from Bio import SeqIO

###########################
## USER ARGUMENT PARSING ##
###########################
parser = argparse.ArgumentParser(description='Read arguments')
parser.add_argument('-a', '--annotation', type=str, help='Input CDS annotation file (gff format)')
parser.add_argument('-f', '--fasta', type=str, help='Input genome sequence file (fasta format)')
parser.add_argument('-u', '--upstreamlength', type=str, help='Length of upstream sequence to annotate')
parser.add_argument('-l', '--label', type=str, help='Label of field in GFF file that contains the gene name')
args = parser.parse_args()
# input annotation file parsing
InputAnnotation = args.annotation
if InputAnnotation is not None:
	print("Input annotation file is " + InputAnnotation)
else:
	print("ERROR: no input annotation file (-a) specified")
	sys.exit(0)
# input fasta file parsing
InputFasta = args.fasta
if InputFasta is not None:
	print("Input genome file is " + InputFasta)
else:
	print("ERROR: no input genome file (-f) specified")
	sys.exit(0)
# upstream length parsing
if str.isdigit(args.upstreamlength):
	Length = int(args.upstreamlength)
	if Length>0:
		print("Upstream length = "+str(Length))
	else:
		print("ERROR: upstream length (-u) must be >0")
else:
	print("ERROR: upstream length (-u) not specified as an integer")
	sys.exit(0)
# gene name label parsing
GeneLabel = args.label
if GeneLabel is not None:
	print("Label used for finding gene name in GFF fields = " + GeneLabel)
else:
	print("ERROR: no gene name field label (-l) specified")
	sys.exit(0)

#################
## ACTUAL CODE ##
#################
# function to find start position of 1st exon for each gene - takes a GFF file, returns a dict of gene:start position pairs
def UpstreamFinder(GFF, Fasta, UpstreamLength, Label):
	Genes = []
	Chromosomes = {}
	ExonStarts = {}
	Strands = {}
	# go through input, recording gene names in original order and the Chromosome, Strand and Exon Starts for each gene
	for line in open(GFF,"r"):
		if not line.startswith("#"):
			temp = line.split("\t")
			if temp[2] == "CDS":
				# extract data from GFF fields
				Chromosome = temp[0]
				Strand = temp[6]
				if Strand == "+":
					Start = temp[3]
				elif Strand == "-":
					Start = temp[4]
				NameMatch = re.search(Label+"\=[^\;]*",line)
				Name = NameMatch.group().replace(Label+"=","")
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
			UpstreamStarts[i]=ExonStarts[i][0]-UpstreamLength
		elif Strands[i] == "-":
			UpstreamStarts[i]=ExonStarts[i][-1]+1
	UpstreamEnds = {}
	for i in Genes:
		if Strands[i] == "+":
			UpstreamEnds[i]=ExonStarts[i][0]-1
		elif Strands[i] == "-":
			UpstreamEnds[i]=ExonStarts[i][-1]+UpstreamLength
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
	# go through genome fasta file, recording the length for each chromosome
	ChromosomeLengths = {}
	for record in SeqIO.parse(Fasta, "fasta"):
		ChromosomeLengths[record.id] = len(str(record.seq))
	# clean up all lists	
	ToDelete = []
	for i in range(len(Genes)):
		# remove entries where upstream start coordinate is < 1
		if UpstreamStarts[Genes[i]]<1:
			ToDelete.append(i)
		# remove entries where upstream end coordinate is > chromosome length
		elif UpstreamEnds[Genes[i]]>ChromosomeLengths[ChromosomesList[i]]:
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
result=UpstreamFinder(GFF = InputAnnotation, Fasta = InputFasta, UpstreamLength=Length, Label=GeneLabel)
OutputString = ""
for i in result:
	OutputString += i[0] + "\tUpstreamFinder\tUpstreamRegion\t" + str(i[1]) + "\t" + str(i[2]) + "\t.\t" + i[3] + "\t.\tID=" + i[4].rstrip("\n") + "\n"
outfile = open("Upstream.gff","wt")
outfile.write(OutputString)
outfile.close()
print("Upstream regions written to Upstream.gff")
