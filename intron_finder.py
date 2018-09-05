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
parser.add_argument('-l', '--label', type=str, help='Label of field in GFF file that contains the gene name')
args = parser.parse_args()
# input annotation file parsing
InputAnnotation = args.annotation
if InputAnnotation is not None:
	print("Input annotation file is " + InputAnnotation)
else:
	print("ERROR: no input annotation file (-a) specified")
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
def exon_finder(GFF, Label):
	Genes = []
	Chromosomes = {}
	ExonStarts = {}
	ExonEnds = {}
	Strands = {}
	# go through input, recording gene names in original order and the Chromosome, Strand and Exon Starts for each gene
	for line in open(GFF,"r"):
		if not line.startswith("#"):
			temp = line.split("\t")
			if temp[2] == "CDS":
				# extract data from GFF fields
				Chromosome = temp[0]
				Strand = temp[6]	
				Start = temp[3]
				End = temp[4]
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
				# add Exon end to dict, either as a new list or appended to the existing list
				if Name not in ExonEnds:
					ExonEnds[Name] = [int(End)]
				else:
					ExonEnds[Name].append(int(End))
				# add strand to dict if it isn't there already
				if Name not in Strands:
					Strands[Name] = Strand
	# create pairs of coordinates for start-end of each exon
	for i in ExonStarts:
		ExonStarts[i].sort()
	for i in ExonEnds:
		ExonEnds[i].sort()
	ExonStartEndPairs = {}
	for i in ExonStarts:
		StartEndPairs = []
		for j in range(len(ExonStarts[i])):
			pair = str(ExonStarts[i][j]) + "-" + str(ExonEnds[i][j])
			StartEndPairs.append(pair)
		ExonStartEndPairs[i] = StartEndPairs	
	# screen out duplicate pairs of coordinates
	ExonStartEndPairsNR = {}
	for i in ExonStartEndPairs:
		ExonStartEndPairsNR[i] = list(set(ExonStartEndPairs[i]))
	# convert these non-redundant pairs into two dicts (Start and End)
	ExonStartsNR = {}
	ExonEndsNR = {}
	for i in ExonStartEndPairsNR:
		starts = []
		ends = []
		for j in ExonStartEndPairsNR[i]:
			starts.append(int(j.split("-")[0]))
			ends.append(int(j.split("-")[1]))
		ExonStartsNR[i] = sorted(starts)
		ExonEndsNR[i] = sorted(ends)

	# find introns
	IntronStarts = {}
	IntronEnds = {}
	# go through each gene, calculating the start and end of each intron based on the starts and ends of each exon
	for gene in ExonStartsNR:
		# go through each exon (apart from the last), calculating the start and end of the corresponding intron
		# defining the range here as 0 to (1 - number of exons) to stop after the penultimate exon
		for exon in range(0,len(ExonStartsNR[gene])-1):
			if gene not in IntronStarts:
				# intron starts 1 base downstream of end of upstream exon
				IntronStarts[gene] = [ExonEndsNR[gene][exon]+1]
				# intron ends 1 base upstream of start of downstream exon
				IntronEnds[gene] = [ExonStartsNR[gene][exon+1]-1]
			else:
				# intron starts 1 base downstream of end of upstream exon
				IntronStarts[gene].append(ExonEndsNR[gene][exon]+1)
				# intron starts 1 base upstream of start of downstream exon
				IntronEnds[gene].append(ExonStartsNR[gene][exon+1]-1)
	# find genes that have introns
	GenesWithIntrons = []
	for gene in Genes:
		if gene in IntronStarts:
			GenesWithIntrons.append(gene)
	# make lists of contents of each dict, to arrange their contents in the same gene order as the input file
	ChromosomesOrdered = []
	IntronStartsOrdered = []
	IntronEndsOrdered = []
	StrandsOrdered = []
	for gene in GenesWithIntrons:
		ChromosomesOrdered.append(Chromosomes[gene])
		IntronStartsOrdered.append(IntronStarts[gene])
		IntronEndsOrdered.append(IntronEnds[gene])
		StrandsOrdered.append(Strands[gene])
	# format output
	intron_gff = ""
	for gene in range(len(GenesWithIntrons)):
		for intron in range(len(IntronStartsOrdered[gene])):
			intron_gff += ChromosomesOrdered[gene] + "\tintron_finder\tIntron\t" + str(IntronStartsOrdered[gene][intron]) + "\t" + str(IntronEndsOrdered[gene][intron]) + "\t.\t" + StrandsOrdered[gene] + "\t.\tID=" + GenesWithIntrons[gene] + "\n"
	return(intron_gff)

# final call
intron_gff_string = exon_finder(GFF = InputAnnotation, Label=GeneLabel)
outfile = open("Introns.gff","wt")
outfile.write(intron_gff_string)
outfile.close()
print("Intron annotations written to Introns.gff")
