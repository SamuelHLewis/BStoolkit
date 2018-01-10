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
				NameMatch = re.search("Parent\=[^\;]*",line)
				Name = NameMatch.group().replace("Parent=","")
				print("Currently on an exon for gene " + Name + " starting at " + str(Start))
				# check whether an exon for this gene has been seen before
				# if it hasn't been seen before, add all of its information to the lists
				if Name not in Genes:
					Chromosomes.append(Chromosome)
					Genes.append(Name)
					Starts.append(Start)
					Strands.append(Strand)
					print("First site of exon from gene " + Name)
				# if it has been seen before, compare the start position of this exon to the last element in the Starts list (which will have come from another exon of this gene)
				else:
					# if the gene is on the + strand, replace the last element with this start position if the start position is lower (i.e. more 5' on the + strand)
					# NB: this should never happen if the exon entries are in ascending numerical order, but keeping the check in here for now just in case
					if Strand == "+":
						if Start < Starts[-1]:
							print("Have already seen exon from gene " + Name + ", but this exon has a more 5' start site (" + str(Start) + " vs " + str(Starts[-1]) + " on the + strand)")
							del Starts[-1]
							Starts.append(Start)
					# if the gene is on the - strand, replace the last element with this start position if the start position is higher (i.e. more 5' on the - strand)
					if Strand == "-":
						if Start > Starts[-1]:
							print("Have already seen exon from gene " + Name + ", but this exon has a more 5' start site (" + str(Start) + " vs " + str(Starts[-1]) + " on the - strand)")
							del Starts[-1]
							Starts.append(Start)
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
