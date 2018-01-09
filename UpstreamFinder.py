#!/usr/bin/env python3

import os
import sys
import subprocess
import shutil
import re
import argparse
from collections import OrderedDict

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
	UpstreamLength = int(args.length)
	if UpstreamLength>0:
		print("Upstream length = "+str(UpstreamLength))
	else:
		print("ERROR: upstream length (-l) must be >0")
else:
	print("ERROR: upstream length (-l) not specified as an integer")
	sys.exit(0)

#################
## ACTUAL CODE ##
#################
# function to find start position of 1st exon for each gene - takes a GFF file, returns a dict of gene:start position pairs
def Exon1Finder(GFF):
	Genes = OrderedDict()
	for line in open(GFF,"r"):
		if not line.startswith("#"):
			temp = line.split("\t")
			if temp[2] == "exon":
				Strand = temp[6]
				Start = temp[3]
				Name = temp[8].split(";")[1].replace("Parent=","")
				if Name not in Genes:
					Genes[Name] = Start
	return(Genes)

# final call
result=Exon1Finder(GFF=InputFile)
for i in result:
	print(i + "\t" + result[i])

