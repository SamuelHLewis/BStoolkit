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
	UpstreamLength = int(args.length)
	if UpstreamLength>0:
		print("Upstream length = "+str(UpstreamLength))
	else:
		print("ERROR: upstream length (-l) must be >0")
else:
	print("ERROR: upstream length (-l) not specified as an integer")
	sys.exit(0)
