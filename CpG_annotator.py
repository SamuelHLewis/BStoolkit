#!/usr/bin/env python3

import pandas as pd
import numpy as np
from os import listdir, remove
from os.path import isfile, join
import argparse
import sys

## USER ARGUMENT PARSING ##
# take input from command line
parser = argparse.ArgumentParser(description='Read arguments')
parser.add_argument('-c', '--CG', type=str, help='CG.output file')
parser.add_argument('-g', '--GFF', type=str, help='GFF file')
args = parser.parse_args()
# check that a CG.output file has been given
if args.CG:
	input_CG = args.CG
else:
	print("ERROR: no CG.output file (--CG) specified")
	sys.exit(0)
# assign the GFF filename
if args.GFF:
	input_GFF = args.GFF
else:
	input_GFF = "No GFF"

def CpG_GFF_annotator(CG, GFF = "No GFF", flank_length = 1000):
	"""
	A function that takes a CG and GFF file - with header rows, as formatted by CpG_GFF_splitter.py - and adds annotation data to each CG entry
	"""
	## CpG PARSING ##
	# read in CpG data
	CpG = pd.read_csv(CG, sep = "\t")
	# add extra columns filled with "Unannotated" (these will only be changed if the CpG site falls into a feature)
	CpG["Feature_Name"] = ["Unannotated" for i in range(len(CpG))]
	CpG["Feature_Type"] = ["Unannotated" for i in range(len(CpG))]
	CpG["Feature_Size"] = ["Unannotated" for i in range(len(CpG))]
	CpG["Feature_Total_Coding_Exons"] = ["Unannotated" for i in range(len(CpG))]
	CpG["Feature_Number"] = ["Unannotated" for i in range(len(CpG))]
	CpG["Feature_Start_Distance"] = ["Unannotated" for i in range(len(CpG))]
	CpG["Gene_Start_Distance"] = ["Unannotated" for i in range(len(CpG))]
	
	# test whether a GFF file exists for this chromosome; if not, all CpGs will be left unannotated
	if GFF != "No GFF":
		## GFF IMPORT AND PREPROCESSING ##
		# read in GFF
		features = pd.read_csv(GFF, sep = "\t")
#		print("After GFF reading:")
#		print(features)
		# make a new dataframe holding the coding span for each gene
		coding_spans = pd.DataFrame(features[features["Feature"] == "CDS"].groupby("Metadata")["Start"].min()).join(pd.DataFrame(features[features["Feature"] == "CDS"].groupby("Metadata")["End"].max())).join(pd.DataFrame(features[features["Feature"] == "CDS"].groupby("Metadata")["Strand"].first()))
#		print("Coding spans:")
#		print(coding_spans.head())

		## CODING EXON COUNTING ##
		# store the total number of coding exons per gene in a dataframe (for use later)
		coding_exon_counts = features[features["Feature"] == "CDS"].groupby("Metadata")["Feature"].count()
#		print(coding_exon_counts)

		## FEATURE NUMBERING ##
		# add a Feature_Number column to the feature dataframe with each entry set to 1 by default
		features["Feature_Number"] = 1
		# for each gene
		for feature in features.itertuples():
			if feature.Feature == "gene":
				# remove the CDS regions for this gene from the features dataframe (temporarily)
				CDS_for_gene = pd.DataFrame(features[(features["Metadata"] == feature.Metadata) & (features["Feature"] == "CDS")])
				features.drop(CDS_for_gene.index.tolist(), inplace = True)
				# add incremental feature numbers to CDS, reversing these numbers if the gene is - strand
				if feature.Strand == "+":
					CDS_for_gene["Feature_Number"] = list(range(1,CDS_for_gene.shape[0]+1))
				elif feature.Strand == "-":
					CDS_for_gene["Feature_Number"] = list(range(CDS_for_gene.shape[0], 0, -1))
				# add the numbered CDS annotations back into the features dataframe
				features = features.append(CDS_for_gene)
				# remove the introns (between CDS) for this gene from the features dataframe (temporarily)
				introns_for_gene = pd.DataFrame(features[(features["Metadata"] == feature.Metadata) & (features["Feature"] == "intron")])
				features.drop(introns_for_gene.index.tolist(), inplace = True)
				# add incremental numbers to introns, reversing these numbers if the gene is - strand
				if feature.Strand == "+":
					introns_for_gene["Feature_Number"] = list(range(1,introns_for_gene.shape[0]+1))
				elif feature.Strand == "-":
					introns_for_gene["Feature_Number"] = list(range(introns_for_gene.shape[0], 0, -1))
				# add the numbered intron annotations back into the features dataframe
				features = features.append(introns_for_gene)
		# sort features dataframe by index to put CDS and introns back with their other annotations
		features.sort_index(inplace = True)

		## FEATURE LENGTH COMPUTATION ##
		# add a Feature_Length column based on the start and end of each feature
		features["Feature_Length"] = features["End"] - features["Start"] + 1
		# for each gene, set the length of each 5'UTR as the total length of all 5'UTRs for this gene, and the length of each 3'UTR as the total length of all 3'UTRs for this gene
		# doing this because we don't care how many UTR exons there are, we just want to treat the UTR as a single (spliced) feature
		for feature in features.itertuples():
			if feature.Feature == "gene":
				# remove 5'UTRs for this gene from the features dataframe (temporarily)
				five_prime_UTR_for_gene = pd.DataFrame(features[(features["Metadata"] == feature.Metadata) & (features["Feature"] == "five_prime_UTR")])
				features.drop(five_prime_UTR_for_gene.index.tolist(), inplace = True)
				# calculate total length, and change each length to this value
				total_five_prime_UTR_length = five_prime_UTR_for_gene["Feature_Length"].sum()
				five_prime_UTR_for_gene["Feature_Length"] = total_five_prime_UTR_length
				# add the 5'UTRs back into the features dataframe
				features = features.append(five_prime_UTR_for_gene)
				# remove 3'UTRs for this gene from the features dataframe (temporarily)
				three_prime_UTR_for_gene = pd.DataFrame(features[(features["Metadata"] == feature.Metadata) & (features["Feature"] == "three_prime_UTR")])
				features.drop(three_prime_UTR_for_gene.index.tolist(), inplace = True)
				# calculate total length, and change each length to this value
				total_three_prime_UTR_length = three_prime_UTR_for_gene["Feature_Length"].sum()
				three_prime_UTR_for_gene["Feature_Length"] = total_three_prime_UTR_length
				# add the 3'UTRs back into the features dataframe
				features = features.append(three_prime_UTR_for_gene)
		# sort features dataframe by index to put UTRs back with their other annotations
		features.sort_index(inplace = True)
#		print("Final processed GFF contents:")
#		print(features)
		
		## CpG ANNOTATION ##
		# compute the total span of each gene
		gene_spans = pd.DataFrame(features[features["Feature"] == "gene"].groupby("Metadata")["Start"].min()).join(pd.DataFrame(features[features["Feature"] == "gene"].groupby("Metadata")["End"].max()))
#		print("Gene spans:")
#		print(gene_spans)
		# compute the left noncoding span for each gene
		# deal with + strand genes first
		left_noncoding_spans = pd.DataFrame(features[(features["Feature"] == "five_prime_UTR") & (features["Strand"] == "+")].groupby("Metadata")["Start"].min()).join(pd.DataFrame(features[(features["Feature"] == "five_prime_UTR") & (features["Strand"] == "+")].groupby("Metadata")["End"].max()))
		# then add the - strand genes
		left_noncoding_spans = left_noncoding_spans.append(pd.DataFrame(features[(features["Feature"] == "three_prime_UTR") & (features["Strand"] == "-")].groupby("Metadata")["Start"].min()).join(pd.DataFrame(features[(features["Feature"] == "three_prime_UTR") & (features["Strand"] == "-")].groupby("Metadata")["End"].max())))
		# compute the right noncoding span for each gene
		# deal with + strand genes first
		right_noncoding_spans = pd.DataFrame(features[(features["Feature"] == "three_prime_UTR") & (features["Strand"] == "+")].groupby("Metadata")["Start"].min()).join(pd.DataFrame(features[(features["Feature"] == "three_prime_UTR") & (features["Strand"] == "+")].groupby("Metadata")["End"].max()))
		# then add the - strand genes
		right_noncoding_spans = right_noncoding_spans.append(pd.DataFrame(features[(features["Feature"] == "five_prime_UTR") & (features["Strand"] == "-")].groupby("Metadata")["Start"].min()).join(pd.DataFrame(features[(features["Feature"] == "five_prime_UTR") & (features["Strand"] == "-")].groupby("Metadata")["End"].max())))
#		print("Left non-coding spans =")
#		print(left_noncoding_spans)
#		print("Right non-coding spans =")
#		print(right_noncoding_spans)
		# for each CpG, test if it falls within each feature on this chromosome (ignoring gene body annotations)
		for site in CpG.itertuples():
			for feature in features[features["Feature"] != "gene"].itertuples():
				if site.Chromosome == feature.Chromosome and site.Position in list(range(int(feature.Start), int(feature.End)+1)):
					# store the feature name, type and number
					temp_name = feature.Metadata
					temp_type = feature.Feature
					temp_number = feature.Feature_Number
					# if the feature is a 5'UTR, we want to calculate the feature distance based on all 5'UTRs (i.e. the non-coding region):
					if feature.Feature == "five_prime_UTR":
						# if the gene is on the + strand
						if feature.Strand == "+":
							# store the feature start as the start of the left noncoding region
							temp_start = left_noncoding_spans.loc[feature.Metadata]["Start"]
						# if the gene is on the - strand
						if feature.Strand == "-":
							# store the feature start as the start of the right noncoding region
							temp_start = right_noncoding_spans.loc[feature.Metadata]["Start"]
					# if the feature is a 3'UTR, we want to calculate the feature distance based on all 3'UTRs (i.e. the non-coding region):
					elif feature.Feature == "three_prime_UTR":
						# if the gene is on the + strand
						if feature.Strand == "+":
							# store the feature start as the start of the right noncoding region
							temp_start = right_noncoding_spans.loc[feature.Metadata]["Start"]
						# if the gene is on the - strand
						if feature.Strand == "-":
							# store the feature start as the start of the left noncoding region
							temp_start = left_noncoding_spans.loc[feature.Metadata]["Start"]
					else:
						# if the feature is on the + strand
						if feature.Strand == "+":
							# store the feature start as the start of this feature
							temp_start = feature.Start
						# if the feature is on the - strand
						elif feature.Strand == "-":
							# store the feature end as the start of this feature
							temp_start = feature.End
					# calculate the distance to the start of the feature
					temp_start_distance = abs(site.Position - temp_start)
					# for features contained in genes, calculate the CDS count and distance to the gene start
					if feature.Feature in ["five_prime_UTR", "intron_five_prime_UTR", "CDS", "intron", "three_prime_UTR", "three_prime_UTR_intron"]:
						temp_coding_exon_count = coding_exon_counts[feature.Metadata]
						if feature.Strand == "+":
							temp_gene_start_distance = abs(site.Position - gene_spans.loc[feature.Metadata]["Start"])
						elif feature.Strand == "-":
							temp_gene_start_distance = abs(site.Position - gene_spans.loc[feature.Metadata]["End"])
					else:
						temp_coding_exon_count = "NA"
						temp_gene_start_distance = "NA"
					# add/append data to the entry for this CpG site
					if CpG.iat[site.Index, 3] == "Unannotated":
						CpG.iat[site.Index, 3] = feature.Metadata
						CpG.iat[site.Index, 4] = feature.Feature
						CpG.iat[site.Index, 5] = str(feature.Feature_Length)
						CpG.iat[site.Index, 6] = str(temp_coding_exon_count)
						CpG.iat[site.Index, 7] = str(feature.Feature_Number)
						CpG.iat[site.Index, 8] = str(temp_start_distance)
						CpG.iat[site.Index, 9] = str(temp_gene_start_distance)
					else:
						CpG.iat[site.Index, 3] = CpG.iat[site.Index, 3] + "|" + feature.Metadata
						CpG.iat[site.Index, 4] = CpG.iat[site.Index, 4] + "|" + feature.Feature
						CpG.iat[site.Index, 5] = CpG.iat[site.Index, 5] + "|" + str(feature.Feature_Length)
						CpG.iat[site.Index, 6] = CpG.iat[site.Index, 6] + "|" + str(temp_coding_exon_count)
						CpG.iat[site.Index, 7] = CpG.iat[site.Index, 7] + "|" + str(feature.Feature_Number)
						CpG.iat[site.Index, 8] = CpG.iat[site.Index, 8] + "|" + str(temp_start_distance)
						CpG.iat[site.Index, 9] = CpG.iat[site.Index, 9] + "|" + str(temp_gene_start_distance)
#		print(CpG)
	## FILE OUTPUT ##
	filename = CG.split(".CG.output")[0] + ".CGanno.txt"
	CpG.to_csv(filename, sep = "\t", index = False, header = False)
	print("Annotated CpG data written to file", filename)
	return()

# for chromosomes that have annotations
if isfile(input_GFF):
	CpG_GFF_annotator(CG = input_CG, GFF = input_GFF)
# for chromosomes with no annotation
else:
	CpG_GFF_annotator(CG = input_CG)


