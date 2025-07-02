#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/11/04
# version ='3.5'
# ---------------------------------------------------------------------------
'''
This script run Genome Wide Association Study using FaST-LMM verion 0.6.4. It
performs association and permutation test to output significant SNPs. 

It takes as input :
	-g --genotype: the genotype matrix (Plink format)
	-k --kinship: the kinship matrix (Plink format)
	-p --phenotype: path to the phenotypes (*.norm.phen), in format : 
											Strain		CondName
											StrainXXX	0.02
											StrainXXX	0.12
	-c --covariance: covariance matrix (optional)
	-o --output: working and output directory
	-n --nbPermutations: number of permutations for the permutation test

Output will be:
outdir_first_assoc.txt : FaST-LMM association results
outdir_condition_threshold.txt: Value of the threshold
outdir_condition_signif_snps: FaST-LMM results for SNPs above the threshold
'''
# ---------------------------------------------------------------------------
import time
import numpy as np
import pandas as pd
import math
import os
import shutil
from sys import argv
from fastlmm.association import single_snp
from pysnptools.snpreader import Bed, SnpData
from pysnptools.util.mapreduce1.runner import LocalMultiProc
import logging
import argparse
import glob
# ---------------------------------------------------------------------------

start_time=time.time()
### random shuffling of matrix columns
def matrix_permut(x):
	ncol = x.shape[1]
	for i in range(ncol):
		np.random.shuffle(x[:,i])

# =====================
# === GWAS FUNCTION ===
# =====================

def runGWAS(bed_file, kinship_file, pheno_file, covar_file, outdir, nperm, threads):

	# Read phenotypes
	df=pd.read_csv(pheno_file,sep="\t",index_col=0)
	# Get phenotype name
	phen_name = df.columns[0]
	# Get phenotype values
	values = np.array(df.values)
	# Get list of strains
	strainslist=np.c_[np.array(df.index)]
	#strainslist_2el=np.append(strainslist,strainslist,axis=1).astype('<U23')
	strainslist_2el = np.c_[np.array([(x, x) if not '_' in x else (x.split('_')[0], x.split('_')[1]) for x in df.index])] # Split strain name with underscore

	# Create working directory
	if not outdir.endswith("/"): outdir += "/"
	outdir += phen_name
	try:
		os.makedirs(outdir) # create a folder with the outdir for each phen
	except:
		raise Exception("The output directory '"+outdir+"' already exists. Please provide a new output directory. ")
	# Go to work directory
	os.chdir(outdir)

	# If there are permutations to run
	if nperm > 0:
		# Create table with all permutated phenotypes
		phen_topermut=values[:, [0]]
		permut_matrix=np.repeat(phen_topermut,repeats=nperm,axis=1)
		matrix_permut(permut_matrix)

		# Merge the phenotype to the permutated phenotypes
		merge_pheno_value=np.hstack((df,permut_matrix))

		# Merge Condition id
		permHeaders = ["perm"+str(x) for x in range(1, nperm + 1)]
		headers = np.array([phen_name]+permHeaders)

		# Init list that will contain the best PValue of each association
		bestPVal = []
	else:
		merge_pheno_value=df
		headers = np.array([phen_name])

	# If more than 100 permutations is asked, multiple associations must be done by batch of 101 (1 first association + 100 permutations)
	for i in range(int(math.ceil((nperm+1)/101))):
		if (i+1)*101 < nperm+1:
			nAsso2run = 101
		else:
			nAsso2run = nperm+1 - (i * 101)

		if nperm > 0: # If there are permutations, take a subset of maximum 101 phenotypes
			pheno2run = np.c_[merge_pheno_value[:,i*101:i*101+nAsso2run]]
		else: # Else take the whole pheno table
			pheno2run = np.c_[merge_pheno_value[:]]
		# Create phenotype object
		phenoSNPread=SnpData(iid=strainslist_2el, sid=headers[i*101:i*101+nAsso2run], val=pheno2run)

		#Run association
		runner = LocalMultiProc(threads)
		start_time=time.time()
		if covar_file == "" and kinship_file == "": # If no covariance matrix or kinship matrix
			results_df = single_snp(test_snps = bed_file, pheno = phenoSNPread, count_A1=True, runner=runner, show_snp_fract_var_exp = True)
		elif covar_file == "": # If only kinship matrix
			results_df = single_snp(test_snps = bed_file, G0 = kinship_file, pheno = phenoSNPread, count_A1=True, runner=runner, show_snp_fract_var_exp = True)
		elif kinship_file == "": # If only covariance matrix
			results_df = single_snp(test_snps = bed_file, pheno = phenoSNPread, covar = covar_file, count_A1=True, runner=runner, show_snp_fract_var_exp = True)
		else: # If kinship and covariance matrix
			results_df = single_snp(test_snps = bed_file, G0 = kinship_file, pheno = phenoSNPread, covar = covar_file, count_A1=True, runner=runner, show_snp_fract_var_exp = True)
		print("TrackTime GWAS {}: Batch{}, {} phenotypes: Association lasted {} seconds ".format(pheno_file.split("/")[-1], i, nAsso2run, time.time() - start_time))

		# If there are permutations, discriminate permutations from the first association
		if nperm > 0:
			# If first association in this batch
			if phen_name in headers[i*101:i*101+nAsso2run]:
				# Split results in first association and permutations
				first_association = results_df[results_df.Pheno == phen_name]
				permutations = results_df[results_df.Pheno != phen_name]
				# Write first association to file
				first_association.to_csv((outdir + "/" + phen_name + ".first_assoc.txt.gz"),index=False,sep="\t")
			else:
				# Batch only contains permutations
				permutations = results_df

			# Gather the best pValue of each permutation
			# Keep only the best PValue of each permutation (=Pheno)
			# It is the first row since ordered by increasing PValue, and drop_duplicates keeps by default the first row
			FirstPValues = permutations.drop_duplicates(subset=['Pheno'])
			FirstPValues = FirstPValues[FirstPValues.Pheno != phen_name]
			bestPVal += list(FirstPValues.PValue)
			del permutations

		else: # If no permutations
			first_association = results_df
			# Write first association to file
			first_association.to_csv((outdir + "/" + phen_name + ".first_assoc.txt.gz"),index=False,sep="\t")

		del results_df

	# If there are permutations, compute pValue threshold, and write threshold and Signif SNPs files
	if nperm > 0:
		# Get PValue threshold
		threshold=np.percentile(bestPVal, 5)

		# Write threshold to file
		thresh = open(outdir + "/" + phen_name +".threshold.txt", "w")
		thresh.write("x\n"+str(threshold)+"\n")
		thresh.close()

		# Get signif SNPs above threshold and write to file
		signif_snps=first_association.loc[first_association.loc[:]['PValue'] < threshold][:]
		signif_snps.to_csv((outdir + "/" + phen_name +".signif_snps.txt"),index=False,sep="\t")

	del first_association

	# Remove temporary files (if some are created)
	try:
		shutil.rmtree(outdir+"/runs")
	except:
		pass
	for F in glob.glob(outdir+"/.work*"):
		shutil.rmtree(F)

# =================
# === MAIN CODE ===
# =================
# FaST-LMM parmeters
logging.basicConfig(level=logging.ERROR)

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genotype", help="genotype matrix in Plink format", required=True)
parser.add_argument("-k", "--kinship", help="kinship matrix in Plink format", type = str, default = "")
parser.add_argument("-p", "--phenotypes", help="phenotype data", required=True, nargs='+')
parser.add_argument("-o", "--outdir", help="Output directory", required=True)
parser.add_argument("-c", "--covariance", help="covariance matrix (optional)", type=str, default="")
parser.add_argument("-n", "--nbPermutations", help="number of permutations for the permutation test (default = 100)", type=int, default=100)
parser.add_argument("-t", "--threads", help="Number of threads to use", type=int, default=1)

# Read arguments from the command line
args = parser.parse_args()

# variant matrix in PLINK format
bed_file = args.genotype
# variant matrix in PLINK format
kinship_file = args.kinship
if kinship_file != "" and not kinship_file.endswith(".bed"):
	kinship_file += ".bed" # FaST-LMM expects a .bed file as kinship, not just the plink prefix
# Phenotypes in format, multiple files allowed
# strain	Cond1	Cond2	Cond3
# StrainXXX	0.02	0.98	0.14
# StrainXXX	0.12	0.52	0.65
pheno_files = args.phenotypes
# Covariance matrix in format
# Strain	Strain	Covariant1	Covariant2 ...
covar_file = args.covariance
# Output directory
outdir = args.outdir
# Number of permutations
nperm = args.nbPermutations
# Number of threads
threads = args.threads

# Run function
for pheno in pheno_files:
	start_time=time.time()
	runGWAS(bed_file, kinship_file, pheno, covar_file, outdir, nperm, threads)
	print("TrackTime GWAS {}: Total time: Lasted {} seconds ".format(pheno.split("/")[-1], time.time() - start_time))




