#!/bin/bash
#SBATCH -p fast               # partition all / fast
#SBATCH -n 1                      # number of cores
#SBATCH -c 32                      # number of cores
#SBATCH --mem-per-cpu=200
#SBATCH -o slurm.%N.%j.out        # STDOUT
#SBATCH -e slurm.%N.%j.err        # STDERR

# =============================================================================
# Project         : 1000 ONT Genome assemblies
# title           : 2_GeneClustering.sh
# description     : This script uses data generated with 
#					1_GeneratePrimaryData.sh to	run gene clustering on samples 
# 					with HQ genome assemblies. 
# author          : vloegler
# date            : 2023/04/27
# version         : 3.0
# usage           : sbatch 2_GeneClustering.sh
# =============================================================================

module load ncbi-blast+

DATADIR=/shared/home/vloegler/Pangenome_Sace1000ONT/03-data
WORKDIR=/shared/home/vloegler/Pangenome_Sace1000ONT/04-analysis
SCRIPTDIR=/shared/home/vloegler/Pangenome_Sace1000ONT/02-scripts

rm -rf $WORKDIR/07_GeneClustering_HQStrains

python $SCRIPTDIR/PythonScripts/GeneClusteringParallel.py \
	--allGenesFasta $WORKDIR/05_BlastNonOrthoAllVSAll_HQStrains/AllGenes_Redundant.min0.1kb.fasta \
	--referenceGenesFasta $WORKDIR/05_BlastNonOrthoAllVSAll_HQStrains/orf_coding_all_R64-4-1_20230830.min0.1kb.fasta \
	--referenceProteinsFasta $DATADIR/orf_trans_all_R64-4-1_20230830.fasta \
	--annotationDir $WORKDIR/01_TransferAnnotation \
	--orthoGenesDir $WORKDIR/02_RefGenes \
	--blastAllVSAll $WORKDIR/05_BlastNonOrthoAllVSAll_HQStrains \
	--output $WORKDIR/07_GeneClustering_HQStrains \
	--nbThreads 128

python $SCRIPTDIR/PythonScripts/RemoveRedundancy.py \
	-pn $WORKDIR/07_GeneClustering_HQStrains/Pangenome.cds.fna \
	-pp $WORKDIR/07_GeneClustering_HQStrains/Pangenome.pep.faa
