#!/bin/bash
#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=200              # Memory pool for all cores (see also --mem-per-cpu) 
#SBATCH -o slurm.%N.%j.out      # File to which STDOUT will be written
#SBATCH -e slurm.%N.%j.err      # File to which STDERR will be written

SCRIPTDIR=/shared/home/vloegler/Pangenome_Sace1000ONT/02-scripts
DATADIR=/shared/home/vloegler/Pangenome_Sace1000ONT/03-data
WORKDIR=/shared/home/vloegler/Pangenome_Sace1000ONT/04-analysis

rm -rf $WORKDIR/12_FinalAnnotTable
mkdir $WORKDIR/12_FinalAnnotTable
cd $WORKDIR/12_FinalAnnotTable

python $SCRIPTDIR/PythonScripts/MergeTables.py \
	-l $WORKDIR/07_GeneClustering_HQStrains/Pangenome.NonRef.Localization.tsv \
	-p $WORKDIR/07_GeneClustering_HQStrains/Pangenome.NoRedundancy.cds.fna \
	-r $WORKDIR/09_OriginInference/Pangenome.origin.tsv \
	-f $WORKDIR/10_FunctionInference/Pangenome.FunctionTransfer.tsv \
	-i $WORKDIR/11_InterProScan/Pangenome.NoRedundancy.NonRef.pep.faa.InterProScan.tsv \
	-o Pangenome.Annotations.tsv

