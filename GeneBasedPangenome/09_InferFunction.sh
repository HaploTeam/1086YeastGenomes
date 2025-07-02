#!/bin/bash
#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory pool for all cores (see also --mem-per-cpu) 
#SBATCH -o slurm.%N.%j.out      # File to which STDOUT will be written
#SBATCH -e slurm.%N.%j.err      # File to which STDERR will be written

SCRIPTDIR=/shared/home/vloegler/Pangenome_Sace1000ONT/02-scripts
DATADIR=/shared/home/vloegler/Pangenome_Sace1000ONT/03-data
WORKDIR=/shared/home/vloegler/Pangenome_Sace1000ONT/04-analysis

rm -rf $WORKDIR/10_FunctionInference
mkdir $WORKDIR/10_FunctionInference
cd $WORKDIR/10_FunctionInference

python $SCRIPTDIR/PythonScripts/InferGeneFunction.py -f $WORKDIR/07_GeneClustering_HQStrains/Pangenome.NoRedundancy.pep.faa -b $WORKDIR/08_BlastOnRefSeq/BlastResults/PangenomeNonRef_on_RefSeq_Shen2018_Yue2017.blastp -o Pangenome.FunctionTransfer.tsv -d $SCRIPTDIR/datasets -g $DATADIR/gene2go.gz
