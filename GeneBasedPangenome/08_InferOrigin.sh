#!/bin/bash
#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory pool for all cores (see also --mem-per-cpu) 
#SBATCH -o slurm.%N.%j.out      # File to which STDOUT will be written
#SBATCH -e slurm.%N.%j.err      # File to which STDERR will be written

module load ncbi-blast+

SCRIPTDIR=/shared/home/vloegler/Pangenome_Sace1000ONT/02-scripts
DATADIR=/shared/home/vloegler/Pangenome_Sace1000ONT/03-data
WORKDIR=/shared/home/vloegler/Pangenome_Sace1000ONT/04-analysis

rm -rf $WORKDIR/09_OriginInference
mkdir $WORKDIR/09_OriginInference
cd $WORKDIR/09_OriginInference

PANGENOME=$WORKDIR/07_GeneClustering_HQStrains/Pangenome.NoRedundancy.pep.faa
BLASTRESULTS=$WORKDIR/08_BlastOnRefSeq/BlastResults/PangenomeNonRef_on_RefSeq_Shen2018_Yue2017.blastp

python $SCRIPTDIR/PythonScripts/InferGeneOrigin.py -f $PANGENOME -b $BLASTRESULTS -o Pangenome.origin.tsv
