#!/bin/bash
#MSUB -r FlyeAssembler
#MSUB -n 64
#MSUB -c 8
#MSUB -Q long
#MSUB -T 259200
#MSUB -m scratch,work
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

echo $(date)

module load extenv/fg
module load python
module load r
module load glost

# ==================================================================
# This script runs the Flye assembler on Nanopore raw sequence data
# ==================================================================

# Software path
FLYE=/ccc/work/cont007/fg0006/loeglerv/Soft/Flye-2.9/bin/flye

# Run commands file
COMMANDS=commandsFlye.cmd
rm -f $COMMANDS

# Path to the raw data, where there is 1 fasta.gz file per strain read set
BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/2_ASSEMBLY_Phased
READSDIR=$BATCHDIR/07_DownsampleLR
WORKDIR=$BATCHDIR/09_FlyeAssembly


# Assemble reads with 40X Downsampling
for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	for HP in HP1 HP2
	do
		if [ -f $READSDIR/$S.$HP.Cov40X.SizePrior.minLength1kb.minQual9.fastq.gz ]
		then
			READS=$READSDIR/$S.$HP.Cov40X.SizePrior.minLength1kb.minQual9.fastq.gz
			PREFIX=$S.$HP.Cov40X.FlyeAssembly

			# Output directory
			OUT=$WORKDIR/$PREFIX

			mkdir $OUT
			echo "$FLYE --nano-raw $READS --out-dir $OUT --threads 4 --iterations 3" >> $COMMANDS
		fi
		for COV in TotalReads Cov30X
		do
			READS=$READSDIR/$S.$HP.$COV.*fastq.gz
			PREFIX=$S.$HP.$COV.FlyeAssembly

			# Output directory
			OUT=$WORKDIR/$PREFIX

			mkdir $OUT
			echo "$FLYE --nano-raw $READS --out-dir $OUT --threads 4 --iterations 3" >> $COMMANDS
		done
	done
done

echo "--- RUNNING COMMANDS ---"
ccc_mprun glost_launch $COMMANDS

rm -f $COMMANDS

echo $(date)


