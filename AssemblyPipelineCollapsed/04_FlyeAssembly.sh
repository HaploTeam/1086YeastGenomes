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
READSDIR=/ccc/scratch/cont007/fg0006/loeglerv/1_ASSEMBLY_Collapsed/01_SeqData
BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/1_ASSEMBLY_Collapsed
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/1_ASSEMBLY_Collapsed/03_FlyeAssemblies


# Assemble reads with 40X Downsampling
for S in $(cat $BATCHDIR/Batch2_Down40X.txt)
do
	READS=$READSDIR/$S.Cov40X.SizePrior.minLength1kb.minQual9.fastq.gz
	PREFIX=$S.Cov40X.FlyeAssembly

	# Output directory
	OUT=$WORKDIR/$PREFIX

	mkdir $OUT
	echo "$FLYE --nano-raw $READS --out-dir $OUT --threads 4 --iterations 3" >> $COMMANDS
done


# Assemble reads with 30X Downsampling & RawData
for S in $(cat $BATCHDIR/Batch2.txt)
do
	for COV in Cov30X TotalReads
	do
		if [ $COV = "Cov30X" ]
		then
			READS=$READSDIR/$S.Cov30X.SizePrior.minLength1kb.minQual9.fastq.gz
		else
			READS=$READSDIR/$S.NanoporeReads.fastq.gz
		fi

		PREFIX=$S.$COV.FlyeAssembly

		# Output directory
		OUT=$WORKDIR/$PREFIX

		mkdir $OUT
		echo "$FLYE --nano-raw $READS --out-dir $OUT --threads 4 --iterations 3" >> $COMMANDS
	done
done




echo "--- RUNNING COMMANDS ---"
ccc_mprun glost_launch $COMMANDS

rm -f $COMMANDS

echo $(date)


