#!/bin/bash
#MSUB -r NecatAssembler
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
# This script runs the NECAT assembler on Nanopore raw sequence data
# ==================================================================

# Software path
NECAT=$CCCWORKDIR/Soft/NECAT/Linux-amd64/bin/necat.pl

# Run commands file
COMMANDS=commandsNECAT.cmd
rm -f $COMMANDS*

# Path to the raw data, where there is 1 fasta.gz file per strain read set
READSDIR=/ccc/scratch/cont007/fg0006/loeglerv/1_ASSEMBLY_Collapsed/01_SeqData
BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/1_ASSEMBLY_Collapsed

# Assemble reads with 40X Downsampling
for S in $(cat $BATCHDIR/Batch2_Down40X.txt)
do
	READS=$READSDIR/$S.Cov40X.SizePrior.minLength1kb.minQual9.fastq.gz
	PREFIX=$S.Cov40X.NecatAssembly
	echo $READS > reads_$PREFIX.txt
	echo "PROJECT=$PREFIX
ONT_READ_LIST=reads_$PREFIX.txt
GENOME_SIZE=12000000
THREADS=4
MIN_READ_LENGTH=1000
PREP_OUTPUT_COVERAGE=40
OVLP_FAST_OPTIONS=-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000
OVLP_SENSITIVE_OPTIONS=-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000
CNS_FAST_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
CNS_SENSITIVE_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
TRIM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400
ASM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400
NUM_ITER=2
CNS_OUTPUT_COVERAGE=30
CLEANUP=1
USE_GRID=false
GRID_NODE=0
GRID_OPTIONS=
SMALL_MEMORY=0
FSA_OL_FILTER_OPTIONS=
FSA_ASSEMBLE_OPTIONS=
FSA_CTG_BRIDGE_OPTIONS=
POLISH_CONTIGS=true" > $PREFIX.config.txt

	echo "$NECAT correct $PREFIX.config.txt" >> $COMMANDS"1"
	echo "$NECAT assemble $PREFIX.config.txt" >> $COMMANDS"2"
	echo "$NECAT bridge $PREFIX.config.txt" >> $COMMANDS"3"
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
	
		PREFIX=$S.$COV.NecatAssembly
		echo $READS > reads_$PREFIX.txt
		echo "PROJECT=$PREFIX
ONT_READ_LIST=reads_$PREFIX.txt
GENOME_SIZE=12000000
THREADS=4
MIN_READ_LENGTH=1000
PREP_OUTPUT_COVERAGE=40
OVLP_FAST_OPTIONS=-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000
OVLP_SENSITIVE_OPTIONS=-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000
CNS_FAST_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
CNS_SENSITIVE_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
TRIM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400
ASM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400
NUM_ITER=2
CNS_OUTPUT_COVERAGE=30
CLEANUP=1
USE_GRID=false
GRID_NODE=0
GRID_OPTIONS=
SMALL_MEMORY=0
FSA_OL_FILTER_OPTIONS=
FSA_ASSEMBLE_OPTIONS=
FSA_CTG_BRIDGE_OPTIONS=
POLISH_CONTIGS=true" > $PREFIX.config.txt

		echo "$NECAT correct $PREFIX.config.txt" >> $COMMANDS"1"
		echo "$NECAT assemble $PREFIX.config.txt" >> $COMMANDS"2"
		echo "$NECAT bridge $PREFIX.config.txt" >> $COMMANDS"3"
	done
done




echo "--- RUNNING COMMANDS ---"
ccc_mprun glost_launch $COMMANDS"1"
ccc_mprun glost_launch $COMMANDS"2"
ccc_mprun glost_launch $COMMANDS"3"

rm -f $COMMANDS*

echo $(date)


