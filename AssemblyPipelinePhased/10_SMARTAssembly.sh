#!/bin/bash
#MSUB -r SMARTAssembler
#MSUB -n 100
#MSUB -c 5
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

# ================================================================
# This script runs the Sdn assembler on Nanopore raw sequence data
# ================================================================

# Software path
SMARTDENOVO=/ccc/work/cont007/fg0006/loeglerv/Soft/smartdenovo-master/smartdenovo.pl

# Run commands file
COMMANDS=commandsSMART.cmd
rm -f $COMMANDS*

# Path to the raw data, where there is 1 fasta.gz file per strain read set
BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/2_ASSEMBLY_Phased
READSDIR=$BATCHDIR/07_DownsampleLR
NECATDIR=$BATCHDIR/08_NecatAssembly
WORKDIR=$BATCHDIR/10_SMARTAssembly


# Assemble reads with 40X Downsampling
for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	for HP in HP1 HP2
	do
		if [ -f $READSDIR/$S.$HP.Cov40X.SizePrior.minLength1kb.minQual9.fastq.gz ]
		then
			PREFIX=$S.$HP.Cov40X
			# NECAT-corrected and trimmed Nanopore reads
			FASTAGZ=$NECATDIR/$PREFIX.NecatAssembly/trimReads.fasta.gz

			# Output directory
			OUT=$WORKDIR/$PREFIX.SMARTAssembly
			mkdir $OUT

			NEWFASTAGZ=$OUT/$(basename $FASTAGZ)
			#cp $FASTAGZ $NEWFASTAGZ
			#gzip -d $NEWFASTAGZ
			echo "cp $FASTAGZ $NEWFASTAGZ && gzip -d $NEWFASTAGZ" >> $COMMANDS"1"

			FASTA=$OUT/$(basename $FASTAGZ .gz)

			# Prefix of the output files
			PREFIX=$OUT/$PREFIX".SMARTdenovoAssembly"

			echo "$SMARTDENOVO -p $PREFIX -t 10 -c 1 $FASTA > $PREFIX.mak && make -f $PREFIX.mak" >> $COMMANDS"2"
		fi
		for COV in TotalReads Cov30X
		do
			PREFIX=$S.$HP.$COV
			# NECAT-corrected and trimmed Nanopore reads
			FASTAGZ=$NECATDIR/$PREFIX.NecatAssembly/trimReads.fasta.gz

			# Output directory
			OUT=$WORKDIR/$PREFIX.SMARTAssembly
			mkdir $OUT

			NEWFASTAGZ=$OUT/$(basename $FASTAGZ)
			#cp $FASTAGZ $NEWFASTAGZ
			#gzip -d $NEWFASTAGZ
			echo "cp $FASTAGZ $NEWFASTAGZ && gzip -d $NEWFASTAGZ" >> $COMMANDS"1"

			FASTA=$OUT/$(basename $FASTAGZ .gz)

			# Prefix of the output files
			PREFIX=$OUT/$PREFIX".SMARTdenovoAssembly"

			echo "$SMARTDENOVO -p $PREFIX -t 10 -c 1 $FASTA > $PREFIX.mak && make -f $PREFIX.mak" >> $COMMANDS"2"
		done
	done
done

echo "--- RUNNING COMMANDS ---"
ccc_mprun glost_launch $COMMANDS"1"
ccc_mprun glost_launch $COMMANDS"2"

rm -f $COMMANDS*

echo $(date)


