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
NECATDIR=/ccc/scratch/cont007/fg0006/loeglerv/1_ASSEMBLY_Collapsed/02_NecatAssemblies
BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/1_ASSEMBLY_Collapsed
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/1_ASSEMBLY_Collapsed/04_SMARTAssemblies


# Assemble reads with 40X Downsampling
for S in $(cat $BATCHDIR/Batch2_Down40X.txt)
do
	PREFIX=$S.Cov40X
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


# Assemble reads with 30X Downsampling & RawData
for S in $(cat $BATCHDIR/Batch2.txt)
do
	for COV in Cov30X TotalReads
	do
		PREFIX=$S.$COV
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




echo "--- RUNNING COMMANDS ---"
ccc_mprun glost_launch $COMMANDS"1"
ccc_mprun glost_launch $COMMANDS"2"

rm -f $COMMANDS*

echo $(date)


