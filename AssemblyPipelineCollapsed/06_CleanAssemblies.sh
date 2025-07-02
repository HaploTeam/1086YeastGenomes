#!/bin/bash
#MSUB -r CleanAssemblies
#MSUB -n 128
#MSUB -c 1
#MSUB -Q long
#MSUB -T 259200
#MSUB -m scratch,work
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

echo $(date)

module load extenv/fg
module load python
module load blast+
module load glost

# ===========================================
# This script runs removeTranslocHaplotigs.py
# ===========================================


# Run commands file
COMMANDS=commandsCLEAN.cmd
rm -f $COMMANDS

# Path to the raw data, where there is 1 fasta.gz file per strain read set
BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/1_ASSEMBLY_Collapsed
WORKDIR=$BATCHDIR/06_CleanedAssemblies
DRAFTDIR=$BATCHDIR/05_RawAssemblies

CLEAN=/ccc/work/cont007/fg0006/loeglerv/GenomeAssemblyTools/removeRedundantContigs.py

for S in $(cat $BATCHDIR/Batch2.txt)
do
	for A in $(ls $DRAFTDIR/$S*)
	do
		O=$WORKDIR/$(basename $A .fasta)
		echo "python $CLEAN -d $A -o $O" >> $COMMANDS
	done
done


echo "--- RUNNING COMMANDS ---"
ccc_mprun glost_launch $COMMANDS

rm -f $COMMANDS

echo $(date)


