#!/bin/bash
#MSUB -r NuclearAssemblies
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


# Run commands file
COMMANDS=commandsNUCLEAR.cmd
rm -f $COMMANDS

BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/1_ASSEMBLY_Collapsed
WORKDIR=$BATCHDIR/07_nuclearAssemblies
DRAFTDIR=$BATCHDIR/06_CleanedAssemblies
NUCLEARDB=/ccc/work/cont007/fg0006/loeglerv/Sace_DataBaseNuclearChr/SaceNuclearChromosomes
GETNUCLEAR="python /ccc/work/cont007/fg0006/loeglerv/GenomeAssemblyTools/getNuclearContigs.py"


for S in $(cat $BATCHDIR/Batch2.txt)
do
	for A in $(ls $DRAFTDIR/$S*.NoRedundantContigs.fasta)
	do
		O=$WORKDIR/$(basename $A .fasta)
		echo "$GETNUCLEAR -d $A -db $NUCLEARDB -o $O" >> $COMMANDS
	done
done


echo "--- RUNNING COMMANDS ---"
ccc_mprun glost_launch $COMMANDS

rm -f $COMMANDS

echo $(date)


