#!/bin/bash
#MSUB -r RenameContigs
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

BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/1_ASSEMBLY_Collapsed
DRAFTDIR=$BATCHDIR/09_BestAssemblies
WORKDIR=$BATCHDIR/10_RenamedAssemblies

CMD="RenameContigs.cmd"
rm -f $CMD

RENAME="python /ccc/work/cont007/fg0006/loeglerv/GenomeAssemblyTools/renameContigs.py"
REF=/ccc/work/cont007/fg0006/loeglerv/Sace_reference/Sace_S288c_reference_FullMatrixID.fna.masked.NoTelo

for S in $(cat $BATCHDIR/Batch2.txt)
do
	DRAFT=$DRAFTDIR/$S.*.fasta
	OUT=$WORKDIR/$(basename $DRAFT .fasta).Renamed.fasta
	echo "$RENAME -r $REF -d $DRAFT -o $OUT -i -t 20" >> $CMD
done
ccc_mprun glost_launch $CMD
rm -f $CMD

