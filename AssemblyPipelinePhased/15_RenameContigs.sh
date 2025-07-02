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

BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/2_ASSEMBLY_Phased
DRAFTDIR=$BATCHDIR/15_BestAssemblies
WORKDIR=$BATCHDIR/16_RenamedAssemblies

CMD="RenameContigs.cmd"
rm -f $CMD

RENAME="python /ccc/work/cont007/fg0006/loeglerv/GenomeAssemblyTools/renameContigs.py"
REF=/ccc/work/cont007/fg0006/loeglerv/Sace_reference/Sace_S288c_reference_FullMatrixID.fna.masked.NoTelo

for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	for HP in HP1 HP2
	do
		DRAFT=$DRAFTDIR/$S.$HP.*.fasta
		OUT=$WORKDIR/$(basename $DRAFT .fasta).Renamed.fasta
		echo "$RENAME -r $REF -d $DRAFT -o $OUT -P ${S}_${HP} -i -t 20" >> $CMD
	done
done
ccc_mprun glost_launch $CMD
rm -f $CMD

for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	cat $WORKDIR/$S.* > $WORKDIR/$S.PhasedGenome.fasta
done

