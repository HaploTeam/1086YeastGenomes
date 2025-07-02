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
SCAFFDIR=$BATCHDIR/20_Scaffolded
POLISHDIR=$BATCHDIR/19_Scaffolding
WORKDIR=$BATCHDIR/21_FINAL

CMD="RenameContigs.cmd"
rm -f $CMD*

CHECKCENTRO="python /ccc/work/cont007/fg0006/loeglerv/GenomeAssemblyTools/findDoubleCentromeres.py"
RENAME="python /ccc/work/cont007/fg0006/loeglerv/GenomeAssemblyTools/renameContigs.py"
REORDER="python /ccc/work/cont007/fg0006/loeglerv/GenomeAssemblyTools/reorderContigs.py"
CENTRO=/ccc/work/cont007/fg0006/loeglerv/Sace_reference/Sace_S288c_Centromeres.fasta
REF=/ccc/work/cont007/fg0006/loeglerv/Sace_reference/Sace_S288c_reference_FullMatrixID.fna
REFMASKED=/ccc/work/cont007/fg0006/loeglerv/Sace_reference/Sace_S288c_reference_FullMatrixID.fna.masked.NoTelo

# CYH name equivalence
NAME_EQ=/ccc/work/cont007/fg0006/loeglerv/Sace_NameEq/CYH_name_eq_ReseqInfo.tsv

for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	for HP in HP1 HP2
	do
		# check if scaffolding merged chromosomes
		DRAFT=$SCAFFDIR/${S}_${HP}*.fasta
		echo "$CHECKCENTRO -d $DRAFT -c $CENTRO -o $WORKDIR/$S.$HP.checkCentro.out" >> ${CMD}1
	done
done

ccc_mprun glost_launch ${CMD}1

for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	for HP in HP1 HP2
	do
		# Get number of merged chromosomes to choose between polished or scaffolded assembly
		nbMerged=$(tail -1 $WORKDIR/$S.$HP.checkCentro.out | cut -f 2)
		if [ $nbMerged -eq 0 ]; then
			DRAFT=$SCAFFDIR/${S}_${HP}*.fasta
		else
			DRAFT=$POLISHDIR/$S.$HP.*.fasta
		fi

		# Get the standardized strain name
		illuS=$(awk -v S=$S '{if ($1 == S) print $4}' $NAME_EQ)

		OUT=$WORKDIR/$illuS.$HP.Renamed.fasta
		echo "$RENAME -r $REFMASKED -d $DRAFT -o $OUT -P ${illuS}_${HP} -i -t 20" >> ${CMD}2
		OUT2=$WORKDIR/$illuS.$HP.Final.fasta
		echo "$REORDER -r $REF -d $OUT -o $OUT2 -t 20" >> ${CMD}3
	done
done
ccc_mprun glost_launch ${CMD}2
ccc_mprun glost_launch ${CMD}3
rm -f $CMD*

rm -f $WORKDIR/*.checkCentro.out
rm -f $WORKDIR/*.Renamed.fasta
