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
SCAFFDIR=$BATCHDIR/14_Scaffolded
POLISHDIR=$BATCHDIR/12_Polishing
WORKDIR=$BATCHDIR/15_FINAL

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


for S in $(cat $BATCHDIR/Batch1.txt $BATCHDIR/Batch2.txt $BATCHDIR/Batch3.txt | sort -u)
do
	# check if the scaffolding step merged chromosomes
	DRAFT=$(ls $SCAFFDIR/${S}_*.fasta | grep -E "${S}_(TotalReads|Cov30X|Cov40X)")
	echo "$CHECKCENTRO -d $DRAFT -c $CENTRO -o $WORKDIR/$S.checkCentro.out" >> ${CMD}1
done

ccc_mprun glost_launch ${CMD}1

for S in $(cat $BATCHDIR/Batch1.txt $BATCHDIR/Batch2.txt $BATCHDIR/Batch3.txt | sort -u)
do
	# Get number of merged chromosomes to choose between polished or scaffolded assembly
	nbMerged=$(tail -1 $WORKDIR/$S.checkCentro.out | cut -f 2)
	if [ $nbMerged -eq 0 ]; then
		DRAFT=$(ls $SCAFFDIR/${S}_*.fasta | grep -E "${S}_(TotalReads|Cov30X|Cov40X)")
	else
		DRAFT=$POLISHDIR/$S.*.hapoGPolished1.fasta
	fi

	# Get the standardized strain name
	illuS=$(awk -v S=$S '{if ($1 == S) print $4}' $NAME_EQ)

	OUT=$WORKDIR/$illuS.Renamed.fasta
	echo "$RENAME -r $REFMASKED -d $DRAFT -o $OUT -P $illuS -i -t 20" >> ${CMD}2
	OUT2=$WORKDIR/$illuS.Final.fasta
	echo "$REORDER -r $REF -d $OUT -o $OUT2 -t 20" >> ${CMD}3
done
ccc_mprun glost_launch ${CMD}2
ccc_mprun glost_launch ${CMD}3
#rm -f $CMD*

rm -f $WORKDIR/*.Renamed.fasta
rm -f $WORKDIR/*.checkCentro.out
