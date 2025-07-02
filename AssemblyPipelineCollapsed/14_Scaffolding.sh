#!/bin/bash
#MSUB -r RagoutScaffolding
#MSUB -n 5
#MSUB -c 1
#MSUB -Q long
#MSUB -T 259200
#MSUB -m scratch,work
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

echo $(date)

module load ragout
module load glost

BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/1_ASSEMBLY_Collapsed
DRAFTDIR=$BATCHDIR/12_Polishing
WORKDIR=$BATCHDIR/13_Scaffolding
REF=/ccc/work/cont007/fg0006/loeglerv/Sace_reference/Sace_S288c_reference_FullMatrixID.fna

CMD="Scaffolding.cmd"
rm -f $CMD*

for S in $(cat $BATCHDIR/Batch2.txt)
do
	DRAFT=$DRAFTDIR/$S.*.hapoGPolished1.fasta
	# Add Draft_ beforedraft sequence names because it cannot be the same ID than the ref chr
	sed -i 's/>/>Draft_/g' $DRAFT
	# Convert . to _ in the prefix
	PREFIX=$(basename $DRAFT .fasta | tr "." "_")

	O=$WORKDIR/$PREFIX
	mkdir $O

	echo "" > $O/recipeFile.txt

	# Create recipe file
	#reference and target genome names (required)
	echo ".references = Sace_S288c" >> $O/recipeFile.txt
	echo ".target = $PREFIX" >> $O/recipeFile.txt

	#paths to genome fasta files (required for Sibelia)
	echo "Sace_S288c.fasta = $REF" >> $O/recipeFile.txt
	echo "$PREFIX.fasta = "$(echo $DRAFT) >> $O/recipeFile.txt

	#reference to use for scaffold naming (optional)
	#echo ".naming_ref = Sace_S288c" >> $O/recipeFile.txt

	echo "ragout -o $O -t 20 --solid-scaffolds $O/recipeFile.txt" >> $CMD
	#echo "1 ragout -o $O -t 20 --solid-scaffolds $O/recipeFile.txt" >> $CMD
done

ccc_mprun glost_launch $CMD
#ccc_mprun -f $CMD
rm -f $CMD*

for S in $(cat $BATCHDIR/Batch2.txt)
do
	DRAFT=$DRAFTDIR/$S.*.hapoGPolished1.fasta
	PREFIX=$(basename $DRAFT .fasta | tr "." "_")
	O=$WORKDIR/$PREFIX
	# Merge scaffolded and unplaced contigs
	cat $O/$PREFIX"_scaffolds.fasta" $O/$PREFIX"_unplaced.fasta" > $O/$PREFIX"_ScaffAndUnplaced.fasta"
done

# Remove sibelia-workdir and Ragout output
rm -rf $WORKDIR/*/sibelia-workdir
rm -f $WORKDIR/*/ragout.log
rm -f $WORKDIR/*/recipeFile.txt


