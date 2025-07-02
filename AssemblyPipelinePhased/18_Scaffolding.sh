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
module load python
module load glost

BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/2_ASSEMBLY_Phased
DRAFTDIR=$BATCHDIR/18_Polishing
WORKDIR=$BATCHDIR/19_Scaffolding
REF=/ccc/work/cont007/fg0006/loeglerv/Sace_reference/Sace_S288c_reference_FullMatrixID.fna
EXTRACT="python /ccc/work/cont007/fg0006/loeglerv/GenomeAssemblyTools/extractFragment.py"

CMD="Scaffolding.cmd"
rm -f $CMD*

for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	DRAFT=$DRAFTDIR/$S.*.hapoGPolished1.fasta
	for HP in HP1 HP2
	do
		DRAFTHP=$WORKDIR/$S.$HP.medakaPolished1.hapoGPolished1.fasta
		# Extract haplotype
		echo "$EXTRACT -f $DRAFT -p ${S}_${HP} -o $DRAFTHP" >> ${CMD}1

		# Add Draft_ beforedraft sequence names because it cannot be the same ID than the ref chr
		echo "sed -i 's/>/>Draft_/g' $DRAFTHP" >> ${CMD}2
		# Convert . to _ in the prefix
		PREFIX=$(basename $DRAFTHP .fasta | tr "." "_")

		O=$WORKDIR/$PREFIX
		mkdir $O

		echo "" > $O/recipeFile.txt

		# Create recipe file
		#reference and target genome names (required)
		echo ".references = Sace_S288c" >> $O/recipeFile.txt
		echo ".target = $PREFIX" >> $O/recipeFile.txt

		#paths to genome fasta files (required for Sibelia)
		echo "Sace_S288c.fasta = $REF" >> $O/recipeFile.txt
		echo "$PREFIX.fasta = "$(echo $DRAFTHP) >> $O/recipeFile.txt

		#reference to use for scaffold naming (optional)
		#echo ".naming_ref = Sace_S288c" >> $O/recipeFile.txt

		echo "ragout -o $O -t 20 --solid-scaffolds $O/recipeFile.txt" >> ${CMD}3
	done
done

ccc_mprun glost_launch ${CMD}1
ccc_mprun glost_launch ${CMD}2
ccc_mprun glost_launch ${CMD}3
rm -f $CMD*

for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	for HP in HP1 HP2
	do
		DRAFTHP=$WORKDIR/$S.$HP.medakaPolished1.hapoGPolished1.fasta
		PREFIX=$(basename $DRAFTHP .fasta | tr "." "_")
		O=$WORKDIR/$PREFIX
		# Merge scaffolded and unplaced contigs
		cat $O/$PREFIX"_scaffolds.fasta" $O/$PREFIX"_unplaced.fasta" > $O/$PREFIX"_ScaffAndUnplaced.fasta"
	done
done

# Remove sibelia-workdir and Ragout output
rm -rf $WORKDIR/*/sibelia-workdir
rm -f $WORKDIR/*/ragout.log
rm -f $WORKDIR/*/recipeFile.txt


