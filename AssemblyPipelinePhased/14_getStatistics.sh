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

BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/2_ASSEMBLY_Phased
WORKDIR=$BATCHDIR/13_nuclearAssemblies

# ===================
# Assembly statistics
# ===================
# To Run in background
#ASSSTAT=/ccc/work/cont007/fg0006/loeglerv/Soft/assembly-stats-1.0.1-docker1/build/assembly-stats
#ASSEMBLIES=""
#for S in $(cat $BATCHDIR/Batch2_Phased.txt)
#do
#ASSEMBLIES=$ASSEMBLIES" $S.*.NoRedundantContigs.Nuclear.fasta"
#done
#$ASSSTAT -t $ASSEMBLIES > Batch2_AssemblyStatistics.tsv &


# ======================
# Find merge chromosomes
# ======================
CMD=findMergedChr.cmd
rm -f $CMD
CENTRO=/ccc/work/cont007/fg0006/loeglerv/Sace_reference/Sace_S288c_Centromeres.fasta
FIND="python /ccc/work/cont007/fg0006/loeglerv/GenomeAssemblyTools/findDoubleCentromeres.py"

OUT=$WORKDIR/Batch2_MergedChromosomes.tsv
for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	for F in $WORKDIR/$S.*.NoRedundantContigs.Nuclear.fasta
	do
		PREFIX=$(basename $F .fasta)
		echo "$FIND -c $CENTRO -d $F -o $WORKDIR/$PREFIX.MergedChromosomes.tsv" >> $CMD
	done
done

ccc_mprun glost_launch $CMD
rm -f $CMD

echo -e "Assembly\tNbMergedContigs" > $OUT
for F in $WORKDIR/*.MergedChromosomes.tsv
do
	tail -1 $F >> $OUT
done
rm -f $WORKDIR/*.MergedChromosomes.tsv

# =============================
# Check presence of chromosomes
# =============================
CMD=checkChrPres.cmd
rm -f $CMD
REF=/ccc/work/cont007/fg0006/loeglerv/Sace_reference/Sace_S288c_reference_FullMatrixID.fna
CHECK="python /ccc/work/cont007/fg0006/loeglerv/GenomeAssemblyTools/checkChromosomePresence.py"

OUT=$WORKDIR/Batch2_ChrPresence.tsv
for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	for F in $WORKDIR/$S.*.NoRedundantContigs.Nuclear.fasta
	do
		PREFIX=$(basename $F .fasta)
		echo "$CHECK -r $REF -d $F -o $WORKDIR/$PREFIX.ChrPresence.tsv -p 80 -pid chromosome1=75 -t 20" >> $CMD
	done
done

ccc_mprun glost_launch $CMD
rm -f $CMD

echo -e "Assembly\tNbChrPresent" > $OUT
for F in $WORKDIR/*.ChrPresence.tsv
do
	tail -1 $F >> $OUT
done
rm -f $WORKDIR/*.ChrPresence.tsv

# ==============================
# Check presence of mitochondria
# ==============================
CMD=checkMitoPres.cmd
rm -f $CMD
REF=/ccc/work/cont007/fg0006/loeglerv/Sace_reference/Sace_S288c_MitoChromosome.fna
CHECK="python /ccc/work/cont007/fg0006/loeglerv/GenomeAssemblyTools/checkChromosomePresence.py"

OUT=$WORKDIR/Batch2_MitoPresence.tsv
for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	for F in $WORKDIR/$S.*.NoRedundantContigs.Nuclear.fasta
	do
		PREFIX=$(basename $F .fasta)
		echo "$CHECK -r $REF -d $F -o $WORKDIR/$PREFIX.MitoPresence.tsv -pid MT=50 -t 20" >> $CMD
	done
done

ccc_mprun glost_launch $CMD
rm -f $CMD

echo -e "Assembly\tNbMitoPresent" > $OUT
for F in $WORKDIR/*.MitoPresence.tsv
do
	tail -1 $F >> $OUT
done
rm -f $WORKDIR/*.MitoPresence.tsv



# ===========================
# Nb contigs to cover 95% ref
# ===========================
CMD=commandsNbChrTo95.cmd
rm -f $CMD

NBCONTIGS="python /ccc/work/cont007/fg0006/loeglerv/GenomeAssemblyTools/NbContigsToXCoverage.py"

REF=/ccc/work/cont007/fg0006/loeglerv/Sace_reference/Sace_S288c_reference_FullMatrixID.fna
OUT=$WORKDIR/Batch2_NbContigsTo95Coverage.tsv
for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	for F in $WORKDIR/$S.*.NoRedundantContigs.Nuclear.fasta
	do
		PREFIX=$(basename $F .fasta)
		echo "$NBCONTIGS -r $REF -d $F -o $WORKDIR/$PREFIX.nbContigsTo95Coverage.tsv -p 95 -t 20" >> $CMD
	done
done

ccc_mprun glost_launch $CMD
rm -f $CMD

echo -e "Assembly\tNbContigsTo95" > $OUT
for F in $WORKDIR/*.nbContigsTo95Coverage.tsv
do
	tail -1 $F >> $OUT
done
rm -f $WORKDIR/*.nbContigsTo95Coverage.tsv



