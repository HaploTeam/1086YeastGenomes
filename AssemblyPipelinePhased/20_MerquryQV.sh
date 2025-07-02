#!/bin/bash
#MSUB -r MerquryQV
#MSUB -n 128
#MSUB -c 2
#MSUB -Q long
#MSUB -T 259200
#MSUB -m scratch,work
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

echo $(date)

module load extenv/fg
module load glost

# Export meryl to path
export PATH=/ccc/work/cont007/fg0006/loeglerv/Soft/meryl-1.4/bin/:$PATH
export MERQURY=/ccc/work/cont007/fg0006/loeglerv/Soft/merqury

BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/2_ASSEMBLY_Phased
FINALDIR=$BATCHDIR/21_FINAL
WORKDIR=$BATCHDIR/22_MerquryQV
READSDIR=$BATCHDIR/17_IlluminaSeqData

CMD=$WORKDIR/Meryl_kmer.cmd
rm -f $CMD

# CYH name equivalence
NAME_EQ=/ccc/work/cont007/fg0006/loeglerv/Sace_NameEq/CYH_name_eq_ReseqInfo.tsv

for S in $(cat $BATCHDIR/Batch1_Phased.txt $BATCHDIR/Batch2_Phased.txt $BATCHDIR/Batch3_Phased.txt | sort -u)
do
	# Get the standardized strain name
	illuS=$(awk -v S=$S '{if ($1 == S) print $4}' $NAME_EQ)

	# Create one directory per assembly
	mkdir $WORKDIR/$illuS.merquryQV
	cd $WORKDIR/$illuS.merquryQV

	# Define kmer size (17 for yeast, obtained from sh /ccc/workflash/cont007/fg0006/loeglerv/Soft/merqury/best_k.sh <genome_size>)
	k=17

	# Make kmer database from illumina reads
	#echo "cd $PWD && meryl memory=2g threads=2 k=$k count $READSDIR/$S.illuminaSeqData_*.fastq.gz output $illuS.meryl" >> $CMD
	# Or links reads database from Collapsed assemblies
	ln -s /ccc/scratch/cont007/fg0006/loeglerv/1_ASSEMBLY_Collapsed/16_MerquryQV/$illuS.merquryQV/$illuS.meryl
done

ccc_mprun glost_launch $CMD
rm -f $CMD

CMD=$WORKDIR/runMerqury.cmd
rm -f $CMD
for S in $(cat $BATCHDIR/Batch1_Phased.txt $BATCHDIR/Batch2_Phased.txt $BATCHDIR/Batch3_Phased.txt | sort -u)
do
	# Get the standardized strain name
	illuS=$(awk -v S=$S '{if ($1 == S) print $4}' $NAME_EQ)
	# Go to strain directory
	cd $WORKDIR/$illuS.merquryQV

	# Link mercury script and genome assembly
	ln -s /ccc/workflash/cont007/fg0006/loeglerv/Soft/merqury/merqury.sh
	ln -s $FINALDIR/$illuS.HP1.Final.fasta
	ln -s $FINALDIR/$illuS.HP2.Final.fasta
	echo "cd $PWD && ./merqury.sh $illuS.meryl $illuS.HP1.Final.fasta $illuS.HP2.Final.fasta $illuS.out" >> $CMD
done

ccc_mprun glost_launch $CMD
rm -f $CMD

cd $WORKDIR
for F in */*.out.qv
do
	S=$(basename $F .out.qv)
	awk -v S=$S '{if ($1 == "Both") print S"\t"$4}' $F
done > Phased.qv
