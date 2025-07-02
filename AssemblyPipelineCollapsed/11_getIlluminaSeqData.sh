#!/bin/bash
#MSUB -r getIlluSeq
#MSUB -n 128
#MSUB -c 1
#MSUB -Q long
#MSUB -T 259200
#MSUB -m scratch,work,store
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

echo $(date)

module load extenv/fg
module load python
module load glost

BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/1_ASSEMBLY_Collapsed
WORKDIR=$BATCHDIR/11_IlluminaSeqData
NAME_EQ=/ccc/work/cont007/fg0006/loeglerv/Sace_NameEq/CYH_name_eq.tsv
ILLUDIR=/ccc/store/cont007/fg0006/fg0006/rawdata/projet_BCM
EXTDIR=/ccc/store/cont007/fg0006/fg0006/rawdata/external
MEXDIR=/ccc/scratch/cont007/fg0006/loeglerv/MexMapping/MexicanAgaveIlluminaReads

CMD="getIlluminaSeqData.cmd"
rm -f $CMD

for S in $(cat $BATCHDIR/Batch2.txt)
do
	# Get standardized Name
	illuS=$(awk -v S=$S '{if ($1 == S) print $4}' $NAME_EQ)
	# Seq Paths
	in1=$ILLUDIR/$illuS/RunsSolexa/*/*_[0-9]_1_*.fastq.gz
	in2=$ILLUDIR/$illuS/RunsSolexa/*/*_[0-9]_2_*.fastq.gz
	gianni1=$EXTDIR/Gianni/$illuS/${illuS}_1.fastq.gz
	gianni2=$EXTDIR/Gianni/$illuS/${illuS}_2.fastq.gz
	maitreya1=$EXTDIR/Maitreya_strains/$illuS/${illuS}_1.fastq.gz
	maitreya2=$EXTDIR/Maitreya_strains/$illuS/${illuS}_2.fastq.gz
	yjm1=$EXTDIR/YJM/$illuS/${illuS}_1.fastq.gz
	yjm2=$EXTDIR/YJM/$illuS/${illuS}_2.fastq.gz
	mex1=$MEXDIR/$illuS.*_R1_clean.fastq.gz
	mex2=$MEXDIR/$illuS.*_R2_clean.fastq.gz

	# Merge files
	echo "cat $in1 $gianni1 $maitreya1 $yjm1 $mex1 > $WORKDIR/$S.illuminaSeqData_1.fastq.gz" >> $CMD
	echo "cat $in2 $gianni2 $maitreya2 $yjm2 $mex2 > $WORKDIR/$S.illuminaSeqData_2.fastq.gz" >> $CMD
done
ccc_mprun glost_launch $CMD
rm -f $CMD

