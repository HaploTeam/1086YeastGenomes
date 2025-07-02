#!/bin/bash
#MSUB -r Phasing
#MSUB -n 92
#MSUB -c 1
#MSUB -T 86400
#MSUB -m scratch,work
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load whatshap
module load samtools
module load glost


BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/2_ASSEMBLY_Phased
WORKDIR=$BATCHDIR/06_LRtagged
INBAM=$BATCHDIR/04_LRmapping
INVCF=$BATCHDIR/05_phasedVCF
REF=/ccc/work/cont007/fg0006/loeglerv/Sace_reference/Sace_S288c_reference_FullMatrixID.fna

CMD=$WORKDIR/LRtagging.cmd
rm -f $CMD

for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	BAM=$INBAM/$S.NanoporeReads.bam
	VCF=$INVCF/$S.LongshotCall.noLOH.phased.vcf.gz
	OUT=$WORKDIR/$S.NanoporeReads.haplotagged.bam
	PREFIX=$WORKDIR/$S.NanoporeReads.haplotagged

	#whatshap haplotag -o $OUT --reference $REF --output-threads 4 $VCF $BAM
	#samtools view -h --tag HP:1 $OUT | samtools fastq -0 $PREFIX.HP1.fastq.gz
	#samtools view -h --tag HP:2 $OUT | samtools fastq -0 $PREFIX.HP2.fastq.gz
	#samtools view -h $OUT | grep -v "HP:i:" | samtools fastq -0 $PREFIX.noHPtag.fastq.gz
	#rm -f $OUT

	echo "whatshap haplotag -o $OUT --reference $REF --output-threads 4 $VCF $BAM && samtools view -h --tag HP:1 $OUT | samtools fastq -0 $PREFIX.HP1.fastq.gz && samtools view -h --tag HP:2 $OUT | samtools fastq -0 $PREFIX.HP2.fastq.gz && samtools view -h $OUT | grep -v "HP:i:" | samtools fastq -0 $PREFIX.noHPtag.fastq.gz && rm -f $OUT" >> $CMD
done

ccc_mprun glost_launch $CMD

rm -f $CMD
