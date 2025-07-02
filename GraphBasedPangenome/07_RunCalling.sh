#!/bin/bash
#MSUB -r MinigraphCactus
#MSUB -n 64
#MSUB -c 8
#MSUB -T 86400
#MSUB -q milan
#MSUB -m scratch,work,store
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load glost

SCRIPTDIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/02-scripts
DATADIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/03-data
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/04-analysis
vg=$CCCWORKDIR/Soft/vg

mkdir $WORKDIR/09_Calling
cd $WORKDIR/09_Calling

CMD=runCalling.cmd
rm -f $CMD

# Call variants
for P in $WORKDIR/08_Packing/*.pack
do
	S=$(basename $P .pack)
	if [ ! -f $S.vcf.gz ]; then
		echo "timeout 1800 $vg call $WORKDIR/06_ReferenceGraph/SacePangenomeGraph.500Haplotypes.giraffe.gbz --genotype-snarls --all-snarls --snarls $WORKDIR/06_ReferenceGraph/SacePangenomeGraph.500Haplotypes.snarls.pb -k $P --sample $S -t 16 --ref-sample S288c | bcftools view -Oz -o $S.vcf.gz && bcftools index $S.vcf.gz" >> $CMD
	fi
done

ccc_mprun glost_launch $CMD
rm -f $CMD
