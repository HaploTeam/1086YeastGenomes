#!/bin/bash
#MSUB -r MinigraphCactus
#MSUB -n 24
#MSUB -c 32
#MSUB -T 86400
#MSUB -q milan
#MSUB -m scratch,work,store
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load glost
module load vg

SCRIPTDIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/02-scripts
DATADIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/03-data
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/04-analysis
#vg=$CCCWORKDIR/Soft/vg

mkdir $WORKDIR/07_Mapping
cd $WORKDIR/07_Mapping

CMD=runMapping.cmd
rm -f $CMD

# Read mapping
for R1 in $CCCSCRATCHDIR/FullSaceMatrix/03-data/SequencingReads/[A-C]*_1.fastq.gz
do
	S=$(basename $R1 _1.fastq.gz)
	if [ ! -f $S.gam ]; then
		R2=$CCCSCRATCHDIR/FullSaceMatrix/03-data/SequencingReads/${S}_2.fastq.gz
		echo "vg giraffe -Z $WORKDIR/06_ReferenceGraph/SacePangenomeGraph.500Haplotypes.giraffe.gbz -d $WORKDIR/06_ReferenceGraph/SacePangenomeGraph.500Haplotypes.dist -m $WORKDIR/06_ReferenceGraph/SacePangenomeGraph.500Haplotypes.min -f $R1 -f $R2 --fragment-mean 350 --fragment-stdev 100 --threads 64 -p -N $S -b fast > $S.gam" >> $CMD
	fi
done

ccc_mprun glost_launch $CMD
rm -f $CMD
