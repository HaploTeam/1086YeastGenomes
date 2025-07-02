#!/bin/bash
#MSUB -r MC
#MSUB -n 1
#MSUB -c 128
#MSUB -T 86400
#MSUB -q milan
#MSUB -m scratch,work
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load cactus

SCRIPTDIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/02-scripts
DATADIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/03-data
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/04-analysis

cd $WORKDIR

# ==========================================
# STEP 1: Generate the total list of strains
ls -1 $DATADIR/Sace_GenomeAssemblies_1011ONTProject_VL_20230927/*/BestAssemblies/*.Final.fasta | rev | cut -d "/" -f 1 | rev | \
	cut -d "." -f 1 | sort -u | \
	grep -v AIS > Strains.txt

# =========================================================
# STEP 2: Split the list in collapsed and phased assemblies
rm -f StrainsPhased.txt
rm -f StrainsCollapsed.txt
for S in $(cat Strains.txt)
do
	if [ -f $DATADIR/Sace_GenomeAssemblies_1011ONTProject_VL_20230927/PhasedAssemblies/BestAssemblies/$S.HP1.Final.fasta ] || [ -f $DATADIR/Sace_GenomeAssemblies_1011ONTProject_VL_20230927/PhenovarAssemblies/BestAssemblies/$S.HP1.Final.fasta ]
	then
		echo $S >> StrainsPhased.txt
	else
		echo $S >> StrainsCollapsed.txt
	fi
done

# Make list of high quality genome assemblies
cat $DATADIR/Sace_GenomeAssemblies_1011ONTProject_VL_20230927/Table_PhasedAssembliesStatistics.tsv | grep HP[1-2] | cut -f 3,16 | sort -u > Phased.qv
cat $DATADIR/Sace_GenomeAssemblies_1011ONTProject_VL_20230927/Table_PhenovarAssembliesStatistics.tsv | grep HP[1-2] | cut -f 3,17 | sort -u >> Phased.qv
cat $DATADIR/Sace_GenomeAssemblies_1011ONTProject_VL_20230927/Table_CollapsedAssembliesStatistics.tsv | grep -v HP[1-2] | grep -v MerquryQV | cut -f 3,15 | sort -u > Collapsed.qv
cat $DATADIR/Sace_GenomeAssemblies_1011ONTProject_VL_20230927/Table_PhenovarAssembliesStatistics.tsv | grep -v HP[1-2] | cut -f 3,17 | sort -u >> Collapsed.qv

for S in $(cat StrainsPhased.txt)
do
	SCORE=$(awk -v S=$S '{if ($1 == S) print $2}' $WORKDIR/Phased.qv | cut -d "." -f 1)
	if (( $SCORE >= 40 )); then
		echo $S
	fi
done > StrainsHighQuality.txt
for S in $(cat StrainsCollapsed.txt)
do
	SCORE=$(awk -v S=$S '{if ($1 == S) print $2}' $WORKDIR/Collapsed.qv | cut -d "." -f 1)
	if (( $SCORE >= 40 )); then
		echo $S
	fi
done >> StrainsHighQuality.txt
rm -f Phased.qv Collapsed.qv

# ===============================================================
# STEP 3: Select genomes that represent the largest amount of SVs

rm -rf $WORKDIR/03_MC_Input
mkdir $WORKDIR/03_MC_Input
cd $WORKDIR/03_MC_Input
mkdir AllGenomes
ln -s $DATADIR/Sace_GenomeAssemblies_1011ONTProject_VL_20230927/*/BestAssemblies/*.Final.fasta AllGenomes/
ln -s $DATADIR/Sace_S288c_reference_FullMatrixID.fna AllGenomes/

VCF=$DATADIR/SV.1087Samples.vcf.gz
python $SCRIPTDIR/SelectHaplotypes.py -v $VCF -s $WORKDIR/Strains.txt -d $WORKDIR/03_MC_Input/AllGenomes -o $WORKDIR/03_MC_Input/SelectedHaplotypes

# ============================
# STEP 4: Run Minigraph Cactus

rm -rf $WORKDIR/04_MC_GraphPangenome
mkdir $WORKDIR/04_MC_GraphPangenome
cd $WORKDIR/04_MC_GraphPangenome

# Create input file
echo -e "S288c\t$WORKDIR/03_MC_Input/SelectedHaplotypes/Sace_S288c_reference_FullMatrixID.fna" > input.txt
for F in $WORKDIR/03_MC_Input/SelectedHaplotypes/*.fasta
do
	S=$(basename $F .Final.fasta | sed 's/\.HP/\./g')
	echo -e "$S\t$F" >> input.txt
done

rm -rf tmpDirCactus
mkdir tmpDirCactus
TMPDIR=$PWD/tmpDirCactus

rm -rf JobStore
cactus-pangenome ./JobStore input.txt --workDir tmpDirCactus --maxDisk 180G --maxCores 20 --outDir . --outName SacePangenomeGraph.500Haplotypes --reference S288c --vcf --gfa --gbz
