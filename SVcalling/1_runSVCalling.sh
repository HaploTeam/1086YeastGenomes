#!/bin/bash
#SBATCH -p fast               # partition all / fast
#SBATCH -n 40                      # number of cores
#SBATCH --mem-per-cpu=200
#SBATCH -o slurm.%N.%j.out        # STDOUT
#SBATCH -e slurm.%N.%j.err        # STDERR

# =============================================================================
# Project         : 1000 ONT Genome assemblies
# title           : SVCalling.sh
# description     : This script callv Structural Variants (SV) using MUM&Co 
#					(https://github.com/SAMtoBAM/MUMandCo) and Jasmine SV
#					(https://github.com/mkirsche/Jasmine). 
#                   It uses phased genome assemblies to infer the zygosity
#					of each SV. 
# author          : vloegler
# date            : 2023/01/30
# version         : 1.0
# usage           : bash SVCalling.sh
# =============================================================================

# Load Mummer, bcftools and samtools
module load mummer/4.0.0beta2
module load bcftools/1.16
module load samtools/1.16.1
module load ncbi-blast+/2.14.0
# Use MUM&Co v3.9 (Custom version similar to 3.8)
# See https://github.com/VLoegler/MUMandCo/blob/patch-1/mumandco_v3.9.sh
MUMANDCo=/home/vloegler/SVCalling_Sace1000ONT/02-scripts/MUMandCo/mumandco_v3.9.sh

DATADIR=/home/vloegler/SVCalling_Sace1000ONT/03-data
ASMRX=$DATADIR/Sace_GenomeAssemblies_1011ONTProject_VL_20230927 # Assembly archive
WORKDIR=/home/vloegler/SVCalling_Sace1000ONT/04-analysis
SCRIPTDIR=/home/vloegler/SVCalling_Sace1000ONT/02-scripts

REF=$DATADIR/Sace_S288c_reference_FullMatrixID.fna

# ==========================================
# STEP 1: Generate the total list of strains
cat $ASMRX/Table_* | cut -f 3 | grep -v "Standardized_name" | sort -u > $WORKDIR/Strains.txt

# =========================================================
# STEP 2: Split the list in collapsed and phased assemblies
cat $ASMRX/Table_PhasedAssembliesStatistics.tsv | cut -f 3 | grep -v "Standardized_name" | sort -u > $WORKDIR/StrainsPhased.txt
cat $ASMRX/Table_PhenovarAssembliesStatistics.tsv | grep ".HP1" | cut -f 3 | grep -v "Standardized_name" | sort -u >> $WORKDIR/StrainsPhased.txt

sed 's/^/^/g' $WORKDIR/StrainsPhased.txt > $WORKDIR/StrainsPhased.txt2
sed -i 's/$/$/g' $WORKDIR/StrainsPhased.txt2
grep -v -f $WORKDIR/StrainsPhased.txt2 $WORKDIR/Strains.txt > $WORKDIR/StrainsCollapsed.txt
rm -f $WORKDIR/StrainsPhased.txt2
echo "Strain lists generated"

# Link phenovar Asm in the phased and collapsed folders
for F in $(cat $ASMRX/Table_PhenovarAssembliesStatistics.tsv | cut -f 9 | grep -v "Filename")
do
	if [ $(echo $F | cut -d "." -f 2) = "HP1" ] || [ $(echo $F | cut -d "." -f 2) = "HP2" ]; then
		ln -f -s $ASMRX/PhenovarAssemblies/BestAssemblies/$F $ASMRX/PhasedAssemblies/BestAssemblies/
	else
		ln -f -s $ASMRX/PhenovarAssemblies/BestAssemblies/$F $ASMRX/CollapsedAssemblies/BestAssemblies/
	fi
done

# =================================
# STEP 3: Call SV for each assembly
#rm -rf $WORKDIR/01_SV_VCF
#mkdir $WORKDIR/01_SV_VCF
cd $WORKDIR/01_SV_VCF


# Call SV for collapsed assemblies
for S in $(cat $WORKDIR/StrainsCollapsed.txt)
do
	ASSEMBLY=$ASMRX/CollapsedAssemblies/BestAssemblies/$S.Final.fasta
	if [ ! -d $WORKDIR/01_SV_VCF/${S}_output ]; then
		mkdir $WORKDIR/01_SV_VCF/$S.MumAndCoWorkdir
		cd $WORKDIR/01_SV_VCF/$S.MumAndCoWorkdir/
		cut -d " " -f 1 $ASSEMBLY > $(basename $ASSEMBLY)
		srun --exclusive -N1 -n1 bash $MUMANDCo -r $REF -q $(basename $ASSEMBLY) -g 12000000 -t 4 -o $S &
	fi
	ASSEMBLY2=$ASMRX/CollapsedAssemblies/SecondBestAssemblies/$S.SecondBest.fasta
	if [ -f $ASSEMBLY2 ]; then
		if [ ! -d $WORKDIR/01_SV_VCF/${S}_2_output ]; then
			mkdir $WORKDIR/01_SV_VCF/${S}_2.MumAndCoWorkdir
			cd $WORKDIR/01_SV_VCF/${S}_2.MumAndCoWorkdir/
			cut -d " " -f 1 $ASSEMBLY2 > $(basename $ASSEMBLY2)
			srun --exclusive -N1 -n1 bash $MUMANDCo -r $REF -q $(basename $ASSEMBLY2) -g 12000000 -t 4 -o ${S}_2 &
		fi
	fi
done
wait

cd $WORKDIR/01_SV_VCF
for S in $(cat $WORKDIR/StrainsCollapsed.txt)
do
	mv $WORKDIR/01_SV_VCF/$S.MumAndCoWorkdir/${S}_output $WORKDIR/01_SV_VCF/
	rm -rf $WORKDIR/01_SV_VCF/$S.MumAndCoWorkdir
	mv $WORKDIR/01_SV_VCF/${S}_2.MumAndCoWorkdir/${S}_2_output $WORKDIR/01_SV_VCF/
	rm -rf $WORKDIR/01_SV_VCF/${S}_2.MumAndCoWorkdir
done

# Call SV for phased assemblies
for S in $(cat $WORKDIR/StrainsPhased.txt)
do
	for HP in HP1 HP2
	do
		ASSEMBLY=$ASMRX/PhasedAssemblies/BestAssemblies/$S.$HP.Final.fasta
		if [ ! -d $WORKDIR/01_SV_VCF/${S}.${HP}_output ]; then
			mkdir $WORKDIR/01_SV_VCF/$S.$HP.MumAndCoWorkdir
			cd $WORKDIR/01_SV_VCF/$S.$HP.MumAndCoWorkdir/
			cut -d " " -f 1 $ASSEMBLY > $(basename $ASSEMBLY)
			srun --exclusive -N1 -n1 bash $MUMANDCo -r $REF -q $(basename $ASSEMBLY) -g 12000000 -t 4 -o $S.$HP &
		fi
		ASSEMBLY2=$ASMRX/PhasedAssemblies/SecondBestAssemblies/$S.$HP.SecondBest.fasta
		if [ -f $ASSEMBLY2 ]; then
			if [ ! -d $WORKDIR/01_SV_VCF/${S}_2.${HP}_output ]; then
				mkdir $WORKDIR/01_SV_VCF/${S}_2.$HP.MumAndCoWorkdir
				cd $WORKDIR/01_SV_VCF/${S}_2.$HP.MumAndCoWorkdir/
				cut -d " " -f 1 $ASSEMBLY2 > $(basename $ASSEMBLY2)
				srun --exclusive -N1 -n1 bash $MUMANDCo -r $REF -q $(basename $ASSEMBLY2) -g 12000000 -t 4 -o ${S}_2.$HP &
			fi
		fi
	done
done
wait

cd $WORKDIR/01_SV_VCF
for S in $(cat $WORKDIR/StrainsPhased.txt)
do
	for HP in HP1 HP2
	do
		mv $WORKDIR/01_SV_VCF/$S.$HP.MumAndCoWorkdir/${S}.${HP}_output $WORKDIR/01_SV_VCF/
		rm -rf $WORKDIR/01_SV_VCF/$S.$HP.MumAndCoWorkdir
		mv $WORKDIR/01_SV_VCF/${S}_2.$HP.MumAndCoWorkdir/${S}_2.${HP}_output $WORKDIR/01_SV_VCF/
		rm -rf $WORKDIR/01_SV_VCF/${S}_2.$HP.MumAndCoWorkdir
	done
done

echo "SV Calling done"

