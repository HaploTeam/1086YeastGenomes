#!/bin/bash
#SBATCH -p fast               # partition all / fast
#SBATCH -n 40                      # number of cores
#SBATCH --mem-per-cpu=200
#SBATCH -o slurm.%N.%j.out        # STDOUT
#SBATCH -e slurm.%N.%j.err        # STDERR

# =============================================================================
# Project         : 1000 ONT Genome assemblies
# title           : 1_GeneratePrimaryData.sh
# description     : This script generates data required for gene clustering, 
#					and generates the list of strains with high quality 
#					assembly (Merqury Quality Value >= 40). 
# author          : vloegler
# date            : 2023/04/26
# version         : 3.3
# usage           : sbatch 1_GeneratePrimaryData.sh
# =============================================================================

# Load blast+
module load ncbi-blast+

DATADIR=/shared/home/vloegler/Pangenome_Sace1000ONT/03-data
WORKDIR=/shared/home/vloegler/Pangenome_Sace1000ONT/04-analysis
SCRIPTDIR=/shared/home/vloegler/Pangenome_Sace1000ONT/02-scripts
ASMRX=$DATADIR/Sace_GenomeAssemblies_1011ONTProject_VL_20230927 # Assembly archive
PHASEDGENOMES=$ASMRX/PhasedAssemblies/Annotations
COLGENOMES=$ASMRX/CollapsedAssemblies/Annotations
PHENOVARGENOMES=$ASMRX/PhenovarAssemblies/Annotations
cd $WORKDIR


if false; then
	# ==========================================
	# STEP 1: Generate the total list of strains
	cat $ASMRX/Table_* | cut -f 3 | grep -v "Standardized_name" | grep -v AIS | sort -u > $WORKDIR/Strains.txt

	# =========================================================
	# STEP 2: Split the list in collapsed and phased assemblies
	cat $ASMRX/Table_PhasedAssembliesStatistics.tsv | cut -f 3 | grep -v "Standardized_name" | sort -u > $WORKDIR/StrainsPhased.txt
	cat $ASMRX/Table_PhenovarAssembliesStatistics.tsv | grep ".HP1" | cut -f 3 | grep -v "Standardized_name" | sort -u | grep -v AIS >> $WORKDIR/StrainsPhased.txt

	sed 's/^/^/g' $WORKDIR/StrainsPhased.txt > $WORKDIR/StrainsPhased.txt2
	sed -i 's/$/$/g' $WORKDIR/StrainsPhased.txt2
	grep -v -f $WORKDIR/StrainsPhased.txt2 $WORKDIR/Strains.txt > $WORKDIR/StrainsCollapsed.txt
	rm -f $WORKDIR/StrainsPhased.txt2
	echo "Strain lists generated"

	# Make list of high quality genome assemblies
	rm -f $WORKDIR/StrainsCollapsedHighQuality.txt
	for S in $(cat $WORKDIR/StrainsCollapsed.txt)
	do
		MERQV=0
		InCollapsed=$(cat $ASMRX/Table_CollapsedAssembliesStatistics.tsv | awk -v S=$S '{if ($3 == S) print $0}' | wc -l)
		if [ $InCollapsed -ge 1 ]; then
			MERQV=$(cat $ASMRX/Table_CollapsedAssembliesStatistics.tsv | awk -F '\t' -v S=$S '{if ($3 == S) print $15}')
		else # Assembly from Phenovar
			MERQV=$(cat $ASMRX/Table_PhenovarAssembliesStatistics.tsv | awk -F '\t' -v S=$S '{if ($3 == S) print $17}')		
		fi
		echo -e "$S\t$MERQV"
		# Test if Merqury QV is higher than 40
		if (( $(echo "$MERQV >= 40" | bc -l) )); then
			echo $S >> $WORKDIR/StrainsCollapsedHighQuality.txt
		fi
	done

	rm -f $WORKDIR/StrainsPhasedHighQuality.txt
	for S in $(cat $WORKDIR/StrainsPhased.txt)
	do
		InPhased=$(cat $ASMRX/Table_PhasedAssembliesStatistics.tsv | awk -v S=$S '{if ($3 == S) print $0}' | wc -l)
		if [ $InPhased -ge 1 ]; then
			MERQV=$(cat $ASMRX/Table_PhasedAssembliesStatistics.tsv | awk -F '\t' -v S=$S '{if ($3 == S) print $16}' | head -1)
		else # Assembly from Phenovar
			MERQV=$(cat $ASMRX/Table_PhenovarAssembliesStatistics.tsv | awk -F '\t' -v S=$S '{if ($3 == S) print $17}' | head -1)
		fi
		# Test if Merqury QV is higher than 40
		if (( $(echo "$MERQV >= 40" | bc -l) )); then
			echo $S >> $WORKDIR/StrainsPhasedHighQuality.txt
		fi
	done

	echo "Strain lists generated"


	# Link phenovar Annotations in the phased and collapsed folders
	for F in $(cat $ASMRX/Table_PhenovarAssembliesStatistics.tsv | cut -f 9 | grep -v "Filename")
	do
		S=$(basename $F .Final.fasta)
		if [ $(echo $F | cut -d "." -f 2) = "HP1" ] || [ $(echo $F | cut -d "." -f 2) = "HP2" ]; then
			ln -f -s $ASMRX/PhenovarAssemblies/Annotations/$S.nuclear_genome.Final.cds.fa $ASMRX/PhasedAssemblies/Annotations/
			ln -f -s $ASMRX/PhenovarAssemblies/Annotations/$S.nuclear_genome.Final.pep.fa $ASMRX/PhasedAssemblies/Annotations/
			ln -f -s $ASMRX/PhenovarAssemblies/Annotations/$S.nuclear_genome.Final.gff3 $ASMRX/PhasedAssemblies/Annotations/
		else
			ln -f -s $ASMRX/PhenovarAssemblies/Annotations/$S.nuclear_genome.Final.cds.fa $ASMRX/CollapsedAssemblies/Annotations/
			ln -f -s $ASMRX/PhenovarAssemblies/Annotations/$S.nuclear_genome.Final.pep.fa $ASMRX/CollapsedAssemblies/Annotations/
			ln -f -s $ASMRX/PhenovarAssemblies/Annotations/$S.nuclear_genome.Final.gff3 $ASMRX/CollapsedAssemblies/Annotations/
		fi
	done

	# ==================================================================
	# STEP 3: Transfer reference annotation to the denovo annotation of the assemblies
	# Do it for all genome assemblies (All collapsed and all phased)
	rm -rf $WORKDIR/01_TransferAnnotation
	mkdir $WORKDIR/01_TransferAnnotation
	cd $WORKDIR/01_TransferAnnotation

	TRANSFER_ANNOT="python $SCRIPTDIR/PythonScripts/TransferAnnotations.py"

	for S in $(cat $WORKDIR/Strains.txt)
	do
		if [ -f $COLGENOMES/$S.nuclear_genome.Final.cds.fa ]; then
			srun --exclusive -N1 -n1 $TRANSFER_ANNOT -r $DATADIR/orf_coding_all_R64-4-1_20230830.fasta -f $COLGENOMES/$S.nuclear_genome.Final -o $S.nuclear_genome.transAnnot &
		fi
	done
	wait

	for S in $(cat $WORKDIR/StrainsPhased.txt)
	do
		for HP in HP1 HP2
		do
			srun --exclusive -N1 -n1 $TRANSFER_ANNOT -r $DATADIR/orf_coding_all_R64-4-1_20230830.fasta -f $PHASEDGENOMES/$S.$HP.nuclear_genome.Final -o $S.$HP.nuclear_genome.transAnnot &
		done
	done
	wait

	rm -rf Workdir_TransferAnnotations_*

	# ==================================================================
	# STEP 4: Separate Genes with and without orthology in the reference
	rm -rf $WORKDIR/02_RefGenes
	mkdir $WORKDIR/02_RefGenes
	rm -rf $WORKDIR/03_NonRefGenes
	mkdir $WORKDIR/03_NonRefGenes
	cd $WORKDIR

	EXTRACT_FASTA="python $SCRIPTDIR/GenomeAssemblyTools/extractFragment.py"
	FILTER_SIZE="python $SCRIPTDIR/GenomeAssemblyTools/filterContigSize.py"

	for S in $(cat $WORKDIR/StrainsCollapsed.txt)
	do
		echo "grep '>' $WORKDIR/01_TransferAnnotation/$S.nuclear_genome.transAnnot.cds.fa | tr -d '>' | cut -d '|' -f 3 > $S.tmp1" > $S.cmd.sh
		echo "grep '>' $WORKDIR/01_TransferAnnotation/$S.nuclear_genome.transAnnot.cds.fa | tr -d '>' | cut -d '|' -f 1 > $S.tmp2" >> $S.cmd.sh
		echo "grep '>' $WORKDIR/01_TransferAnnotation/$S.nuclear_genome.transAnnot.cds.fa | tr -d '>' > $S.tmp3" >> $S.cmd.sh
		echo "paste $S.tmp1 $S.tmp2 $S.tmp3 > $S.tmp" >> $S.cmd.sh
		echo "rm -f $S.tmp[1-3]" >> $S.cmd.sh
		echo "awk '{if (\$1 == \$2) print \$3}' $S.tmp > $S.NonOrtho.txt" >> $S.cmd.sh
		echo "$EXTRACT_FASTA -f $WORKDIR/01_TransferAnnotation/$S.nuclear_genome.transAnnot.cds.fa -v -F $S.NonOrtho.txt -o 02_RefGenes/$S.OrthoGenes.cds.fa" >> $S.cmd.sh
		echo "$EXTRACT_FASTA -f $WORKDIR/01_TransferAnnotation/$S.nuclear_genome.transAnnot.cds.fa -F $S.NonOrtho.txt -o 03_NonRefGenes/$S.NonOrthoGenes.cds.fa.tmp" >> $S.cmd.sh
		echo "cut -d '|' -f 1 03_NonRefGenes/$S.NonOrthoGenes.cds.fa.tmp > 03_NonRefGenes/$S.NonOrthoGenes.cds.fa" >> $S.cmd.sh
		echo "$FILTER_SIZE -f 03_NonRefGenes/$S.NonOrthoGenes.cds.fa -m 0.1 -o 03_NonRefGenes/$S.NonOrthoGenes.cds" >> $S.cmd.sh
		echo "rm -f 03_NonRefGenes/$S.NonOrthoGenes.cds.fa.tmp" >> $S.cmd.sh
		echo "rm -f $S.tmp" >> $S.cmd.sh
		echo "rm -f $S.NonOrtho.txt" >> $S.cmd.sh
	done

	for S in $(cat $WORKDIR/StrainsPhased.txt)
	do
		for HP in HP1 HP2
		do
			echo "grep '>' $WORKDIR/01_TransferAnnotation/$S.$HP.nuclear_genome.transAnnot.cds.fa | tr -d '>' | cut -d '|' -f 3 > $S.$HP.tmp1" > $S.$HP.cmd.sh
			echo "grep '>' $WORKDIR/01_TransferAnnotation/$S.$HP.nuclear_genome.transAnnot.cds.fa | tr -d '>' | cut -d '|' -f 1 > $S.$HP.tmp2" >> $S.$HP.cmd.sh
			echo "grep '>' $WORKDIR/01_TransferAnnotation/$S.$HP.nuclear_genome.transAnnot.cds.fa | tr -d '>' > $S.$HP.tmp3" >> $S.$HP.cmd.sh
			echo "paste $S.$HP.tmp1 $S.$HP.tmp2 $S.$HP.tmp3 > $S.$HP.tmp" >> $S.$HP.cmd.sh
			echo "rm -f $S.$HP.tmp[1-3]" >> $S.$HP.cmd.sh
			echo "awk '{if (\$1 == \$2) print \$3}' $S.$HP.tmp > $S.$HP.NonOrtho.txt" >> $S.$HP.cmd.sh
			echo "$EXTRACT_FASTA -f $WORKDIR/01_TransferAnnotation/$S.$HP.nuclear_genome.transAnnot.cds.fa -v -F $S.$HP.NonOrtho.txt -o 02_RefGenes/$S.$HP.OrthoGenes.cds.fa" >> $S.$HP.cmd.sh
			echo "$EXTRACT_FASTA -f $WORKDIR/01_TransferAnnotation/$S.$HP.nuclear_genome.transAnnot.cds.fa -F $S.$HP.NonOrtho.txt -o 03_NonRefGenes/$S.$HP.NonOrthoGenes.cds.fa.tmp" >> $S.$HP.cmd.sh
			echo "cut -d '|' -f 1 03_NonRefGenes/$S.$HP.NonOrthoGenes.cds.fa.tmp > 03_NonRefGenes/$S.$HP.NonOrthoGenes.cds.fa" >> $S.$HP.cmd.sh
			echo "$FILTER_SIZE -f 03_NonRefGenes/$S.$HP.NonOrthoGenes.cds.fa -m 0.1 -o 03_NonRefGenes/$S.$HP.NonOrthoGenes.cds" >> $S.$HP.cmd.sh
			echo "rm -f 03_NonRefGenes/$S.$HP.NonOrthoGenes.cds.fa.tmp" >> $S.$HP.cmd.sh
			echo "rm -f $S.$HP.tmp" >> $S.$HP.cmd.sh
			echo "rm -f $S.$HP.NonOrtho.txt" >> $S.$HP.cmd.sh
		done
	done

	for F in *.cmd.sh
	do
		srun --exclusive -N1 -n1 bash $F &
	done
	wait
	echo "Genes segregated between ortho and non ortho genes to the ref. "
	rm -f *.cmd.sh
fi

# ================================
# STEP 5: Non ref All VS All blast
rm -rf $WORKDIR/04_BlastNonOrthoAllVSAll
mkdir $WORKDIR/04_BlastNonOrthoAllVSAll
cd $WORKDIR/04_BlastNonOrthoAllVSAll

python $SCRIPTDIR/GenomeAssemblyTools/filterContigSize.py -f $DATADIR/orf_coding_all_R64-4-1_20230830.fasta -m 0.1 -o orf_coding_all_R64-4-1_20230830
cat $WORKDIR/03_NonRefGenes/*.min0.1kb.fasta orf_coding_all_R64-4-1_20230830.min0.1kb.fasta > AllGenes_Redundant.min0.1kb.fasta
makeblastdb -in AllGenes_Redundant.min0.1kb.fasta -dbtype nucl -title AllGenes_Redundant.min0.1kb.fasta

IDENTITY=95

# Run blast against all for reference genes
blastn -query orf_coding_all_R64-4-1_20230830.min0.1kb.fasta -db AllGenes_Redundant.min0.1kb.fasta -num_threads 4 -dust no -perc_identity $IDENTITY -strand plus -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -out Reference.vsAll.blastn

# Run blast against all for genes of each strain
for S in $(cat $WORKDIR/StrainsCollapsed.txt)
do
	srun --exclusive -N1 -n1 blastn -query $WORKDIR/03_NonRefGenes/$S.NonOrthoGenes.cds.min0.1kb.fasta -db AllGenes_Redundant.min0.1kb.fasta -dust no -num_threads 4 -perc_identity $IDENTITY -strand plus -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -out $S.vsAll.blastn &
done
wait

for S in $(cat $WORKDIR/StrainsPhased.txt)
do
	for HP in HP1 HP2
	do
		srun --exclusive -N1 -n1 blastn -query $WORKDIR/03_NonRefGenes/$S.$HP.NonOrthoGenes.cds.min0.1kb.fasta -db AllGenes_Redundant.min0.1kb.fasta -num_threads 4 -dust no -perc_identity $IDENTITY -strand plus -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -out $S.$HP.vsAll.blastn &
	done
done
wait

# =========================================================
# STEP 6: Non ref All VS All blast for High quality genomes
rm -rf $WORKDIR/05_BlastNonOrthoAllVSAll_HQStrains
mkdir $WORKDIR/05_BlastNonOrthoAllVSAll_HQStrains
cd $WORKDIR/05_BlastNonOrthoAllVSAll_HQStrains

python $SCRIPTDIR/GenomeAssemblyTools/filterContigSize.py -f $DATADIR/orf_coding_all_R64-4-1_20230830.fasta -m 0.1 -o orf_coding_all_R64-4-1_20230830
cat orf_coding_all_R64-4-1_20230830.min0.1kb.fasta > AllGenes_Redundant.min0.1kb.fasta
for S in $(cat $WORKDIR/StrainsCollapsedHighQuality.txt $WORKDIR/StrainsPhasedHighQuality.txt)
do
	if [ -f $WORKDIR/03_NonRefGenes/${S}.NonOrthoGenes.cds.min0.1kb.fasta ]; then
		cat $WORKDIR/03_NonRefGenes/${S}.NonOrthoGenes.cds.min0.1kb.fasta >> AllGenes_Redundant.min0.1kb.fasta
	fi
	if [ -f $WORKDIR/03_NonRefGenes/${S}.HP1.NonOrthoGenes.cds.min0.1kb.fasta ]; then
		cat $WORKDIR/03_NonRefGenes/${S}.HP1.NonOrthoGenes.cds.min0.1kb.fasta >> AllGenes_Redundant.min0.1kb.fasta
	fi
	if [ -f $WORKDIR/03_NonRefGenes/${S}.HP2.NonOrthoGenes.cds.min0.1kb.fasta ]; then
		cat $WORKDIR/03_NonRefGenes/${S}.HP2.NonOrthoGenes.cds.min0.1kb.fasta >> AllGenes_Redundant.min0.1kb.fasta
	fi
done
makeblastdb -in AllGenes_Redundant.min0.1kb.fasta -dbtype nucl -title AllGenes_Redundant.min0.1kb.fasta

IDENTITY=95

# Run blast against all for reference genes
blastn -query orf_coding_all_R64-4-1_20230830.min0.1kb.fasta -db AllGenes_Redundant.min0.1kb.fasta -num_threads 4 -dust no -perc_identity $IDENTITY -strand plus -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -out Reference.vsAll.blastn

# Run blast against all for genes of each strain
for S in $(cat $WORKDIR/StrainsCollapsedHighQuality.txt)
do
	srun --exclusive -N1 -n1 blastn -query $WORKDIR/03_NonRefGenes/$S.NonOrthoGenes.cds.min0.1kb.fasta -db AllGenes_Redundant.min0.1kb.fasta -dust no -num_threads 4 -perc_identity $IDENTITY -strand plus -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -out $S.vsAll.blastn &
done
wait

for S in $(cat $WORKDIR/StrainsPhasedHighQuality.txt)
do
	for HP in HP1 HP2
	do
		srun --exclusive -N1 -n1 blastn -query $WORKDIR/03_NonRefGenes/$S.$HP.NonOrthoGenes.cds.min0.1kb.fasta -db AllGenes_Redundant.min0.1kb.fasta -num_threads 4 -dust no -perc_identity $IDENTITY -strand plus -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -out $S.$HP.vsAll.blastn &
	done
done
wait
