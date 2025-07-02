#!/bin/bash
#SBATCH -p fast               # partition all / fast
#SBATCH -n 1                      # number of cores
#SBATCH --mem-per-cpu=2000
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

# ======================
# STEP 4: Annotate SV ID
cd $WORKDIR/01_SV_VCF
for VCF in $WORKDIR/01_SV_VCF/*/*.SVs_all.vcf
do
	PREFIX=$(basename $VCF .SVs_all.vcf)
	bcftools annotate --set-id "%SVTYPE\_%CHROM\_%POS\_%qCHR\_%qSTART" $VCF | bcftools sort --temp-dir $WORKDIR/01_SV_VCF/${PREFIX}_output | awk -v PREFIX=$PREFIX 'BEGIN {OFS = "\t"; n = 0 }; {if ($1 ~ /^chromosome/) {n = ++n; $3 = PREFIX":"n"_"$3}; print $0 }' > $WORKDIR/01_SV_VCF/${PREFIX}_output/${PREFIX}.SVs_all.IDset.vcf
done

echo "ID annotation done"


# ==============================
# STEP 5: Merge SVs with Jasmine

# Load Jasmine
source /home/vloegler/.bashrc
conda activate JasmineEnv

rm -rf $WORKDIR/02_JasmineMerging
mkdir $WORKDIR/02_JasmineMerging

cd $WORKDIR/02_JasmineMerging

for VCF in $WORKDIR/01_SV_VCF/*_output/*.SVs_all.IDset.vcf
do
	ln -s $VCF $(basename $VCF)
	echo $(basename $VCF) >> listVCFs.txt
done

jasmine file_list=listVCFs.txt out_file=SV.JasmineMerged.vcf
rmdir output
conda deactivate


# Remove translocations
grep -v "<TRA>" $WORKDIR/02_JasmineMerging/SV.JasmineMerged.vcf > $WORKDIR/02_JasmineMerging/SV.JasmineMerged.NoTransloc.vcf
# Sort and bgzip
bcftools sort $WORKDIR/02_JasmineMerging/SV.JasmineMerged.NoTransloc.vcf | bcftools view -Oz -o $WORKDIR/02_JasmineMerging/SV.JasmineMerged.NoTransloc.vcf.gz
bcftools index $WORKDIR/02_JasmineMerging/SV.JasmineMerged.NoTransloc.vcf.gz

# ============================================================
# STEP 6: Merge translocations
rm -rf $WORKDIR/03_TranslocMerge
mkdir $WORKDIR/03_TranslocMerge
cd $WORKDIR/03_TranslocMerge

python $SCRIPTDIR/TranslocMerge.py -v $WORKDIR/*/*.IDset.vcf -kb 10 -o TranslocationMerged.vcf
bcftools view TranslocationMerged.vcf -Oz -o TranslocationMerged.vcf.gz
bcftools index TranslocationMerged.vcf.gz
rm -f TranslocationMerged.vcf

# ============================================================
# STEP 7: Convert Jasmine output to classical heterozygous VCF
rm -rf $WORKDIR/04_PopulationVCF
mkdir $WORKDIR/04_PopulationVCF
cd $WORKDIR/04_PopulationVCF

# Combine Jasmine output and translocations
bcftools concat -a $WORKDIR/03_TranslocMerge/TranslocationMerged.vcf.gz $WORKDIR/02_JasmineMerging/SV.JasmineMerged.NoTransloc.vcf.gz -Ov -o SV.JasmineMerged.TranslocMerged.vcf
# Transform Jasmine VCF to multi-sample VCF
NSamples=$(cat $WORKDIR/Strains.txt | wc -l)
python $SCRIPTDIR/JasmineToVCF.py -v SV.JasmineMerged.TranslocMerged.vcf -i $WORKDIR/02_JasmineMerging/listVCFs.txt -d $WORKDIR/02_JasmineMerging -o SV.${NSamples}Samples.vcf
rm -f SV.JasmineMerged.TranslocMerged.vcf
# Get insertion sequences
OUTFILE=SV.${NSamples}Samples.InsertionSequences.fasta
rm -f $OUTFILE
for SV in $(bcftools query -f '%ID\n' SV.${NSamples}Samples.vcf | grep '_INS_')
do
	PREFIX=$(echo $SV | cut -d ":" -f 1 | cut -d "_" -f 2-)
	ID=$(echo $SV | cut -d "_" -f 2-)
	echo ">$SV" >> $OUTFILE
	grep $ID $WORKDIR/01_SV_VCF/${PREFIX}_output/${PREFIX}.SVs_all.IDset.vcf | cut -f 5 >> $OUTFILE
done

# Annotate SVs
bcftools query -f '%ID\n' SV.${NSamples}Samples.vcf > OldIDs.txt
bcftools annotate --set-id "%SVTYPE\_%CHROM\_%POS\_len_%SVLEN" SV.${NSamples}Samples.vcf | awk 'BEGIN {OFS = "\t"; n = 0 }; {if ($1 ~ /^chromosome/) {n = ++n; $3 = n"_"$3}; print $0 }' | bcftools view -o SV.${NSamples}Samples.IDset.vcf
rm -f SV.${NSamples}Samples.vcf
mv SV.${NSamples}Samples.IDset.vcf SV.${NSamples}Samples.vcf
bcftools query -f '%ID\n' SV.${NSamples}Samples.vcf > NewIDs.txt

# Change IDs in insertion fasta file
paste OldIDs.txt NewIDs.txt >> IDsEquivalence.tsv
rm -f OldIDs.txt NewIDs.txt
python $SCRIPTDIR/changeFastaID.py -f SV.${NSamples}Samples.InsertionSequences.fasta -id IDsEquivalence.tsv -o SV.${NSamples}Samples.InsertionSequences.IDset.fasta
rm -f SV.${NSamples}Samples.InsertionSequences.fasta
mv SV.${NSamples}Samples.InsertionSequences.IDset.fasta SV.${NSamples}Samples.InsertionSequences.fasta


# Get deletion sequences
OUTFILE=SV.${NSamples}Samples.DeletionSequences.fasta
rm -f $OUTFILE
for DEL in $(bcftools query -f '%ID\n' SV.${NSamples}Samples.vcf | grep '_DEL_')
do
	CHR=$(echo $DEL | cut -d "_" -f 3)
	S=$(echo $DEL | cut -d "_" -f 4)
	L=$(echo $DEL | cut -d "_" -f 6 | tr -d '-')
	E=$(( $S + $L ))
	python $SCRIPTDIR/GenomeAssemblyTools/extractFragment.py -f $REF -id $CHR -s $S -e $E | sed "s/>/>$DEL /g" >> $OUTFILE
done

# Get duplication sequences
OUTFILE=SV.${NSamples}Samples.DuplicationSequences.fasta
rm -f $OUTFILE
for DUP in $(bcftools query -f '%ID\n' SV.${NSamples}Samples.vcf | grep '_DUP_')
do
	CHR=$(echo $DUP | cut -d "_" -f 3)
	S=$(echo $DUP | cut -d "_" -f 4)
	L=$(echo $DUP | cut -d "_" -f 6)
	E=$(( $S + $L ))
	python $SCRIPTDIR/GenomeAssemblyTools/extractFragment.py -f $REF -id $CHR -s $S -e $E | sed "s/>/>$DUP /g" >> $OUTFILE
done

# Get contraction sequences
OUTFILE=SV.${NSamples}Samples.ContractionSequences.fasta
rm -f $OUTFILE
for CONTR in $(bcftools query -f '%ID\n' SV.${NSamples}Samples.vcf | grep '_CONTR_')
do
	CHR=$(echo $CONTR | cut -d "_" -f 3)
	S=$(echo $CONTR | cut -d "_" -f 4)
	L=$(echo $CONTR | cut -d "_" -f 6)
	E=$(( $S + $L ))
	python $SCRIPTDIR/GenomeAssemblyTools/extractFragment.py -f $REF -id $CHR -s $S -e $E | sed "s/>/>$CONTR /g" >> $OUTFILE
done

# Get inversion sequences
OUTFILE=SV.${NSamples}Samples.InversionSequences.fasta
rm -f $OUTFILE
for INV in $(bcftools query -f '%ID\n' SV.${NSamples}Samples.vcf | grep '_INV_')
do
	CHR=$(echo $INV | cut -d "_" -f 3)
	S=$(echo $INV | cut -d "_" -f 4)
	L=$(echo $INV | cut -d "_" -f 6)
	E=$(( $S + $L ))
	python $SCRIPTDIR/GenomeAssemblyTools/extractFragment.py -f $REF -id $CHR -s $S -e $E | sed "s/>/>$INV /g" >> $OUTFILE
done

# ============================================================
# STEP 8: Annotate Ty related SVs

# Merge Ty fasta and LTR fasta
cat $DATADIR/allty-ltr.fasta $DATADIR/exhaustiveTy.fasta > TyDB.fasta

# Identify Ty-related SVs for each type of SVs
python $SCRIPTDIR/AnnotateTySVs.py -f SV.${NSamples}Samples.InsertionSequences.fasta -ty TyDB.fasta -o TyRelatedInsertions.txt
python $SCRIPTDIR/AnnotateTySVs.py -f SV.${NSamples}Samples.DeletionSequences.fasta -ty TyDB.fasta -o TyRelatedDeletions.txt
python $SCRIPTDIR/AnnotateTySVs.py -f SV.${NSamples}Samples.DuplicationSequences.fasta -ty TyDB.fasta -o TyRelatedDuplications.txt
python $SCRIPTDIR/AnnotateTySVs.py -f SV.${NSamples}Samples.ContractionSequences.fasta -ty TyDB.fasta -o TyRelatedContractions.txt
python $SCRIPTDIR/AnnotateTySVs.py -f SV.${NSamples}Samples.InversionSequences.fasta -ty TyDB.fasta -o TyRelatedInversions.txt
rm -f TyDB.fasta

# Change IDs in fasta files
cat TyRelated*.txt > TyRelatedIDs.txt
sed 's/$/_TyRelated/g' TyRelatedIDs.txt > TyRelatedIDs2.txt
paste TyRelatedIDs.txt TyRelatedIDs2.txt > TyRelatedIDs.tsv
rm -f TyRelated*.txt

python $SCRIPTDIR/changeFastaID.py -f SV.${NSamples}Samples.InsertionSequences.fasta -id TyRelatedIDs.tsv -o SV.${NSamples}Samples.InsertionSequences.IDset.fasta
python $SCRIPTDIR/changeFastaID.py -f SV.${NSamples}Samples.DeletionSequences.fasta -id TyRelatedIDs.tsv -o SV.${NSamples}Samples.DeletionSequences.IDset.fasta
python $SCRIPTDIR/changeFastaID.py -f SV.${NSamples}Samples.DuplicationSequences.fasta -id TyRelatedIDs.tsv -o SV.${NSamples}Samples.DuplicationSequences.IDset.fasta
python $SCRIPTDIR/changeFastaID.py -f SV.${NSamples}Samples.ContractionSequences.fasta -id TyRelatedIDs.tsv -o SV.${NSamples}Samples.ContractionSequences.IDset.fasta
python $SCRIPTDIR/changeFastaID.py -f SV.${NSamples}Samples.InversionSequences.fasta -id TyRelatedIDs.tsv -o SV.${NSamples}Samples.InversionSequences.IDset.fasta

rm -f SV.${NSamples}Samples.*Sequences.fasta
for F in SV.${NSamples}Samples.*Sequences.IDset.fasta
do
	mv $F $(basename $F .IDset.fasta).fasta
done

# Change IDs in VCF file
python $SCRIPTDIR/changeVcfID.py -v SV.${NSamples}Samples.vcf -id TyRelatedIDs.tsv -o SV.${NSamples}Samples.IDset.vcf
rm -f SV.${NSamples}Samples.vcf
mv SV.${NSamples}Samples.IDset.vcf SV.${NSamples}Samples.vcf

# Gzip VCF
bcftools view SV.${NSamples}Samples.vcf -Oz -o SV.${NSamples}Samples.vcf.gz
rm -f SV.${NSamples}Samples.vcf

# Get Common SVs
bcftools view --min-af 0.05:minor SV.${NSamples}Samples.vcf.gz -Oz -o SV.${NSamples}Samples.MAF5.vcf.gz

# ============================================================
# STEP 9: Build neighbor joining tree based on SVs
rm -rf $WORKDIR/05_NJTree
mkdir $WORKDIR/05_NJTree
cd $WORKDIR/05_NJTree

plink --vcf $WORKDIR/04_PopulationVCF/SV.${NSamples}Samples.vcf.gz --out SV.${NSamples}Samples.plink --make-bed --allow-extra-chr
Rscript $SCRIPTDIR/NJTree_SVs.R SV.${NSamples}Samples.plink



