#!/bin/bash
#MSUB -r FilterSNPs
#MSUB -n 1
#MSUB -c 4
#MSUB -T 86400
#MSUB -m scratch,work,store
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load bcftools
module load python

SCRIPTDIR=/ccc/scratch/cont007/fg0006/loeglerv/GWAS_Sace1000ONT/02-scripts
DATADIR=/ccc/scratch/cont007/fg0006/loeglerv/GWAS_Sace1000ONT/03-data
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/GWAS_Sace1000ONT/04-analysis

# Process SNPs and INDELs matrices

#rm -rf $WORKDIR/01_Genotypes
#mkdir $WORKDIR/01_Genotypes
cd $WORKDIR/01_Genotypes

if false; then
    # First remove AIS strain, change XTRA_DGZ to JLL and filter the 966 strains
    # Then Filter for DP10, GQ20, missing genotype, Excess of heterozygosity
    # Finally separate SNPs and INDELs

    bcftools query -l $DATADIR/full1087Matrix.AllPositions.vcf.gz | sed 's/XTRA_DGZ/JLL/g' > NewNames.txt
    bcftools view $DATADIR/full1087Matrix.AllPositions.vcf.gz | bcftools reheader --samples NewNames.txt | \
        bcftools norm -m -any | \
        bcftools view --samples ^AIS | \
        bcftools view --samples-file $DATADIR/966Strains_GWAS.txt | \
        bcftools +fill-tags | \
        bcftools view --min-ac 1 | \
        bcftools norm -m +any | \
        bcftools +setGT -- -t q -n . -e 'FMT/DP>=10' | \
        bcftools +setGT -- -t q -n . -i 'FMT/GQ<20' | \
        bcftools +fill-tags | \
        bcftools view -i 'F_MISSING<0.01' | \
        bcftools view -e 'ExcHet < 0.99' | \
        bcftools +fill-tags | \
        bcftools view --min-ac 1 --threads 8 -O z -o full966Matrix.Var.DP10.GQ20.99pNonMiss.ExcHet99.vcf.gz

    # Filter SNPs
    bcftools view -e 'ALT="*" || (type!="snp" && type!="ref")' full966Matrix.Var.DP10.GQ20.99pNonMiss.ExcHet99.vcf.gz | \
        bcftools annotate --set-id 'rs_%CHROM\_%POS\_%REF\_%FIRST_ALT' | \
        bcftools view -m2 -M2 | \
        bcftools view --threads 8 -O z -o full966Matrix.Var.DP10.GQ20.99pNonMiss.ExcHet99.SNPs.biallelic.vcf.gz

    # Filter INDELs
    bcftools view -e '(type!="indel" && type!="ref")' full966Matrix.Var.DP10.GQ20.99pNonMiss.ExcHet99.vcf.gz | \
        bcftools annotate --set-id 'indel_%CHROM\_%POS' | \
        bcftools view -m2 -M2 | \
        bcftools view --threads 8 -O z -o full966Matrix.Var.DP10.GQ20.99pNonMiss.ExcHet99.INDELs.biallelic.vcf.gz

    # Change chromosome names
    ( echo "chromosome1 1"
      echo "chromosome2 2"
      echo "chromosome3 3"
      echo "chromosome4 4"
      echo "chromosome5 5"
      echo "chromosome6 6"
      echo "chromosome7 7"
      echo "chromosome8 8"
      echo "chromosome9 9"
      echo "chromosome10 10"
      echo "chromosome11 11"
      echo "chromosome12 12"
      echo "chromosome13 13"
      echo "chromosome14 14"
      echo "chromosome15 15"
      echo "chromosome16 16") > ChrNames.txt
    # Change chr names with bcftools
    bcftools annotate --rename-chrs ChrNames.txt full966Matrix.Var.DP10.GQ20.99pNonMiss.ExcHet99.SNPs.biallelic.vcf.gz | bcftools view -Oz -o full966Matrix.Var.DP10.GQ20.99pNonMiss.ExcHet99.SNPs.biallelic.ChrIDSet.vcf.gz --threads 8
    bcftools annotate --rename-chrs ChrNames.txt full966Matrix.Var.DP10.GQ20.99pNonMiss.ExcHet99.INDELs.biallelic.vcf.gz | bcftools view -Oz -o full966Matrix.Var.DP10.GQ20.99pNonMiss.ExcHet99.INDELs.biallelic.ChrIDSet.vcf.gz --threads 8
    bcftools annotate --rename-chrs ChrNames.txt $DATADIR/SV.1087Samples.vcf.gz | bcftools view --samples ^AIS -Oz -o SV.1086Samples.ChrIDSet.vcf.gz --threads 8
    rm -f ChrNames.txt NewNames.txt

    # Create CNV matrix
    python $SCRIPTDIR/convertCNVsToPlink.py -i $DATADIR/full1086Sace_MatrixCNV.tsv.gz --gff $DATADIR/Sace_SGD_R64-4-1_20230830.gff --additionalLocation $DATADIR/Pangenome.NonRef.Localization.tsv -o full1086.CNVs.plink
fi
cd $WORKDIR/01_Genotypes

# Convert to plink
plink --vcf full966Matrix.Var.DP10.GQ20.99pNonMiss.ExcHet99.SNPs.biallelic.ChrIDSet.vcf.gz --out full966.SNPs.plink --make-bed
plink --vcf full966Matrix.Var.DP10.GQ20.99pNonMiss.ExcHet99.INDELs.biallelic.ChrIDSet.vcf.gz --out full966.INDELs.plink --make-bed
plink --vcf SV.1086Samples.ChrIDSet.vcf.gz --keep $DATADIR/966Strains_GWAS.plink.txt --out full966.SVs.plink --make-bed
plink --file full1086.CNVs.plink --keep $DATADIR/966Strains_GWAS.plink.txt --out full966.CNVs.plink --make-bed

# Add Type information to bim files
sed -i 's/rs_/SNP_/g' full966.SNPs.plink.bim
sed -i 's/indel_/INDEL_/g' full966.INDELs.plink.bim
awk -i inplace 'BEGIN{OFS="\t"} {$2 = "SV_"$2; print $0}' full966.SVs.plink.bim
awk -i inplace 'BEGIN{OFS="\t"} {$2 = "CNV_"$2; print $0}' full966.CNVs.plink.bim


# Filter MAF
plink --bfile full966.SNPs.plink --maf 0.05 --out full966.SNPs.MAF5.plink --make-bed
plink --bfile full966.INDELs.plink --maf 0.05 --out full966.INDELs.MAF5.plink --make-bed
plink --bfile full966.SVs.plink --maf 0.05 --out full966.SVs.MAF5.plink --make-bed
plink --bfile full966.CNVs.plink --maf 0.05 --out full966.CNVs.MAF5.plink --make-bed

# Prune LD
plink --bfile full966.SNPs.MAF5.plink --indep-pairwise 50kb 1 0.8 --out full966.SNPs.MAF5.indep50k1k08r2.plink
plink --bfile full966.SNPs.MAF5.plink --extract full966.SNPs.MAF5.indep50k1k08r2.plink.prune.in --out full966.SNPs.MAF5.indep50k1k08r2.plink --make-bed
plink --bfile full966.INDELs.MAF5.plink --indep-pairwise 50kb 1 0.8 --out full966.INDELs.MAF5.indep50k1k08r2.plink
plink --bfile full966.INDELs.MAF5.plink --extract full966.INDELs.MAF5.indep50k1k08r2.plink.prune.in --out full966.INDELs.MAF5.indep50k1k08r2.plink --make-bed
plink --bfile full966.SVs.MAF5.plink --indep-pairwise 50kb 1 0.8 --out full966.SVs.MAF5.indep50k1k08r2.plink
plink --bfile full966.SVs.MAF5.plink --extract full966.SVs.MAF5.indep50k1k08r2.plink.prune.in --out full966.SVs.MAF5.indep50k1k08r2.plink --make-bed
plink --bfile full966.CNVs.MAF5.plink --indep-pairwise 50kb 1 0.8 --out full966.CNVs.MAF5.indep50k1k08r2.plink
plink --bfile full966.CNVs.MAF5.plink --extract full966.CNVs.MAF5.indep50k1k08r2.plink.prune.in --out full966.CNVs.MAF5.indep50k1k08r2.plink --make-bed

# Combine all matrices
plink --bfile full966.SNPs.MAF5.indep50k1k08r2.plink --bmerge full966.INDELs.MAF5.indep50k1k08r2.plink --out full966.SNPs.INDELs.MAF5.indep50k1k08r2.plink --make-bed
plink --bfile full966.SNPs.INDELs.MAF5.indep50k1k08r2.plink --bmerge full966.SVs.MAF5.indep50k1k08r2.plink --out full966.SNPs.INDELs.SVs.MAF5.indep50k1k08r2.plink --make-bed
plink --bfile full966.SNPs.INDELs.SVs.MAF5.indep50k1k08r2.plink --bmerge full966.CNVs.MAF5.indep50k1k08r2.plink --out full966.SNPs.INDELs.SVs.CNVs.MAF5.indep50k1k08r2.plink --make-bed
rm -f full966.SNPs.INDELs.MAF5.indep50k1k08r2.plink.*
rm -f full966.SNPs.INDELs.SVs.MAF5.indep50k1k08r2.plink.*


# Compute LD for each type of variant
plink --bfile full966.SNPs.MAF5.plink --r2 --out full966.SNPs.MAF5.LD_0.8r2_50kb.plink --ld-window 100000 --ld-window-kb 50 --ld-window-r2 0.8
plink --bfile full966.INDELs.MAF5.plink --r2 --out full966.INDELs.MAF5.LD_0.8r2_50kb.plink --ld-window 100000 --ld-window-kb 50 --ld-window-r2 0.8
plink --bfile full966.SVs.MAF5.plink --r2 --out full966.SVs.MAF5.LD_0.8r2_50kb.plink --ld-window 100000 --ld-window-kb 50 --ld-window-r2 0.8
plink --bfile full966.CNVs.MAF5.plink --r2 --out full966.CNVs.MAF5.LD_0.8r2_50kb.plink --ld-window 100000 --ld-window-kb 50 --ld-window-r2 0.8
cat *.MAF5.LD_0.8r2_50kb.plink.ld | gzip > full966.SNPs.INDELs.SVs.CNVs.MAF5.LD_0.8r2_50kb.InependentlyPerTypeOfVariant.plink.ld.gz
