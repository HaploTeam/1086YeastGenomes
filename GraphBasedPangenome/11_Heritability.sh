#!/bin/bash
#MSUB -r MinigraphCactus
#MSUB -n 64
#MSUB -c 2
#MSUB -Q long
#MSUB -T 259200
#MSUB -q milan
#MSUB -m scratch,work,store
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load bcftools
module load r
module load glost

SCRIPTDIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/02-scripts
DATADIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/03-data
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/04-analysis

LDAK=$SCRIPTDIR/ldak5.2.linux
PREFIX=GraphGenotyping.937Samples.Reheader.DP2.TrimAlt.Atomize.VariantsAnnotated

mkdir $WORKDIR/11_HeritabilityEstimates
cd $WORKDIR/11_HeritabilityEstimates

# Filter 966 samples
bcftools view --samples-file $DATADIR/966Strains_GWAS.txt --force-samples $WORKDIR/10_PopulationVCF/GraphGenotyping.2874Samples.Reheader.DP2.TrimAlt.Atomize.VariantsAnnotated.vcf.gz -Oz -o $PREFIX.vcf.gz

# Convert to plink
plink --vcf $PREFIX.vcf.gz --out $PREFIX.plink --make-bed -aec --vcf-half-call 'missing'

# Remove duplicated variants
cut -f 2 $PREFIX.plink.bim | sort | uniq -d > DuplicatedVariants.txt
plink --bfile $PREFIX.plink --out $PREFIX.NoDuplicate.plink --make-bed -aec --exclude DuplicatedVariants.txt

# Change genotypes to 1 and 2 to prevent problematic multicharacter genotypes
cut -f 1-4 $PREFIX.NoDuplicate.plink.bim | sed 's/$/\t1\t2/g' > $PREFIX.NoDuplicate.plink.bim2
rm -f $PREFIX.NoDuplicate.plink.bim
mv $PREFIX.NoDuplicate.plink.bim2 $PREFIX.NoDuplicate.plink.bim
# Change chromosome names
sed -i 's/chromosome//g' $PREFIX.NoDuplicate.plink.bim

# Compute variant weights with LDAK-Thin model
$LDAK --thin thin.GraphGenotyping --bfile $PREFIX.NoDuplicate.plink --window-prune .98 --window-kb 20
awk < thin.GraphGenotyping.in '{print $1, 1}' > weights.GraphGenotyping.thin

# Compute kinship
$LDAK --calc-kins-direct LDAK-Thin.GraphGenotyping --bfile $PREFIX.NoDuplicate.plink --weights weights.GraphGenotyping.thin --power -.25

# Normalize phenotypes
mkdir NormPheno
cp $DATADIR/PhenotypesPlink/*.phen NormPheno/


CMD=NormalizePheno.cmd

for PHENO in NormPheno/*.plink.phen
do
	echo "Rscript $SCRIPTDIR/rank_based_Inverse_Normal_Transformation.PlinkFormat.R $PHENO"
done > $CMD
ccc_mprun glost_launch $CMD
rm -f $CMD
rm -f $WORKDIR/11_Heritability/NormPheno/*.plink.phen

# Compute heritability
CMD=ComputeHeritability.cmd
mkdir HeritabilityEstimates
for PHENO in NormPheno/*.norm.phen
do
	COND=$(basename $PHENO .plink.norm.phen)
	echo "$LDAK --reml HeritabilityEstimates/$COND --pheno $PHENO --covar $DATADIR/1086Samples_Ploidy.covar.tsv --grm LDAK-Thin.GraphGenotyping --constrain YES && rm -f HeritabilityEstimates/$COND.coeff HeritabilityEstimates/$COND.cross HeritabilityEstimates/$COND.indi.* HeritabilityEstimates/$COND.progress HeritabilityEstimates/$COND.reg.* HeritabilityEstimates/$COND.share HeritabilityEstimates/$COND.vars"
done > $CMD

ccc_mprun glost_launch $CMD
rm -f $CMD

echo -e 'Cond\tHerAll' > Heritability_GraphPangenome.tsv
for PHENO in NormPheno/*.norm.phen
do
	COND=$(basename $PHENO .plink.norm.phen)
	H=$(grep "Her_All" HeritabilityEstimates/$COND.reml | cut -d ' ' -f 2)
	echo -e "$COND\t$H"
done >> Heritability_GraphPangenome.tsv


