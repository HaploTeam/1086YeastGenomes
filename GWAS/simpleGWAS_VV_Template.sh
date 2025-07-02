#!/bin/bash
#MSUB -r simple_GWAS
#MSUB -n 1
#MSUB -c 8
#MSUB -m work,scratch
#MSUB -Q long
#MSUB -T 259200
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

echo $(date)

module purge
module load extenv/fg
module load fastlmm
module load r

# ================================================================
# This script aims to run a Genome-Wide Association (GWA) using
# the FaST-LMM software. 
# ================================================================
# ==============
# = INPUT DATA =
# ==============
WORKDIR=

# Genotype in plink format (basename of the corresponding .bed .bim .fam)
GENO=

# Phenotype directory 
PHENODIR=

# Kinship matrix in plink format (basename of the corresponding .bed .bim .fam)
KINSHIP= # Leave "" if Genotypes as kinship

# Covariance matrix, plink format with (columns Strain Strain Cov with no header)
COVAR= # Leave "" if no covariance matrix

# Number of permutation
NB_PERM=

# GFF file to add CDS information to the output
#GFF="" # Leave "" if not wanted

# Phenotype directory must contains .phen files to be normalized or .norm.phen files already normalized
# with columnes Strain PhenName (with header present)

# ==============================
# = PATH OF THE SOURCE SCRIPTS =
# ==============================
SCRIPT_DIR=/ccc/work/cont007/fg0006/loeglerv/Scripts/GWAS_Scripts/src
RANKBASED_INT=$SCRIPT_DIR/rank_based_Inverse_Normal_Transformation.R # Phenotype normalization
FASTLMM=$SCRIPT_DIR/runFaSTLMMv6.4.py # run FaST-LMM
CALC_GIF=$SCRIPT_DIR/calc_GIF.R # Compute the Genomic Inflation Factor
QQMAN=$SCRIPT_DIR/qqman_script.R # Create Manhattan and QQ plots
LINKAGE=$SCRIPT_DIR/addLinkageGroups.py # Add linkage group information

# ==========
# == MAIN ==
# ==========
echo -e "--- RUNNING FaST-LMM GWAS ---\n"
echo "Genotypes: "$GENO
echo "Phenotypes: "$PHENODIR
if [ ! -z $KINSHIP ]; then
    echo "Kinship: "$KINSHIP
fi
if [ ! -z $COVAR ]; then
    echo "Covariance: "$COVAR
fi
echo "Nb permutations: "$NB_PERM
nbPheno=$(ls -l $PHENODIR | grep phen | grep -cv .norm.phen)
echo -e "\nRunning GWAS on $nbPheno phenotypes"

cd $WORKDIR
# Create output directory
OUTDIR=$WORKDIR/$(date +%Y%m%d_%H%M%S)_fastlmm_results
mkdir -p $OUTDIR

# STEP1: PHENOTYPE NORMALIZATION
# ==============================
nbNormalisation=0 # Pheno to normalize counter
phenoToNormalize="" # List of pheno to normalize
for pheno in ${PHENODIR}/*.phen # For each phenotype
do
    # If phenotype is already normalized, or the norm file exists
    if [[ $pheno == *.norm.phen ]] || [[ -f $PHENODIR/$(basename $pheno .phen).norm.phen ]]
    then
        : # do nothing
    else # If phenotype is not normalized
        # Pheno has to be normalized
        Rscript $RANKBASED_INT $pheno
        let nbNormalisation++
        phenoToNormalize=$phenoToNormalize" "$(basename $pheno .phen)
    fi
done
echo "$nbNormalisation phenotypes to normalize"
if [ $nbNormalisation != 0 ]; then
    echo $phenoToNormalize
fi
# /!\ Phenotype files must not contain NA values. The INT normalisation script should remove NA values
# but be careful when using another normalisation script. FaST-LMM cannot handle NA values.  

# STEP2: GWAS
# ===========
for normPheno in $PHENODIR/*.norm.phen
do
    if [ ! -z $KINSHIP ] && [ ! -z $COVAR ]; then
        ccc_mprun python $FASTLMM -g $GENO -p $normPheno -o $OUTDIR -k $KINSHIP -c $COVAR -n $NB_PERM -t $BRIDGE_MSUB_NCORE
    elif [ ! -z $KINSHIP ]; then
        ccc_mprun python $FASTLMM -g $GENO -p $normPheno -o $OUTDIR -k $KINSHIP -n $NB_PERM -t $BRIDGE_MSUB_NCORE
    elif [ ! -z $COVAR ]; then
        ccc_mprun python $FASTLMM -g $GENO -p $normPheno -o $OUTDIR -c $COVAR -n $NB_PERM -t $BRIDGE_MSUB_NCORE
    else
        ccc_mprun python $FASTLMM -g $GENO -p $normPheno -o $OUTDIR -n $NB_PERM -t $BRIDGE_MSUB_NCORE
    fi
done

# STEP3: GIF, LINKAGE GROUPS and PLOTS
# ====================================
for normPheno in $PHENODIR/*.norm.phen
do
    COND=$(head -1 $normPheno | cut -f 2)
    #Rscript $QQMAN $OUTDIR/$COND/$COND.first_assoc.txt $OUTDIR/$COND/$COND.threshold.txt fastlmm207
    Rscript $CALC_GIF $OUTDIR/$COND/$COND.first_assoc.txt
    python $LINKAGE --genotype $GENO --results $OUTDIR/$COND/$COND.signif_snps.txt --window_size 50 -r2 0.5 --prefix_type SNP INDEL SV CNV
done


echo $(date)

