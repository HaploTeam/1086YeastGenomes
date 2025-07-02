#!/bin/bash

NB_PARA=5

WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/GWAS_RNASeq
GENO=$WORKDIR/genotypes/full969Matrix.SNPs.Biallelic.Var.DP10.99percNonMissing.MAF5.NoSubTelo20kb.plink
PHENODIR=$WORKDIR/phenoTest
KINSHIP=""
COVAR=""
NB_PERM=100

# Create new Pheno dir
for K in $(seq 1 $NB_PARA)
do
	mkdir $WORKDIR/phenotypes_Para$K
done

# Split phenotypes among the K pheno directory
K=1
for F in $(ls $PHENODIR | grep .phen | grep -v .norm.phen)
do
	cp $PHENODIR/$(basename $F .phen).* $WORKDIR/phenotypes_Para$K/
	let K++
	if [ $K -eq $((NB_PARA + 1)) ]; then
		K=1
	fi
done

# Create GWAS script and run
for K in $(seq 1 $NB_PARA)
do
	cp -f /ccc/work/cont007/fg0006/loeglerv/Scripts/GWAS_Scripts/simpleGWAS_VV_Template.sh $WORKDIR/simpleGWAS_VV_Para$K.sh
	sed -i "s#%N.%j.out#Para${K}.%N.%j.out#g" $WORKDIR/simpleGWAS_VV_Para$K.sh
	sed -i "s#%N.%j.err#Para${K}.%N.%j.err#g" $WORKDIR/simpleGWAS_VV_Para$K.sh
	sed -i "s#simple_GWAS#simple_GWAS_P${K}#g" $WORKDIR/simpleGWAS_VV_Para$K.sh
	sed -i "s#WORKDIR=#WORKDIR=${WORKDIR}#g" $WORKDIR/simpleGWAS_VV_Para$K.sh
	sed -i "s#GENO=#GENO=${GENO}#g" $WORKDIR/simpleGWAS_VV_Para$K.sh
	sed -i "s#PHENODIR=#PHENODIR=${WORKDIR}/phenotypes_Para${K}#g" $WORKDIR/simpleGWAS_VV_Para$K.sh
	sed -i "s#KINSHIP=#KINSHIP=\"${KINSHIP}\"#g" $WORKDIR/simpleGWAS_VV_Para$K.sh
	sed -i "s#COVAR=#COVAR=\"${COVAR}\"#g" $WORKDIR/simpleGWAS_VV_Para$K.sh
	sed -i "s#NB_PERM=#NB_PERM=${NB_PERM}#g" $WORKDIR/simpleGWAS_VV_Para$K.sh
	sed -i "s#_fastlmm_results#_fastlmm_results_Para${K}#g" $WORKDIR/simpleGWAS_VV_Para$K.sh
	ccc_msub $WORKDIR/simpleGWAS_VV_Para$K.sh
done


