library(calibrate, lib.loc = "/ccc/work/cont007/fg0006/loeglerv/Soft/Rlibraries/")
library(qqman, lib.loc = "/ccc/work/cont007/fg0006/loeglerv/Soft/Rlibraries/")

library(qqman)
library(viridis)

directory = getwd()

args <- commandArgs(trailingOnly = TRUE)

fileIn <- paste(args[1])
threshold <- paste(args[2])
assoc_df <- read.table(fileIn, header=TRUE, sep='\t')
if (threshold!='0') {
  thresh_df <- read.table(threshold, header=TRUE)
  thresh=thresh_df$x[1]
} else {
  thresh=0
}
#print(thresh)
options(warn=1)
zz <- file(paste(fileIn, "_R.log", sep=""), open="wt")
sink(zz, type="message")

toHilight <- F

if(length(args)>3) {
    causal_snps <- paste(args[3])
    if(file.exists(causal_snps)){
        toHilight <-T
        snps_to_hilight <- read.table(causal_snps)

    } else {
        print("causal SNP file doesn't exist")
    }
    gwas_soft <- paste(args[4])
    
} else {
    gwas_soft <- paste(args[3])
}


if (gwas_soft=='gemma') {

    if(toHilight){
        library(pROC)
        ### ROC curves
        postscript(paste(fileIn, "_roc.eps", sep=""))
        predictor <- assoc_df$rs %in% snps_to_hilight$V1
        response_wald <- assoc_df$p_wald
        response_lrt <- assoc_df$p_lrt
        response_score <- assoc_df$p_score
        rocobj1 <- plot.roc(predictor, response_wald, algorithm=2, main="ROC curves", percent=T, print.auc=T, print.auc.y=20, print.auc.x=95, col="#1c61b6")
        rocobj2 <- plot.roc(predictor, response_lrt, add=T, algorithm=2, percent=T, print.auc=T, print.auc.y=15, print.auc.x=95, col="#008600")
        rocobj3 <- plot.roc(predictor, response_score, add=T, algorithm=2, percent=T, print.auc=T, print.auc.y=10, print.auc.x=95, col="#840000")
        legend("bottomright", legend=c("wald", "lrt", "score"), col=c("#1c61b6", "#008600", "#840000"), lwd=2)
        warnings()
        dev.off()
    
        ### Manhattan and qq-plots
        # WALD 
        colnames(assoc_df) <- c("CHR", "SNP", "BP", "N_MISS", "ALLELE1", "ALLELE0", "AF", "BETA", "SE", "L_REMLE", "L_MLE", "P", "P_LRT", "P_SCORE")
        #colnames(assoc_df) <- c("CHR", "RS", "PS", "N_MISS", "ALLELE1", "ALLELE0", "AF", "BETA_1", "BETA_2", "BETA_3", "BETA_4", "BETA_5", "BETA_6", "BETA_7", "VBETA_1_1", "VBETA_1_2", "VBETA_1_3", "VBETA_1_4", "VBETA_1_5", "VBETA_1_6", "VBETA_1_7", "VBETA_2_2", "VBETA_2_3", "VBETA_2_4", "VBETA_2_5", "VBETA_2_6", "VBETA_2_7", "VBETA_3_3", "VBETA_3_4", "VBETA_3_5", "VBETA_3_6", "VBETA_3_7", "VBETA_4_4", "VBETA_4_5", "VBETA_4_6", "VBETA_4_7", "VBETA_5_5", "VBETA_5_6", "VBETA_5_7", "VBETA_6_6", "VBETA_6_7", "VBETA_7_7", "P", "P_LRT", "P_SCORE")
        postscript(paste(fileIn, "_wald.eps", sep=""))
        manhattan(assoc_df, highlight=snps_to_hilight$V1 , suggestiveline = F, genomewideline = F)
        warnings()
        dev.off()
        postscript(paste(fileIn, "_wald_qqplot.eps", sep=""))
        qq(assoc_df$P, main="Q-Q plot of GWAS p-values (wald)")
        warnings()
        dev.off()

        # LRT
        colnames(assoc_df) <- c("CHR", "SNP", "BP", "N_MISS", "ALLELE1", "ALLELE0", "AF", "BETA", "SE", "L_REMLE", "L_MLE", "P_WALD", "P", "P_SCORE")
        #colnames(assoc_df) <- c("CHR", "RS", "PS", "N_MISS", "ALLELE1", "ALLELE0", "AF", "BETA_1", "BETA_2", "BETA_3", "BETA_4", "BETA_5", "BETA_6", "BETA_7", "VBETA_1_1", "VBETA_1_2", "VBETA_1_3", "VBETA_1_4", "VBETA_1_5", "VBETA_1_6", "VBETA_1_7", "VBETA_2_2", "VBETA_2_3", "VBETA_2_4", "VBETA_2_5", "VBETA_2_6", "VBETA_2_7", "VBETA_3_3", "VBETA_3_4", "VBETA_3_5", "VBETA_3_6", "VBETA_3_7", "VBETA_4_4", "VBETA_4_5", "VBETA_4_6", "VBETA_4_7", "VBETA_5_5", "VBETA_5_6", "VBETA_5_7", "VBETA_6_6", "VBETA_6_7", "VBETA_7_7", "P_WALD", "P", "P_SCORE")
        postscript(paste(fileIn, "_lrt.eps", sep=""))
        manhattan(assoc_df, highlight=snps_to_hilight$V1 , suggestiveline = F, genomewideline = F)
        warnings()
        dev.off()
        postscript(paste(fileIn, "_lrt_qqplot.eps", sep=""))
        qq(assoc_df$P, main="Q-Q plot of GWAS p-values (lrt)")
        warnings()
        dev.off()

        # SCORE
        colnames(assoc_df) <- c("CHR", "SNP", "BP", "N_MISS", "ALLELE1", "ALLELE0", "AF", "BETA", "SE", "L_REMLE", "L_MLE", "P_WALD", "P_LRT", "P")
        #colnames(assoc_df) <- c("CHR", "RS", "PS", "N_MISS", "ALLELE1", "ALLELE0", "AF", "BETA_1", "BETA_2", "BETA_3", "BETA_4", "BETA_5", "BETA_6", "BETA_7", "VBETA_1_1", "VBETA_1_2", "VBETA_1_3", "VBETA_1_4", "VBETA_1_5", "VBETA_1_6", "VBETA_1_7", "VBETA_2_2", "VBETA_2_3", "VBETA_2_4", "VBETA_2_5", "VBETA_2_6", "VBETA_2_7", "VBETA_3_3", "VBETA_3_4", "VBETA_3_5", "VBETA_3_6", "VBETA_3_7", "VBETA_4_4", "VBETA_4_5", "VBETA_4_6", "VBETA_4_7", "VBETA_5_5", "VBETA_5_6", "VBETA_5_7", "VBETA_6_6", "VBETA_6_7", "VBETA_7_7", "P_WALD", "P_LRT", "P")
        postscript(paste(fileIn, "_score.eps", sep=""))
        manhattan(assoc_df, highlight=snps_to_hilight$V1, suggestiveline = F, genomewideline = F )
        warnings()
        dev.off()
        postscript(paste(fileIn, "_score_qqplot.eps", sep=""))
        qq(assoc_df$P, main="Q-Q plot of GWAS p-values (score)")
        warnings()    
        dev.off()
    
    } else {
        ### Manhattan and qq-plots
        # WALD
        colnames(assoc_df) <- c("CHR", "SNP", "BP", "N_MISS", "ALLELE1", "ALLELE0", "AF", "BETA", "SE", "L_REMLE", "L_MLE", "P", "P_LRT", "P_SCORE")
        #colnames(assoc_df) <- c("CHR", "RS", "BP", "N_MISS", "ALLELE1", "ALLELE0", "AF", "BETA_1", "BETA_2", "BETA_3", "BETA_4", "BETA_5", "BETA_6", "BETA_7", "VBETA_1_1", "VBETA_1_2", "VBETA_1_3", "VBETA_1_4", "VBETA_1_5", "VBETA_1_6", "VBETA_1_7", "VBETA_2_2", "VBETA_2_3", "VBETA_2_4", "VBETA_2_5", "VBETA_2_6", "VBETA_2_7", "VBETA_3_3", "VBETA_3_4", "VBETA_3_5", "VBETA_3_6", "VBETA_3_7", "VBETA_4_4", "VBETA_4_5", "VBETA_4_6", "VBETA_4_7", "VBETA_5_5", "VBETA_5_6", "VBETA_5_7", "VBETA_6_6", "VBETA_6_7", "VBETA_7_7", "P", "P_LRT", "P_SCORE")
        postscript(paste(fileIn, "_wald.eps", sep=""))
        jpeg(paste(fileIn, "_wald.jpeg", sep=""))
        manhattan(assoc_df, col=c("cyan3", "darkorange", "red", "forestgreen", "slateblue3", "honeydew4", "brown","blue", "darkgoldenrod1", "darkmagenta", "firebrick", "darkorchid2", "hotpink4", "mediumaquamarine", "rosybrown3", "royalblue"), suggestiveline = F, genomewideline = F )
        warnings()
        dev.off()
        postscript(paste(fileIn, "_wald_qqplot.eps", sep=""))
        qq(assoc_df$P, main="Q-Q plot of GWAS p-values (wald)")
        warnings()
        dev.off()

        # LRT
        colnames(assoc_df) <- c("CHR", "SNP", "BP", "N_MISS", "ALLELE1", "ALLELE0", "AF", "BETA", "SE", "L_REMLE", "L_MLE", "P_WALD", "P", "P_SCORE")
        #colnames(assoc_df) <- c("CHR", "RS", "BP", "N_MISS", "ALLELE1", "ALLELE0", "AF", "BETA_1", "BETA_2", "BETA_3", "BETA_4", "BETA_5", "BETA_6", "BETA_7", "VBETA_1_1", "VBETA_1_2", "VBETA_1_3", "VBETA_1_4", "VBETA_1_5", "VBETA_1_6", "VBETA_1_7", "VBETA_2_2", "VBETA_2_3", "VBETA_2_4", "VBETA_2_5", "VBETA_2_6", "VBETA_2_7", "VBETA_3_3", "VBETA_3_4", "VBETA_3_5", "VBETA_3_6", "VBETA_3_7", "VBETA_4_4", "VBETA_4_5", "VBETA_4_6", "VBETA_4_7", "VBETA_5_5", "VBETA_5_6", "VBETA_5_7", "VBETA_6_6", "VBETA_6_7", "VBETA_7_7", "P_WALD", "P", "P_SCORE")
        postscript(paste(fileIn, "_lrt.eps", sep=""))
        jpeg(paste(fileIn, "_lrt.jpeg", sep=""))
        manhattan(assoc_df, col=c("cyan3", "darkorange", "red", "forestgreen", "slateblue3", "honeydew4", "brown","blue", "darkgoldenrod1", "darkmagenta", "firebrick", "darkorchid2", "hotpink4", "mediumaquamarine", "rosybrown3", "royalblue"), suggestiveline = F, genomewideline = F )
        warnings()
        dev.off()
        postscript(paste(fileIn, "_lrt_qqplot.eps", sep=""))
        qq(assoc_df$P, main="Q-Q plot of GWAS p-values (lrt)")
        warnings()
        dev.off()

        # SCORE
        colnames(assoc_df) <- c("CHR", "SNP", "BP", "N_MISS", "ALLELE1", "ALLELE0", "AF", "BETA", "SE", "L_REMLE", "L_MLE", "P_WALD", "P_LRT", "P")
        #colnames(assoc_df) <- c("CHR", "RS", "BP", "N_MISS", "ALLELE1", "ALLELE0", "AF", "BETA_1", "BETA_2", "BETA_3", "BETA_4", "BETA_5", "BETA_6", "BETA_7", "VBETA_1_1", "VBETA_1_2", "VBETA_1_3", "VBETA_1_4", "VBETA_1_5", "VBETA_1_6", "VBETA_1_7", "VBETA_2_2", "VBETA_2_3", "VBETA_2_4", "VBETA_2_5", "VBETA_2_6", "VBETA_2_7", "VBETA_3_3", "VBETA_3_4", "VBETA_3_5", "VBETA_3_6", "VBETA_3_7", "VBETA_4_4", "VBETA_4_5", "VBETA_4_6", "VBETA_4_7", "VBETA_5_5", "VBETA_5_6", "VBETA_5_7", "VBETA_6_6", "VBETA_6_7", "VBETA_7_7", "P_WALD", "P_LRT", "P")
        postscript(paste(fileIn, "_score.jpeg", sep=""))
        jpeg(paste(fileIn, "_score.eps", sep=""))
        manhattan(assoc_df, col=c("cyan3", "darkorange", "red", "forestgreen", "slateblue3", "honeydew4", "brown","blue", "darkgoldenrod1", "darkmagenta", "firebrick", "darkorchid2", "hotpink4", "mediumaquamarine", "rosybrown3", "royalblue"), suggestiveline = F, genomewideline = F )
        warnings()
        dev.off()
        postscript(paste(fileIn, "_score_qqplot.eps", sep=""))
        qq(assoc_df$P, main="Q-Q plot of GWAS p-values (score)")
        warnings()
        dev.off()
    }

} else if (gwas_soft=='fastlmm') {
    colnames(assoc_df) <- c("SNP", "CHR", "GENETIC_DIST", "BP", "P", "Q", "N", "NB_EXCLUDED", "EXCLUSION_START", "NULL_LOGLIKE", "ALT_LOGLIKE", "SNP_WEIGHT", "SNP_WEIGHT_SE", "ODDSRATIO", "WALD", "NULL_LOG_D", "NULL_GENETIC_VAR", "NULL_RESIDUAL_VAR", "NULL_BIAS", "NULL_COV01WEIGHT")
    
    postscript(paste(fileIn, "_p.eps", sep=""))
    manhattan(assoc_df, col=c("cyan3", "darkorange", "red", "forestgreen", "slateblue3", "honeydew4", "brown","blue", "darkgoldenrod1", "darkmagenta", "firebrick", "darkorchid2", "hotpink4", "mediumaquamarine", "rosybrown3", "royalblue"), suggestiveline = F, genomewideline =F, annotatePval= thresh )
    warnings()
    dev.off()
    jpeg(paste(fileIn, "_qq.jpeg", sep=""))
    qq(assoc_df$P, main="Q-Q plot of GWAS p-values")
    warnings()
    dev.off()
    # colnames(assoc_df) <- c("SNP", "CHR", "GENETIC_DIST", "BP", "P_VAL", "Q", "N", "NB_EXCLUDED", "EXCLUSION_START", "NULL_LOGLIKE", "ALT_LOGLIKE", "NSP_WEIGHT", "SNP_WEIGHT_SE", "ODDSRATIO", "P", "NULL_LOG_D", "NULL_GENETIC_VAR", "NULL_RESIDUAL_VAR", "NULL_BIAS", "NULL_COV01WEIGHT")
    # postscript(paste(fileIn, "_wald.eps", sep=""))
    # manhattan(assoc_df, col=c("cyan3", "darkorange", "red", "forestgreen", "slateblue3", "honeydew4", "brown","blue", "darkgoldenrod1", "darkmagenta", "firebrick", "darkorchid2", "hotpink4", "mediumaquamarine", "rosybrown3", "royalblue") )
    # warnings()
    # dev.off()

} else if (gwas_soft=='fastlmm207') {
    assoc_df <- assoc_df[1:11]
    assoc_df <- assoc_df[complete.cases(assoc_df[,6]),]
    colnames(assoc_df) <- c("sid_index","SNP","CHR","GenDist","BP","P","SnpWeight","SnpWeightSE","SnpFractVarExpl","Mixing","Nullh2")
    jpeg(paste(fileIn, "_p.jpg", sep=""))
    if(toHilight){
      library(pROC)
      # OLD COLORS
      #manhattan(assoc_df, col=c("cyan3", "darkorange", "red", "forestgreen", "slateblue3", "honeydew4", "brown","blue", "darkgoldenrod1", "darkmagenta", "firebrick", "darkorchid2", "hotpink4", "mediumaquamarine", "rosybrown3", "royalblue"), highlight=snps_to_hilight$V1 , suggestiveline = F, genomewideline = -log10(thresh))
      # NEW VIRIDIS COL
      viridisCol = viridis(7)[c((seq(0, 100, by = 3) %% 7) + 1)]
      # Bicolor
      biCol = rep(c("cornflowerblue", "gray40"), 40)
      manhattan(assoc_df, col=biCol, highlight=snps_to_hilight$V1 , suggestiveline = F, genomewideline = -log10(thresh))
    } else {
      # OLD COLORS
      #manhattan(assoc_df, col=c("cyan3", "darkorange", "red", "forestgreen", "slateblue3", "honeydew4", "brown","blue", "darkgoldenrod1", "darkmagenta", "firebrick", "darkorchid2", "hotpink4", "mediumaquamarine", "rosybrown3", "royalblue"), suggestiveline = F, genomewideline = -log10(thresh))
      # NEW VIRIDIS COL
      viridisCol = viridis(7)[c((seq(0, 100, by = 3) %% 7) + 1)]
      # Bicolor
      biCol = rep(c("cornflowerblue", "gray40"), 40)
      manhattan(assoc_df, col=biCol, suggestiveline = F, genomewideline = -log10(thresh))
    }
    warnings()
    dev.off()
    jpeg(paste(fileIn, "_qq.jpg", sep=""))
    qq(assoc_df$P, main="Q-Q plot of GWAS p-values")
    warnings()
    dev.off()

} else {
    print("The software used for gwas is not recognized : ")
    print(gwas_soft)
}

