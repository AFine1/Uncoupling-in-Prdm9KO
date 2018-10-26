#   The purpose of this file is to identify differentially expressed genes between conditions


##################################################################################################################
############ Load packages #######################################################################################
##################################################################################################################

library("sva")
library("contrast")
library("Biobase")
library("biomaRt")
library("gplots")
library("VennDiagram")


##################################################################################################################
############ Aquire Data #########################################################################################
###### P9DL_1 ####################################################################################################

#Set working directory
setwd("/Users/s-fine/Desktop/Carter/Projects/P9/3 DifferentialExpression/")

#Load ComBat adjusted dataset
setwd("Data")
load("P9DL_gene.rdt")
setwd("..")

#Make an all positive version
#P9DL_pos = P9DL_1 + abs(min(P9DL_1))


##################################################################################################################
############ Run linear model for DE genes #######################################################################
###### P9DL_1  ->  P9K_coeftable + P9K_serrtable + P9K_pvaltable #################################################

#Run Linear Regression for Genotype and Age Effect
genotype_all = as.factor(c("A","A","A","B","B","C","C","C","C","C","A","A","A","A","B","B","B","C","C","C","C","A","A","A","A","A","A","A","B","B","B","C","C","C","C","C","A","A","A","A","A"))
age_all = as.factor(c("A","A","A","A","A","A","A","A","A","A","A","A","A","B","B","B","B","B","B","B","B","B","B","B","B","B","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C"))
LM_all = apply(P9DL_1,1,function(x) lm(x ~ genotype_all + age_all + genotype_all*age_all))
#GLM_all = apply(P9DL_pos,1,function(x) glm(x ~ genotype_all + age_all + genotype_all*age_all,family="poisson"))


#make output tables
coeftable_DL = matrix(NA, nrow = nrow(P9DL_1), ncol =12 )
colnames(coeftable_DL) = c("WT_08_HT_08","WT_08_MT_08","HT_08_MT_08","WT_12_HT_12","WT_12_MT_12","HT_12_MT_12","WT_16_HT_16","WT_16_MT_16","HT_16_MT_16","WT_08_WT_12","WT_08_WT_16","WT_12_WT_16")
rownames(coeftable_DL) = rownames(P9DL_1)

pvaltable_DL = matrix(NA, nrow = nrow(P9DL_1), ncol =12 )
colnames(pvaltable_DL) = c("WT_08_HT_08","WT_08_MT_08","HT_08_MT_08","WT_12_HT_12","WT_12_MT_12","HT_12_MT_12","WT_16_HT_16","WT_16_MT_16","HT_16_MT_16","WT_08_WT_12","WT_08_WT_16","WT_12_WT_16")
rownames(pvaltable_DL) = rownames(P9DL_1)

SEtable_DL = matrix(NA, nrow = nrow(P9DL_1), ncol =12 )
colnames(SEtable_DL) = c("WT_08_HT_08","WT_08_MT_08","HT_08_MT_08","WT_12_HT_12","WT_12_MT_12","HT_12_MT_12","WT_16_HT_16","WT_16_MT_16","HT_16_MT_16","WT_08_WT_12","WT_08_WT_16","WT_12_WT_16")
rownames(SEtable_DL) = rownames(P9DL_1)

#fill output tables
for (i in 1:nrow(P9DL_1)){
  pvaltable_DL[i,1] = summary(LM_all[[i]])$coefficients[2,"Pr(>|t|)"]
  SEtable_DL[i,1] = summary(LM_all[[i]])$coefficients[2,"Std. Error"]
  coeftable_DL[i,1] = summary(LM_all[[i]])$coefficients[2,"Estimate"]
  pvaltable_DL[i,2] = summary(LM_all[[i]])$coefficients[3,"Pr(>|t|)"]
  SEtable_DL[i,2] = summary(LM_all[[i]])$coefficients[3,"Std. Error"]
  coeftable_DL[i,2] = summary(LM_all[[i]])$coefficients[3,"Estimate"]
  pvaltable_DL[i,10] = summary(LM_all[[i]])$coefficients[4,"Pr(>|t|)"]
  SEtable_DL[i,10] = summary(LM_all[[i]])$coefficients[4,"Std. Error"]
  coeftable_DL[i,10] = summary(LM_all[[i]])$coefficients[4,"Estimate"]
  pvaltable_DL[i,11] = summary(LM_all[[i]])$coefficients[5,"Pr(>|t|)"]
  SEtable_DL[i,11] = summary(LM_all[[i]])$coefficients[5,"Std. Error"]
  coeftable_DL[i,11] = summary(LM_all[[i]])$coefficients[5,"Estimate"]
  AA_GBC=contrast(LM_all[[i]], list(genotype_all = "C",age_all = "A"),list(genotype_all = "B",age_all = "A"))
  pvaltable_DL[i,3] = AA_GBC$Pvalue
  SEtable_DL[i,3] = AA_GBC$SE
  coeftable_DL[i,3] = AA_GBC$Contrast
  AB_GAB=contrast(LM_all[[i]], list(genotype_all = "B",age_all = "B"),list(genotype_all = "A",age_all = "B"))
  pvaltable_DL[i,4] = AB_GAB$Pvalue
  SEtable_DL[i,4] = AB_GAB$SE
  coeftable_DL[i,4] = AB_GAB$Contrast
  AB_GAC=contrast(LM_all[[i]], list(genotype_all = "C",age_all = "B"),list(genotype_all = "A",age_all = "B"))
  pvaltable_DL[i,5] = AB_GAC$Pvalue
  SEtable_DL[i,5] = AB_GAC$SE
  coeftable_DL[i,5] = AB_GAC$Contrast
  AB_GBC=contrast(LM_all[[i]], list(genotype_all = "C",age_all = "B"),list(genotype_all = "B",age_all = "B"))
  pvaltable_DL[i,6] = AB_GBC$Pvalue
  SEtable_DL[i,6] = AB_GBC$SE
  coeftable_DL[i,6] = AB_GBC$Contrast
  AC_GAB=contrast(LM_all[[i]], list(genotype_all = "B",age_all = "C"),list(genotype_all = "A",age_all = "C"))
  pvaltable_DL[i,7] = AC_GAB$Pvalue
  SEtable_DL[i,7] = AC_GAB$SE
  coeftable_DL[i,7] = AC_GAB$Contrast 
  AC_GAC=contrast(LM_all[[i]], list(genotype_all = "C",age_all = "C"),list(genotype_all = "A",age_all = "C"))
  pvaltable_DL[i,8] = AC_GAC$Pvalue
  SEtable_DL[i,8] = AC_GAC$SE
  coeftable_DL[i,8] = AC_GAC$Contrast
  AC_GBC=contrast(LM_all[[i]], list(genotype_all = "C",age_all = "C"),list(genotype_all = "B",age_all = "C"))
  pvaltable_DL[i,9] = AC_GBC$Pvalue
  SEtable_DL[i,9] = AC_GBC$SE
  coeftable_DL[i,9] = AC_GBC$Contrast 
  ABC_GA=contrast(LM_all[[i]], list(genotype_all = "A",age_all = "C"),list(genotype_all = "A",age_all = "B"))
  pvaltable_DL[i,12] = ABC_GA$Pvalue
  SEtable_DL[i,12] = ABC_GA$SE
  coeftable_DL[i,12] = ABC_GA$Contrast 
}

#save data
setwd("Data")
save(pvaltable_DL,file="pvaltable_DL.rdt")
save(SEtable_DL,file="SEtable_DL.rdt")
save(coeftable_DL,file="coeftable_DL.rdt")
load("pvaltable_DL.rdt")
load("SEtable_DL.rdt")
load("coeftable_DL.rdt")
setwd("..")




#Rename coefficient table
P9K_coeftable = coeftable_DL
colnames(P9K_coeftable) = c("S_P9K_COEF_08D.W_08D.H","S_P9K_COEF_08D.W_08D.M","S_P9K_COEF_08D.H_08D.M",
                            "S_P9K_COEF_12D.W_12D.H","S_P9K_COEF_12D.W_12D.M","S_P9K_COEF_12D.H_12D.M",
                            "S_P9K_COEF_16D.W_16D.H","S_P9K_COEF_16D.W_16D.M","S_P9K_COEF_16D.H_16D.M",
                            "S_P9K_COEF_08D.W_12D.W","S_P9K_COEF_08D.W_16D.W","S_P9K_COEF_12D.W_16D.W")

#Rename Standard Error table
P9K_serrtable = SEtable_DL
colnames(P9K_serrtable) = c("S_P9K_SERR_08D.W_08D.H","S_P9K_SERR_08D.W_08D.M","S_P9K_SERR_08D.H_08D.M",
                            "S_P9K_SERR_12D.W_12D.H","S_P9K_SERR_12D.W_12D.M","S_P9K_SERR_12D.H_12D.M",
                            "S_P9K_SERR_16D.W_16D.H","S_P9K_SERR_16D.W_16D.M","S_P9K_SERR_16D.H_16D.M",
                            "S_P9K_SERR_08D.W_12D.W","S_P9K_SERR_08D.W_16D.W","S_P9K_SERR_12D.W_16D.W")

#Rename pvalue table
P9K_pvaltable = pvaltable_DL
colnames(P9K_pvaltable) = c("S_P9K_PVAL_08D.W_08D.H","S_P9K_PVAL_08D.W_08D.M","S_P9K_PVAL_08D.H_08D.M",
                            "S_P9K_PVAL_12D.W_12D.H","S_P9K_PVAL_12D.W_12D.M","S_P9K_PVAL_12D.H_12D.M",
                            "S_P9K_PVAL_16D.W_16D.H","S_P9K_PVAL_16D.W_16D.M","S_P9K_PVAL_16D.H_16D.M",
                            "S_P9K_PVAL_08D.W_12D.W","S_P9K_PVAL_08D.W_16D.W","S_P9K_PVAL_12D.W_16D.W")


##################################################################################################################
############ Add other important information #####################################################################
###### P9DL_1  ->  P9DL_isw  #####################################################################################

#run add info function
P9DL_i = addinfo(P9DL_1,bygene=TRUE)

#add specific substages
P9DL_is = AddSpecificSubstages(P9DL_i)

#add WT patterns
P9DL_isw = AddWTPatterns(P9DL_is)

#rename columns
colnames(P9DL_isw) = c("D_P9K_W_08D_03_02","D_P9K_W_08D_08_05","D_P9K_W_08D_08_06",
                       "D_P9K_H_08D_10_07","D_P9K_H_08D_14_09",
                       "D_P9K_M_08D_03_01","D_P9K_M_08D_08_03","D_P9K_M_08D_08_04","D_P9K_M_08D_10_08","D_P9K_M_08D_14_10",
                       "D_RYT_W_08D_ZZ_01","D_RYT_W_08D_ZZ_02","D_RYT_W_08D_ZZ_03",
                       
                       "D_P9K_W_12D_11_15",
                       "D_P9K_H_12D_01_11","D_P9K_H_12D_09_13","D_P9K_H_12D_12_18",
                       "D_P9K_M_12D_01_12","D_P9K_M_12D_09_14","D_P9K_M_12D_11_16","D_P9K_M_12D_12_17",
                       "D_RYT_W_12D_ZZ_01","D_RYT_W_12D_ZZ_02","D_RYT_W_12D_ZZ_03","D_RYT_W_12D_ZZ_04","D_RYT_W_12D_ZZ_05",
                       
                       "D_P9K_W_16D_04_24","D_P9K_W_16D_05_25",
                       "D_P9K_H_16D_02_21","D_P9K_H_16D_06_28","D_P9K_H_16D_07_30",
                       "D_P9K_M_16D_02_22","D_P9K_M_16D_04_23","D_P9K_M_16D_05_26","D_P9K_M_16D_06_27","D_P9K_M_16D_07_29",
                       "D_RYT_W_16D_ZZ_01","D_RYT_W_16D_ZZ_02","D_RYT_W_16D_ZZ_03","D_RYT_W_16D_ZZ_04","D_RYT_W_16D_ZZ_05",
                       
                       "I_CNUM","I_GEID","I_GETY","I_GSRT","I_GEND","I_ENTZ","S_RYT_SpecSub","S_RYT_WTPat")

#rename data
P9K_Data = P9DL_isw
P9DL_isw = P9K_Data

#save data
setwd("Results")
save(P9K_Data,file="P9K_Data.rdt")
save(P9DL_isw,file="P9DL_isw.rdt")
setwd("..")

##################################################################################################################
############ Combine datasets ####################################################################################
###### P9K_Data + P9K_coeftable + P9K_serrtable + P9K_pvaltable  -> P9K_Analysis #################################

#group data
P9K_Analysis = cbind(P9K_Data,P9K_coeftable,P9K_pvaltable,P9K_serrtable)

#save data
setwd("Results")
save(P9K_Analysis,file="P9K_Analysis.rdt")
write.csv(P9K_Analysis,file="/Users/s-fine/Desktop/Carter/Projects/P9/4 GO Enrichment/Data/Background.csv")
setwd("..")

##################################################################################################################
############ Identify DE genes ###################################################################################
###### P9K_Analysis -> P9K_Analysis ##############################################################################

#Load data
setwd("Results")
load("P9K_Analysis.rdt")
setwd("..")

#Run multiple testing correction
S_P9K_PVHA_08D.W_08D.M = p.adjust(P9K_Analysis$S_P9K_PVAL_08D.W_08D.M)
S_P9K_PVFA_08D.W_08D.M = p.adjust(P9K_Analysis$S_P9K_PVAL_08D.W_08D.M,method="fdr")

S_P9K_PVHA_12D.W_12D.M = p.adjust(P9K_Analysis$S_P9K_PVAL_12D.W_12D.M)
S_P9K_PVFA_12D.W_12D.M = p.adjust(P9K_Analysis$S_P9K_PVAL_12D.W_12D.M,method="fdr")

S_P9K_PVHA_16D.W_16D.M = p.adjust(P9K_Analysis$S_P9K_PVAL_16D.W_16D.M)
S_P9K_PVFA_16D.W_16D.M = p.adjust(P9K_Analysis$S_P9K_PVAL_16D.W_16D.M,method="fdr")

#Add new pvalues to table
P9K_Analysis$S_P9K_PVHA_08D.W_08D.M = S_P9K_PVHA_08D.W_08D.M
P9K_Analysis$S_P9K_PVFA_08D.W_08D.M = S_P9K_PVFA_08D.W_08D.M
P9K_Analysis$S_P9K_PVHA_12D.W_12D.M = S_P9K_PVHA_12D.W_12D.M
P9K_Analysis$S_P9K_PVFA_12D.W_12D.M = S_P9K_PVFA_12D.W_12D.M
P9K_Analysis$S_P9K_PVHA_16D.W_16D.M = S_P9K_PVHA_16D.W_16D.M
P9K_Analysis$S_P9K_PVFA_16D.W_16D.M = S_P9K_PVFA_16D.W_16D.M

#save data
setwd("Results")
save(P9K_Analysis,file="P9K_Analysis.rdt")
load("P9K_Analysis.rdt")
setwd("..")


##################################################################################################################
############ Make DE lists for each age for various cutoffs ######################################################
###### P9K_Analysis -> DE gene lists #############################################################################

#Make DE lists for each age for various cutoffs

#################
#Day 8

P9K_WT08_MT08_P0.05 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVAL_08D.W_08D.M)<0.05),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_08D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_08D.W_08D.M",colnames(P9K_Analysis)))]
  P9K_WT08_MT08_P0.05_up  = P9K_WT08_MT08_P0.05[which(P9K_WT08_MT08_P0.05$S_P9K_COEF_08D.W_08D.M>0),] 
    P9K_WT08_MT08_P0.05_LFC0.5_up  = P9K_WT08_MT08_P0.05_up[which(P9K_WT08_MT08_P0.05_up$S_P9K_COEF_08D.W_08D.M>0.5),] 
    P9K_WT08_MT08_P0.05_LFC1.0_up  = P9K_WT08_MT08_P0.05_up[which(P9K_WT08_MT08_P0.05_up$S_P9K_COEF_08D.W_08D.M>1.0),]
  P9K_WT08_MT08_P0.05_dn  = P9K_WT08_MT08_P0.05[which(P9K_WT08_MT08_P0.05$S_P9K_COEF_08D.W_08D.M<0),] 
    P9K_WT08_MT08_P0.05_LFC0.5_dn  = P9K_WT08_MT08_P0.05_dn[which(P9K_WT08_MT08_P0.05_dn$S_P9K_COEF_08D.W_08D.M<(-0.5)),] 
    P9K_WT08_MT08_P0.05_LFC1.0_dn  = P9K_WT08_MT08_P0.05_dn[which(P9K_WT08_MT08_P0.05_dn$S_P9K_COEF_08D.W_08D.M<(-1.0)),]

P9K_WT08_MT08_P0.01 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVAL_08D.W_08D.M)<0.01),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_08D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_08D.W_08D.M",colnames(P9K_Analysis)))]
  P9K_WT08_MT08_P0.01_up  = P9K_WT08_MT08_P0.01[which(P9K_WT08_MT08_P0.01$S_P9K_COEF_08D.W_08D.M>0),] 
    P9K_WT08_MT08_P0.01_LFC0.5_up  = P9K_WT08_MT08_P0.01_up[which(P9K_WT08_MT08_P0.01_up$S_P9K_COEF_08D.W_08D.M>0.5),] 
    P9K_WT08_MT08_P0.01_LFC1.0_up  = P9K_WT08_MT08_P0.01_up[which(P9K_WT08_MT08_P0.01_up$S_P9K_COEF_08D.W_08D.M>1.0),]
  P9K_WT08_MT08_P0.01_dn  = P9K_WT08_MT08_P0.01[which(P9K_WT08_MT08_P0.01$S_P9K_COEF_08D.W_08D.M<0),] 
    P9K_WT08_MT08_P0.01_LFC0.5_dn  = P9K_WT08_MT08_P0.01_dn[which(P9K_WT08_MT08_P0.01_dn$S_P9K_COEF_08D.W_08D.M<(-0.5)),] 
    P9K_WT08_MT08_P0.01_LFC1.0_dn  = P9K_WT08_MT08_P0.01_dn[which(P9K_WT08_MT08_P0.01_dn$S_P9K_COEF_08D.W_08D.M<(-1.0)),]

P9K_WT08_MT08_F0.05 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVFA_08D.W_08D.M)<0.05),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_08D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_08D.W_08D.M",colnames(P9K_Analysis)))]
  P9K_WT08_MT08_F0.05_up  = P9K_WT08_MT08_F0.05[which(P9K_WT08_MT08_F0.05$S_P9K_COEF_08D.W_08D.M>0),] 
    P9K_WT08_MT08_F0.05_LFC0.5_up  = P9K_WT08_MT08_F0.05_up[which(P9K_WT08_MT08_F0.05_up$S_P9K_COEF_08D.W_08D.M>0.5),] 
    P9K_WT08_MT08_F0.05_LFC1.0_up  = P9K_WT08_MT08_F0.05_up[which(P9K_WT08_MT08_F0.05_up$S_P9K_COEF_08D.W_08D.M>1.0),]
  P9K_WT08_MT08_F0.05_dn  = P9K_WT08_MT08_F0.05[which(P9K_WT08_MT08_F0.05$S_P9K_COEF_08D.W_08D.M<0),] 
    P9K_WT08_MT08_F0.05_LFC0.5_dn  = P9K_WT08_MT08_F0.05_dn[which(P9K_WT08_MT08_F0.05_dn$S_P9K_COEF_08D.W_08D.M<(-0.5)),] 
    P9K_WT08_MT08_F0.05_LFC1.0_dn  = P9K_WT08_MT08_F0.05_dn[which(P9K_WT08_MT08_F0.05_dn$S_P9K_COEF_08D.W_08D.M<(-1.0)),]

P9K_WT08_MT08_F0.01 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVFA_08D.W_08D.M)<0.01),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_08D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_08D.W_08D.M",colnames(P9K_Analysis)))]
  P9K_WT08_MT08_F0.01_up  = P9K_WT08_MT08_F0.01[which(P9K_WT08_MT08_F0.01$S_P9K_COEF_08D.W_08D.M>0),] 
    P9K_WT08_MT08_F0.01_LFC0.5_up  = P9K_WT08_MT08_F0.01_up[which(P9K_WT08_MT08_F0.01_up$S_P9K_COEF_08D.W_08D.M>0.5),]
    P9K_WT08_MT08_F0.01_LFC1.0_up  = P9K_WT08_MT08_F0.01_up[which(P9K_WT08_MT08_F0.01_up$S_P9K_COEF_08D.W_08D.M>1.0),]
  P9K_WT08_MT08_F0.01_dn  = P9K_WT08_MT08_F0.01[which(P9K_WT08_MT08_F0.01$S_P9K_COEF_08D.W_08D.M<0),] 
    P9K_WT08_MT08_F0.01_LFC0.5_dn  = P9K_WT08_MT08_F0.01_dn[which(P9K_WT08_MT08_F0.01_dn$S_P9K_COEF_08D.W_08D.M<(-0.5)),] 
    P9K_WT08_MT08_F0.01_LFC1.0_dn  = P9K_WT08_MT08_F0.01_dn[which(P9K_WT08_MT08_F0.01_dn$S_P9K_COEF_08D.W_08D.M<(-1.0)),]

P9K_WT08_MT08_H0.05 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVHA_08D.W_08D.M)<0.05),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_08D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_08D.W_08D.M",colnames(P9K_Analysis)))]
  P9K_WT08_MT08_H0.05_up  = P9K_WT08_MT08_H0.05[which(P9K_WT08_MT08_H0.05$S_P9K_COEF_08D.W_08D.M>0),] 
    P9K_WT08_MT08_H0.05_LFC0.5_up  = P9K_WT08_MT08_H0.05_up[which(P9K_WT08_MT08_H0.05_up$S_P9K_COEF_08D.W_08D.M>0.5),] 
    P9K_WT08_MT08_H0.05_LFC1.0_up  = P9K_WT08_MT08_H0.05_up[which(P9K_WT08_MT08_H0.05_up$S_P9K_COEF_08D.W_08D.M>1.0),]
  P9K_WT08_MT08_H0.05_dn  = P9K_WT08_MT08_H0.05[which(P9K_WT08_MT08_H0.05$S_P9K_COEF_08D.W_08D.M<0),] 
    P9K_WT08_MT08_H0.05_LFC0.5_dn  = P9K_WT08_MT08_H0.05_dn[which(P9K_WT08_MT08_H0.05_dn$S_P9K_COEF_08D.W_08D.M<(-0.5)),] 
    P9K_WT08_MT08_H0.05_LFC1.0_dn  = P9K_WT08_MT08_H0.05_dn[which(P9K_WT08_MT08_H0.05_dn$S_P9K_COEF_08D.W_08D.M<(-1.0)),]

P9K_WT08_MT08_H0.01 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVHA_08D.W_08D.M)<0.01),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_08D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_08D.W_08D.M",colnames(P9K_Analysis)))]
  P9K_WT08_MT08_H0.01_up  = P9K_WT08_MT08_H0.01[which(P9K_WT08_MT08_H0.01$S_P9K_COEF_08D.W_08D.M>0),] 
    P9K_WT08_MT08_H0.01_LFC0.5_up  = P9K_WT08_MT08_H0.01_up[which(P9K_WT08_MT08_H0.01_up$S_P9K_COEF_08D.W_08D.M>0.5),]
    P9K_WT08_MT08_H0.01_LFC1.0_up  = P9K_WT08_MT08_H0.01_up[which(P9K_WT08_MT08_H0.01_up$S_P9K_COEF_08D.W_08D.M>1.0),]
  P9K_WT08_MT08_H0.01_dn  = P9K_WT08_MT08_H0.01[which(P9K_WT08_MT08_H0.01$S_P9K_COEF_08D.W_08D.M<0),] 
    P9K_WT08_MT08_H0.01_LFC0.5_dn  = P9K_WT08_MT08_H0.01_dn[which(P9K_WT08_MT08_H0.01_dn$S_P9K_COEF_08D.W_08D.M<(-0.5)),] 
    P9K_WT08_MT08_H0.01_LFC1.0_dn  = P9K_WT08_MT08_H0.01_dn[which(P9K_WT08_MT08_H0.01_dn$S_P9K_COEF_08D.W_08D.M<(-1.0)),]

setwd("Results/DE_Gene_Lists/Day 8/")

write.csv(P9K_WT08_MT08_F0.01,file="P9K_WT08_MT08_F0.01.csv")
write.csv(P9K_WT08_MT08_H0.01,file="P9K_WT08_MT08_H0.01.csv")

write.csv(P9K_WT08_MT08_P0.05_LFC0.5_up,file="Pvalue/P05_LFC05_up.csv")
write.csv(P9K_WT08_MT08_P0.05_LFC0.5_dn,file="Pvalue/P05_LFC05_dn.csv")
write.csv(P9K_WT08_MT08_P0.05_LFC1.0_up,file="Pvalue/P05_LFC10_up.csv")
write.csv(P9K_WT08_MT08_P0.05_LFC1.0_dn,file="Pvalue/P05_LFC10_dn.csv")

write.csv(P9K_WT08_MT08_P0.01_LFC0.5_up,file="Pvalue/P01_LFC05_up.csv")
write.csv(P9K_WT08_MT08_P0.01_LFC0.5_dn,file="Pvalue/P01_LFC05_dn.csv")
write.csv(P9K_WT08_MT08_P0.01_LFC1.0_up,file="Pvalue/P01_LFC10_up.csv")
write.csv(P9K_WT08_MT08_P0.01_LFC1.0_dn,file="Pvalue/P01_LFC10_dn.csv")


write.csv(P9K_WT08_MT08_F0.05_LFC0.5_up,file="FDR/F05_LFC05_up.csv")
write.csv(P9K_WT08_MT08_F0.05_LFC0.5_dn,file="FDR/F05_LFC05_dn.csv")
write.csv(P9K_WT08_MT08_F0.05_LFC1.0_up,file="FDR/F05_LFC10_up.csv")
write.csv(P9K_WT08_MT08_F0.05_LFC1.0_dn,file="FDR/F05_LFC10_dn.csv")

write.csv(P9K_WT08_MT08_F0.01_LFC0.5_up,file="FDR/F01_LFC05_up.csv")
write.csv(P9K_WT08_MT08_F0.01_LFC0.5_dn,file="FDR/F01_LFC05_dn.csv")
write.csv(P9K_WT08_MT08_F0.01_LFC1.0_up,file="FDR/F01_LFC10_up.csv")
write.csv(P9K_WT08_MT08_F0.01_LFC1.0_dn,file="FDR/F01_LFC10_dn.csv")


write.csv(P9K_WT08_MT08_H0.05_LFC0.5_up,file="Holm/H05_LFC05_up.csv")
write.csv(P9K_WT08_MT08_H0.05_LFC0.5_dn,file="Holm/H05_LFC05_dn.csv")
write.csv(P9K_WT08_MT08_H0.05_LFC1.0_up,file="Holm/H05_LFC10_up.csv")
write.csv(P9K_WT08_MT08_H0.05_LFC1.0_dn,file="Holm/H05_LFC10_dn.csv")

write.csv(P9K_WT08_MT08_H0.01_LFC0.5_up,file="Holm/H01_LFC05_up.csv")
write.csv(P9K_WT08_MT08_H0.01_LFC0.5_dn,file="Holm/H01_LFC05_dn.csv")
write.csv(P9K_WT08_MT08_H0.01_LFC1.0_up,file="Holm/H01_LFC10_up.csv")
write.csv(P9K_WT08_MT08_H0.01_LFC1.0_dn,file="Holm/H01_LFC10_dn.csv")
setwd("..")


#################
#Day 12

P9K_WT12_MT12_P0.05 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVAL_12D.W_12D.M)<0.05),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_12D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_12D.W_12D.M",colnames(P9K_Analysis)))]
  P9K_WT12_MT12_P0.05_up  = P9K_WT12_MT12_P0.05[which(P9K_WT12_MT12_P0.05$S_P9K_COEF_12D.W_12D.M>0),] 
    P9K_WT12_MT12_P0.05_LFC0.5_up  = P9K_WT12_MT12_P0.05_up[which(P9K_WT12_MT12_P0.05_up$S_P9K_COEF_12D.W_12D.M>0.5),] 
    P9K_WT12_MT12_P0.05_LFC1.0_up  = P9K_WT12_MT12_P0.05_up[which(P9K_WT12_MT12_P0.05_up$S_P9K_COEF_12D.W_12D.M>1.0),]
  P9K_WT12_MT12_P0.05_dn  = P9K_WT12_MT12_P0.05[which(P9K_WT12_MT12_P0.05$S_P9K_COEF_12D.W_12D.M<0),] 
    P9K_WT12_MT12_P0.05_LFC0.5_dn  = P9K_WT12_MT12_P0.05_dn[which(P9K_WT12_MT12_P0.05_dn$S_P9K_COEF_12D.W_12D.M<(-0.5)),] 
    P9K_WT12_MT12_P0.05_LFC1.0_dn  = P9K_WT12_MT12_P0.05_dn[which(P9K_WT12_MT12_P0.05_dn$S_P9K_COEF_12D.W_12D.M<(-1.0)),]

P9K_WT12_MT12_P0.01 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVAL_12D.W_12D.M)<0.01),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_12D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_12D.W_12D.M",colnames(P9K_Analysis)))]
  P9K_WT12_MT12_P0.01_up  = P9K_WT12_MT12_P0.01[which(P9K_WT12_MT12_P0.01$S_P9K_COEF_12D.W_12D.M>0),] 
    P9K_WT12_MT12_P0.01_LFC0.5_up  = P9K_WT12_MT12_P0.01_up[which(P9K_WT12_MT12_P0.01_up$S_P9K_COEF_12D.W_12D.M>0.5),] 
    P9K_WT12_MT12_P0.01_LFC1.0_up  = P9K_WT12_MT12_P0.01_up[which(P9K_WT12_MT12_P0.01_up$S_P9K_COEF_12D.W_12D.M>1.0),]
  P9K_WT12_MT12_P0.01_dn  = P9K_WT12_MT12_P0.01[which(P9K_WT12_MT12_P0.01$S_P9K_COEF_12D.W_12D.M<0),] 
    P9K_WT12_MT12_P0.01_LFC0.5_dn  = P9K_WT12_MT12_P0.01_dn[which(P9K_WT12_MT12_P0.01_dn$S_P9K_COEF_12D.W_12D.M<(-0.5)),] 
    P9K_WT12_MT12_P0.01_LFC1.0_dn  = P9K_WT12_MT12_P0.01_dn[which(P9K_WT12_MT12_P0.01_dn$S_P9K_COEF_12D.W_12D.M<(-1.0)),]

P9K_WT12_MT12_F0.05 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVFA_12D.W_12D.M)<0.05),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_12D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_12D.W_12D.M",colnames(P9K_Analysis)))]
  P9K_WT12_MT12_F0.05_up  = P9K_WT12_MT12_F0.05[which(P9K_WT12_MT12_F0.05$S_P9K_COEF_12D.W_12D.M>0),] 
    P9K_WT12_MT12_F0.05_LFC0.5_up  = P9K_WT12_MT12_F0.05_up[which(P9K_WT12_MT12_F0.05_up$S_P9K_COEF_12D.W_12D.M>0.5),] 
    P9K_WT12_MT12_F0.05_LFC1.0_up  = P9K_WT12_MT12_F0.05_up[which(P9K_WT12_MT12_F0.05_up$S_P9K_COEF_12D.W_12D.M>1.0),]
  P9K_WT12_MT12_F0.05_dn  = P9K_WT12_MT12_F0.05[which(P9K_WT12_MT12_F0.05$S_P9K_COEF_12D.W_12D.M<0),] 
    P9K_WT12_MT12_F0.05_LFC0.5_dn  = P9K_WT12_MT12_F0.05_dn[which(P9K_WT12_MT12_F0.05_dn$S_P9K_COEF_12D.W_12D.M<(-0.5)),] 
    P9K_WT12_MT12_F0.05_LFC1.0_dn  = P9K_WT12_MT12_F0.05_dn[which(P9K_WT12_MT12_F0.05_dn$S_P9K_COEF_12D.W_12D.M<(-1.0)),]

P9K_WT12_MT12_F0.01 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVFA_12D.W_12D.M)<0.01),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_12D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_12D.W_12D.M",colnames(P9K_Analysis)))]
  P9K_WT12_MT12_F0.01_up  = P9K_WT12_MT12_F0.01[which(P9K_WT12_MT12_F0.01$S_P9K_COEF_12D.W_12D.M>0),] 
    P9K_WT12_MT12_F0.01_LFC0.5_up  = P9K_WT12_MT12_F0.01_up[which(P9K_WT12_MT12_F0.01_up$S_P9K_COEF_12D.W_12D.M>0.5),]
    P9K_WT12_MT12_F0.01_LFC1.0_up  = P9K_WT12_MT12_F0.01_up[which(P9K_WT12_MT12_F0.01_up$S_P9K_COEF_12D.W_12D.M>1.0),]
  P9K_WT12_MT12_F0.01_dn  = P9K_WT12_MT12_F0.01[which(P9K_WT12_MT12_F0.01$S_P9K_COEF_12D.W_12D.M<0),] 
    P9K_WT12_MT12_F0.01_LFC0.5_dn  = P9K_WT12_MT12_F0.01_dn[which(P9K_WT12_MT12_F0.01_dn$S_P9K_COEF_12D.W_12D.M<(-0.5)),] 
    P9K_WT12_MT12_F0.01_LFC1.0_dn  = P9K_WT12_MT12_F0.01_dn[which(P9K_WT12_MT12_F0.01_dn$S_P9K_COEF_12D.W_12D.M<(-1.0)),]

P9K_WT12_MT12_H0.05 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVHA_12D.W_12D.M)<0.05),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_12D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_12D.W_12D.M",colnames(P9K_Analysis)))]
  P9K_WT12_MT12_H0.05_up  = P9K_WT12_MT12_H0.05[which(P9K_WT12_MT12_H0.05$S_P9K_COEF_12D.W_12D.M>0),] 
    P9K_WT12_MT12_H0.05_LFC0.5_up  = P9K_WT12_MT12_H0.05_up[which(P9K_WT12_MT12_H0.05_up$S_P9K_COEF_12D.W_12D.M>0.5),] 
    P9K_WT12_MT12_H0.05_LFC1.0_up  = P9K_WT12_MT12_H0.05_up[which(P9K_WT12_MT12_H0.05_up$S_P9K_COEF_12D.W_12D.M>1.0),]
  P9K_WT12_MT12_H0.05_dn  = P9K_WT12_MT12_H0.05[which(P9K_WT12_MT12_H0.05$S_P9K_COEF_12D.W_12D.M<0),] 
    P9K_WT12_MT12_H0.05_LFC0.5_dn  = P9K_WT12_MT12_H0.05_dn[which(P9K_WT12_MT12_H0.05_dn$S_P9K_COEF_12D.W_12D.M<(-0.5)),] 
    P9K_WT12_MT12_H0.05_LFC1.0_dn  = P9K_WT12_MT12_H0.05_dn[which(P9K_WT12_MT12_H0.05_dn$S_P9K_COEF_12D.W_12D.M<(-1.0)),]

P9K_WT12_MT12_H0.01 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVHA_12D.W_12D.M)<0.01),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_12D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_12D.W_12D.M",colnames(P9K_Analysis)))]
  P9K_WT12_MT12_H0.01_up  = P9K_WT12_MT12_H0.01[which(P9K_WT12_MT12_H0.01$S_P9K_COEF_12D.W_12D.M>0),] 
    P9K_WT12_MT12_H0.01_LFC0.5_up  = P9K_WT12_MT12_H0.01_up[which(P9K_WT12_MT12_H0.01_up$S_P9K_COEF_12D.W_12D.M>0.5),]
    P9K_WT12_MT12_H0.01_LFC1.0_up  = P9K_WT12_MT12_H0.01_up[which(P9K_WT12_MT12_H0.01_up$S_P9K_COEF_12D.W_12D.M>1.0),]
  P9K_WT12_MT12_H0.01_dn  = P9K_WT12_MT12_H0.01[which(P9K_WT12_MT12_H0.01$S_P9K_COEF_12D.W_12D.M<0),] 
    P9K_WT12_MT12_H0.01_LFC0.5_dn  = P9K_WT12_MT12_H0.01_dn[which(P9K_WT12_MT12_H0.01_dn$S_P9K_COEF_12D.W_12D.M<(-0.5)),] 
    P9K_WT12_MT12_H0.01_LFC1.0_dn  = P9K_WT12_MT12_H0.01_dn[which(P9K_WT12_MT12_H0.01_dn$S_P9K_COEF_12D.W_12D.M<(-1.0)),]

setwd("Day 12/")

write.csv(P9K_WT12_MT12_F0.05,file="P9K_WT12_MT12_F0.05.csv")
write.csv(P9K_WT12_MT12_F0.01,file="P9K_WT12_MT12_F0.01.csv")
write.csv(P9K_WT12_MT12_H0.01,file="P9K_WT12_MT12_H0.01.csv")

write.csv(P9K_WT12_MT12_P0.05_LFC0.5_up,file="Pvalue/P05_LFC05_up.csv")
write.csv(P9K_WT12_MT12_P0.05_LFC0.5_dn,file="Pvalue/P05_LFC05_dn.csv")
write.csv(P9K_WT12_MT12_P0.05_LFC1.0_up,file="Pvalue/P05_LFC10_up.csv")
write.csv(P9K_WT12_MT12_P0.05_LFC1.0_dn,file="Pvalue/P05_LFC10_dn.csv")

write.csv(P9K_WT12_MT12_P0.01_LFC0.5_up,file="Pvalue/P01_LFC05_up.csv")
write.csv(P9K_WT12_MT12_P0.01_LFC0.5_dn,file="Pvalue/P01_LFC05_dn.csv")
write.csv(P9K_WT12_MT12_P0.01_LFC1.0_up,file="Pvalue/P01_LFC10_up.csv")
write.csv(P9K_WT12_MT12_P0.01_LFC1.0_dn,file="Pvalue/P01_LFC10_dn.csv")


write.csv(P9K_WT12_MT12_F0.05_LFC0.5_up,file="FDR/F05_LFC05_up.csv")
write.csv(P9K_WT12_MT12_F0.05_LFC0.5_dn,file="FDR/F05_LFC05_dn.csv")
write.csv(P9K_WT12_MT12_F0.05_LFC1.0_up,file="FDR/F05_LFC10_up.csv")
write.csv(P9K_WT12_MT12_F0.05_LFC1.0_dn,file="FDR/F05_LFC10_dn.csv")

write.csv(P9K_WT12_MT12_F0.01_LFC0.5_up,file="FDR/F01_LFC05_up.csv")
write.csv(P9K_WT12_MT12_F0.01_LFC0.5_dn,file="FDR/F01_LFC05_dn.csv")
write.csv(P9K_WT12_MT12_F0.01_LFC1.0_up,file="FDR/F01_LFC10_up.csv")
write.csv(P9K_WT12_MT12_F0.01_LFC1.0_dn,file="FDR/F01_LFC10_dn.csv")


write.csv(P9K_WT12_MT12_H0.05_LFC0.5_up,file="Holm/H05_LFC05_up.csv")
write.csv(P9K_WT12_MT12_H0.05_LFC0.5_dn,file="Holm/H05_LFC05_dn.csv")
write.csv(P9K_WT12_MT12_H0.05_LFC1.0_up,file="Holm/H05_LFC10_up.csv")
write.csv(P9K_WT12_MT12_H0.05_LFC1.0_dn,file="Holm/H05_LFC10_dn.csv")

write.csv(P9K_WT12_MT12_H0.01_LFC0.5_up,file="Holm/H01_LFC05_up.csv")
write.csv(P9K_WT12_MT12_H0.01_LFC0.5_dn,file="Holm/H01_LFC05_dn.csv")
write.csv(P9K_WT12_MT12_H0.01_LFC1.0_up,file="Holm/H01_LFC10_up.csv")
write.csv(P9K_WT12_MT12_H0.01_LFC1.0_dn,file="Holm/H01_LFC10_dn.csv")
setwd("..")


#################
#Day 16

P9K_WT16_MT16_P0.05 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVAL_16D.W_16D.M)<0.05),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_16D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_16D.W_16D.M",colnames(P9K_Analysis)))]
  P9K_WT16_MT16_P0.05_up  = P9K_WT16_MT16_P0.05[which(P9K_WT16_MT16_P0.05$S_P9K_COEF_16D.W_16D.M>0),] 
    P9K_WT16_MT16_P0.05_LFC0.5_up  = P9K_WT16_MT16_P0.05_up[which(P9K_WT16_MT16_P0.05_up$S_P9K_COEF_16D.W_16D.M>0.5),] 
    P9K_WT16_MT16_P0.05_LFC1.0_up  = P9K_WT16_MT16_P0.05_up[which(P9K_WT16_MT16_P0.05_up$S_P9K_COEF_16D.W_16D.M>1.0),]
  P9K_WT16_MT16_P0.05_dn  = P9K_WT16_MT16_P0.05[which(P9K_WT16_MT16_P0.05$S_P9K_COEF_16D.W_16D.M<0),] 
    P9K_WT16_MT16_P0.05_LFC0.5_dn  = P9K_WT16_MT16_P0.05_dn[which(P9K_WT16_MT16_P0.05_dn$S_P9K_COEF_16D.W_16D.M<(-0.5)),] 
    P9K_WT16_MT16_P0.05_LFC1.0_dn  = P9K_WT16_MT16_P0.05_dn[which(P9K_WT16_MT16_P0.05_dn$S_P9K_COEF_16D.W_16D.M<(-1.0)),]

P9K_WT16_MT16_P0.01 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVAL_16D.W_16D.M)<0.01),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_16D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_16D.W_16D.M",colnames(P9K_Analysis)))]
  P9K_WT16_MT16_P0.01_up  = P9K_WT16_MT16_P0.01[which(P9K_WT16_MT16_P0.01$S_P9K_COEF_16D.W_16D.M>0),] 
    P9K_WT16_MT16_P0.01_LFC0.5_up  = P9K_WT16_MT16_P0.01_up[which(P9K_WT16_MT16_P0.01_up$S_P9K_COEF_16D.W_16D.M>0.5),] 
    P9K_WT16_MT16_P0.01_LFC1.0_up  = P9K_WT16_MT16_P0.01_up[which(P9K_WT16_MT16_P0.01_up$S_P9K_COEF_16D.W_16D.M>1.0),]
  P9K_WT16_MT16_P0.01_dn  = P9K_WT16_MT16_P0.01[which(P9K_WT16_MT16_P0.01$S_P9K_COEF_16D.W_16D.M<0),] 
    P9K_WT16_MT16_P0.01_LFC0.5_dn  = P9K_WT16_MT16_P0.01_dn[which(P9K_WT16_MT16_P0.01_dn$S_P9K_COEF_16D.W_16D.M<(-0.5)),] 
    P9K_WT16_MT16_P0.01_LFC1.0_dn  = P9K_WT16_MT16_P0.01_dn[which(P9K_WT16_MT16_P0.01_dn$S_P9K_COEF_16D.W_16D.M<(-1.0)),]

P9K_WT16_MT16_F0.05 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVFA_16D.W_16D.M)<0.05),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_16D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_16D.W_16D.M",colnames(P9K_Analysis)))]
  P9K_WT16_MT16_F0.05_up  = P9K_WT16_MT16_F0.05[which(P9K_WT16_MT16_F0.05$S_P9K_COEF_16D.W_16D.M>0),] 
    P9K_WT16_MT16_F0.05_LFC0.5_up  = P9K_WT16_MT16_F0.05_up[which(P9K_WT16_MT16_F0.05_up$S_P9K_COEF_16D.W_16D.M>0.5),] 
    P9K_WT16_MT16_F0.05_LFC1.0_up  = P9K_WT16_MT16_F0.05_up[which(P9K_WT16_MT16_F0.05_up$S_P9K_COEF_16D.W_16D.M>1.0),]
  P9K_WT16_MT16_F0.05_dn  = P9K_WT16_MT16_F0.05[which(P9K_WT16_MT16_F0.05$S_P9K_COEF_16D.W_16D.M<0),] 
    P9K_WT16_MT16_F0.05_LFC0.5_dn  = P9K_WT16_MT16_F0.05_dn[which(P9K_WT16_MT16_F0.05_dn$S_P9K_COEF_16D.W_16D.M<(-0.5)),] 
    P9K_WT16_MT16_F0.05_LFC1.0_dn  = P9K_WT16_MT16_F0.05_dn[which(P9K_WT16_MT16_F0.05_dn$S_P9K_COEF_16D.W_16D.M<(-1.0)),]

P9K_WT16_MT16_F0.01 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVFA_16D.W_16D.M)<0.01),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_16D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_16D.W_16D.M",colnames(P9K_Analysis)))]
  P9K_WT16_MT16_F0.01_up  = P9K_WT16_MT16_F0.01[which(P9K_WT16_MT16_F0.01$S_P9K_COEF_16D.W_16D.M>0),] 
    P9K_WT16_MT16_F0.01_LFC0.5_up  = P9K_WT16_MT16_F0.01_up[which(P9K_WT16_MT16_F0.01_up$S_P9K_COEF_16D.W_16D.M>0.5),]
    P9K_WT16_MT16_F0.01_LFC1.0_up  = P9K_WT16_MT16_F0.01_up[which(P9K_WT16_MT16_F0.01_up$S_P9K_COEF_16D.W_16D.M>1.0),]
  P9K_WT16_MT16_F0.01_dn  = P9K_WT16_MT16_F0.01[which(P9K_WT16_MT16_F0.01$S_P9K_COEF_16D.W_16D.M<0),] 
    P9K_WT16_MT16_F0.01_LFC0.5_dn  = P9K_WT16_MT16_F0.01_dn[which(P9K_WT16_MT16_F0.01_dn$S_P9K_COEF_16D.W_16D.M<(-0.5)),] 
    P9K_WT16_MT16_F0.01_LFC1.0_dn  = P9K_WT16_MT16_F0.01_dn[which(P9K_WT16_MT16_F0.01_dn$S_P9K_COEF_16D.W_16D.M<(-1.0)),]
    P9K_WT16_MT16_F0.01_LFC0.5 = P9K_WT16_MT16_F0.01[which(abs(P9K_WT16_MT16_F0.01$S_P9K_COEF_16D.W_16D.M)>0.5),]
    
P9K_WT16_MT16_H0.05 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVHA_16D.W_16D.M)<0.05),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_16D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_16D.W_16D.M",colnames(P9K_Analysis)))]
  P9K_WT16_MT16_H0.05_up  = P9K_WT16_MT16_H0.05[which(P9K_WT16_MT16_H0.05$S_P9K_COEF_16D.W_16D.M>0),] 
    P9K_WT16_MT16_H0.05_LFC0.5_up  = P9K_WT16_MT16_H0.05_up[which(P9K_WT16_MT16_H0.05_up$S_P9K_COEF_16D.W_16D.M>0.5),] 
    P9K_WT16_MT16_H0.05_LFC1.0_up  = P9K_WT16_MT16_H0.05_up[which(P9K_WT16_MT16_H0.05_up$S_P9K_COEF_16D.W_16D.M>1.0),]
  P9K_WT16_MT16_H0.05_dn  = P9K_WT16_MT16_H0.05[which(P9K_WT16_MT16_H0.05$S_P9K_COEF_16D.W_16D.M<0),] 
    P9K_WT16_MT16_H0.05_LFC0.5_dn  = P9K_WT16_MT16_H0.05_dn[which(P9K_WT16_MT16_H0.05_dn$S_P9K_COEF_16D.W_16D.M<(-0.5)),] 
    P9K_WT16_MT16_H0.05_LFC1.0_dn  = P9K_WT16_MT16_H0.05_dn[which(P9K_WT16_MT16_H0.05_dn$S_P9K_COEF_16D.W_16D.M<(-1.0)),]

P9K_WT16_MT16_H0.01 = P9K_Analysis[which(as.numeric(P9K_Analysis$S_P9K_PVHA_16D.W_16D.M)<0.01),c(grep("D_[A-Z][A-Z_1:9][A-Z]_[A-Z]_16D",colnames(P9K_Analysis)),grep("I_",colnames(P9K_Analysis)),grep("S_RYT_",colnames(P9K_Analysis)),grep("S_P9K_[A-Z][A-Z][A-Z][A-Z]_16D.W_16D.M",colnames(P9K_Analysis)))]
  P9K_WT16_MT16_H0.01_up  = P9K_WT16_MT16_H0.01[which(P9K_WT16_MT16_H0.01$S_P9K_COEF_16D.W_16D.M>0),] 
    P9K_WT16_MT16_H0.01_LFC0.5_up  = P9K_WT16_MT16_H0.01_up[which(P9K_WT16_MT16_H0.01_up$S_P9K_COEF_16D.W_16D.M>0.5),]
    P9K_WT16_MT16_H0.01_LFC1.0_up  = P9K_WT16_MT16_H0.01_up[which(P9K_WT16_MT16_H0.01_up$S_P9K_COEF_16D.W_16D.M>1.0),]
  P9K_WT16_MT16_H0.01_dn  = P9K_WT16_MT16_H0.01[which(P9K_WT16_MT16_H0.01$S_P9K_COEF_16D.W_16D.M<0),] 
    P9K_WT16_MT16_H0.01_LFC0.5_dn  = P9K_WT16_MT16_H0.01_dn[which(P9K_WT16_MT16_H0.01_dn$S_P9K_COEF_16D.W_16D.M<(-0.5)),] 
    P9K_WT16_MT16_H0.01_LFC1.0_dn  = P9K_WT16_MT16_H0.01_dn[which(P9K_WT16_MT16_H0.01_dn$S_P9K_COEF_16D.W_16D.M<(-1.0)),]

setwd("Day 16/")

write.csv(P9K_WT16_MT16_F0.01,file="P9K_WT16_MT16_F0.01.csv")

write.csv(P9K_WT16_MT16_F0.05,file="Results/DE_Gene_Lists/Day 16/FDR/P9K_WT16_MT16_F0.05.csv")

write.csv(P9K_WT16_MT16_F0.01_LFC0.5,file="P9K_WT16_MT16_F0.01_LFC0.5.csv")
write.csv(P9K_WT16_MT16_H0.01,file="P9K_WT16_MT16_H0.01.csv")

write.csv(P9K_WT16_MT16_P0.05_LFC0.5_up,file="Pvalue/P05_LFC05_up.csv")
write.csv(P9K_WT16_MT16_P0.05_LFC0.5_dn,file="Pvalue/P05_LFC05_dn.csv")
write.csv(P9K_WT16_MT16_P0.05_LFC1.0_up,file="Pvalue/P05_LFC10_up.csv")
write.csv(P9K_WT16_MT16_P0.05_LFC1.0_dn,file="Pvalue/P05_LFC10_dn.csv")

write.csv(P9K_WT16_MT16_P0.01_LFC0.5_up,file="Pvalue/P01_LFC05_up.csv")
write.csv(P9K_WT16_MT16_P0.01_LFC0.5_dn,file="Pvalue/P01_LFC05_dn.csv")
write.csv(P9K_WT16_MT16_P0.01_LFC1.0_up,file="Pvalue/P01_LFC10_up.csv")
write.csv(P9K_WT16_MT16_P0.01_LFC1.0_dn,file="Pvalue/P01_LFC10_dn.csv")


write.csv(P9K_WT16_MT16_F0.05_LFC0.5_up,file="FDR/F05_LFC05_up.csv")
write.csv(P9K_WT16_MT16_F0.05_LFC0.5_dn,file="FDR/F05_LFC05_dn.csv")
write.csv(P9K_WT16_MT16_F0.05_LFC1.0_up,file="FDR/F05_LFC10_up.csv")
write.csv(P9K_WT16_MT16_F0.05_LFC1.0_dn,file="FDR/F05_LFC10_dn.csv")

write.csv(P9K_WT16_MT16_F0.01_LFC0.5_up,file="FDR/F01_LFC05_up.csv")
write.csv(P9K_WT16_MT16_F0.01_LFC0.5_dn,file="FDR/F01_LFC05_dn.csv")
write.csv(P9K_WT16_MT16_F0.01_LFC1.0_up,file="FDR/F01_LFC10_up.csv")
write.csv(P9K_WT16_MT16_F0.01_LFC1.0_dn,file="FDR/F01_LFC10_dn.csv")


write.csv(P9K_WT16_MT16_H0.05_LFC0.5_up,file="Holm/H05_LFC05_up.csv")
write.csv(P9K_WT16_MT16_H0.05_LFC0.5_dn,file="Holm/H05_LFC05_dn.csv")
write.csv(P9K_WT16_MT16_H0.05_LFC1.0_up,file="Holm/H05_LFC10_up.csv")
write.csv(P9K_WT16_MT16_H0.05_LFC1.0_dn,file="Holm/H05_LFC10_dn.csv")

write.csv(P9K_WT16_MT16_H0.01_LFC0.5_up,file="Holm/H01_LFC05_up.csv")
write.csv(P9K_WT16_MT16_H0.01_LFC0.5_dn,file="Holm/H01_LFC05_dn.csv")
write.csv(P9K_WT16_MT16_H0.01_LFC1.0_up,file="Holm/H01_LFC10_up.csv")
write.csv(P9K_WT16_MT16_H0.01_LFC1.0_dn,file="Holm/H01_LFC10_dn.csv")
setwd("..")


##################################################################################################################
############ Make volcano plots ##################################################################################
###### DE gene lists -> plots ####################################################################################

#Made matrix/vectors for the colors by which you want your volacno plots

#Day 8, commented lines are to select shape by p-value too
color8 = matrix(NA, ncol=2, nrow=nrow(P9K_Analysis))
for (i in 1:nrow(color8)){
  if (P9K_Analysis[i,"S_P9K_PVHA_08D.W_08D.M"] < 0.05){
    color8[i,1]= "darkcyan"
    color8[i,2]=19
  }
  else{
    color8[i,1]= "black"
    color8[i,2]=1
  }
}

#Plot Volcano plot
plot(as.numeric(P9K_Analysis$S_P9K_COEF_08D.W_08D.M),-log10(as.numeric(P9K_Analysis$S_P9K_PVAL_08D.W_08D.M)),col=color8[,1])
#plot(as.numeric(P9K_Analysis$S_P9K_COEF_08D.W_08D.M),-log10(as.numeric(P9K_Analysis$S_P9K_PVAL_08D.W_08D.M)),col=color8[,1],pch=as.numeric(color8[,2]))

#Day 12, commented lines are to select shape by p-value too
color12 = matrix(NA, ncol=2, nrow=nrow(P9K_Analysis))
for (i in 1:nrow(color12)){
  if (P9K_Analysis[i,"S_P9K_PVHA_12D.W_12D.M"] < 0.05){
    color12[i,1]= "darkcyan"
    color12[i,2]=19
  }
  else{
    color12[i,1]= "black"
    color12[i,2]=1
  }
}

#Plot Volcano plot
plot(as.numeric(P9K_Analysis$S_P9K_COEF_12D.W_12D.M),-log10(as.numeric(P9K_Analysis$S_P9K_PVAL_12D.W_12D.M)),col=color12[,1])
#plot(as.numeric(P9K_Analysis$S_P9K_COEF_12D.W_12D.M),-log10(as.numeric(P9K_Analysis$S_P9K_PVAL_12D.W_12D.M)),col=color12[,1],pch=as.numeric(color12[,2]))

#Day 16, commented lines are to select shape by p-value too
color16 = matrix(NA, ncol=2, nrow=nrow(P9K_Analysis))
for (i in 1:nrow(color16)){
  if (P9K_Analysis[i,"S_P9K_PVHA_16D.W_16D.M"] < 0.05){
    color16[i,1]= "darkcyan"
    color16[i,2]=19
  }
  else{
    color16[i,1]= "black"
    color16[i,2]=1
  }
}


#Plot Volcano plot
plot(as.numeric(P9K_Analysis$S_P9K_COEF_16D.W_16D.M),-log10(as.numeric(P9K_Analysis$S_P9K_PVAL_16D.W_16D.M)),col=color16[,1])
#plot(as.numeric(P9K_Analysis$S_P9K_COEF_12D.W_12D.M),-log10(as.numeric(P9K_Analysis$S_P9K_PVAL_12D.W_12D.M)),col=color12[,1],pch=as.numeric(color12[,2]))
plot(-as.numeric(P9K_Analysis$S_P9K_COEF_12D.W_16D.W),-log10(as.numeric(P9K_Analysis$S_P9K_PVAL_12D.W_16D.W)),col=color16[,1])


##################################################################################################################
############ Make overlap plots ##################################################################################
###### DE gene lists -> plots ####################################################################################

D08_up = rownames(P9K_WT08_MT08_H0.05_up)
D12_up = rownames(P9K_WT12_MT12_H0.05_up)
D16_up = rownames(P9K_WT16_MT16_H0.05_up)
length(D08_up)
length(D12_up)
length(D16_up)
length(which(D08_up %in% D12_up))
length(which(D08_up %in% D16_up))
length(which(D12_up %in% D16_up))
length(which(D08_up[which(D08_up %in% D12_up)] %in% D16_up))

dev.off()
draw.triple.venn(area1 = 52,
                 area2 = 240,
                 area3 = 1195,
                 n12 = 9,
                 n13 = 9,
                 n23 = 34,
                 n123 = 4,
                 #fill = c("#AA5585","#D4C26A","#407F7F"),
                 #col = "dark grey",
                 rotation.degree = 60)

D08_dn = rownames(P9K_WT08_MT08_H0.05_dn)
D12_dn = rownames(P9K_WT12_MT12_H0.05_dn)
D16_dn = rownames(P9K_WT16_MT16_H0.05_dn)
length(D08_dn)
length(D12_dn)
length(D16_dn)
length(which(D08_dn %in% D12_dn))
length(which(D08_dn %in% D16_dn))
length(which(D12_dn %in% D16_dn))
length(which(D08_dn[which(D08_dn %in% D12_dn)] %in% D16_dn))

dev.off()
draw.triple.venn(area1 = 47,
                 area2 = 236,
                 area3 = 1242,
                 n12 = 7,
                 n13 = 4,
                 n23 = 8,
                 n123 = 2,
                 #fill = c("#AA5585","#D4C26A","#407F7F"),
                 #col = "dark grey",
                 rotation.degree = 60)


P9K_Analysis[D12_up[which(D12_up %in% D16_up)],"S_RYT_SpecSub"]
D08_dn[which(D08_dn[which(D08_dn %in% D12_dn)] %in% D16_dn)]


##################################################################################################################
############ Make heatmap ########################################################################################
###### DE gene lists -> plots ####################################################################################

P9K_Matrix = as.matrix(P9DL_1)
LPD_Matrix = as.matrix(P9K_LPD[,1:41])
colnames(LPD_Matrix) = substr(colnames(LPD_Matrix),7,10)
LPD_Matrix_16 = LPD_Matrix[,grep("16",colnames(LPD_Matrix))]
LPD_Matrix_16_WM = LPD_Matrix_16[,-grep("H",colnames(LPD_Matrix_16))]


LPD_to_LLZ_Matrix = as.matrix(P9K_LPD_lt_to_LLZ[,1:41])
colnames(LPD_to_LLZ_Matrix) = substr(colnames(LPD_to_LLZ_Matrix),7,10)
LPD_to_LLZ_Matrix_16 = LPD_to_LLZ_Matrix[,grep("16",colnames(LPD_to_LLZ_Matrix))]
LPD_to_LLZ_Matrix_16_WM = LPD_to_LLZ_Matrix_16[,-grep("H",colnames(LPD_to_LLZ_Matrix_16))]




heatmap.2(P9K_Matrix,density.info="none",symm=F,symkey=F,symbreaks=F,trace="none")
heatmap.2(LPD_Matrix,density.info="none",symm=F,symkey=F,symbreaks=F,trace="none",srtCol=90)
heatmap.2(LPD_Matrix_16,density.info="none",symm=F,symkey=F,symbreaks=F,trace="none",srtCol=90)
heatmap.2(LPD_Matrix_16_WM,density.info="none",symm=F,symkey=F,symbreaks=F,trace="none",srtCol=90)
heatmap.2(LPD_to_LLZ_Matrix,density.info="none",symm=F,symkey=F,symbreaks=F,trace="none",srtCol=90)
heatmap.2(LPD_to_LLZ_Matrix_16,density.info="none",symm=F,symkey=F,symbreaks=F,trace="none",srtCol=90)
heatmap.2(LPD_to_LLZ_Matrix_16_WM,density.info="none",symm=F,symkey=F,symbreaks=F,trace="none",srtCol=80)



#,dendrogram="none",col=palette,density.info="none",
 #         Colv="NA",margins=c(12,9),notecol="black",Rowv=FALSE,vline=F,trace="none",
  #        xlab="",breaks=pal.breaks,
   #       symm=F,symkey=F,symbreaks=F,cexCol=2,
    #      labRow="",key=T,keysize=1)



# Dont bother reading this


# D12 = P9K_Analysis[which(P9K_Analysis$S_P9K_PVFA_12D.W_12D.M<0.01),]
# D12 = D12[which(D12$I_ENTZ!="NA"&D12$I_ENTZ!="Z"),]
# write.csv(D12,file="/Users/s-fine/Desktop/Carter/Projects/P9/3 DifferentialExpression/Results/DE_Gene_Lists/D12.csv")
# 
# sex_P9K_Analysis = P9K_Analysis[which(P9K_Analysis$I_CNUM=="X" | P9K_Analysis$I_CNUM=="Y"),]
# sex_P9K_Analysis_N = P9K_Analysis[-which(P9K_Analysis$I_CNUM=="X" | P9K_Analysis$I_CNUM=="Y" | P9K_Analysis$I_CNUM=="Z"),]
# 
# colorsex = matrix(NA, ncol=2, nrow=nrow(P9K_Analysis))
# for (i in 1:nrow(colorsex)){
#   if (P9K_Analysis[i,"I_CNUM"] == "X" | P9K_Analysis[i,"I_CNUM"] == "Y"){
#     colorsex[i,1]= "darkcyan"
#     colorsex[i,2]=19
#   }
#   else{
#     colorsex[i,1]= "black"
#     colorsex[i,2]=1
#   }
# }
# plot(as.numeric(P9K_Analysis$S_P9K_COEF_16D.W_16D.M),-log10(as.numeric(P9K_Analysis$S_P9K_PVAL_16D.W_16D.M)),col=colorsex[,1])
# vioplot(sex_P9K_Analysis_N$S_P9K_COEF_16D.W_16D.M,sex_P9K_Analysis$S_P9K_COEF_16D.W_16D.M)
# 
# 
# 
# 
# 
# 
# 
# 
# plot(sex_P9K_Analysis_N$S_P9K_COEF_16D.W_16D.M,-log10(sex_P9K_Analysis_N$S_P9K_PVAL_16D.W_16D.M),ylab=NA,xlab=NA,xlim=c(-7,3))
# plot(sex_P9K_Analysis$S_P9K_COEF_16D.W_16D.M,-log10(sex_P9K_Analysis$S_P9K_PVAL_16D.W_16D.M),ylab=NA,xlab=NA,xlim=c(-7,3))
# 
# 
# 
# plot(sex_P9K_Analysis_N$S_P9K_COEF_12D.W_12D.M,-log10(sex_P9K_Analysis_N$S_P9K_PVAL_12D.W_12D.M),ylab=NA,xlab=NA,ylim=c(0,17),xlim=c(-5,5))
# plot(sex_P9K_Analysis$S_P9K_COEF_12D.W_12D.M,-log10(sex_P9K_Analysis$S_P9K_PVAL_12D.W_12D.M),ylab=NA,xlab=NA,ylim=c(0,17),xlim=c(-5,5))
# 
# plot(sex_P9K_Analysis_N$S_P9K_COEF_08D.W_08D.M,-log10(sex_P9K_Analysis_N$S_P9K_PVAL_08D.W_08D.M),ylab=NA,xlab=NA,ylim=c(0,22),xlim=c(-3,4))
# plot(sex_P9K_Analysis$S_P9K_COEF_08D.W_08D.M,-log10(sex_P9K_Analysis$S_P9K_PVAL_08D.W_08D.M),ylab=NA,xlab=NA,ylim=c(0,22),xlim=c(-3,4))
# 
# 
# which(sex_P9K_Analysis$S_P9K_COEF_16D.W_16D.M<(-3))
# sex_P9K_Analysis[672,]
# 
# head(P9K_Analysis)
# #Upregulated genes at 12
# plotgeneP9("Rpl3") #ribosome, translation
# plotgeneP9("Tagap") #GTPase, G-protein
# plotgeneP9("Tcte3") #apoptosis in sperm/testis
# plotgeneP9("Anxa5") #linkted to Tcte3, related to apoptosis
# plotgeneP9("Rnaset2a") #ribonuclease, wiki says its linked to angiogenesis, linked to apoptosis
# plotgeneP9("Rab42") #GTPase
# plotgeneP9("Msx2") #RAS signaling (RAS/P21 are a complex), RAS are small GTPase, apoptosis regulatgor in limb
# plotgeneP9("Ang4") #angiogenesis
# plotgeneP9("Eif4ebp3") #translation
# plotgeneP9("Npw") #g-proteins
# plotgeneP9("Otud6a") #ubiquitination
# 
# 
# 
# #Downregulated genes at 12
# plotgeneP9("Atp6v0c") #ATPase
# plotgeneP9("Rnaset2b") #Weird that Rnaset2a is in the up...
# plotgeneP9("Rpl29") #ribosome, translation
# plotgeneP9("Tagap1") #TAGAP?!?
# plotgeneP9("Nuggc") #GTPase activity
# 
# plotgeneP9("Decr2") #GTPase activity
# 
# 
# D16 = P9K_Analysis[which(P9K_Analysis$S_P9K_PVHA_16D.W_16D.M<0.01),]
# D16 = D16[which(D16$I_ENTZ!="NA"&D16$I_ENTZ!="Z"),]
# write.csv(D16,file="/Users/s-fine/Desktop/Carter/Projects/P9/3 DifferentialExpression/Results/DE_Gene_Lists/D16.csv")
# 
# #Day16 Up
# 
# plotgeneP9("Srgap1") #GTPase activity
# plotgeneP9("Defb15") #completely gone to there, defense??
# plotgeneP9("Hoxd13") #completely gone to there
# plotgeneP9("Ang") #angiogenesis, ribonuclease?
# plotgeneP9("Insl3") #Related to fertlity
# plotgeneP9("Lbx2") #transcription factor
# plotgeneP9("Ecm1") #angiogenesis
# plotgeneP9("Lhcgr") #g-protein
# plotgeneP9("Rapgef3") #GTPase
# 
# plotgeneP9("Ldhc") #GTPase
# 
# 
# 
# plotgeneP9("Gm25767")
# 
# 
# 
# 
# plotgeneP9("Cdkn1a")
# 
# 
# 
# 
# 
