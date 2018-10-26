library(gplots)
library(ggplot2)
library(reshape2)


#Set up directory
setwd("/Users/s-fine/Desktop/Carter/Projects/P9/6 PMCA/")

#Load data
setwd("Results/")
load(file = "WtSub.rdt")
load(file = "UnAdj_MutSub.rdt")
load("Sp_KOGenes.rdt")
load("Pl_KOGenes.rdt")
load("EL_KOGenes.rdt")
load("LLZ_KOGenes.rdt")
load("PLike_KOGenes.rdt")
load("EP_KOGenes.rdt")
setwd("../")
load("/Users/s-fine/Desktop/Carter/Projects/P9/3 DifferentialExpression/Results/P9K_Analysis.rdt")

#Make WT lists
Sp_WTGenes = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g2[[1]]),]
Pl_WTGenes = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g2[[2]]),]
EL_WTGenes = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g2[[3]]),]
LLZ_WTGenes = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g2[[4]]),]
EP_WTGenes = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g2[[6]]),]
LPD_WTGenes = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g2[[8]]),]
PiRNA = P9K_Analysis[1:55,]

#Make plotting functions
substage_vioplot_wm = function(data_fine){
  ggplot(data_fine,aes(x=Age,y=value)) +
  geom_violin(aes(fill=factor(Genotype2)),position="dodge") +
  scale_fill_manual(values=c("blue","red","blue","red","blue","red")) +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
}

substage_boxplot_wm = function(data_fine){
  ggplot(data_fine,aes(x=Age,y=value)) +
    geom_boxplot(aes(position="dodge",fill=factor(Genotype2))) +
    scale_fill_manual(values=c("blue","red","blue","red","blue","red")) +
    theme_bw() +
    theme(legend.position="none") +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_line(colour="black"),
          panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
}

#MAke plotable data for piRNAs
PiRNAmeans = PiRNA[,grep("D_",colnames(PiRNA))]
PiRNAmeans$Wt_08 = rowMeans(PiRNAmeans[,grep("_W_08",colnames(PiRNAmeans))])
PiRNAmeans$Mt_08 = rowMeans(PiRNAmeans[,grep("_M_08",colnames(PiRNAmeans))])
PiRNAmeans$Wt_12 = rowMeans(PiRNAmeans[,grep("_W_12",colnames(PiRNAmeans))])
PiRNAmeans$Mt_12 = rowMeans(PiRNAmeans[,grep("_M_12",colnames(PiRNAmeans))])
PiRNAmeans$Wt_16 = rowMeans(PiRNAmeans[,grep("_W_16",colnames(PiRNAmeans))])
PiRNAmeans$Mt_16 = rowMeans(PiRNAmeans[,grep("_M_16",colnames(PiRNAmeans))])
PiRNA_z = t(scale(t(PiRNAmeans[,c(42:47)]),center = TRUE,scale = TRUE))
PiRNA_zm = melt(PiRNA_z)
PiRNA_zm$Age = substr(PiRNA_zm$Var2,4,5)
PiRNA_zm$Genotype = substr(PiRNA_zm$Var2,1,2)
PiRNA_zm$Genotype2 <- factor(PiRNA_zm$Genotype, c("Wt", "Mt"))


#Makde plotable data for SP
Sp_KOmeans = Sp_KOGenes[,grep("D_",colnames(Sp_KOGenes))]
Sp_KOmeans$Wt_08 = rowMeans(Sp_KOmeans[,grep("_W_08",colnames(Sp_KOmeans))])
Sp_KOmeans$Mt_08 = rowMeans(Sp_KOmeans[,grep("_M_08",colnames(Sp_KOmeans))])
Sp_KOmeans$Wt_12 = rowMeans(Sp_KOmeans[,grep("_W_12",colnames(Sp_KOmeans))])
Sp_KOmeans$Mt_12 = rowMeans(Sp_KOmeans[,grep("_M_12",colnames(Sp_KOmeans))])
Sp_KOmeans$Wt_16 = rowMeans(Sp_KOmeans[,grep("_W_16",colnames(Sp_KOmeans))])
Sp_KOmeans$Mt_16 = rowMeans(Sp_KOmeans[,grep("_M_16",colnames(Sp_KOmeans))])
Sp_KO_z = t(scale(t(Sp_KOmeans[,c(42:47)]),center = TRUE,scale = TRUE))
Sp_KO_zm = melt(Sp_KO_z)
Sp_KO_zm$Age = substr(Sp_KO_zm$Var2,4,5)
Sp_KO_zm$Genotype = substr(Sp_KO_zm$Var2,1,2)
Sp_KO_zm$Genotype2 <- factor(Sp_KO_zm$Genotype, c("Wt", "Mt"))

Sp_WTmeans = Sp_WTGenes[,grep("D_",colnames(Sp_WTGenes))]
Sp_WTmeans$Wt_08 = rowMeans(Sp_WTmeans[,grep("_W_08",colnames(Sp_WTmeans))])
Sp_WTmeans$Mt_08 = rowMeans(Sp_WTmeans[,grep("_M_08",colnames(Sp_WTmeans))])
Sp_WTmeans$Wt_12 = rowMeans(Sp_WTmeans[,grep("_W_12",colnames(Sp_WTmeans))])
Sp_WTmeans$Mt_12 = rowMeans(Sp_WTmeans[,grep("_M_12",colnames(Sp_WTmeans))])
Sp_WTmeans$Wt_16 = rowMeans(Sp_WTmeans[,grep("_W_16",colnames(Sp_WTmeans))])
Sp_WTmeans$Mt_16 = rowMeans(Sp_WTmeans[,grep("_M_16",colnames(Sp_WTmeans))])
Sp_WT_z = t(scale(t(Sp_WTmeans[,c(42:47)]),center = TRUE,scale = TRUE))
Sp_WT_zm = melt(Sp_WT_z)
Sp_WT_zm$Age = substr(Sp_WT_zm$Var2,4,5)
Sp_WT_zm$Genotype = substr(Sp_WT_zm$Var2,1,2)
Sp_WT_zm$Genotype2 <- factor(Sp_WT_zm$Genotype, c("Wt", "Mt"))

substage_vioplot_wm(Sp_KO_zm)
substage_vioplot_wm(Sp_WT_zm)

substage_boxplot_wm(Sp_KO_zm)
substage_boxplot_wm(Sp_WT_zm)



#Makde plotable data for Pl

Pl_KOmeans = Pl_KOGenes[,grep("D_",colnames(Pl_KOGenes))]
Pl_KOmeans$Wt_08 = rowMeans(Pl_KOmeans[,grep("_W_08",colnames(Pl_KOmeans))])
Pl_KOmeans$Mt_08 = rowMeans(Pl_KOmeans[,grep("_M_08",colnames(Pl_KOmeans))])
Pl_KOmeans$Wt_12 = rowMeans(Pl_KOmeans[,grep("_W_12",colnames(Pl_KOmeans))])
Pl_KOmeans$Mt_12 = rowMeans(Pl_KOmeans[,grep("_M_12",colnames(Pl_KOmeans))])
Pl_KOmeans$Wt_16 = rowMeans(Pl_KOmeans[,grep("_W_16",colnames(Pl_KOmeans))])
Pl_KOmeans$Mt_16 = rowMeans(Pl_KOmeans[,grep("_M_16",colnames(Pl_KOmeans))])
Pl_KO_z = t(scale(t(Pl_KOmeans[,c(42:47)]),center = TRUE,scale = TRUE))
Pl_KO_zm = melt(Pl_KO_z)
Pl_KO_zm$Age = substr(Pl_KO_zm$Var2,4,5)
Pl_KO_zm$Genotype = substr(Pl_KO_zm$Var2,1,2)
Pl_KO_zm$Genotype2 <- factor(Pl_KO_zm$Genotype, c("Wt", "Mt"))


Pl_WTmeans = Pl_WTGenes[,grep("D_",colnames(Pl_WTGenes))]
Pl_WTmeans$Wt_08 = rowMeans(Pl_WTmeans[,grep("_W_08",colnames(Pl_WTmeans))])
Pl_WTmeans$Mt_08 = rowMeans(Pl_WTmeans[,grep("_M_08",colnames(Pl_WTmeans))])
Pl_WTmeans$Wt_12 = rowMeans(Pl_WTmeans[,grep("_W_12",colnames(Pl_WTmeans))])
Pl_WTmeans$Mt_12 = rowMeans(Pl_WTmeans[,grep("_M_12",colnames(Pl_WTmeans))])
Pl_WTmeans$Wt_16 = rowMeans(Pl_WTmeans[,grep("_W_16",colnames(Pl_WTmeans))])
Pl_WTmeans$Mt_16 = rowMeans(Pl_WTmeans[,grep("_M_16",colnames(Pl_WTmeans))])
Pl_WT_z = t(scale(t(Pl_WTmeans[,c(42:47)]),center = TRUE,scale = TRUE))
Pl_WT_zm = melt(Pl_WT_z)
Pl_WT_zm$Age = substr(Pl_WT_zm$Var2,4,5)
Pl_WT_zm$Genotype = substr(Pl_WT_zm$Var2,1,2)
Pl_WT_zm$Genotype2 <- factor(Pl_WT_zm$Genotype, c("Wt", "Mt"))




#Makde plotable data for EL

EL_KOmeans = EL_KOGenes[,grep("D_",colnames(EL_KOGenes))]
EL_KOmeans$Wt_08 = rowMeans(EL_KOmeans[,grep("_W_08",colnames(EL_KOmeans))])
EL_KOmeans$Mt_08 = rowMeans(EL_KOmeans[,grep("_M_08",colnames(EL_KOmeans))])
EL_KOmeans$Wt_12 = rowMeans(EL_KOmeans[,grep("_W_12",colnames(EL_KOmeans))])
EL_KOmeans$Mt_12 = rowMeans(EL_KOmeans[,grep("_M_12",colnames(EL_KOmeans))])
EL_KOmeans$Wt_16 = rowMeans(EL_KOmeans[,grep("_W_16",colnames(EL_KOmeans))])
EL_KOmeans$Mt_16 = rowMeans(EL_KOmeans[,grep("_M_16",colnames(EL_KOmeans))])
EL_KO_z = t(scale(t(EL_KOmeans[,c(42:47)]),center = TRUE,scale = TRUE))
EL_KO_zm = melt(EL_KO_z)
EL_KO_zm$Age = substr(EL_KO_zm$Var2,4,5)
EL_KO_zm$Genotype = substr(EL_KO_zm$Var2,1,2)
EL_KO_zm$Genotype2 <- factor(EL_KO_zm$Genotype, c("Wt", "Mt"))


EL_WTmeans = EL_WTGenes[,grep("D_",colnames(EL_WTGenes))]
EL_WTmeans$Wt_08 = rowMeans(EL_WTmeans[,grep("_W_08",colnames(EL_WTmeans))])
EL_WTmeans$Mt_08 = rowMeans(EL_WTmeans[,grep("_M_08",colnames(EL_WTmeans))])
EL_WTmeans$Wt_12 = rowMeans(EL_WTmeans[,grep("_W_12",colnames(EL_WTmeans))])
EL_WTmeans$Mt_12 = rowMeans(EL_WTmeans[,grep("_M_12",colnames(EL_WTmeans))])
EL_WTmeans$Wt_16 = rowMeans(EL_WTmeans[,grep("_W_16",colnames(EL_WTmeans))])
EL_WTmeans$Mt_16 = rowMeans(EL_WTmeans[,grep("_M_16",colnames(EL_WTmeans))])
EL_WT_z = t(scale(t(EL_WTmeans[,c(42:47)]),center = TRUE,scale = TRUE))
EL_WT_zm = melt(EL_WT_z)
EL_WT_zm$Age = substr(EL_WT_zm$Var2,4,5)
EL_WT_zm$Genotype = substr(EL_WT_zm$Var2,1,2)
EL_WT_zm$Genotype2 <- factor(EL_WT_zm$Genotype, c("Wt", "Mt"))


#Makde plotable data for LLZ

LLZ_KOmeans = LLZ_KOGenes[,grep("D_",colnames(LLZ_KOGenes))]
LLZ_KOmeans$Wt_08 = rowMeans(LLZ_KOmeans[,grep("_W_08",colnames(LLZ_KOmeans))])
LLZ_KOmeans$Mt_08 = rowMeans(LLZ_KOmeans[,grep("_M_08",colnames(LLZ_KOmeans))])
LLZ_KOmeans$Wt_12 = rowMeans(LLZ_KOmeans[,grep("_W_12",colnames(LLZ_KOmeans))])
LLZ_KOmeans$Mt_12 = rowMeans(LLZ_KOmeans[,grep("_M_12",colnames(LLZ_KOmeans))])
LLZ_KOmeans$Wt_16 = rowMeans(LLZ_KOmeans[,grep("_W_16",colnames(LLZ_KOmeans))])
LLZ_KOmeans$Mt_16 = rowMeans(LLZ_KOmeans[,grep("_M_16",colnames(LLZ_KOmeans))])
LLZ_KO_z = t(scale(t(LLZ_KOmeans[,c(42:47)]),center = TRUE,scale = TRUE))
LLZ_KO_zm = melt(LLZ_KO_z)
LLZ_KO_zm$Age = substr(LLZ_KO_zm$Var2,4,5)
LLZ_KO_zm$Genotype = substr(LLZ_KO_zm$Var2,1,2)
LLZ_KO_zm$Genotype2 <- factor(LLZ_KO_zm$Genotype, c("Wt", "Mt"))


LLZ_WTmeans = LLZ_WTGenes[,grep("D_",colnames(LLZ_WTGenes))]
LLZ_WTmeans$Wt_08 = rowMeans(LLZ_WTmeans[,grep("_W_08",colnames(LLZ_WTmeans))])
LLZ_WTmeans$Mt_08 = rowMeans(LLZ_WTmeans[,grep("_M_08",colnames(LLZ_WTmeans))])
LLZ_WTmeans$Wt_12 = rowMeans(LLZ_WTmeans[,grep("_W_12",colnames(LLZ_WTmeans))])
LLZ_WTmeans$Mt_12 = rowMeans(LLZ_WTmeans[,grep("_M_12",colnames(LLZ_WTmeans))])
LLZ_WTmeans$Wt_16 = rowMeans(LLZ_WTmeans[,grep("_W_16",colnames(LLZ_WTmeans))])
LLZ_WTmeans$Mt_16 = rowMeans(LLZ_WTmeans[,grep("_M_16",colnames(LLZ_WTmeans))])
LLZ_WT_z = t(scale(t(LLZ_WTmeans[,c(42:47)]),center = TRUE,scale = TRUE))
LLZ_WT_zm = melt(LLZ_WT_z)
LLZ_WT_zm$Age = substr(LLZ_WT_zm$Var2,4,5)
LLZ_WT_zm$Genotype = substr(LLZ_WT_zm$Var2,1,2)
LLZ_WT_zm$Genotype2 <- factor(LLZ_WT_zm$Genotype, c("Wt", "Mt"))


#Makde plotable data for EP

EP_KOmeans = EP_KOGenes[,grep("D_",colnames(EP_KOGenes))]
EP_KOmeans$Wt_08 = rowMeans(EP_KOmeans[,grep("_W_08",colnames(EP_KOmeans))])
EP_KOmeans$Mt_08 = rowMeans(EP_KOmeans[,grep("_M_08",colnames(EP_KOmeans))])
EP_KOmeans$Wt_12 = rowMeans(EP_KOmeans[,grep("_W_12",colnames(EP_KOmeans))])
EP_KOmeans$Mt_12 = rowMeans(EP_KOmeans[,grep("_M_12",colnames(EP_KOmeans))])
EP_KOmeans$Wt_16 = rowMeans(EP_KOmeans[,grep("_W_16",colnames(EP_KOmeans))])
EP_KOmeans$Mt_16 = rowMeans(EP_KOmeans[,grep("_M_16",colnames(EP_KOmeans))])
EP_KO_z = t(scale(t(EP_KOmeans[,c(42:47)]),center = TRUE,scale = TRUE))
EP_KO_zm = melt(EP_KO_z)
EP_KO_zm$Age = substr(EP_KO_zm$Var2,4,5)
EP_KO_zm$Genotype = substr(EP_KO_zm$Var2,1,2)
EP_KO_zm$Genotype2 <- factor(EP_KO_zm$Genotype, c("Wt", "Mt"))


EP_WTmeans = EP_WTGenes[,grep("D_",colnames(EP_WTGenes))]
EP_WTmeans$Wt_08 = rowMeans(EP_WTmeans[,grep("_W_08",colnames(EP_WTmeans))])
EP_WTmeans$Mt_08 = rowMeans(EP_WTmeans[,grep("_M_08",colnames(EP_WTmeans))])
EP_WTmeans$Wt_12 = rowMeans(EP_WTmeans[,grep("_W_12",colnames(EP_WTmeans))])
EP_WTmeans$Mt_12 = rowMeans(EP_WTmeans[,grep("_M_12",colnames(EP_WTmeans))])
EP_WTmeans$Wt_16 = rowMeans(EP_WTmeans[,grep("_W_16",colnames(EP_WTmeans))])
EP_WTmeans$Mt_16 = rowMeans(EP_WTmeans[,grep("_M_16",colnames(EP_WTmeans))])
EP_WT_z = t(scale(t(EP_WTmeans[,c(42:47)]),center = TRUE,scale = TRUE))
EP_WT_zm = melt(EP_WT_z)
EP_WT_zm$Age = substr(EP_WT_zm$Var2,4,5)
EP_WT_zm$Genotype = substr(EP_WT_zm$Var2,1,2)
EP_WT_zm$Genotype2 <- factor(EP_WT_zm$Genotype, c("Wt", "Mt"))

#Makde plotable data for PLike

PLike_KOmeans = PLike_KOGenes[,grep("D_",colnames(PLike_KOGenes))]
PLike_KOmeans$Wt_08 = rowMeans(PLike_KOmeans[,grep("_W_08",colnames(PLike_KOmeans))])
PLike_KOmeans$Mt_08 = rowMeans(PLike_KOmeans[,grep("_M_08",colnames(PLike_KOmeans))])
PLike_KOmeans$Wt_12 = rowMeans(PLike_KOmeans[,grep("_W_12",colnames(PLike_KOmeans))])
PLike_KOmeans$Mt_12 = rowMeans(PLike_KOmeans[,grep("_M_12",colnames(PLike_KOmeans))])
PLike_KOmeans$Wt_16 = rowMeans(PLike_KOmeans[,grep("_W_16",colnames(PLike_KOmeans))])
PLike_KOmeans$Mt_16 = rowMeans(PLike_KOmeans[,grep("_M_16",colnames(PLike_KOmeans))])
PLike_KO_z = t(scale(t(PLike_KOmeans[,c(42:47)]),center = TRUE,scale = TRUE))
PLike_KO_zm = melt(PLike_KO_z)
PLike_KO_zm$Age = substr(PLike_KO_zm$Var2,4,5)
PLike_KO_zm$Genotype = substr(PLike_KO_zm$Var2,1,2)
PLike_KO_zm$Genotype2 <- factor(PLike_KO_zm$Genotype, c("Wt", "Mt"))

#Makde plotable data for LPD

LPD_WTmeans = LPD_WTGenes[,grep("D_",colnames(LPD_WTGenes))]
LPD_WTmeans$Wt_08 = rowMeans(LPD_WTmeans[,grep("_W_08",colnames(LPD_WTmeans))])
LPD_WTmeans$Mt_08 = rowMeans(LPD_WTmeans[,grep("_M_08",colnames(LPD_WTmeans))])
LPD_WTmeans$Wt_12 = rowMeans(LPD_WTmeans[,grep("_W_12",colnames(LPD_WTmeans))])
LPD_WTmeans$Mt_12 = rowMeans(LPD_WTmeans[,grep("_M_12",colnames(LPD_WTmeans))])
LPD_WTmeans$Wt_16 = rowMeans(LPD_WTmeans[,grep("_W_16",colnames(LPD_WTmeans))])
LPD_WTmeans$Mt_16 = rowMeans(LPD_WTmeans[,grep("_M_16",colnames(LPD_WTmeans))])
LPD_WT_z = t(scale(t(LPD_WTmeans[,c(42:47)]),center = TRUE,scale = TRUE))
LPD_WT_zm = melt(LPD_WT_z)
LPD_WT_zm$Age = substr(LPD_WT_zm$Var2,4,5)
LPD_WT_zm$Genotype = substr(LPD_WT_zm$Var2,1,2)
LPD_WT_zm$Genotype2 <- factor(LPD_WT_zm$Genotype, c("Wt", "Mt"))

#Makde plotable data for LPD

Shift_means = P9K_LPD_lt_to_LLZ[,grep("D_",colnames(P9K_LPD_lt_to_LLZ))]
Shift_means$Wt_08 = rowMeans(Shift_means[,grep("_W_08",colnames(Shift_means))])
Shift_means$Mt_08 = rowMeans(Shift_means[,grep("_M_08",colnames(Shift_means))])
Shift_means$Wt_12 = rowMeans(Shift_means[,grep("_W_12",colnames(Shift_means))])
Shift_means$Mt_12 = rowMeans(Shift_means[,grep("_M_12",colnames(Shift_means))])
Shift_means$Wt_16 = rowMeans(Shift_means[,grep("_W_16",colnames(Shift_means))])
Shift_means$Mt_16 = rowMeans(Shift_means[,grep("_M_16",colnames(Shift_means))])
Shift_z = t(scale(t(Shift_means[,c(42:47)]),center = TRUE,scale = TRUE))
Shift_zm = melt(Shift_z)
Shift_zm$Age = substr(Shift_zm$Var2,4,5)
Shift_zm$Genotype = substr(Shift_zm$Var2,1,2)
Shift_zm$Genotype2 <- factor(Shift_zm$Genotype, c("Wt", "Mt"))




substage_vioplot_wm(Shift_zm)
substage_vioplot_wm(LPD_WT_zm)






LPD_WTmeans = LPD_WTGenes[,grep("D_",colnames(LPD_WTGenes))]
LPD_WTmeans$Wt_08 = rowMeans(LPD_WTmeans[,grep("_W_08",colnames(LPD_WTmeans))])
LPD_WTmeans$Mt_08 = rowMeans(LPD_WTmeans[,grep("_M_08",colnames(LPD_WTmeans))])
LPD_WTmeans$Wt_12 = rowMeans(LPD_WTmeans[,grep("_W_12",colnames(LPD_WTmeans))])
LPD_WTmeans$Mt_12 = rowMeans(LPD_WTmeans[,grep("_M_12",colnames(LPD_WTmeans))])
LPD_WTmeans$Wt_16 = rowMeans(LPD_WTmeans[,grep("_W_16",colnames(LPD_WTmeans))])
LPD_WTmeans$Mt_16 = rowMeans(LPD_WTmeans[,grep("_M_16",colnames(LPD_WTmeans))])
LPD_WT_z = t(scale(t(LPD_WTmeans[,c(42:47)]),center = TRUE,scale = TRUE))
LPD_WT_zm = melt(LPD_WT_z)



ggplot(LPD_WT_zm,aes(x=Var2,y=value))+
  geom_boxplot()

heatmap.2(as.matrix(LPD_WTmeans[,c(grep("Wt_",colnames(LPD_WTmeans)),
                                  grep("Mt_",colnames(LPD_WTmeans)))]),
          trace = 'none',
           scale='row',
          # key=FALSE,
          Rowv = FALSE,
          Colv = FALSE,
          dendrogram = 'none',
          col = colorRampPalette(c('white', '#a6761d'))(20),
          labCol = FALSE,
          labRow = FALSE)


heatmap.2(as.matrix(Sp_KOmeans[,c(grep("Wt_",colnames(Sp_KOmeans)),
                                  grep("Mt_",colnames(Sp_KOmeans)))]),
          trace = 'none',
         # scale='row',
         # key=FALSE,
          Rowv = FALSE,
          Colv = FALSE,
          dendrogram = 'none',
          col = colorRampPalette(c('blue','white', '#666666'))(20),
          labCol = FALSE,
          labRow = FALSE)

heatmap.2(as.matrix(Sp_KOmeans[,c(grep("Wt_",colnames(Sp_KOmeans)),
                                  grep("Mt_",colnames(Sp_KOmeans)))]),
          trace = 'none',
          scale='row',
          # key=FALSE,
          Rowv = FALSE,
          Colv = FALSE,
          dendrogram = 'none',
          col = colorRampPalette(c('white', '#666666'))(20),
          #labCol = FALSE,
          labRow = FALSE)

heatmap.2(as.matrix(Sp_KOmeans[,c(grep("Wt_",colnames(Sp_KOmeans)),
                                  grep("Mt_",colnames(Sp_KOmeans)))]),
          trace = 'none',
          # scale='row',
          # key=FALSE,
          Rowv = FALSE,
          Colv = FALSE,
          dendrogram = 'none',
          col = colorRampPalette(c('blue','white', '#666666'))(20),
          labCol = FALSE,
          labRow = FALSE)

heatmap.2(as.matrix(Sp_KOmeans[,c(grep("Wt_",colnames(Sp_KOmeans)),
                                  grep("Mt_",colnames(Sp_KOmeans)))]),
          trace = 'none',
          scale='row',
          # key=FALSE,
          Rowv = FALSE,
          Colv = FALSE,
          dendrogram = 'none',
          col = colorRampPalette(c('white', '#666666'))(20),
          #labCol = FALSE,
          labRow = FALSE)


heatmap.2(as.matrix(Pl_KOmeans[,c(grep("Wt_",colnames(Pl_KOmeans)),
                                  grep("Mt_",colnames(Pl_KOmeans)))]),
          trace = 'none',
          scale='row',
          key=FALSE,
          Rowv = FALSE,
          Colv = FALSE,
          dendrogram = 'none',
          col = colorRampPalette(c('white', '#1b9e77'))(20),
          labCol = FALSE,
          labRow = FALSE)





heatmap.2(as.matrix(Sp_KOGenes[,c(grep("_W_",colnames(Sp_KOGenes)),
                                  grep("_M_",colnames(Sp_KOGenes)))]),
          trace = 'none',
          scale='row',
          key=FALSE,
          Rowv = FALSE,
          Colv = FALSE,
          dendrogram = 'none',
          col = colorRampPalette(c('white', '#666666'))(20),
          labCol = FALSE,
          labRow = FALSE)
heatmap(as.matrix(Sp_KOGenes[,grep("_M_",colnames(Sp_KOGenes))]))

heatmap.2(as.matrix(Sp_KOGenes[,grep("_W_",colnames(Sp_KOGenes))]),trace = 'none',scale='row',Rowv = FALSE,Colv = FALSE)


Sp_KOGenes[1:5,grep("_M_",colnames(Sp_KOGenes))]






#subset substage specific genes
P9K_LPD = P9K_Analysis[which(P9K_Analysis$S_RYT_SpecSub=="LPD"),]
P9K_EP = P9K_Analysis[which(P9K_Analysis$S_RYT_SpecSub=="EP"),]
P9K_LPD_lt = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g2[[8]]),]
P9K_EP_lt = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g2[[6]]),]
P9K_LLZ = P9K_Analysis[which(P9K_Analysis$S_RYT_SpecSub=="LLZ"),]
P9K_LLZ_lt = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g2[[4]]),]

which(g_H[[1]][[3]] %in% g2[[1]])
which(g_H[[4]][[3]] %in% g2[[4]])
which(g_H[[4]][[3]] %in% g2[[6]])
which(g_H[[4]][[3]] %in% g2[[8]])

g_H[[4]][[3]]

Sp_KOG_H 
Pl_KOG_H
EL_KOG_H 
LLZ_KOG_H
PLike_KOG_H
EP_KOG_H





#subset substage specific gene changes
P9K_LPD_to_LLZ = P9K_LPD[which(P9K_LPD$I_GEID %in% g[[4]][[4]]),]
P9K_LPD_lt_to_LLZ = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g[[4]][[4]][which(g[[4]][[4]] %in% g2[[8]])]),]
P9K_EP_to_LLZ = P9K_EP[which(P9K_EP$I_GEID %in% g[[4]][[4]]),]
P9K_EP_lt_to_LLZ = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g[[4]][[4]][which(g[[4]][[4]] %in% g2[[6]])]),]
P9K_LLZ_to_LLZ = P9K_LLZ[which(P9K_LLZ$I_GEID %in% g[[4]][[4]]),]
P9K_LLZ_lt_to_LLZ = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g[[4]][[4]][which(g[[4]][[4]] %in% g2[[4]])]),]

#Subset data by significance in change and direction
P9K_LPD_lt_16NC = P9K_LPD_lt[which(P9K_LPD_lt$S_P9K_PVFA_16D.W_16D.M>0.05),]
P9K_LPD_lt_16UP = P9K_LPD_lt[which(P9K_LPD_lt$S_P9K_PVFA_16D.W_16D.M<0.05 & P9K_LPD_lt$S_P9K_COEF_16D.W_16D.M > 0),]
P9K_LPD_lt_16DN = P9K_LPD_lt[which(P9K_LPD_lt$S_P9K_PVHA_16D.W_16D.M<0.05 & P9K_LPD_lt$S_P9K_COEF_16D.W_16D.M < 0),]


#Save data
setwd("Results/")
setwd("Lists_For_GO//")
write.csv(P9K_Analysis,file="P9K_Analysis.csv")
write.csv(P9K_LPD,file="P9K_LPD.csv")
write.csv(P9K_EP,file="P9K_EP.csv")
write.csv(P9K_LLZ,file="P9K_LLZ.csv")
write.csv(P9K_LPD_lt,file="P9K_LPD_lt.csv")
write.csv(P9K_EP_lt,file="P9K_EP_lt.csv")
write.csv(P9K_LLZ_lt,file="P9K_LLZ_lt.csv")

write.csv(P9K_LPD_to_LLZ,file="P9K_LPD_to_LLZ.csv")
write.csv(P9K_LPD_lt_to_LLZ,file="P9K_LPD_lt_to_LLZ.csv")
write.csv(P9K_EP_to_LLZ,file="P9K_EP_to_LLZ.csv")
write.csv(P9K_EP_lt_to_LLZ,file="P9K_EP_lt_to_LLZ.csv")
write.csv(P9K_LLZ_to_LLZ,file="P9K_LLZ_to_LLZ.csv")
write.csv(P9K_LLZ_lt_to_LLZ,file="P9K_LLZ_lt_to_LLZ.csv")

write.csv(P9K_LPD_lt_16NC,file="P9K_LPD_lt_16NC.csv")
write.csv(P9K_LPD_lt_16UP,file="P9K_LPD_lt_16UP.csv")
write.csv(P9K_LPD_lt_16DN,file="P9K_LPD_lt_16DN.csv")


setwd("../")
setwd("../")







plot(P9K_LPD$S_P9K_COEF_16D.W_16D.M,-log(P9K_LPD$S_P9K_PVAL_16D.W_16D.M))

plot(P9K_EP$S_P9K_COEF_16D.W_16D.M,-log(P9K_EP$S_P9K_PVAL_16D.W_16D.M))


P9K_LPD_Up = P9K_LPD[which(P9K_LPD$S_P9K_COEF_16D.W_16D.M>=0),]

plot(P9K_LPD_Up$S_P9K_COEF_16D.W_16D.M,-log(P9K_LPD_Up$S_P9K_PVAL_16D.W_16D.M))

g[[4]][[4]][which(g[[4]][[4]] %in% g2[[8]])]
g[[4]][[4]][which(g[[4]][[4]] %in% g2[[6]])]

P9K_LPD_to_LLZ = P9K_LPD[which(P9K_LPD$I_GEID %in% g[[4]][[4]]),]
P9K_LPD_to_LLZ = P9K_LPD[which(P9K_LPD$I_GEID %in% g[[4]][[4]]),]
P9K_EP_to_LLZ = P9K_EP[which(P9K_EP$I_GEID %in% g[[4]][[4]]),]


P9K_EP_to_LLZ_2 = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g[[4]][[4]][which(g[[4]][[4]] %in% g2[[6]])]),]

P9K_LPD_to_LLZ_2 = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g[[4]][[4]][which(g[[4]][[4]] %in% g2[[8]])]),]


save(P9K_LPD_to_LLZ_2,file="/Users/s-fine/Desktop/Carter/Projects/P9/Substage_Specificity/Results/P9K_LPD_to_LLZ_2.rdt")







plot(P9K_LPD$S_P9K_COEF_16D.W_16D.M,-log(P9K_LPD$S_P9K_PVAL_16D.W_16D.M))
points(P9K_LPD_to_LLZ$S_P9K_COEF_16D.W_16D.M,-log(P9K_LPD_to_LLZ$S_P9K_PVAL_16D.W_16D.M),col="red")


plot(P9K_EP_to_LLZ_2$S_P9K_COEF_16D.W_16D.M,-log(P9K_EP_to_LLZ_2$S_P9K_PVAL_16D.W_16D.M))
plot(P9K_LPD_to_LLZ_2$S_P9K_COEF_16D.W_16D.M,-log(P9K_LPD_to_LLZ_2$S_P9K_PVAL_16D.W_16D.M))


plot(P9K_LPD_to_LLZ$S_P9K_COEF_16D.W_16D.M,-log(P9K_LPD_to_LLZ$S_P9K_PVAL_16D.W_16D.M))

P9K_LPD_to_LLZ["C2cd4b",]


P9K_Analysis[which(P9K_Analysis$S_P9K_COEF_12D.W_16D.W>1 & P9K_Analysis$S_P9K_PVFA_16D.W_16D.M>(0.01) & P9K_Analysis$I_GEID %in% g2[[8]]),]
P9K_Analysis[which(P9K_Analysis$S_P9K_COEF_16D.W_16D.M>0),]

plotgeneP9("Prdm9","Holm")


length(which(P9K_LPD$S_P9K_PVHA_16D.W_16D.M<0.05 & P9K_LPD$S_P9K_COEF_16D.W_16D.M<0))*100/nrow(P9K_LPD)

length(which(P9K_EP$S_P9K_PVHA_16D.W_16D.M<0.05 & P9K_EP$S_P9K_COEF_16D.W_16D.M<0))*100/nrow(P9K_EP)

length(which(P9K_EP_to_LLZ_2$S_P9K_PVHA_16D.W_16D.M<0.05 & P9K_EP_to_LLZ_2$S_P9K_COEF_16D.W_16D.M<0))*100/nrow(P9K_EP_to_LLZ_2)

length(which(P9K_LPD_to_LLZ_2$S_P9K_PVHA_16D.W_16D.M<0.05 & P9K_LPD_to_LLZ_2$S_P9K_COEF_16D.W_16D.M<0))*100/nrow(P9K_LPD_to_LLZ_2)





