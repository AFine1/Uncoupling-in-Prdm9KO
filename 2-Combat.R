#   The purpose of this file is to use ComBat to normalize systematic variance between samples

#   It takes data files that contain RNAseq TPM from Baseline and P9 datasets, merges them, and 
#   corrects between datasets and litters. PCA plots are saved for sample groups and overall 
#   changes, and one output file is created: P9DL_1, locally saved at 
#   Desktop/Carter/Projects/P9/Normalization/Results/P9DL_gene.rdt; however, that file feeds into
#   Differential_Expression.R to identify differentially expressed genes and create the main data file.



##################################################################################################################
############ Load packages #######################################################################################
##################################################################################################################

library("sva")
library("contrast")
library("Biobase")
library("biomaRt")
library(scatterplot3d)
attach(mtcars)
library(Rcmdr)
library(rgl)


##################################################################################################################
############ Aquire Data #########################################################################################
###### AllData ###################################################################################################

#Set working directory
setwd("/Users/s-fine/Desktop/Carter/Projects/P9/2 Normalization/")

# #Load Robyn's Dataset
# setwd("Data")
# RobynTable=read.table("all_samples_gene_expression_RSEM_TPM_GRCm38-75_lnc-piRNA.txt", header = T)
# setwd("..")
# 
# #Change Column Names
# colnames(RobynTable)=c("gene_id","R_08_1","R_08_2","R_08_3","R_08_4","R_08_5","R_10_1","R_10_2","R_10_3","R_10_4","R_10_5","R_12_1","R_12_2","R_12_3","R_12_4","R_12_5","R_14_1","R_14_2","R_14_3","R_14_6","R_14_7","R_16_1","R_16_2","R_16_3","R_16_4","R_16_5","R_18_1","R_18_2","R_18_3","R_18_4","R_18_5")
# RobynTable = RobynTable[,-4]
# rownames(RobynTable)=RobynTable[,1]
# 
# #Change piRNA nomenclature to be the same
# rownames(MutantData[grep("pi-",rownames(MutantData)),])
# 
# #Load P9 and WT Data
# setwd("Data")
# load("P9MData.rdt")
# load("WTData.rdt")
# setwd("..")
# 
# #Identify the overlapping genes between the two
# MutantDataOverlap=MutantData[c(1:85430),]
# 
# #Change piRNA nomenclature to be the same
# for (i in grep("pi-",rownames(MutantDataOverlap))){
#   rownames(MutantDataOverlap)[i]=rownames(RobynTable)[i]
# }
# 
# #Save file of WT data
# save(RobynTable,file="/Users/s-fine/Desktop/Carter/Projects/P9/Normalization/Data/WTData.rdt")
# 
# #Take important columns from WT data
# WTData=RobynTable[,c(1,2,3,9,10,11,12,13,19,20,21,22,23)]
# 
# #Merge important columns from both
# AllData=cbind(MutantDataOverlap,WTData)
# 
# #Save data
# save(AllData,file="/Users/s-fine/Desktop/Carter/Projects/P9/Normalization/Data/AllData.rdt")


# Load Data
setwd("Data")
load("AllData.rdt")
setwd("..")


##################################################################################################################
############ Re-Format Data ######################################################################################
###### AllData  -> P9_L2 #########################################################################################

#log transform Data
P9_L2 = log2(AllData+1)

#Limit to important data
P9_L2i = P9_L2[,c(1:14,16:19,21:43)]


##################################################################################################################
############ Display data with PCA ###############################################################################
###### P9_L2i  -> plot ###########################################################################################

#Old PCA with label
P9_L2i_var=P9_L2i[which(apply(P9_L2i, 1, var) != 0),]
original_model = prcomp(t(P9_L2i_var),scale.=T)
pca_o <- predict(original_model)
pca_o_1 <- predict(original_model)[,1]
pca_o_2 <- predict(original_model)[,2]
pca_o_3 <- predict(original_model)[,3]
colors_o = c("blue","blue","blue","purple","purple","red","red","red","red","red","blue","purple","purple","purple","red","red","red","red","blue","blue","purple","purple","purple","red","red","red","red","red","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue")
symbols_o = as.numeric(c(16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,15,15,15,15,15,15,15,15,15,15,1,1,1,2,2,2,2,2,0,0,0,0,0))
plot(pca_o_1,pca_o_2,col=colors_o,pch=symbols_o,main="Original RNAseq PCA",xlab="PC1",ylab="PC2",xlim=c(-225,290))
legend(200,-55,title = "Age",c("Day 8","Day 12","Day 16"),pch=c(15,17,16),col=1)
legend(200,10,title = "Genotypes ",c("WT     ","Het","Mut"),pch=18,col=c("blue","purple","red"))
legend(199.5,63,title = "Dataset ",c("Baseline","Prdm9"),pch=c(23,18),col=c("black"))

#Plot PCA without label
plot(pca_o_1,pca_o_2,col=colors_o,pch=symbols_o,main="Original RNAseq PCA",xlab="PC1",ylab="PC2",xlim=c(-225,225))

#Plot PCA without axis adjustment
plot(pca_o_1,pca_o_2,col=colors_o,pch=symbols_o,main="Original RNAseq PCA",xlab="PC1",ylab="PC2")


##################################################################################################################
############ Run COMBAT ##########################################################################################
###### P9_L2i  -> P9_all_adj_a ###################################################################################

#Correct for dataset, accounting for genotype and age
P9_L2i_var=P9_L2i[which(apply(P9_L2i, 1, var) != 0),]
P9_L2i_var_m = as.matrix(P9_L2i_var)
batch_all = c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","2","2","2","2","2","2","2","2","2","2","2","2","2")
#batch_alli = c("3","8","8","10","14","3","8","8","10","14","11","1","9","12","1","9","11","12","4","5","2","6","7","2","4","5","6","7","20","20","20","30","30","30","30","30","40","40","40","40","40")
mod_all = matrix(NA, ncol=2,nrow=41)
colnames(mod_all)=c("genotype","age")
rownames(mod_all)=colnames(P9_L2i_var)
mod_all[c(1,2,3,11,19,20,29:41),1]="A"
mod_all[c(4,5,12,13,14,21,22,23),1]="B"
mod_all[c(6,7,8,9,10,15,16,17,18,24,25,26,27,28),1]="C"
mod_all[c(1:10,29,30,31),2]="D"
mod_all[c(11:18,32,33,34,35,36),2]="E"
mod_all[c(19:28,37,38,39,40,41),2]="F"
P9_L2i_var_adj = ComBat(P9_L2i_var,as.factor(batch_all),mod=as.factor(mod_all),par.prior = TRUE)
#P9_L2i_var_adji = ComBat(P9_L2i_var,as.factor(batch_alli),mod=as.factor(mod_all),par.prior = TRUE)

#Plot PCA Result
P9_L2i_var_adj_var=P9_L2i_var_adj[which(apply(P9_L2i_var_adj, 1, var) != 0),]
dataset_model = prcomp(t(P9_L2i_var_adj_var),scale.=T)
pca_d <- predict(dataset_model)
pca_d_1 <- predict(dataset_model)[,1]
pca_d_2 <- predict(dataset_model)[,2]
pca_d_3 <- predict(dataset_model)[,3]
colors_d = c("blue","blue","blue","purple","purple","red","red","red","red","red","blue","purple","purple","purple","red","red","red","red","blue","blue","purple","purple","purple","red","red","red","red","red","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue")
symbols_d = as.numeric(c(16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,15,15,15,15,15,15,15,15,15,15,1,1,1,2,2,2,2,2,0,0,0,0,0))
plot(pca_d_1,pca_d_2,col=colors_d,pch=symbols_d,main="New RNAseq PCA",xlab="PC1",ylab="PC2")
plot(pca_d_1,pca_d_2,col=colors_d,pch=symbols_d,main="New RNAseq PCA",xlab="PC1",ylab="PC2",xlim=c(-225,225))

#################
#Day 8 
P9_L2i_var_adj_8 = P9_L2i_var_adj[,c(1:10,29,30,31)]
P9_L2i_var_adj_8 = as.matrix(P9_L2i_var_adj_8)

#Plot Day 8 PCA
P9_L2i_var_adj_8_var=P9_L2i_var_adj_8[which(apply(P9_L2i_var_adj_8, 1, var) != 0),]
D8_model = prcomp(t(P9_L2i_var_adj_8_var),scale.=T)
pca_D8 <- predict(D8_model)
pca_D8_1 <- predict(D8_model)[,1]
pca_D8_2 <- predict(D8_model)[,2]
colors_D8 = c(4,9,9,11,14,4,9,9,11,14,"red","red","red")
colors_D8 = c("blue","blue","blue","purple","purple","red","red","red","red","red","blue","blue","blue")
plot(pca_D8_1,pca_D8_2,col=colors_D8,pch=as.numeric(c(16,16,16,16,16,16,16,16,16,16,21,21,21)),main="Day 8 All RNAseq PCA",xlab="PC1",ylab="PC2")
text(pca_D8_1,pca_D8_2,pos=4,c(1,2,2,3,4,1,2,2,3,4,"B","B","B"))

#Adjust Day 8 for Litter
P9_8_all_n = P9_L2i_var_adj_8[which(rowMax(P9_L2i_var_adj_8)!=0),]
batch_8_all = c("3","8","8","10","14","3","8","8","10","14","30","30","30")
mod_8_all = matrix(NA, ncol=1,nrow=13)
colnames(mod_8_all)=c("genotype")
rownames(mod_8_all)=colnames(P9_8_all_n)
mod_8_all[c(1:3,11,12,13),1]="A"
mod_8_all[c(4,5),1]="B"
mod_8_all[c(6:10),1]="C"
P9_8_all_adj = ComBat(P9_8_all_n,batch_8_all,mod=as.factor(mod_8_all),par.prior = TRUE)
P9_8_all_adj_a = rbind(P9_8_all_adj,P9_8_all_adj[which(rowMax(P9_8_all_adj)==0),])

#Plot PCA Result
P9_8_all_adj_var=P9_8_all_adj_a[which(apply(P9_8_all_adj_a, 1, var) != 0),]
D8_all_model = prcomp(t(P9_8_all_adj_var),scale.=T)
pca_D8_all <- predict(D8_all_model)
pca_D8_all_1 <- predict(D8_all_model)[,1]
pca_D8_all_2 <- predict(D8_all_model)[,2]
colors_D8_all = c("blue","blue","blue","purple","purple","red","red","red","red","red","blue","blue","blue")
plot(pca_D8_all_1,pca_D8_all_2,col=colors_D8_all,pch=as.numeric(c(16,16,16,16,16,16,16,16,16,16,21,21,21)),main="Day 8 All RNAseq PCA",xlab="PC1",ylab="PC2")
text(pca_D8_all_1,pca_D8_all_2,pos=4,c(1,2,2,3,4,1,2,2,3,4,"B","B","B"))

################
#Day 12 
P9_L2i_var_adj_12 = P9_L2i_var_adj[,c(11:18,32,33,34,35,36)]
P9_L2i_var_adj_12 = as.matrix(P9_L2i_var_adj_12)

#Plot PCA Result
P9_L2i_var_adj_12_var=P9_L2i_var_adj_12[which(apply(P9_L2i_var_adj_12, 1, var) != 0),]
D12_model = prcomp(t(P9_L2i_var_adj_12_var),scale.=T)
pca_D12 <- predict(D12_model)
pca_D12_1 <- predict(D12_model)[,1]
pca_D12_2 <- predict(D12_model)[,2]
pca_D12_3 <- predict(D12_model)[,3]
colors_D12 = c("blue","purple","purple","purple","red","red","red","red","blue","blue","blue","blue","blue")
plot(pca_D12_1,pca_D12_2,col=colors_D12,pch=(c(17,17,17,17,17,17,17,17,24,24,24,24,24)),main="Day 12 All RNAseq PCA",xlab="PC1",ylab="PC2")
text(pca_D12_1,pca_D12_2,pos=4,c(1,2,3,4,2,3,1,4,"B","B","B","B","B"))

#Adjust Day 12 for Litter
P9_12_all_n = P9_L2i_var_adj_12[which(rowMax(P9_L2i_var_adj_12)!=0),]
batch_12_all = c("11","1","9","12","1","9","11","12","30","30","30","30","30")
mod_12_all = matrix(NA, ncol=1,nrow=13)
colnames(mod_12_all)=c("genotype")
rownames(mod_12_all)=colnames(P9_12_a)
mod_12_all[c(1,9,10,11,12,13),1]="A"
mod_12_all[c(2,3,4),1]="B"
mod_12_all[c(5,6,7,8),1]="C"
P9_12_all_adj = ComBat(P9_12_all_n,batch_12_all,mod=as.factor(mod_12_all),par.prior = TRUE)
P9_12_all_adj_a = rbind(P9_12_all_adj,P9_12_all_n[which(apply(P9_12_all_n, 1, var) ==0),])

#Plot PCA Result
P9_12_all_adj_var=P9_12_all_adj_a[which(apply(P9_12_all_adj_a, 1, var) != 0),]
D12_all_model = prcomp(t(P9_12_all_adj_var),scale.=T)
pca_D12_all <- predict(D12_all_model)
pca_D12_all_1 <- predict(D12_all_model)[,1]
pca_D12_all_2 <- predict(D12_all_model)[,2]
pca_D12_all_3 <- predict(D12_all_model)[,3]
colors_D12_all = c("blue","purple","purple","purple","red","red","red","red","blue","blue","blue","blue","blue")
plot(pca_D12_all_1,pca_D12_all_2,col=colors_D12_all,pch=(c(17,17,17,17,17,17,17,17,24,24,24,24,24)),main="Day 12 All RNAseq PCA",xlab="PC1",ylab="PC2")
text(pca_D12_all_1,pca_D12_all_2,pos=4,c(1,2,3,4,2,3,1,4,"B","B","B","B","B"))

###########
#Day 16
P9_L2i_var_adj_16 = P9_L2i_var_adj[,c(19:28,37,38,39,40,41)]
P9_L2i_var_adj_16 = as.matrix(P9_L2i_var_adj_16)

#Plot PCA before
P9_L2i_var_adj_16_var=P9_L2i_var_adj_16[which(apply(P9_L2i_var_adj_16, 1, var) != 0),]
D16_model = prcomp(t(P9_L2i_var_adj_16_var),scale.=T)
pca_D16 <- predict(D16_model)
pca_D16_1 <- predict(D16_model)[,1]
pca_D16_2 <- predict(D16_model)[,2]
pca_D16_3 <- predict(D16_model)[,3]
colors_D16 = c("blue","blue","purple","purple","purple","red","red","red","red","red","blue","blue","blue","blue","blue")
plot(pca_D16_1,pca_D16_2,col=colors_D16,pch=(c(15,15,15,15,15,15,15,15,15,15,22,22,22,22,22)),main="New RNAseq PCA",xlab="PC1",ylab="PC2")
text(pca_D16_1,pca_D16_2,pos=4,c(1,2,3,4,5,3,1,2,4,5,"B","B","B","B","B"))

#Adjust Day 16 for Litter
P9_16_all_n = P9_L2i_var_adj_16[which(rowMax(P9_L2i_var_adj_16)!=0),]
batch_16_all = c("4","5","2","6","7","2","4","5","6","7","30","30","30","30","30")
mod_16_all = matrix(NA, ncol=1,nrow=15)
colnames(mod_16_all)=c("genotype")
rownames(mod_16_all)=colnames(P9_16_a)
mod_16_all[c(1,2,11,12,13,14,15),1]="A"
mod_16_all[c(3,4,5),1]="B"
mod_16_all[c(6,7,8,9,10),1]="C"
P9_16_all_adj = ComBat(P9_16_all_n,batch_16_all,mod=as.factor(mod_16_all),par.prior = TRUE)
P9_16_all_adj_a = rbind(P9_16_all_adj,P9_L2i_var_adj_16[which(rowMax(P9_L2i_var_adj_16)==0),])

#Plot PCA Result
P9_16_all_adj_a_var=P9_16_all_adj_a[which(apply(P9_16_all_adj_a, 1, var) != 0),]
D16_all_model = prcomp(t(P9_16_all_adj_a_var),scale.=T)
pca_D16_all <- predict(D16_all_model)
pca_D16_all_1 <- predict(D16_all_model)[,1]
pca_D16_all_2 <- predict(D16_all_model)[,2]
pca_D16_all_3 <- predict(D16_all_model)[,3]
colors_D16_all = c("blue","blue","purple","purple","purple","red","red","red","red","red","blue","blue","blue","blue","blue")
plot(pca_D16_all_1,pca_D16_all_2,col=colors_D16_all,pch=(c(15,15,15,15,15,15,15,15,15,15,22,22,22,22,22)),main="New RNAseq PCA",xlab="PC1",ylab="PC2")
text(pca_D16_all_1,pca_D16_all_2,pos=4,c(1,2,3,4,5,3,1,2,4,5,"B","B","B","B","B"))

#################
#All

#Group All Days
P9_all_adj_a = cbind(P9_8_all_adj_a,P9_12_all_adj_a,P9_16_all_adj_a)

#Alex Start

#Singular value decomposition (log transformed and mean centered)
#svd(matrix)
#output is x$v small, diagonal, PCs
#x$u big one gene by PCs, 
#x$w small this is what you plot
#

#Plot PCA Result
P9_all_adj_a_var=P9_all_adj_a[which(apply(P9_all_adj_a, 1, var) != 0),]
all_model = prcomp(t(P9_all_adj_a_var),scale.=T)
pca_a <- predict(all_model)
barplot(100*all_model$sdev^2 / sum(all_model$sdev^2))
(100*all_model$sdev^2 / sum(all_model$sdev^2))[1]
(100*all_model$sdev^2 / sum(all_model$sdev^2))[2]
(100*all_model$sdev^2 / sum(all_model$sdev^2))[3]
barplot(100*original_model$sdev^2 / sum(original_model$sdev^2))
(100*original_model$sdev^2 / sum(original_model$sdev^2))[1]
(100*original_model$sdev^2 / sum(original_model$sdev^2))[2]
(100*original_model$sdev^2 / sum(original_model$sdev^2))[3]
pca_a_1 <- predict(all_model)[,1]
pca_a_2 <- predict(all_model)[,2]
pca_a_3 <- predict(all_model)[,3]
pca_a_4 <- predict(all_model)[,4]
pca_a_5 <- predict(all_model)[,5]
pca_a_6 <- predict(all_model)[,6]
pca_a_7 <- predict(all_model)[,7]
pca_a_8 <- predict(all_model)[,8]
pca_a_9 <- predict(all_model)[,9]
pca_a_10 <- predict(all_model)[,10]
pca_a_11 <- predict(all_model)[,11]
pca_a_12 <- predict(all_model)[,12]
pca_a_13 <- predict(all_model)[,13]
pca_a_14 <- predict(all_model)[,14]
colors_a = c("blue","blue","blue","purple","purple","red","red","red","red","red","blue","blue","blue","blue","purple","purple","purple","red","red","red","red","blue","blue","blue","blue","blue","blue","blue","purple","purple","purple","red","red","red","red","red","blue","blue","blue","blue","blue")
symbols_a = as.numeric(c(16,16,16,16,16,16,16,16,16,16,1,1,1,17,17,17,17,17,17,17,17,2,2,2,2,2,15,15,15,15,15,15,15,15,15,15,0,0,0,0,0))
plot(pca_a_1,pca_a_2,col=colors_a,pch=symbols_a,main="New RNAseq PCA",xlab="PC1",ylab="PC2",axes=FALSE)
axis(1,seq(-200,300,50))
axis(2,seq(-200,300,50))


plot(hclust(dist(pca_a_1)))
plot(hclust(dist(pca_a_2)))

pca_a_1
#try with ggplot
PCA_Main = as.data.frame(cbind(pca_a_1,pca_a_2))
ggplot(PCA_Main,aes(x = pca_a_1, y = pca_a_2,color=colors_a))+geom_point(color=colors_a)

#Other plots
plot(pca_a_1,pca_a_2,col=colors_a,pch=symbols_a,main="New RNAseq PCA",xlab="PC1",ylab="PC2")
plot(pca_a_1,pca_a_3,col=colors_a,pch=symbols_a,main="New RNAseq PCA",xlab="PC1",ylab="PC3")
plot(pca_a_1,pca_a_4,col=colors_a,pch=symbols_a,main="New RNAseq PCA",xlab="PC1",ylab="PC4")
plot(pca_a_1,pca_a_5,col=colors_a,pch=symbols_a,main="New RNAseq PCA",xlab="PC1",ylab="PC5")
plot(pca_a_1,pca_a_6,col=colors_a,pch=symbols_a,main="New RNAseq PCA",xlab="PC1",ylab="PC6")
plot(pca_a_1,pca_a_7,col=colors_a,pch=symbols_a,main="New RNAseq PCA",xlab="PC1",ylab="PC7")
plot(pca_a_1,pca_a_8,col=colors_a,pch=symbols_a,main="New RNAseq PCA",xlab="PC1",ylab="PC8")
plot(pca_a_1,pca_a_9,col=colors_a,pch=symbols_a,main="New RNAseq PCA",xlab="PC1",ylab="PC9")
plot(pca_a_1,pca_a_10,col=colors_a,pch=symbols_a,main="New RNAseq PCA",xlab="PC1",ylab="PC10")
plot(pca_a_1,pca_a_11,col=colors_a,pch=symbols_a,main="New RNAseq PCA",xlab="PC1",ylab="PC11")
plot(pca_a_1,pca_a_12,col=colors_a,pch=symbols_a,main="New RNAseq PCA",xlab="PC1",ylab="PC12")
plot(pca_a_1,pca_a_13,col=colors_a,pch=symbols_a,main="New RNAseq PCA",xlab="PC1",ylab="PC13")
plot(pca_a_1,pca_a_14,col=colors_a,pch=symbols_a,main="New RNAseq PCA",xlab="PC1",ylab="PC14")
scatterplot3d(pca_a_1,pca_a_2,pca_a_3,color = colors_a,pch = symbols_a)
scatter3d(pca_a_1,pca_a_2,pca_a_3)
scatter3d(pca_a_2,pca_a_3,pca_a_5)

#Load Data
setwd("Data")
load("P9_all_adj_a.rdt")
setwd("..")


#Try again with SVD
P9_all_adj_a_var=P9_all_adj_a[which(apply(P9_all_adj_a, 1, var) != 0),]
P9_all_adj_a_var_c = t(scale(t(P9_all_adj_a_var),center=TRUE,scale=TRUE))
P9_SVD = svd(P9_all_adj_a_var_c)
svd_sym = c(16,16,16,16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,17,17,17,17,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15)
colors_a = c("blue","blue","blue","purple","purple","red","red","red","red","red","blue","blue","blue","blue","purple","purple","purple","red","red","red","red","blue","blue","blue","blue","blue","blue","blue","purple","purple","purple","red","red","red","red","red","blue","blue","blue","blue","blue")
plot(-P9_SVD$v[,1],P9_SVD$v[,2],col=colors_a,pch=svd_sym)

#Set up dataframe for ggplot
P9_PCA = as.data.frame(P9_SVD$v)
colnames(P9_PCA) = 1:ncol(P9_PCA)
colnames(P9_PCA) = paste("PC_",colnames(P9_PCA),sep="")

#Run this
P9_all_adj_a_var=P9_all_adj_a[which(apply(P9_all_adj_a, 1, var) != 0),]
all_model = prcomp(t(P9_all_adj_a_var),scale.=T)
pca_a <- predict(all_model)
barplot(100*all_model$sdev^2 / sum(all_model$sdev^2))
pca_a_1 <- predict(all_model)[,1]
pca_a_2 <- predict(all_model)[,2]

#Add genotype data
PCA_Main = as.data.frame(cbind(pca_a_1,pca_a_2))
P9_PCA[,"geno_data"] = substr(rownames(PCA_Main),1,1)
P9_PCA[,"genotype"] = "NA"
for (i in 1:nrow(P9_PCA)){
  if (P9_PCA[i,"geno_data"]=="R"){
    P9_PCA[i,"genotype"] = "W"
  }
  else {
    P9_PCA[i,"genotype"] = P9_PCA[i,"geno_data"]
  }
}
P9_PCA[,"age"] = as.factor(substr(rownames(PCA_Main),3,4))

#Set up dataframe pre-ComBat for ggplot
P9_L2i_var_c = t(scale(t(P9_L2i_var),center=TRUE,scale=TRUE))
Pre_SVD = svd(P9_L2i_var_c)
Pre_PCA = as.data.frame(Pre_SVD$v)
colnames(Pre_PCA) = 1:ncol(Pre_PCA)
colnames(Pre_PCA) = paste("PC_",colnames(Pre_PCA),sep="")
plot(-Pre_SVD$v[,1],Pre_SVD$v[,2])

#Add genotype data
Pre_PCA[,"geno_data"] = substr(colnames(P9_L2i_var),1,1)
Pre_PCA[,"genotype"] = "NA"
for (i in 1:nrow(Pre_PCA)){
  if (Pre_PCA[i,"geno_data"]=="R"){
    Pre_PCA[i,"genotype"] = "W"
  }
  else {
    Pre_PCA[i,"genotype"] = Pre_PCA[i,"geno_data"]
  }
}
Pre_PCA[,"age"] = as.factor(substr(colnames(P9_L2i_var),3,4))



colors_bw = c("#252525","#252525","#252525","#636363","#636363","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#252525","#252525","#252525","#252525","#636363","#636363","#636363","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#252525","#252525","#252525","#252525","#252525","#252525","#252525","#636363","#636363","#636363","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#252525","#252525","#252525","#252525","#252525")
sym_bw = c(21,21,21,21,21,21,21,21,21,21,21,21,21,24,24,24,24,24,24,24,24,24,24,24,24,24,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22)

colors_c = c("blue","blue","blue","orchid","orchid","red","red","red","red","red","blue","blue","blue","blue","orchid","orchid","orchid","red","red","red","red","blue","blue","blue","blue","blue","blue","blue","orchid","orchid","orchid","red","red","red","red","red","blue","blue","blue","blue","blue")
colors_d = c("blue","blue","blue","orchid","orchid","red","red","red","red","red","blue","orchid","orchid","orchid","red","red","red","red","blue","blue","orchid","orchid","orchid","red","red","red","red","red","white","white","white","white","white","white","white","white","white","white","white","white","white")

colors_e = c("#a50f15","#a50f15","#a50f15","#fb6a4a","#fb6a4a","#fee5d9","#fee5d9","#fee5d9","#fee5d9","#fee5d9","#a50f15","#a50f15","#a50f15","#a50f15","#fb6a4a","#fb6a4a","#fb6a4a","#fee5d9","#fee5d9","#fee5d9","#fee5d9","#a50f15","#a50f15","#a50f15","#a50f15","#a50f15","#a50f15","#a50f15","#fb6a4a","#fb6a4a","#fb6a4a","#fee5d9","#fee5d9","#fee5d9","#fee5d9","#fee5d9","#a50f15","#a50f15","#a50f15","#a50f15","#a50f15")


ggplot(P9_PCA,aes(x = -PC_1, y = PC_2,colour=colors_bw))+
  geom_point(color = "black",fill=colors_bw,pch=sym_bw,size=4,stroke=1) +
  theme(panel.border = element_rect(colour = "grey",fill=NA,size=1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

svd_sym2 = c(1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
svd_sym3 = c(21,21,21,21,21,21,21,21,21,21,21,21,21,24,24,24,24,24,24,24,24,24,24,24,24,24,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22)
svd_sym4 = c(21,21,21,21,21,21,21,21,21,21,24,24,24,24,24,24,24,24,22,22,22,22,22,22,22,22,22,22,21,21,21,24,24,24,24,24,22,22,22,22,22)





ggplot(P9_PCA,aes(x = -PC_1, y = PC_2,colour=colors_a))+
  geom_point(color="black",fill=colors_a,pch=svd_sym3,size=4,stroke=1) +
  theme(panel.border = element_rect(colour = "grey",fill=NA,size=1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())


ggplot(P9_PCA,aes(x = -PC_1, y = PC_2,colour=colors_b))+
  geom_point(color="black",fill=colors_b,pch=svd_sym3,size=4,stroke=1) +
  theme(panel.border = element_rect(colour = "grey",fill=NA,size=1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())


ggplot(P9_PCA,aes(x = -PC_1, y = PC_2,colour=colors_c))+
  geom_point(color="black",fill=colors_c,pch=svd_sym3,size=4,stroke=1) +
  theme(panel.border = element_rect(colour = "grey",fill=NA,size=1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) #+ 
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggplot(P9_PCA,aes(x = -PC_1, y = PC_3,colour=colors_c))+
  geom_point(color="black",fill=colors_c,pch=svd_sym3,size=4,stroke=1) +
  theme(panel.border = element_rect(colour = "grey",fill=NA,size=1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggplot(P9_PCA,aes(x = -PC_1, y = PC_2,colour=colors_e))+
  geom_point(color="black",fill=colors_e,pch=svd_sym3,size=4,stroke=1) +
  theme(panel.border = element_rect(colour = "grey",fill=NA,size=1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggplot(Pre_PCA,aes(x = -PC_1, y = PC_2,colour=colors_d))+
  geom_point(color="black",fill=colors_d,pch=svd_sym4,size=4,stroke=1) +
  theme(panel.border = element_rect(colour = "grey",fill=NA,size=1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) #+ 
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggplot(P9_PCA,aes(x = -PC_1, y = PC_3,colour=colors_a))+geom_point(color=colors_a,pch=svd_sym,size=3)
ggplot(P9_PCA,aes(x = -PC_1, y = PC_4,colour=colors_a))+geom_point(color=colors_a,pch=svd_sym,size=3)
ggplot(P9_PCA,aes(x = -PC_1, y = PC_5,colour=colors_a))+geom_point(color=colors_a,pch=svd_sym,size=3)


barplot(100*testc$d^2 / sum(testc$d^2))
hist(testc$u[,2],breaks=100)
abline(v=2*sd(testc$u[,2]))
abline(v=-2*sd(testc$u[,2]))


#Day by day

Pre_d08 = P9_L2i_var[,grep("_08_",colnames(P9_L2i_var))]
Pre_d12 = P9_L2i_var[,grep("_12_",colnames(P9_L2i_var))]
Pre_d16 = P9_L2i_var[,grep("_16_",colnames(P9_L2i_var))]

Pre_d08_var = Pre_d08[which(apply(Pre_d08, 1, var) != 0),]
Pre_d12_var = Pre_d12[which(apply(Pre_d12, 1, var) != 0),]
Pre_d16_var = Pre_d16[which(apply(Pre_d16, 1, var) != 0),]

Pre_d08_c = t(scale(t(Pre_d08_var),center=TRUE,scale=TRUE))
Pre_d12_c = t(scale(t(Pre_d12_var),center=TRUE,scale=TRUE))
Pre_d16_c = t(scale(t(Pre_d16_var),center=TRUE,scale=TRUE))

Pre_SVD_08 = svd(Pre_d08_c)
Pre_PCA_08 = as.data.frame(Pre_SVD_08$v)
colnames(Pre_PCA_08) = 1:ncol(Pre_PCA_08)
colnames(Pre_PCA_08) = paste("PC_",colnames(Pre_PCA_08),sep="")
plot(Pre_SVD_08$v[,1],Pre_SVD_08$v[,2])

(100*Pre_SVD_08$d^2 / sum(Pre_SVD_08$d^2))[2]
(100*P9_SVD$d^2 / sum(P9_SVD$d^2))[3]



Pre_SVD_12 = svd(Pre_d12_c)
Pre_PCA_12 = as.data.frame(Pre_SVD_12$v)
colnames(Pre_PCA_12) = 1:ncol(Pre_PCA_12)
colnames(Pre_PCA_12) = paste("PC_",colnames(Pre_PCA_12),sep="")
plot(Pre_SVD_12$v[,1],Pre_SVD_12$v[,2])

(100*Pre_SVD_12$d^2 / sum(Pre_SVD_12$d^2))[2]


Pre_SVD_16 = svd(Pre_d16_c)
Pre_PCA_16 = as.data.frame(Pre_SVD_16$v)
colnames(Pre_PCA_16) = 1:ncol(Pre_PCA_16)
colnames(Pre_PCA_16) = paste("PC_",colnames(Pre_PCA_16),sep="")
plot(Pre_SVD_16$v[,1],Pre_SVD_16$v[,2])

(100*Pre_SVD_16$d^2 / sum(Pre_SVD_16$d^2))[2]

#Add genotype data
Pre_PCA_08[,"geno_data"] = substr(colnames(Pre_d08_c),1,1)
Pre_PCA_08[,"genotype"] = "NA"
for (i in 1:nrow(Pre_PCA_08)){
  if (Pre_PCA_08[i,"geno_data"]=="R"){
    Pre_PCA_08[i,"genotype"] = "W"
  }
  else {
    Pre_PCA_08[i,"genotype"] = Pre_PCA_08[i,"geno_data"]
  }
}
Pre_PCA_08[,"age"] = as.factor(substr(colnames(Pre_d08_c),3,4))
Pre_PCA_08[,"litter"] = c("1","2","2","3","4","1","2","2","3","4","B","B","B")

Pre_PCA_12[,"geno_data"] = substr(colnames(Pre_d12_c),1,1)
Pre_PCA_12[,"genotype"] = "NA"
for (i in 1:nrow(Pre_PCA_12)){
  if (Pre_PCA_12[i,"geno_data"]=="R"){
    Pre_PCA_12[i,"genotype"] = "W"
  }
  else {
    Pre_PCA_12[i,"genotype"] = Pre_PCA_12[i,"geno_data"]
  }
}
Pre_PCA_12[,"age"] = as.factor(substr(colnames(Pre_d12_c),3,4))
Pre_PCA_12[,"litter"] = c("1","2","3","4","2","3","1","4","B","B","B","B","B")


Pre_PCA_16[,"geno_data"] = substr(colnames(Pre_d16_c),1,1)
Pre_PCA_16[,"genotype"] = "NA"
for (i in 1:nrow(Pre_PCA_16)){
  if (Pre_PCA_16[i,"geno_data"]=="R"){
    Pre_PCA_16[i,"genotype"] = "W"
  }
  else {
    Pre_PCA_16[i,"genotype"] = Pre_PCA_16[i,"geno_data"]
  }
}
Pre_PCA_16[,"age"] = as.factor(substr(colnames(Pre_d16_c),3,4))
Pre_PCA_16[,"litter"] = c("1","2","3","4","5","3","1","2","4","5","B","B","B","B","B")

colors_Pre8 = c("Blue","Blue","Blue","orchid","orchid","red","red","red","red","red","white","white","white")
sym_Pre8 = 21
ggplot(Pre_PCA_08,aes(x = PC_1, y = PC_2,colour=colors_Pre8,label=litter))+
  geom_point(color="black",fill=colors_Pre8,pch=sym_Pre8,size=4,stroke=1) +
  geom_text(hjust = 0, nudge_x = 0.025,color="black")+
  theme(panel.border = element_rect(colour = "grey",fill=NA,size=1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) #+ 
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

colors_Pre12 = c("Blue","orchid","orchid","orchid","red","red","red","red","white","white","white","white","white")
sym_Pre12 = 24

ggplot(Pre_PCA_12,aes(x = -PC_1, y = PC_2,colour=colors_Pre12,label=litter))+
  geom_point(color="black",fill=colors_Pre12,pch=sym_Pre12,size=4,stroke=1) +
  geom_text(hjust = 0, nudge_x = 0.025,color="black")+
  theme(panel.border = element_rect(colour = "grey",fill=NA,size=1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) #+ 
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

colors_Pre16 = c("Blue","blue","orchid","orchid","orchid","red","red","red","red","red","white","white","white","white","white")
sym_Pre16 = 22

ggplot(Pre_PCA_16,aes(x = -PC_1, y = PC_2,colour=colors_Pre16,label=litter))+
  geom_point(color="black",fill=colors_Pre16,pch=sym_Pre16,size=4,stroke=1) +
  geom_text(hjust = 0, nudge_x = 0.025,color="black")+
  theme(panel.border = element_rect(colour = "grey",fill=NA,size=1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) #+ 
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())




Post_d08 = P9_all_adj_a[,grep("_08_",colnames(P9_all_adj_a))]
Post_d12 = P9_all_adj_a[,grep("_12_",colnames(P9_all_adj_a))]
Post_d16 = P9_all_adj_a[,grep("_16_",colnames(P9_all_adj_a))]

Post_d08_var = Post_d08[which(apply(Post_d08, 1, var) != 0),]
Post_d12_var = Post_d12[which(apply(Post_d12, 1, var) != 0),]
Post_d16_var = Post_d16[which(apply(Post_d16, 1, var) != 0),]

Post_d08_c = t(scale(t(Post_d08_var),center=TRUE,scale=TRUE))
Post_d12_c = t(scale(t(Post_d12_var),center=TRUE,scale=TRUE))
Post_d16_c = t(scale(t(Post_d16_var),center=TRUE,scale=TRUE))

Post_SVD_08 = svd(Post_d08_c)
Post_PCA_08 = as.data.frame(Post_SVD_08$v)
colnames(Post_PCA_08) = 1:ncol(Post_PCA_08)
colnames(Post_PCA_08) = paste("PC_",colnames(Post_PCA_08),sep="")
plot(Post_SVD_08$v[,1],Post_SVD_08$v[,2])

(100*Post_SVD_08$d^2 / sum(Post_SVD_08$d^2))[2]


Post_SVD_12 = svd(Post_d12_c)
Post_PCA_12 = as.data.frame(Post_SVD_12$v)
colnames(Post_PCA_12) = 1:ncol(Post_PCA_12)
colnames(Post_PCA_12) = paste("PC_",colnames(Post_PCA_12),sep="")
plot(Post_SVD_12$v[,1],Post_SVD_12$v[,2])

(100*Post_SVD_12$d^2 / sum(Post_SVD_12$d^2))[2]

Post_SVD_16 = svd(Post_d16_c)
Post_PCA_16 = as.data.frame(Post_SVD_16$v)
colnames(Post_PCA_16) = 1:ncol(Post_PCA_16)
colnames(Post_PCA_16) = paste("PC_",colnames(Post_PCA_16),sep="")
plot(Post_SVD_16$v[,1],Post_SVD_16$v[,2])

(100*Post_SVD_16$d^2 / sum(Post_SVD_16$d^2))[2]

#Add genotype data
Post_PCA_08[,"geno_data"] = substr(colnames(Post_d08_c),1,1)
Post_PCA_08[,"genotype"] = "NA"
for (i in 1:nrow(Post_PCA_08)){
  if (Post_PCA_08[i,"geno_data"]=="R"){
    Post_PCA_08[i,"genotype"] = "W"
  }
  else {
    Post_PCA_08[i,"genotype"] = Post_PCA_08[i,"geno_data"]
  }
}
Post_PCA_08[,"age"] = as.factor(substr(colnames(Post_d08_c),3,4))
Post_PCA_08[,"litter"] = c("1","2","2","3","4","1","2","2","3","4","B","B","B")

Post_PCA_12[,"geno_data"] = substr(colnames(Post_d12_c),1,1)
Post_PCA_12[,"genotype"] = "NA"
for (i in 1:nrow(Post_PCA_12)){
  if (Post_PCA_12[i,"geno_data"]=="R"){
    Post_PCA_12[i,"genotype"] = "W"
  }
  else {
    Post_PCA_12[i,"genotype"] = Post_PCA_12[i,"geno_data"]
  }
}
Post_PCA_12[,"age"] = as.factor(substr(colnames(Post_d12_c),3,4))
Post_PCA_12[,"litter"] = c("1","2","3","4","2","3","1","4","B","B","B","B","B")


Post_PCA_16[,"geno_data"] = substr(colnames(Post_d16_c),1,1)
Post_PCA_16[,"genotype"] = "NA"
for (i in 1:nrow(Post_PCA_16)){
  if (Post_PCA_16[i,"geno_data"]=="R"){
    Post_PCA_16[i,"genotype"] = "W"
  }
  else {
    Post_PCA_16[i,"genotype"] = Post_PCA_16[i,"geno_data"]
  }
}
Post_PCA_16[,"age"] = as.factor(substr(colnames(Post_d16_c),3,4))
Post_PCA_16[,"litter"] = c("1","2","3","4","5","3","1","2","4","5","B","B","B","B","B")

colors_Post8 = c("Blue","Blue","Blue","orchid","orchid","red","red","red","red","red","white","white","white")
sym_Post8 = 21
ggplot(Post_PCA_08,aes(x = PC_1, y = PC_2,colour=colors_Post8,label=litter))+
  geom_point(color="black",fill=colors_Post8,pch=sym_Post8,size=4,stroke=1) +
  geom_text(hjust = 0, nudge_x = 0.025,color="black")+
  theme(panel.border = element_rect(colour = "grey",fill=NA,size=1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

colors_Post12 = c("Blue","orchid","orchid","orchid","red","red","red","red","white","white","white","white","white")
sym_Post12 = 24

ggplot(Post_PCA_12,aes(x = -PC_1, y = PC_2,colour=colors_Post12,label=litter))+
  geom_point(color="black",fill=colors_Post12,pch=sym_Post12,size=4,stroke=1) +
  geom_text(hjust = 0, nudge_x = 0.025,color="black")+
  theme(panel.border = element_rect(colour = "grey",fill=NA,size=1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

colors_Post16 = c("Blue","blue","orchid","orchid","orchid","red","red","red","red","red","white","white","white","white","white")
sym_Post16 = 22

ggplot(Post_PCA_16,aes(x = PC_1, y = PC_2,colour=colors_Post16,label=litter))+
  geom_point(color="black",fill=colors_Post16,pch=sym_Post16,size=4,stroke=1) +
  geom_text(hjust = 0, nudge_x = 0.025,color="black")+
  theme(panel.border = element_rect(colour = "grey",fill=NA,size=1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())






#Sava data
setwd("Data")
save(P9_all_adj_a,file="P9_all_adj_a.rdt")
setwd("..")


##################################################################################################################
############ Rename by gene ######################################################################################
###### P9_all_adj_a  -> P9_DL ####################################################################################

#Load Data
setwd("Data")
load("P9_all_adj_a.rdt")
setwd("..")

#Rename them to be by gene name
m_ensembl=useMart("ensembl",dataset="mmusculus_gene_ensembl")
gene_info=getBM(attributes=c("external_gene_name","ensembl_gene_id"),mart=m_ensembl)

#Create new datatable
P9_DL = P9_all_adj_a_var
P9_DL = as.data.frame(P9_DL)
P9_DL$gene_name = "NA"

#Fill gene name column
for (i in 1:nrow(P9_DL)){
  EID=rownames(P9_DL[i,])
  if (! (EID %in% as.character(gene_info$ensembl_gene_id))) {
    P9_DL[i,"gene_name"] = EID
  }
  
  else{
    P9_DL[i,"gene_name"]=as.character(gene_info[which(gene_info$ensembl_gene_id==EID),]$external_gene_name)[1]
  }
}

#Find duplicates
P9DL = P9_DL[-which(duplicated(P9_DL$gene_name)),]
P9Dup = P9_DL[which(duplicated(P9_DL$gene_name)),]

#Make Gene name the row name
rownames(P9DL) = P9DL[,"gene_name"]
P9DL["A530058N18Rik",] = P9DL["A530058N18Rik",1:41] + P9Dup[1,1:41]
P9DL["4930556M19Rik",] = P9DL["4930556M19Rik",1:41] + P9Dup[2,1:41]

#Remove Gene name column
P9DL = P9DL[,1:41]

##################################################################################################################
############ Limit to those with at least 1 adj log(TPM+1) #######################################################
###### P9_DL  -> P9DL_1 ##########################################################################################

#Remove those that dont reach 1 adjusted Log(TPM+1)
P9DL_1 = P9DL[which(rowMax(as.matrix(P9DL))>1),]

#Save data
setwd("Results")
save(P9DL_1,file="P9DL_gene.rdt")
setwd("..")