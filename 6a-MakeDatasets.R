#   The purpose of this file is to make datasets appropriate for PMCA

#   It takes csv files that contain information about cell type proportions, as counted by Yasu via chromasome spreads,
#   as well as pre-processed RNAseq transcript abundances, and outputs re-organized datasets that could be run in PMCA.
#   The data is saved in Desktop/Carter/Projects/P9/6 PMCA/Data/ 


##################################################################################################################
############ Load Packages #######################################################################################
##################################################################################################################

#load packages
require(Biobase)


##################################################################################################################
############ Make gene expression table for unadjusted data ######################################################
###### Y #########################################################################################################

#Set working directory
setwd("/Users/s-fine/Desktop/Carter/Projects/P9/2 Normalization/")

#Load datasets
setwd("Data")
load("AllData.rdt")
setwd("..")

#Grap mutant data
MutantData = AllData[,grep("M_",colnames(AllData))]
HetData = AllData[,grep("H_",colnames(AllData))]

#Log transform data
LogData = log2(MutantData+1)
LogHData = log2(HetData+1)

#Select dat of interest
NewLog = LogData[,c(1:9,11:15)]
NewHLog = LogHData[,c(1:5,7:9)]

#Make dataset a matrix
NewLog = as.matrix(NewLog)
NewHLog = as.matrix(NewHLog)

#Remove genes that are never expressed
NewLog_0 = NewLog[which(rowMax(NewLog)!=0),]
NewHLog_0 = NewHLog[which(rowMax(NewHLog)!=0),]

#Remove genes that never go above Log2(TPM+1) = 2
NewLog_2 = NewLog[which(rowMax(NewLog)>=2),]
NewHLog_2 = NewHLog[which(rowMax(NewHLog)>=2),]

#Save data as Y
Y = NewLog_2
YH = NewHLog_2

#Save Y
setwd("/Users/s-fine/Desktop/Carter/Projects/P9/6 PMCA/Data")
save(Y,file="GeneExpression_NonAdjusted_Y.rdt")
save(YH,file="GeneExpression_NonAdjusted_YH.rdt")
setwd("..")


##################################################################################################################
############ Make gene expression table for adjusted data ########################################################
###### Y1 ########################################################################################################

#Load adjusted data
setwd("/Users/s-fine/Desktop/Carter/Projects/P9/3 DifferentialExpression/Results/")
load(file="P9K_Analysis.rdt")
setwd("/Users/s-fine/Desktop/Carter/Projects/P9/6 PMCA/")

#Just take the data
P9K_PMCA = P9K_Analysis[,1:41]

#Rename rows to be Gene IDs
rownames(P9K_PMCA) = P9K_Analysis$I_GEID

#Identify the lowest gene's expression
min(rowMin(as.matrix(P9K_PMCA))) #result = -2.375718

#Add lowest value to all genes
P9K_PMCA = P9K_PMCA + 2.375718

#Select mutant samples
P9K_PMCA_M = P9K_PMCA[,grep("D_P9K_M",colnames(P9K_PMCA))]

#Remove genes that never go above Log2(TPM+1) = 2
P9K_PMCA_2 = P9K_PMCA_M[which(rowMax(as.matrix(P9K_PMCA_M))>=2),]

#Make dataset a matrix
P9K_PMCA_2 = as.matrix(P9K_PMCA_2)

#Rename column names to match Y
colnames(P9K_PMCA_2) = colnames(NewLog_2)

#Save data as Y1
Y1 = P9K_PMCA_2

#Save Y1
setwd("Data")
save(Y1,file="GeneExpression_Adjusted_Y1.rdt")
setwd("..")


##################################################################################################################
############ Make sample info table ##############################################################################
###### X #########################################################################################################

#Load Sample info
setwd("/Users/s-fine/Desktop/Carter/Projects/P9/1 Cytology/Data")
Prdm9_Sample_Info=read.csv("Prdm9_Sample_info.csv")
setwd("/Users/s-fine/Desktop/Carter/Projects/P9/6 PMCA/")

#Select important rows and columns
P9_Sample_Data = Prdm9_Sample_Info[1:30,6:17]

#Select just the mutant samples
P9_Sample_Data = P9_Sample_Data[grep("MUT",P9_Sample_Data$Genotype),]
P9_Sample_Data_H = P9_Sample_Data[grep("HET",P9_Sample_Data$Genotype),]

#Remove the sample that looked like a het by expression
P9_Sample_Data_N = P9_Sample_Data[-which(rownames(P9_Sample_Data)=="19"),]
P9_Sample_Data_HN = P9_Sample_Data_H[-which(rownames(P9_Sample_Data_H)=="20"),]

#Make the rownames match Y column names
rownames(P9_Sample_Data_N) = colnames(NewLog)
rownames(P9_Sample_Data_HN) = colnames(NewHLog)

#Remove extra info
P9_Sample_Data_Nx = P9_Sample_Data_N[,4:12]
P9_Sample_Data_NHx = P9_Sample_Data_HN[,4:12]

#Transform dataset
P9_Sample_Data_t = t(P9_Sample_Data_Nx)
P9_Sample_Data_tH = t(P9_Sample_Data_NHx)

#Remove substages that are 0
P9_Sample_Data_m = P9_Sample_Data_t[1:7,]
P9_Sample_Data_mH = P9_Sample_Data_tH[1:7,]

#Makde a dataset that summarizes proportions of sample data
P9_Cells = matrix(NA,nrow=nrow(P9_Sample_Data_m),ncol=ncol(P9_Sample_Data_m))
P9_Cells[,1] = P9_Sample_Data_m[,1] / sum(P9_Sample_Data_m[,1])
P9_Cells[,2] = P9_Sample_Data_m[,2] / sum(P9_Sample_Data_m[,2])
P9_Cells[,3] = P9_Sample_Data_m[,3] / sum(P9_Sample_Data_m[,3])
P9_Cells[,4] = P9_Sample_Data_m[,4] / sum(P9_Sample_Data_m[,4])
P9_Cells[,5] = P9_Sample_Data_m[,5] / sum(P9_Sample_Data_m[,5])
P9_Cells[,6] = P9_Sample_Data_m[,6] / sum(P9_Sample_Data_m[,6])
P9_Cells[,7] = P9_Sample_Data_m[,7] / sum(P9_Sample_Data_m[,7])
P9_Cells[,8] = P9_Sample_Data_m[,8] / sum(P9_Sample_Data_m[,8])
P9_Cells[,9] = P9_Sample_Data_m[,9] / sum(P9_Sample_Data_m[,9])
P9_Cells[,10] = P9_Sample_Data_m[,10] / sum(P9_Sample_Data_m[,10])
P9_Cells[,11] = P9_Sample_Data_m[,11] / sum(P9_Sample_Data_m[,11])
P9_Cells[,12] = P9_Sample_Data_m[,12] / sum(P9_Sample_Data_m[,12])
P9_Cells[,13] = P9_Sample_Data_m[,13] / sum(P9_Sample_Data_m[,13])
P9_Cells[,14] = P9_Sample_Data_m[,14] / sum(P9_Sample_Data_m[,14])
rownames(P9_Cells) = rownames(P9_Sample_Data_m)
colnames(P9_Cells) = colnames(P9_Sample_Data_m)


P9_HCells = matrix(NA,nrow=nrow(P9_Sample_Data_mH),ncol=ncol(P9_Sample_Data_mH))
P9_HCells[,1] = P9_Sample_Data_mH[,1] / sum(P9_Sample_Data_mH[,1])
P9_HCells[,2] = P9_Sample_Data_mH[,2] / sum(P9_Sample_Data_mH[,2])
P9_HCells[,3] = P9_Sample_Data_mH[,3] / sum(P9_Sample_Data_mH[,3])
P9_HCells[,4] = P9_Sample_Data_mH[,4] / sum(P9_Sample_Data_mH[,4])
P9_HCells[,5] = P9_Sample_Data_mH[,5] / sum(P9_Sample_Data_mH[,5])
P9_HCells[,6] = P9_Sample_Data_mH[,6] / sum(P9_Sample_Data_mH[,6])
P9_HCells[,7] = P9_Sample_Data_mH[,7] / sum(P9_Sample_Data_mH[,7])
P9_HCells[,8] = P9_Sample_Data_mH[,8] / sum(P9_Sample_Data_mH[,8])

rownames(P9_HCells) = rownames(P9_Sample_Data_mH)
colnames(P9_HCells) = colnames(P9_Sample_Data_mH)

#Save as X
X = P9_Cells
XH = P9_HCells

#Combine Late leptotene and Zygotene
X["L.leptotene",] = X["L.leptotene",] + X["Zygotene",]
XH["L.leptotene",] = XH["L.leptotene",] + XH["Zygotene",]

X = X[-5,]
XH = XH[-5,]

#Change rownames to match
rownames(X)=c("Sp","PreL","EL","LLZ","P-L","EP")
rownames(XH)=c("Sp","PreL","EL","LLZ","P-L","EP")

#Save X
setwd("Data")
save(X,file="SampleInfo_X.rdt")
save(XH,file="SampleInfo_XH.rdt")
setwd("..")


##################################################################################################################
############ Plot PCA for all sample info ########################################################################
###### pca_Cell ##################################################################################################

#Limit sample info to useful rows
P9_Cell_Data = Prdm9_Sample_Info[1:30,6:16]

#Remove swapped samples
P9_Cell_Data_N = P9_Cell_Data[-which(rownames(P9_Cell_Data)=="19" | rownames(P9_Cell_Data)=="20"),]

#Set informative rownames
rownames(P9_Cell_Data_N) = paste(P9_Cell_Data_N[,"Age..dpp."],P9_Cell_Data_N[,"Genotype"],P9_Cell_Data_N[,"Litter"],rownames(P9_Cell_Data_N),sep="_")

#Remove unnecessary data
P9_Cell_Data_U = P9_Cell_Data_N[,4:11]

#Transform dataset
P9_Cell_Data_T = t(P9_Cell_Data_U)

#Make proportion table from data
P9_Pro = matrix(NA,nrow=nrow(P9_Cell_Data_T),ncol=ncol(P9_Cell_Data_T))
for (i in 1:ncol(P9_Cell_Data_T)){
  P9_Pro[,i] = P9_Cell_Data_T[,i] / sum(P9_Cell_Data_T[,i])
}

#Change rownames and column names
rownames(P9_Pro) = rownames(P9_Cell_Data_T)
colnames(P9_Pro) = colnames(P9_Cell_Data_T)

#Run PCA
P9_Cell_model = prcomp(t(P9_Pro),scale.=T)
pca_Cell <- predict(P9_Cell_model)
pca_Cell_1 <- predict(P9_Cell_model)[,1]
pca_Cell_2 <- predict(P9_Cell_model)[,2]
pca_Cell_3 <- predict(P9_Cell_model)[,3]
pca_Cell_4 <- predict(P9_Cell_model)[,4]
pca_Cell_5 <- predict(P9_Cell_model)[,5]
colors_Cell2 = as.factor(P9_Cell_Data_N$Litter)
colors_Cell = c("red","blue","red","red","blue","blue","purple","red","purple","red","purple","red","purple","red","blue","red","red","purple","purple","red","red","blue","blue","red","red","purple","red","purple")
symbols_Cell = as.numeric(c(16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,15,15,15,15,15,15,15,15,15,15))

#Plot PCA
plot(-pca_Cell_1,pca_Cell_2,col=colors_Cell,pch = symbols_Cell)
plot(-pca_Cell_1,pca_Cell_3,col=colors_Cell2,pch = as.numeric(as.factor(colors_Cell)))
plot(-pca_Cell_1,-pca_Cell_3,col=colors_Cell,pch = symbols_Cell)
plot(pca_Cell_1,pca_Cell_4,col=colors_Cell,pch = symbols_Cell)
plot(pca_Cell_1,pca_Cell_5,col=colors_Cell,pch = symbols_Cell)
