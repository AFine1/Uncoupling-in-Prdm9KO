#   The purpose of this file is to take cytological data and create a visual model of variance between sample groups

#   It takes csv files that contain information about cell type proportions, as counted by Yasu via chromasome spreads
#   and outputs one barplot visual for average proportions of cell types across sample groups. The output file is saved
#   in Desktop/Carter/Projects/P9/Cytology/Results/ and Desktop/Carter/Projects/P9/Figures as Cytological_Proportions.png


##################################################################################################################
############ Load Packages #######################################################################################
##################################################################################################################

#load packages
library(ggplot2)
library(reshape2)
library(scales)
library(stringr)


##################################################################################################################
############ Aquire Data #########################################################################################
###### Prdm9_Sample_Info & WT_Sample_Info ########################################################################

#Set working directory
setwd("/Users/s-fine/Desktop/Carter/Projects/P9/1 Cytology/")

#Load datasets
setwd("Data")
Prdm9_Sample_Info=read.csv("Prdm9_Sample_info.csv")
WT_Sample_Info=read.csv("Handel13Cytology2014Feb27.csv")
setwd("..")

#Select columns of interest
Prdm9_Sample_Info = Prdm9_Sample_Info[c(1:18,21:30),c(7:17)]
WT_Sample_Info = WT_Sample_Info[c(1:3,11:15,21:25),c(1,6:14)]


##################################################################################################################
############ Convert to percentages ##############################################################################
###### Prdm9_Sample_Info & WT_Sample_Info -> Prdm9_Sample_Percent & WT_Sample_Percent ############################

#Convert P9 dataset to percent
Prdm9_Sample_Percent = Prdm9_Sample_Info
for (i in 1:nrow(Prdm9_Sample_Percent)){
  total = sum(Prdm9_Sample_Percent[i,3:ncol(Prdm9_Sample_Percent)])
  for (n in 1:(ncol(Prdm9_Sample_Percent)-2)){
    rownum = n+2
    Prdm9_Sample_Percent[i,rownum] = Prdm9_Sample_Info[i,rownum] * 100 / total
  }
}

#Convert Baseline dataset to percent
WT_Sample_Percent = WT_Sample_Info
for (i in 1:nrow(WT_Sample_Percent)){
  total = sum(WT_Sample_Percent[i,2:ncol(WT_Sample_Percent)])
  for (n in 1:(ncol(WT_Sample_Percent)-1)){
    rownum = n+1
    WT_Sample_Percent[i,rownum] = WT_Sample_Info[i,rownum] * 100 / total
  }
}

#Add late & mid pachytene together
WT_Sample_Percent$Late.pachytene = WT_Sample_Percent$Late.pachytene + WT_Sample_Percent$Mid.pachytene
WT_Sample_Percent = WT_Sample_Percent[,-8]


##################################################################################################################
############ Create datasets by genotype #########################################################################
###### Prdm9_Sample_Info + WT_Sample_Info -> P9_Hets & P9_Muts & All_WTs #########################################

#Create het dataset
P9_Hets = Prdm9_Sample_Percent[grep("HET",Prdm9_Sample_Percent$Genotype),]

#Create mutant dataset
P9_Muts = Prdm9_Sample_Percent[grep("MUT",Prdm9_Sample_Percent$Genotype),]

#Find WT samples from P9 dataset
P9_WTs = Prdm9_Sample_Percent[grep("WT",Prdm9_Sample_Percent$Genotype),]
P9_WTs = P9_WTs[,-2]

#Add pachytene-like to Baseline dataset
WT_Sample_Percent$Pachytene.LIKE = 0
WT_Sample_Percent= WT_Sample_Percent[,c(1,2,3,4,5,6,10,7,8,9)]
colnames(WT_Sample_Percent) = colnames(P9_WTs)

#Group P9 WTs with Baseline WTs
All_WTs = rbind(P9_WTs,WT_Sample_Percent)
rownames(All_WTs) = c(1:nrow(All_WTs))


##################################################################################################################
############ Create datasets for proportions of celltypes ########################################################
###### P9_Hets & P9_Muts & All_WTs -> Sample_Proportions #########################################################

#Create matrix
Sample_Proportions = matrix(NA,nrow=7,ncol=9)
rownames(Sample_Proportions) = c("Sp","PL","EL","LLZ","EP","LPD","Plike")
colnames(Sample_Proportions) = c("W_08","W_12","W_16","H_08","H_12","H_16","M_08","M_12","M_16")

#Fill matrix with values
Sample_Proportions[1,"M_08"] = mean(P9_Muts[1:5,3])
Sample_Proportions[2,"M_08"] = mean(P9_Muts[1:5,4])
Sample_Proportions[3,"M_08"] = mean(P9_Muts[1:5,5])
Sample_Proportions[4,"M_08"] = mean(P9_Muts[1:5,6] + P9_Muts[1:5,7])
Sample_Proportions[5,"M_08"] = mean(P9_Muts[1:5,9])
Sample_Proportions[6,"M_08"] = mean(P9_Muts[1:5,10] + P9_Muts[1:5,11])
Sample_Proportions[7,"M_08"] = mean(P9_Muts[1:5,8])
Sample_Proportions[1,"M_12"] = mean(P9_Muts[6:9,3])
Sample_Proportions[2,"M_12"] = mean(P9_Muts[6:9,4])
Sample_Proportions[3,"M_12"] = mean(P9_Muts[6:9,5])
Sample_Proportions[4,"M_12"] = mean(P9_Muts[6:9,6] + P9_Muts[6:9,7])
Sample_Proportions[5,"M_12"] = mean(P9_Muts[6:9,9])
Sample_Proportions[6,"M_12"] = mean(P9_Muts[6:9,10] + P9_Muts[6:9,11])
Sample_Proportions[7,"M_12"] = mean(P9_Muts[6:9,8])
Sample_Proportions[1,"M_16"] = mean(P9_Muts[10:14,3])
Sample_Proportions[2,"M_16"] = mean(P9_Muts[10:14,4])
Sample_Proportions[3,"M_16"] = mean(P9_Muts[10:14,5])
Sample_Proportions[4,"M_16"] = mean(P9_Muts[10:14,6] + P9_Muts[10:14,7])
Sample_Proportions[5,"M_16"] = mean(P9_Muts[10:14,9])
Sample_Proportions[6,"M_16"] = mean(P9_Muts[10:14,10] + P9_Muts[10:14,11])
Sample_Proportions[7,"M_16"] = mean(P9_Muts[10:14,8])

Sample_Proportions[1,"H_08"] = mean(P9_Hets[1:2,3])
Sample_Proportions[2,"H_08"] = mean(P9_Hets[1:2,4])
Sample_Proportions[3,"H_08"] = mean(P9_Hets[1:2,5])
Sample_Proportions[4,"H_08"] = mean(P9_Hets[1:2,6] + P9_Hets[1:2,7])
Sample_Proportions[5,"H_08"] = mean(P9_Hets[1:2,9])
Sample_Proportions[6,"H_08"] = mean(P9_Hets[1:2,10] + P9_Hets[1:2,11])
Sample_Proportions[7,"H_08"] = mean(P9_Hets[1:2,8])
Sample_Proportions[1,"H_12"] = mean(P9_Hets[3:5,3])
Sample_Proportions[2,"H_12"] = mean(P9_Hets[3:5,4])
Sample_Proportions[3,"H_12"] = mean(P9_Hets[3:5,5])
Sample_Proportions[4,"H_12"] = mean(P9_Hets[3:5,6] + P9_Hets[3:5,7])
Sample_Proportions[5,"H_12"] = mean(P9_Hets[3:5,9])
Sample_Proportions[6,"H_12"] = mean(P9_Hets[3:5,10] + P9_Hets[3:5,11])
Sample_Proportions[7,"H_12"] = mean(P9_Hets[3:5,8])
Sample_Proportions[1,"H_16"] = mean(P9_Hets[6:8,3])
Sample_Proportions[2,"H_16"] = mean(P9_Hets[6:8,4])
Sample_Proportions[3,"H_16"] = mean(P9_Hets[6:8,5])
Sample_Proportions[4,"H_16"] = mean(P9_Hets[6:8,6] + P9_Hets[6:8,7])
Sample_Proportions[5,"H_16"] = mean(P9_Hets[6:8,9])
Sample_Proportions[6,"H_16"] = mean(P9_Hets[6:8,10] + P9_Hets[6:8,11])
Sample_Proportions[7,"H_16"] = mean(P9_Hets[6:8,8])

Sample_Proportions[1,"W_08"] = mean(All_WTs[c(1:3,7:9),2])
Sample_Proportions[2,"W_08"] = mean(All_WTs[c(1:3,7:9),3])
Sample_Proportions[3,"W_08"] = mean(All_WTs[c(1:3,7:9),4])
Sample_Proportions[4,"W_08"] = mean(All_WTs[c(1:3,7:9),5] + All_WTs[c(1:3,7:9),6])
Sample_Proportions[5,"W_08"] = mean(All_WTs[c(1:3,7:9),8])
Sample_Proportions[6,"W_08"] = mean(All_WTs[c(1:3,7:9),9] + All_WTs[c(1:3,7:9),10])
Sample_Proportions[7,"W_08"] = mean(All_WTs[c(1:3,7:9),7])
Sample_Proportions[1,"W_12"] = mean(All_WTs[c(4,10:14),2])
Sample_Proportions[2,"W_12"] = mean(All_WTs[c(4,10:14),3])
Sample_Proportions[3,"W_12"] = mean(All_WTs[c(4,10:14),4])
Sample_Proportions[4,"W_12"] = mean(All_WTs[c(4,10:14),5] + All_WTs[c(4,10:14),6])
Sample_Proportions[5,"W_12"] = mean(All_WTs[c(4,10:14),8])
Sample_Proportions[6,"W_12"] = mean(All_WTs[c(4,10:14),9] + All_WTs[c(4,10:14),10])
Sample_Proportions[7,"W_12"] = mean(All_WTs[c(4,10:14),7])
Sample_Proportions[1,"W_16"] = mean(All_WTs[c(5,6,15:19),2])
Sample_Proportions[2,"W_16"] = mean(All_WTs[c(5,6,15:19),3])
Sample_Proportions[3,"W_16"] = mean(All_WTs[c(5,6,15:19),4])
Sample_Proportions[4,"W_16"] = mean(All_WTs[c(5,6,15:19),5] + All_WTs[c(5,6,15:19),6])
Sample_Proportions[5,"W_16"] = mean(All_WTs[c(5,6,15:19),8])
Sample_Proportions[6,"W_16"] = mean(All_WTs[c(5,6,15:19),9] + All_WTs[c(5,6,15:19),10])
Sample_Proportions[7,"W_16"] = mean(All_WTs[c(5,6,15:19),7])

#Redo for just p9 dataset
#Create matrix
Prdm9_Proportions = matrix(NA,nrow=7,ncol=9)
rownames(Prdm9_Proportions) = c("Sp","PL","EL","LLZ","EP","LPD","Plike")
colnames(Prdm9_Proportions) = c("W_08","W_12","W_16","H_08","H_12","H_16","M_08","M_12","M_16")

#Fill matrix with values
Prdm9_Proportions[1,"M_08"] = mean(P9_Muts[1:5,3])
Prdm9_Proportions[2,"M_08"] = mean(P9_Muts[1:5,4])
Prdm9_Proportions[3,"M_08"] = mean(P9_Muts[1:5,5])
Prdm9_Proportions[4,"M_08"] = mean(P9_Muts[1:5,6] + P9_Muts[1:5,7])
Prdm9_Proportions[5,"M_08"] = mean(P9_Muts[1:5,9])
Prdm9_Proportions[6,"M_08"] = mean(P9_Muts[1:5,10] + P9_Muts[1:5,11])
Prdm9_Proportions[7,"M_08"] = mean(P9_Muts[1:5,8])
Prdm9_Proportions[1,"M_12"] = mean(P9_Muts[6:9,3])
Prdm9_Proportions[2,"M_12"] = mean(P9_Muts[6:9,4])
Prdm9_Proportions[3,"M_12"] = mean(P9_Muts[6:9,5])
Prdm9_Proportions[4,"M_12"] = mean(P9_Muts[6:9,6] + P9_Muts[6:9,7])
Prdm9_Proportions[5,"M_12"] = mean(P9_Muts[6:9,9])
Prdm9_Proportions[6,"M_12"] = mean(P9_Muts[6:9,10] + P9_Muts[6:9,11])
Prdm9_Proportions[7,"M_12"] = mean(P9_Muts[6:9,8])
Prdm9_Proportions[1,"M_16"] = mean(P9_Muts[10:14,3])
Prdm9_Proportions[2,"M_16"] = mean(P9_Muts[10:14,4])
Prdm9_Proportions[3,"M_16"] = mean(P9_Muts[10:14,5])
Prdm9_Proportions[4,"M_16"] = mean(P9_Muts[10:14,6] + P9_Muts[10:14,7])
Prdm9_Proportions[5,"M_16"] = mean(P9_Muts[10:14,9])
Prdm9_Proportions[6,"M_16"] = mean(P9_Muts[10:14,10] + P9_Muts[10:14,11])
Prdm9_Proportions[7,"M_16"] = mean(P9_Muts[10:14,8])

Prdm9_Proportions[1,"H_08"] = mean(P9_Hets[1:2,3])
Prdm9_Proportions[2,"H_08"] = mean(P9_Hets[1:2,4])
Prdm9_Proportions[3,"H_08"] = mean(P9_Hets[1:2,5])
Prdm9_Proportions[4,"H_08"] = mean(P9_Hets[1:2,6] + P9_Hets[1:2,7])
Prdm9_Proportions[5,"H_08"] = mean(P9_Hets[1:2,9])
Prdm9_Proportions[6,"H_08"] = mean(P9_Hets[1:2,10] + P9_Hets[1:2,11])
Prdm9_Proportions[7,"H_08"] = mean(P9_Hets[1:2,8])
Prdm9_Proportions[1,"H_12"] = mean(P9_Hets[3:5,3])
Prdm9_Proportions[2,"H_12"] = mean(P9_Hets[3:5,4])
Prdm9_Proportions[3,"H_12"] = mean(P9_Hets[3:5,5])
Prdm9_Proportions[4,"H_12"] = mean(P9_Hets[3:5,6] + P9_Hets[3:5,7])
Prdm9_Proportions[5,"H_12"] = mean(P9_Hets[3:5,9])
Prdm9_Proportions[6,"H_12"] = mean(P9_Hets[3:5,10] + P9_Hets[3:5,11])
Prdm9_Proportions[7,"H_12"] = mean(P9_Hets[3:5,8])
Prdm9_Proportions[1,"H_16"] = mean(P9_Hets[6:8,3])
Prdm9_Proportions[2,"H_16"] = mean(P9_Hets[6:8,4])
Prdm9_Proportions[3,"H_16"] = mean(P9_Hets[6:8,5])
Prdm9_Proportions[4,"H_16"] = mean(P9_Hets[6:8,6] + P9_Hets[6:8,7])
Prdm9_Proportions[5,"H_16"] = mean(P9_Hets[6:8,9])
Prdm9_Proportions[6,"H_16"] = mean(P9_Hets[6:8,10] + P9_Hets[6:8,11])
Prdm9_Proportions[7,"H_16"] = mean(P9_Hets[6:8,8])

Prdm9_Proportions[1,"W_08"] = mean(P9_WTs[c(1:3),2])
Prdm9_Proportions[2,"W_08"] = mean(P9_WTs[c(1:3),3])
Prdm9_Proportions[3,"W_08"] = mean(P9_WTs[c(1:3),4])
Prdm9_Proportions[4,"W_08"] = mean(P9_WTs[c(1:3),5] + P9_WTs[c(1:3),6])
Prdm9_Proportions[5,"W_08"] = mean(P9_WTs[c(1:3),8])
Prdm9_Proportions[6,"W_08"] = mean(P9_WTs[c(1:3),9] + P9_WTs[c(1:3),10])
Prdm9_Proportions[7,"W_08"] = mean(P9_WTs[c(1:3),7])
Prdm9_Proportions[1,"W_12"] = mean(P9_WTs[c(4),2])
Prdm9_Proportions[2,"W_12"] = mean(P9_WTs[c(4),3])
Prdm9_Proportions[3,"W_12"] = mean(P9_WTs[c(4),4])
Prdm9_Proportions[4,"W_12"] = mean(P9_WTs[c(4),5] + P9_WTs[c(4),6])
Prdm9_Proportions[5,"W_12"] = mean(P9_WTs[c(4),8])
Prdm9_Proportions[6,"W_12"] = mean(P9_WTs[c(4),9] + P9_WTs[c(4),10])
Prdm9_Proportions[7,"W_12"] = mean(P9_WTs[c(4),7])
Prdm9_Proportions[1,"W_16"] = mean(P9_WTs[c(5,6),2])
Prdm9_Proportions[2,"W_16"] = mean(P9_WTs[c(5,6),3])
Prdm9_Proportions[3,"W_16"] = mean(P9_WTs[c(5,6),4])
Prdm9_Proportions[4,"W_16"] = mean(P9_WTs[c(5,6),5] + P9_WTs[c(5,6),6])
Prdm9_Proportions[5,"W_16"] = mean(P9_WTs[c(5,6),8])
Prdm9_Proportions[6,"W_16"] = mean(P9_WTs[c(5,6),9] + P9_WTs[c(5,6),10])
Prdm9_Proportions[7,"W_16"] = mean(P9_WTs[c(5,6),7])

##################################################################################################################
############ Create plots of proportions #########################################################################
###### Sample_Proportions -> Plots ###############################################################################

#Create vector of cell types..will not be used
names = c("Sp'gonia","Prelep","Early Lep","Late Lep/Zyg","Early Pach","Late Pach/Dip","Pach-Like")

#plot with old colors
# barplot(Sample_Proportions[,c(1:3,7:9)],
#         names.arg=c("WT 8dpp","WT 12dpp","WT 16dpp","Mut 8dpp","Mut 12dpp","Mut 16dpp"),
#         col=c("grey45","mediumseagreen","darkorange2","darkslateblue","limegreen","tan4","white"),
#         ylab = "Average Proportion")

#define colors
col=rep(c("#666666","#1b9e77","#d95f02","#7570b3","#66a61e","#a6761d","#000000"),6)
col_het=rep(c("#666666","#1b9e77","#d95f02","#7570b3","#66a61e","#a6761d","#000000"),3)

#plot with correct colors
barplot(Sample_Proportions[,c(1:3,7:9)],
        names.arg=c("WT 8dpp","WT 12dpp","WT 16dpp","Mut 8dpp","Mut 12dpp","Mut 16dpp"),
        col=c("#666666","#1b9e77","#d95f02","#7570b3","#66a61e","#a6761d","white"),
        ylab = "Average Proportion")

#plot hets
barplot(Sample_Proportions[,c(4:6)],
        names.arg=c("HT 8dpp","HT 12dpp","HT 16dpp"),
        col=c("#666666","#1b9e77","#d95f02","#7570b3","#66a61e","#a6761d","white"),
        ylab = "Average Proportion")

#plot with correct colors
barplot(Prdm9_Proportions[,c(1:3,7:9)],
        names.arg=c("WT 8dpp","WT 12dpp","WT 16dpp","Mut 8dpp","Mut 12dpp","Mut 16dpp"),
        col=c("#666666","#1b9e77","#d95f02","#7570b3","#66a61e","#a6761d","white"),
        ylab = "Average Proportion")

#Convert to a data frame
Sample_Proportions_df = as.data.frame(Sample_Proportions[,c(1:3,7:9)])
colnames(Sample_Proportions_df) = c("WT 8dpp","WT 12dpp","WT 16dpp","Mut 8dpp","Mut 12dpp","Mut 16dpp")

Het_Proportions_df = as.data.frame(Sample_Proportions[,c(4:6)])
colnames(Het_Proportions_df) = c("HT 8dpp","HT 12dpp","HT 16dpp")


#Convert to a data frame
Prdm9_Proportions_df = as.data.frame(Prdm9_Proportions[,c(1:3,7:9)])
colnames(Prdm9_Proportions_df) = c("WT 8dpp","WT 12dpp","WT 16dpp","Mut 8dpp","Mut 12dpp","Mut 16dpp")


#Make melted data table
Sample_Proportions_m <- melt(cbind(Sample_Proportions_df, ind = rownames(Sample_Proportions_df)), id.vars = c('ind'))
colnames(Sample_Proportions_m) = c("Stage","Age","Percent")

Het_Proportions_m <- melt(cbind(Het_Proportions_df, ind = rownames(Het_Proportions_df)), id.vars = c('ind'))
colnames(Het_Proportions_m) = c("Stage","Age","Percent")


#Make melted data table
Prdm9_Proportions_m <- melt(cbind(Prdm9_Proportions_df, ind = rownames(Prdm9_Proportions_df)), id.vars = c('ind'))
colnames(Prdm9_Proportions_m) = c("Stage","Age","Percent")

#Plot proportions
##With plot labels and background grid
ggplot(Sample_Proportions_m,aes(x = Age, y = Percent,fill=col)) + geom_bar(stat="identity",fill=col)+ scale_fill_manual(name="Substage",breaks=c("SP", "PL", "EL","LLZ","EP","EPD","Plike"),labels=c("SP", "PL", "EL","LLZ","EP","EPD","Plike"))

ggplot(Het_Proportions_m,aes(x = Age, y = Percent,fill=col_het)) + geom_bar(stat="identity",fill=col_het)+ scale_fill_manual(name="Substage",breaks=c("SP", "PL", "EL","LLZ","EP","EPD","Plike"),labels=c("SP", "PL", "EL","LLZ","EP","EPD","Plike"))


##With xy only
ggplot(Sample_Proportions_m,aes(x = Age, y = Percent,fill=col)) + geom_bar(stat="identity",fill=col)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title.y=element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(),axis.text.x=element_blank())

##Alone
ggplot(Sample_Proportions_m,aes(x = Age, y = Percent,fill=col)) + geom_bar(stat="identity",fill=col)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+theme(axis.title.y=element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(),axis.text.x=element_blank(),axis.ticks.y=element_blank(),axis.ticks.x=element_blank())

#Sepparate into WT and MT
WT_Sample_Proportions_m = Sample_Proportions_m[grep("WT",Sample_Proportions_m$Age),]
MT_Sample_Proportions_m = Sample_Proportions_m[grep("Mut",Sample_Proportions_m$Age),]

#Define colors for sub-divided  sections
col_half=rep(c("#666666","#1b9e77","#d95f02","#7570b3","#66a61e","#a6761d","#000000"),3)

#Plot WT half with box and tick marks
ggplot(WT_Sample_Proportions_m,aes(x = Age, y = Percent,fill=col_half)) + 
  geom_bar(stat="identity",fill=col_half) + 
  theme(panel.border = element_rect(colour = "grey",fill=NA,size=1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#Plot MT half with box but not tick marks
ggplot(MT_Sample_Proportions_m,aes(x = Age, y = Percent,fill=col_half)) + geom_bar(stat="identity",fill=col_half)+ theme(panel.border = element_rect(colour = "grey",fill=NA,size=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+theme(axis.title.y=element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(),axis.text.x=element_blank(),axis.ticks.y=element_blank(),axis.ticks.x=element_blank())

ggplot(Het_Proportions_m,aes(x = Age, y = Percent,fill=col_half)) + 
  geom_bar(stat="identity",fill=col_half) + 
  theme(panel.border = element_rect(colour = "grey",fill=NA,size=1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


Sample_Proportions_bd = Sample_Proportions_m[c(1:7,22:28,8:14,29:35,15:21,36:42),]
Sample_Proportions_bd[,"DaysOld"] = as.factor(str_sub(Sample_Proportions_bd[,"Age"],-5,-1))

Prdm9_Proportions_bd = Prdm9_Proportions_m[c(1:7,22:28,8:14,29:35,15:21,36:42),]
Prdm9_Proportions_bd[,"DaysOld"] = as.factor(str_sub(Prdm9_Proportions_bd[,"Age"],-5,-1))

x_spaced <- c(1,1.75,2.65,3.4,4.35,5.1)
Sample_Proportions_SP = Sample_Proportions_bd[which(Sample_Proportions_bd$Stage=="Sp"),]
ggplot(transform(Sample_Proportions_SP, x=x_spaced),aes(x = x, y = Percent,width = 0.75)) +
  geom_bar(stat="identity",aes(fill=factor(DaysOld)),color="black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="none") +
  #ylim(0,30) + 
  scale_fill_manual(values=c("#666666", "#666666","#666666")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),       
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
Prdm9_Proportions_SP = Prdm9_Proportions_bd[which(Prdm9_Proportions_bd$Stage=="Sp"),]
ggplot(transform(Prdm9_Proportions_SP, x=x_spaced),aes(x = x, y = Percent,width = 0.75)) +
  geom_bar(stat="identity",aes(fill=factor(DaysOld)),color="black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="none") +
  #ylim(0,30) + 
  scale_fill_manual(values=c("#666666", "#666666","#666666")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),       
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 


Sample_Proportions_PL = Sample_Proportions_bd[which(Sample_Proportions_bd$Stage=="PL"),]
ggplot(transform(Sample_Proportions_PL, x=x_spaced),aes(x = x, y = Percent,width = 0.75)) +
  geom_bar(stat="identity",aes(fill=factor(DaysOld)),color="black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="none") +
  #ylim(0,30) + 
  scale_fill_manual(values=c("#1b9e77", "#1b9e77","#1b9e77")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),       
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
Prdm9_Proportions_PL = Prdm9_Proportions_bd[which(Prdm9_Proportions_bd$Stage=="PL"),]
ggplot(transform(Prdm9_Proportions_PL, x=x_spaced),aes(x = x, y = Percent,width = 0.75)) +
  geom_bar(stat="identity",aes(fill=factor(DaysOld)),color="black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="none") +
  #ylim(0,30) + 
  scale_fill_manual(values=c("#1b9e77", "#1b9e77","#1b9e77")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),       
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 


Sample_Proportions_EL = Sample_Proportions_bd[which(Sample_Proportions_bd$Stage=="EL"),]
ggplot(transform(Sample_Proportions_EL, x=x_spaced),aes(x = x, y = Percent,width = 0.75)) +
  geom_bar(stat="identity",aes(fill=factor(DaysOld)),color="black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="none") +
  #ylim(0,30) + 
  scale_fill_manual(values=c("#d95f02", "#d95f02","#d95f02")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),       
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
Prdm9_Proportions_EL = Prdm9_Proportions_bd[which(Prdm9_Proportions_bd$Stage=="EL"),]
ggplot(transform(Prdm9_Proportions_EL, x=x_spaced),aes(x = x, y = Percent,width = 0.75)) +
  geom_bar(stat="identity",aes(fill=factor(DaysOld)),color="black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="none") +
  #ylim(0,30) + 
  scale_fill_manual(values=c("#d95f02", "#d95f02","#d95f02")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),       
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

Sample_Proportions_LLZ = Sample_Proportions_bd[which(Sample_Proportions_bd$Stage=="LLZ"),]
ggplot(transform(Sample_Proportions_LLZ, x=x_spaced),aes(x = x, y = Percent,width = 0.75)) +
  geom_bar(stat="identity",aes(fill=factor(DaysOld)),color="black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="none") +
  #ylim(0,30) + 
  scale_fill_manual(values=c("#7570b3", "#7570b3","#7570b3")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),       
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
Prdm9_Proportions_LLZ = Prdm9_Proportions_bd[which(Prdm9_Proportions_bd$Stage=="LLZ"),]
ggplot(transform(Prdm9_Proportions_LLZ, x=x_spaced),aes(x = x, y = Percent,width = 0.75)) +
  geom_bar(stat="identity",aes(fill=factor(DaysOld)),color="black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="none") +
  #ylim(0,30) + 
  scale_fill_manual(values=c("#7570b3", "#7570b3","#7570b3")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),       
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

Sample_Proportions_EP = Sample_Proportions_bd[which(Sample_Proportions_bd$Stage=="EP"),]
ggplot(transform(Sample_Proportions_EP, x=x_spaced),aes(x = x, y = Percent,width = 0.75)) +
  geom_bar(stat="identity",aes(fill=factor(DaysOld)),color="black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="none") +
  #ylim(0,30) + 
  scale_fill_manual(values=c("#66a61e", "#66a61e","#66a61e")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),       
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
Prdm9_Proportions_EP = Prdm9_Proportions_bd[which(Prdm9_Proportions_bd$Stage=="EP"),]
ggplot(transform(Prdm9_Proportions_EP, x=x_spaced),aes(x = x, y = Percent,width = 0.75)) +
  geom_bar(stat="identity",aes(fill=factor(DaysOld)),color="black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="none") +
  #ylim(0,30) + 
  scale_fill_manual(values=c("#66a61e", "#66a61e","#66a61e")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),       
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

Sample_Proportions_LPD = Sample_Proportions_bd[which(Sample_Proportions_bd$Stage=="LPD"),]
ggplot(transform(Sample_Proportions_LPD, x=x_spaced),aes(x = x, y = Percent,width = 0.75)) +
  geom_bar(stat="identity",aes(fill=factor(DaysOld)),color="black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="none") +
  #ylim(0,30) + 
  scale_fill_manual(values=c("#a6761d", "#a6761d","#a6761d")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),       
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
Prdm9_Proportions_LPD = Prdm9_Proportions_bd[which(Prdm9_Proportions_bd$Stage=="LPD"),]
ggplot(transform(Prdm9_Proportions_LPD, x=x_spaced),aes(x = x, y = Percent,width = 0.75)) +
  geom_bar(stat="identity",aes(fill=factor(DaysOld)),color="black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="none") +
  #ylim(0,30) + 
  scale_fill_manual(values=c("#a6761d", "#a6761d","#a6761d")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),       
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

Sample_Proportions_PLike = Sample_Proportions_bd[which(Sample_Proportions_bd$Stage=="Plike"),]
ggplot(transform(Sample_Proportions_PLike, x=x_spaced),aes(x = x, y = Percent,width = 0.75)) +
  geom_bar(stat="identity",aes(fill=factor(DaysOld)),color="black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="none") +
  #ylim(0,30) + 
  scale_fill_manual(values=c("#000000", "#000000","#000000")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),       
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
Prdm9_Proportions_PLike = Prdm9_Proportions_bd[which(Prdm9_Proportions_bd$Stage=="Plike"),]
ggplot(transform(Prdm9_Proportions_PLike, x=x_spaced),aes(x = x, y = Percent,width = 0.75)) +
  geom_bar(stat="identity",aes(fill=factor(DaysOld)),color="black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.position="none") +
  #ylim(0,30) + 
  scale_fill_manual(values=c("#000000", "#000000","#000000")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),       
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 


Sample_Proportions_SP$Age <- factor(Sample_Proportions_SP$Age, levels = Sample_Proportions_SP$Age[order(Sample_Proportions_SP$DaysOld)])

Sample_Proportions_d8 = Sample_Proportions_bd[which(Sample_Proportions_bd$DaysOld ==" 8dpp"),]
Sample_Proportions_d12 = Sample_Proportions_bd[which(Sample_Proportions_bd$DaysOld =="12dpp"),]
Sample_Proportions_d16 = Sample_Proportions_bd[which(Sample_Proportions_bd$DaysOld =="16dpp"),]

Prdm9_Proportions_d8 = Prdm9_Proportions_bd[which(Prdm9_Proportions_bd$DaysOld ==" 8dpp"),]
Prdm9_Proportions_d12 = Prdm9_Proportions_bd[which(Prdm9_Proportions_bd$DaysOld =="12dpp"),]
Prdm9_Proportions_d16 = Prdm9_Proportions_bd[which(Prdm9_Proportions_bd$DaysOld =="16dpp"),]


col_two=rep(c("#666666","#1b9e77","#d95f02","#7570b3","#66a61e","#a6761d","#000000"),2)

ggplot(Sample_Proportions_d8,aes(x = Age, y = Percent,fill=col_two)) +
  geom_bar(stat="identity",fill=col_two) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())
ggplot(Prdm9_Proportions_d8,aes(x = Age, y = Percent,fill=col_two)) +
  geom_bar(stat="identity",fill=col_two) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

ggplot(Sample_Proportions_d12,aes(x = Age, y = Percent,fill=col_two)) +
  geom_bar(stat="identity",fill=col_two) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())
ggplot(Prdm9_Proportions_d12,aes(x = Age, y = Percent,fill=col_two)) +
  geom_bar(stat="identity",fill=col_two) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

ggplot(Sample_Proportions_d16,aes(x = Age, y = Percent,fill=col_two)) +
  geom_bar(stat="identity",fill=col_two) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())
ggplot(Prdm9_Proportions_d16,aes(x = Age, y = Percent,fill=col_two)) +
  geom_bar(stat="identity",fill=col_two) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())
