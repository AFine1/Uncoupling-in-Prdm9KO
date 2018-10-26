#   The purpose of this file is to find changing DSBs and H3K4me3 in the Prdm9-/- samples

#   It takes a data file that contains H3K4me3 peaks and DSB locations in WT and Prdm9-/- germ
#   cell samples and identifies overlapping peaks, changes in peaks, and identifies genes 
#   nearby peaks. The resulting data file is saved in /Users/s-fine/Desktop/Carter/Projects/P9/KO Hotspot/Results

##################################################################################################################
############ Load packages #######################################################################################
##################################################################################################################

library(ChIPseeker)
library(IRanges)
library(GenomicRanges)
library(biomaRt)
library(ggplot2)

##################################################################################################################
############ Aquire Data #########################################################################################
###### HotspotData ###############################################################################################

#Set working directory
setwd("/Users/s-fine/Desktop/Carter/Projects/P9/KO Hotspot/")

#Get data
setwd("Data/")
AffinityData = read.csv("Affinity_seq.csv",header=T)
P9Binding = read.csv("Affinity_seq_info.csv",header=T)
B6_Hotspots = read.csv("B6_Hotspots.csv")
B6_Hotspots$type = "Hotspot"
B6_H3K4me3 = read.csv("B6_H3K4me3.csv")
B6_H3K4me3$type = "H3K4me3"
P9_Hotspots = read.csv("P9_Hotspots.csv")
P9_Hotspots$type = "Hotspot"
P9_H3K4me3 = read.csv("P9_H3K4me3.csv")
P9_H3K4me3$type = "H3K4me3"
PolII = read.csv("testes.polII.peak.csv",header=F)
colnames(PolII) = c("chromosome","location")
setwd("..")
P9_mm10_Affin = read.csv("/Users/s-fine/Desktop/Carter/Projects/Stag3_TB/H3K4me3/Hotspots/B6_affinseq_total_e-2_reads_mm10.csv")
load("/Users/s-fine/Desktop/Carter/Projects/P9/3 DifferentialExpression/Results/P9K_Analysis.rdt")

##################################################################################################################
############ Find Overlapping Peaks ##############################################################################
###### HotspotData ###############################################################################################

#Remove duplicates
B6_Hotspots_c = B6_Hotspots[rownames(unique(B6_Hotspots[,1:3])),]
B6_H3K4me3_c = B6_H3K4me3[rownames(unique(B6_H3K4me3[,1:3])),]
P9_Hotspots_c = P9_Hotspots[rownames(unique(P9_Hotspots[,1:3])),]
P9_H3K4me3_c = P9_H3K4me3[rownames(unique(P9_H3K4me3[,1:3])),]

#Make a range column
B6_Hotspots_c$range = 0
B6_H3K4me3_c$range = 0
P9_Hotspots_c$range = 0
P9_H3K4me3_c$range = 0

#rename columns
colnames(B6_Hotspots_c) = c("chromosome","peak_start","peak_end","peak","strain","type","sample","range")
colnames(B6_H3K4me3_c) = c("chromosome","peak_start","peak_end","peak","strain","type","sample","range")
colnames(P9_Hotspots_c) = c("chromosome","peak_start","peak_end","peak","strain","type","sample","range")
colnames(P9_H3K4me3_c) = c("chromosome","peak_start","peak_end","peak","strain","type","sample","range")

#Rename rows
rownames(B6_Hotspots_c) = 1:nrow(B6_Hotspots_c)
rownames(B6_H3K4me3_c) = 1:nrow(B6_H3K4me3_c)
rownames(P9_Hotspots_c) = 1:nrow(P9_Hotspots_c)
rownames(P9_H3K4me3_c) = 1:nrow(P9_H3K4me3_c)

#Add range values
B6_H3K4me3_c$range = B6_H3K4me3_c$peak_end - B6_H3K4me3_c$peak_start
B6_Hotspots_c$range = B6_Hotspots_c$peak_end - B6_Hotspots_c$peak_start
P9_H3K4me3_c$range = P9_H3K4me3_c$peak_end - P9_H3K4me3_c$peak_start
P9_Hotspots_c$range = P9_Hotspots_c$peak_end - P9_Hotspots_c$peak_start

#make 1kb ranges
B6_Hotspots_1K = B6_Hotspots_c
B6_Hotspots_1K$peak_start = B6_Hotspots_c$peak_start + 500
B6_Hotspots_1K$peak_end = B6_Hotspots_c$peak_end - 500
B6_Hotspots_1K$range = B6_Hotspots_1K$peak_end - B6_Hotspots_1K$peak_start

P9_Hotspots_1K = P9_Hotspots_c
P9_Hotspots_1K$peak_start = P9_Hotspots_c$peak_start + 500
P9_Hotspots_1K$peak_end = P9_Hotspots_c$peak_end - 500
P9_Hotspots_1K$range = P9_Hotspots_1K$peak_end - P9_Hotspots_1K$peak_start

B6_H3K4me3_1K = B6_H3K4me3_c
for (i in 1:nrow(B6_H3K4me3_1K)){
  if (B6_H3K4me3_1K[i,"range"]>1000){
    size = round((B6_H3K4me3_1K[i,"range"] - 1000)/2)
    B6_H3K4me3_1K[i,"peak_start"] = B6_H3K4me3_1K[i,"peak_start"] + size
    B6_H3K4me3_1K[i,"peak_end"] = B6_H3K4me3_1K[i,"peak_end"] - size
  }
  else{
  }
}
B6_H3K4me3_1K$range = B6_H3K4me3_1K$peak_end - B6_H3K4me3_1K$peak_start

P9_H3K4me3_1K = P9_H3K4me3_c
for (i in 1:nrow(P9_H3K4me3_1K)){
  if (P9_H3K4me3_1K[i,"range"]>1000){
    size = round((P9_H3K4me3_1K[i,"range"] - 1000)/2)
    P9_H3K4me3_1K[i,"peak_start"] = P9_H3K4me3_1K[i,"peak_start"] + size
    P9_H3K4me3_1K[i,"peak_end"] = P9_H3K4me3_1K[i,"peak_end"] - size
  }
  else{
  }
}
P9_H3K4me3_1K$range = P9_H3K4me3_1K$peak_end - P9_H3K4me3_1K$peak_start

#Make table of all P9 sites we have listed as getting H3K4me3
P9H3K4me3_is = P9Binding[which(P9Binding$is_P9_H3K4me3=="Y"),]
rownames(P9H3K4me3_is) = 1:nrow(P9H3K4me3_is)

#Make a list of the known P9 binding sites
grP9Binding = with(P9Binding, GRanges(chr, IRanges(start=start, end=end)))

#Make a list of the known P9 affinity seq sites (maybe)
grP9_mm10_Affin = with(P9_mm10_Affin, GRanges(chr, IRanges(start=start, end=end)))

#Make a list of the known P9 H3K4me3 sites
grP9H3K4me3= with(P9H3K4me3_is, GRanges(chr, IRanges(start=start, end=end)))


#Make lists for the WT H3K4me3 and DSB sites
grB6_H3K4me3_c = with(B6_H3K4me3_c, GRanges(chromosome, IRanges(start=peak_start, end=peak_end)))
grB6_Hotspots_c = with(B6_Hotspots_c, GRanges(chromosome, IRanges(start=peak_start, end=peak_end)))
grB6_H3K4me3_1K = with(B6_H3K4me3_1K, GRanges(chromosome, IRanges(start=peak_start, end=peak_end)))
grB6_Hotspots_1K = with(B6_Hotspots_1K, GRanges(chromosome, IRanges(start=peak_start, end=peak_end)))

#Make lists for the P9KO H3K4me3 and DSB sites
grP9_H3K4me3_c = with(P9_H3K4me3_c, GRanges(chromosome, IRanges(start=peak_start, end=peak_end)))
grP9_Hotspots_c = with(P9_Hotspots_c, GRanges(chromosome, IRanges(start=peak_start, end=peak_end)))
grP9_H3K4me3_1K = with(P9_H3K4me3_1K, GRanges(chromosome, IRanges(start=peak_start, end=peak_end)))
grP9_Hotspots_1K = with(P9_Hotspots_1K, GRanges(chromosome, IRanges(start=peak_start, end=peak_end)))

#Find overlaps between WT DSBs and WT H3K4me3
SharedgrB6_Hotspots_c.grB6_H3K4me3_c = findOverlaps(grB6_Hotspots_c, grB6_H3K4me3_c) #Full length
SharedgrB6_Hotspots_c.grB6_H3K4me3_1K = findOverlaps(grB6_Hotspots_c, grB6_H3K4me3_1K) #1KB region
#The percent of WT DSBs that happen at WT H3K4me3 sites
100*length(unique(SharedgrB6_Hotspots_c.grB6_H3K4me3_c@queryHits)) / length(grB6_Hotspots_c)
100*length(unique(SharedgrB6_Hotspots_c.grB6_H3K4me3_1K@queryHits)) / length(grB6_Hotspots_c)

#Find overlaps between P9KO DSBs and P9KO H3K4me3
SharedgrP9_Hotspots_c.grP9_H3K4me3_c = findOverlaps(grP9_Hotspots_c, grP9_H3K4me3_c)
#The percent of P9KO DSBs that happen at P9KO H3K4me3 sites
100*length(unique(SharedgrP9_Hotspots_c.grP9_H3K4me3_c@queryHits)) / length(grP9_Hotspots_c)

#Find overlaps between WT DSBs and P9 Binding Sites
SharedgrB6_Hotspots_c.grP9_mm10_Affin = findOverlaps(grB6_Hotspots_c, grP9_mm10_Affin)
SharedgrB6_Hotspots_1K.grP9_mm10_Affin = findOverlaps(grB6_Hotspots_1K, grP9_mm10_Affin)
#The percent of WT DSBs that happen at P9 Binding Sites
100*length(unique(SharedgrB6_Hotspots_c.grP9_mm10_Affin@queryHits)) / length(grB6_Hotspots_c)
100*length(unique(SharedgrB6_Hotspots_1K.grP9_mm10_Affin@queryHits)) / length(grB6_Hotspots_1K)

#Find overlaps between P9KO DSBs and P9 Binding Sites
SharedgrP9_Hotspots_c.grP9_mm10_Affin = findOverlaps(grP9_Hotspots_c, grP9_mm10_Affin)
SharedgrP9_Hotspots_1K.grP9_mm10_Affin = findOverlaps(grP9_Hotspots_1K, grP9_mm10_Affin)
#The percent of P9KO DSBs that happen at P9 Binding Sites
100*length(unique(SharedgrP9_Hotspots_c.grP9_mm10_Affin@queryHits)) / length(grP9_Hotspots_c)
100*length(unique(SharedgrP9_Hotspots_1K.grP9_mm10_Affin@queryHits)) / length(grP9_Hotspots_1K)

#Define P9 Activated Hotspots in WT
grB6_Hotspots_A = grB6_Hotspots_1K[SharedgrB6_Hotspots_1K.grP9_mm10_Affin@queryHits]

#Define P9 Activated Hotspots in P9KO
grP9_Hotspots_A = grP9_Hotspots_1K[SharedgrP9_Hotspots_1K.grP9_mm10_Affin@queryHits]



#Load gene info
m_ensembl=useMart("ensembl",dataset="mmusculus_gene_ensembl")
#listAttributes(m_ensembl)[1:30,]
gene.TSS_info=getBM(attributes=c("external_gene_name","ensembl_gene_id","chromosome_name","transcription_start_site"),mart=m_ensembl)

#Make list of chromosomes
chromlist = c("1","2","3","4","5","6","7","8","9","10","11","12",
              "13","14","15","16","17","18","19","X","Y")

#Pull out only the 'real' chromosomes
TSS_info = gene.TSS_info[which(gene.TSS_info$chromosome_name %in% chromlist),]

#MAke it into a new file with ranges
tss_table = matrix(NA,nrow=nrow(TSS_info),ncol=4)
rownames(tss_table) = 1:nrow(tss_table)
colnames(tss_table) = c("Gene_Name","chr","start","end")
tss_table[,1] = as.character(TSS_info$external_gene_name)
tss_table[,2] = paste("chr",TSS_info$chromosome_name,sep="")
tss_table[,3] = as.numeric(TSS_info$transcription_start_site) - 1000
tss_table[,4] = as.numeric(TSS_info$transcription_start_site) + 1000

#Make them the right data inputs
tss_table = as.data.frame(tss_table)
tss_table$Gene_Name = as.character(tss_table$Gene_Name)
tss_table$chr = as.factor(tss_table$chr)
tss_table$start = as.numeric(as.character(tss_table$start))
tss_table$end = as.numeric(as.character(tss_table$end))

#Define ranges
grTSS = with(tss_table, GRanges(chr, IRanges(start=start, end=end)))

#Find overlaps between WT DSBs and Promoters (TSS +-1KB)
SharedgrB6_Hotspots_1K.grTSS = findOverlaps(grB6_Hotspots_1K, grTSS)
#Find the percent of WT DSBs at Promoters
100*length(unique(SharedgrB6_Hotspots_1K.grTSS@queryHits))/length(grB6_Hotspots_1K)

#Find overlaps between P9KO DSBs and Promoters (TSS +-1KB)
SharedgrP9_Hotspots_1K.grTSS = findOverlaps(grP9_Hotspots_1K, grTSS)
#Find the percent of P9KO DSBs at Promoters
100*length(unique(SharedgrP9_Hotspots_1K.grTSS@queryHits))/length(grP9_Hotspots_1K)

#Identify genes that are at DSB loci in both mice
P9_DSB_table = unique(tss_table[SharedgrP9_Hotspots_1K.grTSS@subjectHits,]$Gene_Name)
B6_DSB_table = unique(tss_table[SharedgrB6_Hotspots_1K.grTSS@subjectHits,]$Gene_Name)

#Find which genes are in the P9 Analysis
P9_DSB_A = P9_DSB_table[which(P9_DSB_table %in% rownames(P9K_Analysis))]
B6_DSB_A = B6_DSB_table[which(B6_DSB_table %in% rownames(P9K_Analysis))]

#Make a table with peak and gene name info for P9KO
P9gene_sites = P9_Hotspots_1K[SharedgrP9_Hotspots_1K.grTSS@queryHits,]
P9genes = tss_table[SharedgrP9_Hotspots_1K.grTSS@subjectHits,]$Gene_Name
P9gene_sites$P9Genes = P9genes

#Make a table with peak and gene name info for B6
B6gene_sites = B6_Hotspots_1K[SharedgrB6_Hotspots_1K.grTSS@queryHits,]
B6genes = tss_table[SharedgrB6_Hotspots_1K.grTSS@subjectHits,]$Gene_Name
B6gene_sites$B6Genes = B6genes

#Make a file for max peak by gene for B6
B6_Gene_peak = matrix(NA, ncol=1, nrow = length(unique(B6gene_sites$B6Genes)))
colnames(B6_Gene_peak) = c("Peak")
rownames(B6_Gene_peak) <- unique(B6gene_sites$B6Genes)
for (i in 1:nrow(B6_Gene_peak)){
  gene_id = rownames(B6_Gene_peak)[i]
  B6_Gene_peak[i,1] = max(B6gene_sites[which(B6gene_sites$B6Genes == gene_id),"peak"])
}

#Make a file for max peak by gene for P9KO
P9_Gene_peak = matrix(NA, ncol=1, nrow = length(unique(P9gene_sites$P9Genes)))
colnames(P9_Gene_peak) = c("Peak")
rownames(P9_Gene_peak) <- unique(P9gene_sites$P9Genes)
for (i in 1:nrow(P9_Gene_peak)){
  gene_id = rownames(P9_Gene_peak)[i]
  P9_Gene_peak[i,1] = max(P9gene_sites[which(P9gene_sites$P9Genes == gene_id),"peak"])
}

#Add Day12 and Day16 Coef and pval and fdr
P9K_DSB = P9K_Analysis
P9K_DSB$B6peak = 0
P9K_DSB$P9peak = 0

for (i in 1:nrow(P9K_DSB)){
  gene_id = rownames(P9K_DSB)[i]
  if (gene_id %in% rownames(P9_Gene_peak)){
    P9K_DSB[,"P9peak"][i] = P9_Gene_peak[gene_id,]
  }
  if (gene_id %in% rownames(B6_Gene_peak)){
    P9K_DSB[,"B6peak"][i] = B6_Gene_peak[gene_id,]
  }
  else{
  }
  print(i)
}

plot(P9K_DSB$B6peak,-log10(P9K_DSB$S_P9K_PVFA_12D.W_12D.M))
plot(P9K_DSB$P9peak,-log10(P9K_DSB$S_P9K_PVFA_12D.W_12D.M))

cor.test(P9K_DSB[which(P9K_DSB$P9peak!=0 & P9K_DSB$S_P9K_COEF_16D.W_16D.M > 0.5 & P9K_DSB$S_P9K_PVFA_16D.W_16D.M < 0.01),"P9peak"],
         abs(P9K_DSB[which(P9K_DSB$P9peak!=0 & P9K_DSB$S_P9K_COEF_16D.W_16D.M > 0.5& P9K_DSB$S_P9K_PVFA_16D.W_16D.M < 0.01),"S_P9K_COEF_16D.W_16D.M"]),
         method = "spearman")



ggplot(data=P9K_DSB[which(P9K_DSB$P9peak!=0),],aes(x=P9K_DSB[which(P9K_DSB$P9peak!=0),"P9peak"],y=abs(P9K_DSB[which(P9K_DSB$P9peak!=0),"S_P9K_COEF_16D.W_16D.M"]))) +
  geom_point() +
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

P9K_DSB_ch = P9K_DSB[which(P9K_DSB$S_P9K_PVFA_16D.W_16D.M<0.01),]
P9K_DSB_nch = P9K_DSB[-which(P9K_DSB$S_P9K_PVFA_16D.W_16D.M<0.01),]
P9K_DSB_p = P9K_DSB[which(P9K_DSB$P9peak!=0),]

ggplot(data=P9K_DSB_p,aes(x=P9peak,y=-log10(P9K_DSB_p$S_P9K_PVFA_16D.W_16D.M))) +
    geom_point() +
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
ggplot(data=P9K_DSB_p,aes(x=P9peak,y=-log10(P9K_DSB_p$S_P9K_PVFA_16D.W_16D.M))) +
  geom_density2d() +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
ggplot(data=P9K_DSB_p,aes(x=P9peak,y=-log10(P9K_DSB_p$S_P9K_PVFA_16D.W_16D.M))) +
  ylim(-0.4,5.5) + 
  xlim(-0.4,1500) + 
  theme_bw() +
  stat_density2d(aes(fill=..level..),geom="polygon") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggplot(data=P9K_DSB_p,aes(x=P9peak,y=-log10(P9K_DSB_p$S_P9K_PVFA_16D.W_16D.M))) +
  ylim(-4,12) + 
  xlim(-4,8000) + 
  geom_point() +
  theme_bw() +
  stat_density2d(aes(fill=..level..),geom="polygon") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggplot(data=P9K_DSB_p,aes(x=P9peak,y=-P9K_DSB_p$S_P9K_COEF_16D.W_16D.M)) +
  ylim(-1,1) + 
  xlim(-0.4,2000) + 
  theme_bw() +
  stat_density2d(aes(fill=..level..),geom="polygon") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggplot(data=P9K_DSB_p,aes(x=P9peak,y=-log10(P9K_DSB_p$S_P9K_PVFA_16D.W_16D.M))) +
  ylim(-0.4,4.5) + 
  xlim(-0.4,1500) + 
  theme_bw() +
  stat_density2d(aes(fill=..level..),geom="polygon") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
  
plot(P9K_DSB$B6peak,-log10(P9K_DSB$S_P9K_PVFA_16D.W_16D.M))
plot(P9K_DSB$P9peak,-log10(P9K_DSB$S_P9K_PVFA_16D.W_16D.M))

plot(P9K_DSB$B6peak,abs(P9K_DSB$S_P9K_COEF_12D.W_12D.M))
plot(P9K_DSB$P9peak,abs(P9K_DSB$S_P9K_COEF_12D.W_12D.M))

plot(P9K_DSB$B6peak,abs(P9K_DSB$S_P9K_COEF_16D.W_16D.M))
plot(P9K_DSB$B6peak,P9K_DSB$S_P9K_COEF_16D.W_16D.M)
plot(P9K_DSB$P9peak,abs(P9K_DSB$S_P9K_COEF_16D.W_16D.M))
plot(P9K_DSB$P9peak,P9K_DSB$S_P9K_COEF_16D.W_16D.M)

DSB_FT = matrix(NA,nrow=2,ncol=2)



DSB_FT[1,1] = length(which(P9K_DSB$P9peak!="0" & P9K_DSB$S_P9K_PVFA_16D.W_16D.M<0.01 & abs(P9K_DSB$S_P9K_COEF_16D.W_16D.M) > 0.5))
DSB_FT[2,1] = length(which(P9K_DSB$P9peak!="0" & (P9K_DSB$S_P9K_PVFA_16D.W_16D.M>=0.01 | abs(P9K_DSB$S_P9K_COEF_16D.W_16D.M) <= 0.5)))

DSB_FT[1,2] = length(which(P9K_DSB$P9peak=="0" & P9K_DSB$S_P9K_PVFA_16D.W_16D.M<0.01 & abs(P9K_DSB$S_P9K_COEF_16D.W_16D.M) > 0.5))
DSB_FT[2,2] = length(which(P9K_DSB$P9peak=="0" & (P9K_DSB$S_P9K_PVFA_16D.W_16D.M>=0.01 | abs(P9K_DSB$S_P9K_COEF_16D.W_16D.M) <= 0.5)))

fisher.test(DSB_FT,alternative = "greater")
