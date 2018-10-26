#This file plots the expression of key meiotic genes, as well as sex-linked genes.

#Set up directory
setwd("/Users/s-fine/Desktop/Carter/Projects/P9/2.5 BiologicalValidation/Sex Expression/")
library(ggplot2)

#Load data
setwd("Data/")
load(file = "P9K_Analysis.rdt")
setwd("../")


head(P9K_Analysis)
genes = rownames(P9K_Analysis)
Chrom = P9K_Analysis$I_CNUM
GGTest = data.frame(genes,Chrom)

numbers = c(1:21)

sex = c("X","Y","x","y")
test = P9K_Analysis
for (i in 1:nrow(P9K_Analysis)){
  if (P9K_Analysis[i,"I_CNUM"] %in% numbers){
    test[i,"chrom_char"] = "Autosomes"
  }
  else if (P9K_Analysis[i,"I_CNUM"] %in% sex){
    test[i,"chrom_char"] = "XY"
  }
  else {
    test[i,"chrom_char"] = "Unknown"
  }
}

for (i in 1:nrow(test)){
  test[i,"mean_8wt"] = rowMeans(test[i,c(1,2,3,11,12,13)])
}
rowMeans(test[1,c(1,2,3,11,12,13)])

Analysis = melt(cbind(cbind(cbind(P9K_Analysis[,c(1:42,48,51,54,57)],gene=rownames(P9K_Analysis)),chrom_char=test$chrom_char),mean_8wt = test$mean_8wt),id.vars=c('gene','I_CNUM','S_RYT_SpecSub','chrom_char','mean_8wt','S_P9K_COEF_08D.W_08D.M','S_P9K_COEF_12D.W_12D.M','S_P9K_COEF_16D.W_16D.M'))

Analysis$genotype = substr(Analysis$variable,7,7)
Analysis$age = substr(Analysis$variable,9,10)

Analysis$genoage = substr(Analysis$variable,7,10)

Analysis$diff_8wt = Analysis$value - Analysis$mean_8wt

for (i in 1:nrow(Analysis)){
  if (Analysis[i,"age"] == "08"){
    Analysis[i,"coefficient"] = Analysis[i,"S_P9K_COEF_08D.W_08D.M"]
  }
  else if (Analysis[i,"age"] == "12"){
    Analysis[i,"coefficient"] = Analysis[i,"S_P9K_COEF_12D.W_12D.M"]
  }
  else {
    Analysis[i,"coefficient"] = Analysis[i,"S_P9K_COEF_16D.W_16D.M"]
  }
  print(paste(i,"1060219"))
}


setwd("Results/")
save(Analysis,file="Analysis.rdt")
load("Analysis.rdt")
setwd("../")

head(Analysis)

Prdm9_A = Analysis[which(Analysis$gene=="Prdm9"),]
Prdm9_A$genotype2 <- factor(Prdm9_A$genotype, c("W", "H", "M"))
Prdm9plot <- ggplot(Prdm9_A, aes(factor(age), value))


Prdm9plot + geom_dotplot(binaxis = "y",
                         stackdir = "center",
                         aes(fill = factor(genotype2)),
                         position="dodge",
                         binwidth = 0.09) +
  scale_fill_manual(values=c("blue", "orchid","red"))+
  ylim(2.5,6) + 
  theme_bw()+
  theme(legend.position="none") #+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1))



Crem_A = Analysis[which(Analysis$gene=="Crem"),]
Crem_A$genotype2 <- factor(Crem_A$genotype, c("W", "H", "M"))
Cremplot <- ggplot(Crem_A, aes(factor(age), value))


Cremplot + geom_dotplot(binaxis = "y",
                         stackdir = "center",
                         aes(fill = factor(genotype2)),
                         position="dodge",
                         binwidth = 0.09) +
  scale_fill_manual(values=c("blue", "orchid","red"))+
  ylim(3,7) + 
  theme_bw()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1))

Nonm_A = Analysis[which(Analysis$gene=="NONMMUG017504"),]
Nonm_A$genotype2 <- factor(Nonm_A$genotype, c("W", "H", "M"))
Nonm_Aplot <- ggplot(Nonm_A, aes(factor(age), value))
Nonm_Aplot + geom_dotplot(binaxis = "y",
                        stackdir = "center",
                        aes(fill = factor(genotype2)),
                        position="dodge",
                        binwidth = 0.07) +
  scale_fill_manual(values=c("blue", "orchid","red"))+
  ylim(0,3) + 
  theme_bw()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plotgenedotP9 <- function(gene_name,withlabels=FALSE){
  load("/Users/s-fine/Desktop/Carter/Projects/P9/2.5 BiologicalValidation/Sex Expression/Results/Analysis.rdt")
  gene_data = Analysis[which(Analysis$gene==gene_name),]
  gene_data$genotype2 <- factor(gene_data$genotype, c("W", "H", "M"))
  minimum = min(gene_data$value) - sd(gene_data$value)
  maximum = max(gene_data$value) + sd(gene_data$value)
  gene_plot <- ggplot(gene_data, aes(factor(age), value))
  if(withlabels==FALSE){
    gene_plot + geom_dotplot(binaxis = "y",
                             stackdir = "center",
                             aes(fill = factor(genotype2)),
                             position="dodge",
                             binwidth = as.numeric(sd(gene_data$value))/10) +
      scale_fill_manual(values=c("blue", "orchid","red"))+
      ylim(minimum,maximum) + 
      theme_bw()+
      theme(legend.position="none") +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_line(colour="black"),
            panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
  }
  else{
    gene_plot + geom_dotplot(binaxis = "y",
                             stackdir = "center",
                             aes(fill = factor(genotype2)),
                             position="dodge",
                             binwidth = as.numeric(sd(gene_data$value))/10) +
      scale_fill_manual(values=c("blue", "orchid","red"))+
      ylim(minimum,maximum) + 
      theme_bw()+
      theme(legend.position="none") +
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks=element_line(colour="black"),
            panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
  }
}

setwd("~/Desktop/Carter/Projects/P9/Paper/Figures/Figure 5")

dev.off()
pdf("Dmc1_wl.pdf",width = 4,height = 3)
plotgenedotP9("Dmc1",withlabels = TRUE)
dev.off()
pdf("Dmc1_nl.pdf",width = 4,height = 3)
plotgenedotP9("Dmc1",withlabels = FALSE)
dev.off()

Nonm_A = Analysis[which(Analysis$gene=="NONMMUG017504"),]
Nonm_A$genotype2 <- factor(Nonm_A$genotype, c("W", "H", "M"))
Nonm_Aplot <- ggplot(Nonm_A, aes(factor(age), value))
Nonm_Aplot + geom_dotplot(binaxis = "y",
                          stackdir = "center",
                          aes(fill = factor(genotype2)),
                          position="dodge",
                          binwidth = as.numeric(sd(Nonm_A$value))/10) +
  scale_fill_manual(values=c("blue", "orchid","red"))+
  ylim(0,3) + 
  theme_bw()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())




Zfy2_A = Analysis[which(Analysis$gene=="Zfy2"),]
Zfy2_A$genotype2 <- factor(Zfy2_A$genotype, c("W", "H", "M"))
Zfy2plot <- ggplot(Zfy2_A, aes(factor(age), value))



Zfy2plot + geom_dotplot(binaxis = "y",
                         stackdir = "center",
                         aes(fill = factor(genotype2)),
                         position="dodge",
                         binwidth = 0.09) +
  scale_fill_manual(values=c("blue", "orchid","red"))+
  ylim(0,5) + 
  theme_bw()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1))

YY1_A = Analysis[which(Analysis$gene=="Yy1"),]
YY1_A$genotype2 <- factor(YY1_A$genotype, c("W", "H", "M"))
YY1plot <- ggplot(YY1_A, aes(factor(age), value))

YY1plot + geom_dotplot(binaxis = "y",
                          stackdir = "center",
                          aes(fill = factor(genotype2)),
                          position="dodge",
                          binwidth = 0.06) +
  scale_fill_manual(values=c("blue", "orchid","red"))+
  ylim(5,7) + 
  theme_bw()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1))

Morc2bplot + geom_dotplot(binaxis = "y",
                         stackdir = "center",
                         aes(fill = factor(genotype2)),
                         position="dodge",
                         binwidth = 0.15) +
  scale_fill_manual(values=c("blue", "orchid","red"))+
  ylim(-1,6) + 
  theme_bw()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1))

Prdm9plot + geom_boxplot(aes(fill = factor(genotype2))) + 
  ylim(2.5,6) + 
  theme_bw() +  
  scale_fill_manual(values=c("white", "white","white")) +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1)) +
  geom_dotplot(binaxis = "y", stackdir = "center",aes(fill = factor(genotype2)),position="dodge",binwidth = 0.055)

Morc2b_A = Analysis[which(Analysis$gene=="Morc2b"),]
Morc2b_A$genotype2 <- factor(Morc2b_A$genotype, c("W", "H", "M"))
Morc2bplot <- ggplot(Morc2b_A, aes(factor(age), value))
Morc2bplot + geom_boxplot(aes(fill = factor(genotype2))) + 
  ylim(-1,6) + 
  theme_bw() +  
  scale_fill_manual(values=c("blue", "purple","red")) +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1))

Meioc_A = Analysis[which(Analysis$gene=="Meioc"),]
Meioc_A$genotype2 <- factor(Meioc_A$genotype, c("W", "H", "M"))
Meiocplot <- ggplot(Meioc_A, aes(factor(age), value))
Meiocplot + geom_boxplot(aes(fill = factor(genotype2))) + 
 # ylim(0,6) + 
  theme_bw() +  
  scale_fill_manual(values=c("blue", "purple","red")) +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1))

Ccna2_A = Analysis[which(Analysis$gene=="Ccna2"),]
Ccna2_A$genotype2 <- factor(Ccna2_A$genotype, c("W", "H", "M"))
Ccna2plot <- ggplot(Ccna2_A, aes(factor(age), value))
Ccna2plot + geom_boxplot(aes(fill = factor(genotype2))) + 
  # ylim(0,6) + 
  theme_bw() +  
  scale_fill_manual(values=c("blue", "purple","red")) +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1))

Ythdc2_A = Analysis[which(Analysis$gene=="Ythdc2"),]
Ythdc2_A$genotype2 <- factor(Ythdc2_A$genotype, c("W", "H", "M"))
Ythdc2plot <- ggplot(Ythdc2_A, aes(factor(age), value))
Ythdc2plot + geom_boxplot(aes(fill = factor(genotype2))) + 
  # ylim(0,6) + 
  theme_bw() +  
  scale_fill_manual(values=c("blue", "purple","red")) +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1))




Analysis_AXY = Analysis[which(Analysis$chrom_char!="Unknown"),]
Analysis_WM = Analysis_AXY[which(Analysis_AXY$genotype!="H"),]
Analysis_WM$genochrom = paste(Analysis_WM$genotype,Analysis_WM$chrom_char,sep = " ")


Coef_Plot <- ggplot(Analysis_WM, aes(factor(age), coefficient))
Coef_Plot + geom_boxplot(aes(fill = factor(chrom_char))) + 
  ylim(-1,1) + 
  geom_hline(yintercept = 0) + 
  scale_fill_manual(values=c("darkorchid3", "gold1")) + 
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="grey"),
        panel.border = element_rect(colour = "light grey",fill=NA,size=1))



Analysis_W = Analysis_AXY[which(Analysis_AXY$genotype=="W"),]
Analysis_M = Analysis_AXY[which(Analysis_AXY$genotype=="M"),]


chromplotW <- ggplot(Analysis_W, aes(factor(age), diff_8wt))
chromplotM <- ggplot(Analysis_M, aes(factor(age), diff_8wt))
chromplotW + geom_boxplot(aes(fill = factor(chrom_char))) + ylim(-2,2) + theme_bw() + geom_hline(yintercept = 0) + scale_fill_manual(values=c("blue", "light blue"))
chromplotM + geom_boxplot(aes(fill = factor(chrom_char))) + ylim(-2,2) + theme_bw() + geom_hline(yintercept = 0) + scale_fill_manual(values=c("red", "pink"))
chromplotM + geom_boxplot(aes(fill = factor(chrom_char))) + theme_bw() + geom_hline(yintercept = 0)

Analysis_WM$genochrom2 <- factor(Analysis_WM$genochrom, c("W Autosomes", "M Autosomes", "W XY", "M XY"))
chromplotWM <- ggplot(Analysis_WM, aes(factor(age), diff_8wt))
chromplotWM + geom_boxplot(aes(fill = factor(genochrom2))) + 
  ylim(-2,2) + 
  theme_bw() + 
  geom_hline(yintercept = 0) + 
  scale_fill_manual(values=c("blue", "red","light blue", "pink")) + 
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="grey"),
        panel.border = element_rect(colour = "light grey",fill=NA,size=1))


Analysis_WM_S = Analysis_WM[which(Analysis_WM$chrom_char=="XY"),]
Analysis_WM_A = Analysis_WM[which(Analysis_WM$chrom_char=="Autosomes"),]

chromplotWM_S <- ggplot(Analysis_WM_S, aes(factor(age), diff_8wt))
chromplotWM_S + geom_boxplot(aes(fill = factor(genochrom2))) + 
  ylim(-2,2) + 
  theme_bw() + 
  geom_hline(yintercept = 0,size=1.25) + 
  scale_fill_manual(values=c("blue", "red")) + 
  theme(legend.position="none") #+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="grey"),
        panel.border = element_rect(colour = "light grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank())

chromplotWM_A <- ggplot(Analysis_WM_A, aes(factor(age), diff_8wt))
chromplotWM_A + geom_boxplot(aes(fill = factor(genochrom2))) + 
  ylim(-2,2) + 
  theme_bw() + 
  geom_hline(yintercept = 0,size=1.25) + 
  scale_fill_manual(values=c("blue", "red")) + 
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="grey"),
        panel.border = element_rect(colour = "light grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank())

chromplotWM_S + geom_boxplot(aes(fill = factor(genochrom2))) + 
  ylim(-2,2) + 
  theme_bw() + 
  geom_hline(yintercept = 0) + 
  scale_fill_manual(values=c("#252525", "#d9d9d9")) + 
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="grey"),
        panel.border = element_rect(colour = "light grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank())

P9K_Autosomes = P9K_Analysis[which(P9K_Analysis$I_CNUM %in% numbers),]
P9K_Sex = P9K_Analysis[which(P9K_Analysis$I_CNUM %in% sex),]

Autosomes = data.frame()

library("ggplot2")
data("diamonds")
View(diamonds)
data("mtcars")
p <- ggplot(mtcars, aes(factor(cyl), mpg))

p + geom_boxplot()




Prdm9_A = Analysis[which(Analysis$gene=="Prdm9"),]
Prdm9_A$genotype2 <- factor(Prdm9_A$genotype, c("W", "H", "M"))
Prdm9_Aplot <- ggplot(Prdm9_A, aes(factor(age), value)) + geom_dotplot(binaxis = "y",
                         stackdir = "center",
                         aes(fill = factor(genotype2)),
                         position="dodge",
                         binwidth = 0.08) +
  scale_fill_manual(values=c("blue", "orchid","red"))+
  ylim(2.5,6) + 
  theme_bw()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

Piwil1_A = Analysis[which(Analysis$gene=="Piwil1"),]
Piwil1_A$genotype2 <- factor(Piwil1_A$genotype, c("W", "H", "M"))
Piwil1_Aplot <- ggplot(Piwil1_A, aes(factor(age), value)) + geom_dotplot(binaxis = "y",
                         stackdir = "center",
                         aes(fill = factor(genotype2)),
                         position="dodge",
                         binwidth = 0.2) +
  scale_fill_manual(values=c("blue", "orchid","red"))+
  ylim(0,9) + 
  theme_bw()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

Morc2b_A = Analysis[which(Analysis$gene=="Morc2b"),]
Morc2b_A$genotype2 <- factor(Morc2b_A$genotype, c("W", "H", "M"))
Morc2b_Aplot <- ggplot(Morc2b_A, aes(factor(age), value)) + 
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               aes(fill = factor(genotype2)),
               position = position = "dodge",
               binwidth = 0.2) +
  scale_fill_manual(values=c("blue", "orchid","red"))+
  ylim(-1,7) + 
  theme_bw()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
Morc2b_Aplot <- ggplot(Morc2b_A, aes(factor(age), value)) + 
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               aes(fill = factor(genotype2)),
               position="dodge",
               binwidth = 0.2) +
  scale_fill_manual(values=c("blue", "orchid","red"))+
  ylim(-1,7) + 
  theme_bw()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


Crem_A = Analysis[which(Analysis$gene=="Crem"),]
Crem_A$genotype2 <- factor(Crem_A$genotype, c("W", "H", "M"))
Crem_Aplot <- ggplot(Crem_A, aes(factor(age), value)) + 
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               aes(fill = factor(genotype2)),
               position = position_dodge(0.7), #position = "dodge",
               binwidth = 0.06) +
  scale_fill_manual(values=c("blue", "orchid","red"))+
  ylim(3.5,6) + 
  theme_bw()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


Stra8_A = Analysis[which(Analysis$gene=="Stra8"),]
Stra8_A$genotype2 <- factor(Stra8_A$genotype, c("W", "H", "M"))
Stra8_Aplot <- ggplot(Stra8_A, aes(factor(age), value)) + 
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               aes(fill = factor(genotype2)),
               position = position_dodge(0.7), #position = "dodge",
               binwidth = 0.06) +
  scale_fill_manual(values=c("blue", "orchid","red"))+
  ylim(5,8) + 
  theme_bw()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



Spo11_A = Analysis[which(Analysis$gene=="Spo11"),]
Spo11_A$genotype2 <- factor(Spo11_A$genotype, c("W", "H", "M"))
Spo11_Aplot <- ggplot(Spo11_A, aes(factor(age), value)) + 
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               aes(fill = factor(genotype2)),
               position = position_dodge(0.7), #position = "dodge",
               binwidth = 0.15) +
  scale_fill_manual(values=c("blue", "orchid","red"))+
  ylim(0,6) + 
  theme_bw()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

Rec8_A = Analysis[which(Analysis$gene=="Rec8"),]
Rec8_A$genotype2 <- factor(Rec8_A$genotype, c("W", "H", "M"))
Rec8_Aplot <- ggplot(Rec8_A, aes(factor(age), value)) + 
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               aes(fill = factor(genotype2)),
               position = position_dodge(0.7), #position = "dodge",
               binwidth = 0.1) +
  scale_fill_manual(values=c("blue", "orchid","red"))+
  ylim(3,7) + 
  theme_bw()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

Hspa2_A = Analysis[which(Analysis$gene=="Hspa2"),]
Hspa2_A$genotype2 <- factor(Hspa2_A$genotype, c("W", "H", "M"))
Hspa2_Aplot <- ggplot(Hspa2_A, aes(factor(age), value)) + 
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               aes(fill = factor(genotype2)),
               position = position_dodge(0.7), #position = "dodge",
               binwidth = 0.09) +
  scale_fill_manual(values=c("blue", "orchid","red"))+
  ylim(4.5,9) + 
  theme_bw()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_line(colour="black"),
        panel.border = element_rect(colour = "dark grey",fill=NA,size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

Spo11_Aplot #4x3
Prdm9_Aplot

#Run enrichment test
length(which(P9K_Analysis$S_P9K_PVFA_16D.W_16D.M<0.01))
nrow(P9K_Analysis)

P9K_Analysis_XY = P9K_Analysis[which(P9K_Analysis$I_CNUM=="X" | P9K_Analysis$I_CNUM=="Y"),]

length(which(P9K_Analysis_XY$S_P9K_PVFA_16D.W_16D.M<0.01))
nrow(P9K_Analysis_XY)

XY_test = matrix(NA,ncol=2,nrow=2)

XY_test[1,2] = length(which(P9K_Analysis_XY$S_P9K_PVFA_16D.W_16D.M<0.01 & P9K_Analysis_XY$S_P9K_COEF_16D.W_16D.M > 0))
XY_test[2,2] = nrow(P9K_Analysis_XY) - XY_test[1,2]

XY_test[1,1] = length(which(P9K_Analysis$S_P9K_PVFA_16D.W_16D.M<0.01 & P9K_Analysis$S_P9K_COEF_16D.W_16D.M > 0))
XY_test[2,1] = nrow(P9K_Analysis) - XY_test[1,1]

fisher.test(XY_test)

test = matrix(NA,nrow=2,ncol=2)
test[1,2] = 900
test[2,2] = 400
test[1,1] = 300
test[2,1] = 900
fisher.test(test,alternative = "less")

