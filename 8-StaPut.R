#This file compares our LP/D and EP lists to previously published gene lists for transcripts found in STAPUT pachytene samples.

setwd("/Users/s-fine/Desktop/Carter/Projects/P9/STAPUT/")

#load data
setwd("Data/")
STP_Data_TPM = read.table("GSE72833_STA-PUT_TPM.txt",header=T)
setwd("..")
load("/Users/s-fine/Desktop/Carter/Projects/P9/3 DifferentialExpression/Results/P9K_Analysis.rdt")

#log transform data
STP_Data = log2(STP_Data_TPM+1)
STP_Expression = STP_Data[which(rowMaxs(as.matrix(STP_Data))>=3),]

#Convert gene ID to gene name
m_ensembl=useMart("ensembl",dataset="mmusculus_gene_ensembl")
gene_info=getBM(attributes=c("external_gene_name","ensembl_gene_id"),mart=m_ensembl)
STP_Expression$gene_name = "NA"
for (i in 1:nrow(STP_Expression)){
  geneid=rownames(STP_Expression[i,])
  if (geneid %in% as.character(gene_info$ensembl_gene_id)) {
    STP_Expression[i,"gene_name"]=as.character(gene_info[which(gene_info$ensembl_gene_id==geneid),]$external_gene_name[1])
  }
  
  else{
    STP_Expression[i,"gene_name"]=geneid
  }
}
STP_Expression$gene_id = rownames(STP_Expression)
STP_Expression = STP_Expression[-which(STP_Expression$gene_id=="ENSMUSG00000091562"),]
rownames(STP_Expression) = STP_Expression[,"gene_name"]
STP_Exp = STP_Expression[,c(1,2,3,4,6)]

setwd("Data/")
save(STP_Exp,file="STP_Exp.rdt")
load("STP_Exp.rdt")
setwd("..")

length(which(rownames(P9K_Analysis[which(P9K_Analysis$S_RYT_SpecSub=="LLZ"),]) %in% rownames(STP_Exp)))
nrow(P9K_Analysis[which(P9K_Analysis$S_RYT_SpecSub=="LPD"),])
nrow(P9K_Analysis[which(P9K_Analysis$S_RYT_SpecSub=="LLZ"),])

STP_Exp["Spata17",]
STP_Exp["Yy1",]


load("/Users/s-fine/Desktop/Carter/Projects/P9/Substage_Specificity/Results/P9K_LPD_to_LLZ_2.rdt")

length(which(rownames(P9K_LPD_to_LLZ_2) %in% rownames(STP_Exp)))

length(which(rownames(P9K_EP_lt_to_LLZ) %in% rownames(STP_Exp)))

