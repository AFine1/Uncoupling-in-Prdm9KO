#   The purpose of this file is to use PMCA to assign transcripts to specific substages

#    
#    
#    


##################################################################################################################
############ Load packages #######################################################################################
##################################################################################################################



##################################################################################################################
############ Aquire Data #########################################################################################
###### X, Y, and wtdata ##########################################################################################

#set working directory
setwd("/Users/s-fine/Desktop/Carter/Projects/P9/6 PMCA/")

#load scripts
setwd("Code")
source("get.mca.R")
source("get.inter.R")
source("get.scores.R")
source("iterative.proc.R")
source("permutation.proc.R")
source("match.patterns.R")
setwd("..")

#Load data
setwd("Data")
options(stringsAsFactors=F)
load("GeneExpression_NonAdjusted_Y.rdt")
load("GeneExpression_Adjusted_Y1.rdt")
load("SampleInfo_X.rdt")
load("GeneExpression_NonAdjusted_YH.rdt")
load("SampleInfo_XH.rdt")
wt <- read.delim("GSE72833_GeneExpressionTPM.txt")
setwd("..")

##################################################################################################################
############ Run MCA #############################################################################################
###### X, Y, and Y1 ->  mca.real & mca.real_1 ####################################################################

#Run MCA
mca.real <- get.mca(X,Y) # get Zx, Zy, sigma, etc.
mca.real_1 <- get.mca(X,Y1) # get Zx, Zy, sigma, etc.
mca.real_H <- get.mca(XH,YH) # get Zx, Zy, sigma, etc.

#Plot result
plot(mca.real$sigma) # plot of singular values of the covariance matrix
plot(mca.real_1$sigma) # plot of singular values of the covariance matrix
plot(mca.real_H$sigma) # plot of singular values of the covariance matrix


##################################################################################################################
############ Run permutation #####################################################################################
###### mca.real & mca.real_1 -> w & w_1 ##########################################################################

#Set parameters
by=3; method="each"; B=1000; alpha=.05; plot=TRUE;
# by: by which component do you want the FPR <= alpha
# method: "overall" if you want an overall FPR, "each" if you want a FPR for each rowterm of X
# B: number of permutations
# alpha: want the FPR <= alpha
# plot: TRUE will plot one of the FPR distributions across all B permutations.

#Get scores from data for unadjusted
# get scores for the real data. Lower scores --> closer pattern match
scores.real <- get.scores(mca.real$Zx,mca.real$Zy)
scores.neg <- get.scores(-mca.real$Zx,mca.real$Zy) #for anti-associated lists

#Run permuation for unadjusted data
set.seed(B) # for reproducibility
scores.rand <- permutation.proc(X,Y,method=method,B=B) # get scores for all B permutations
w <- apply(mca.real$Zx,2,sd) # get starting window vector
set.seed(B) # for reproducibilty

#Get scores from data for adjusted
scores.real_1 <- get.scores(mca.real_1$Zx,mca.real_1$Zy)
scores.neg_1 <- get.scores(-mca.real_1$Zx,mca.real_1$Zy) #for anti-associated lists

#Run permuation for adjusted data
set.seed(B) # for reproducibility
scores.rand_1 <- permutation.proc(X,Y1,method=method,B=B) # get scores for all B permutations
w_1 <- apply(mca.real_1$Zx,2,sd) # get starting window vector
set.seed(B) # for reproducibilty

#Get scores from data for adjusted
scores.real_H <- get.scores(mca.real_H$Zx,mca.real_H$Zy)
scores.neg_H <- get.scores(-mca.real_H$Zx,mca.real_H$Zy) #for anti-associated lists

#Run permuation for adjusted data
set.seed(B) # for reproducibility
scores.rand_H <- permutation.proc(XH,YH,method=method,B=B) # get scores for all B permutations
w_H <- apply(mca.real_H$Zx,2,sd) # get starting window vector
set.seed(B) # for reproducibilty

##################################################################################################################
############ Get data out ########################################################################################
###### w & w_1 -> it.result & it.result_1 ########################################################################

# tau: controls the width of the window (w/tau) 
#     if you are not getting good results, try making tau smaller or larger (larger = more strict)
it.result <- iterative.proc(scores.rand,alpha=alpha,w=w,method=method,by=by,plot=plot,tau=1) 
it.result_1 <- iterative.proc(scores.rand_1,alpha=alpha,w=w_1,method=method,by=by,plot=plot,tau=1) 
it.result_H <- iterative.proc(scores.rand_H,alpha=alpha,w=w_H,method=method,by=by,plot=plot,tau=1) 

it.result$tau
it.result$FPR
## The way to read the FPR is that for row 1, the FPR for column 3 is .047 and for column 4 is 0.012.
## You will get a finer list as you use more columns so that the genes that map if you use columns 1-3 have a FPR < 0.05


##################################################################################################################
############ Get substage specific genes #########################################################################
###### it.result & it.result_1 -> g & g_1 ########################################################################

# get wopt (optimal window vector)
# g = list of rownames(Y) that match patterns of rownames(X)
#   g[[i]][[j]]: is all the rownames(Y) that match the pattern of rownames(X)[i] when the component = j
#   example: If you want to know what genes map to "LP_Dip" (row 6) 
#         when you use 4 components, g[[6]][[4]]
g <- match.patterns(scores.real,w=it.result$wopt)  

g_1 <- match.patterns(scores.real_1,w=it.result_1$wopt)  

g_H <- match.patterns(scores.real_H,w=it.result_H$wopt)  

##################################################################################################################
############ Get gene lists ######################################################################################
###### g & g_1 -> gene list ######################################################################################

#Make gene lists
Sp_KOG = g[[1]][[4]]
Pl_KOG = g[[2]][[4]]
EL_KOG = g[[3]][[4]]
LLZ_KOG = g[[4]][[4]]
PLike_KOG = g[[5]][[4]]
EP_KOG = g[[6]][[4]]

Sp_KOG_1 = g_1[[1]][[4]]
Pl_KOG_1 = g_1[[2]][[4]]
EL_KOG_1 = g_1[[3]][[4]]
LLZ_KOG_1 = g_1[[4]][[4]]
PLike_KOG_1 = g_1[[5]][[4]]
EP_KOG_1 = g_1[[6]][[4]]


Sp_KOG_H = g_H[[1]][[4]]
Pl_KOG_H = g_H[[2]][[4]]
EL_KOG_H = g_H[[3]][[4]]
LLZ_KOG_H = g_H[[4]][[4]]
PLike_KOG_H = g_H[[5]][[4]]
EP_KOG_H = g_H[[6]][[4]]


Sp_KOG_5 = g_1[[1]][[5]]
Pl_KOG_5 = g_1[[2]][[5]]
EL_KOG_5 = g_1[[3]][[5]]
LLZ_KOG_5 = g_1[[4]][[5]]
PLike_KOG_5 = g_1[[5]][[5]]
EP_KOG_5 = g_1[[6]][[5]]


##################################################################################################################
############ Make wt g ###########################################################################################
###### wt -> g2 ##################################################################################################

#Make g2
g2 <- list()
sub <- c("Sp'gonia","Prelep","Early Lep","Late Lep/Zyg","Late Lep/Zyg/Early Pach",
         "Early Pach","Early Pach/Late Pach/Dip","Late Pach/Dip" )
for (i in 1:length(sub)) {
  g2[[i]] <- wt$gene_id[which(wt$substage==sub[i])]
}


##################################################################################################################
############ Make overlap lists ##################################################################################
##################################################################################################################

Sp_Sp_overlap = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g2[[1]][which(g2[[1]] %in% Sp_KOG)]),]
Pl_Pl_overlap = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g2[[2]][which(g2[[2]] %in% Pl_KOG)]),]
EL_EL_overlap = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g2[[3]][which(g2[[3]] %in% EL_KOG)]),]
LLZ_LLZ_overlap = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g2[[4]][which(g2[[4]] %in% LLZ_KOG)]),]
EP_EP_overlap = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g2[[6]][which(g2[[6]] %in% EP_KOG)]),]
PLike_EP_overlap = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g2[[6]][which(g2[[6]] %in% PLike_KOG)]),]
LLZ_EP_overlap = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g2[[6]][which(g2[[6]] %in% LLZ_KOG)]),]
LLZ_LPD_overlap = P9K_Analysis[which(P9K_Analysis$I_GEID %in% g2[[8]][which(g2[[8]] %in% LLZ_KOG)]),]


##
#Make gene lists
##

Sp_KOGenes = P9K_Analysis[which(P9K_Analysis$I_GEID %in% Sp_KOG),]
Pl_KOGenes = P9K_Analysis[which(P9K_Analysis$I_GEID %in% Pl_KOG),]
EL_KOGenes = P9K_Analysis[which(P9K_Analysis$I_GEID %in% EL_KOG),]
LLZ_KOGenes = P9K_Analysis[which(P9K_Analysis$I_GEID %in% LLZ_KOG),]
PLike_KOGenes = P9K_Analysis[which(P9K_Analysis$I_GEID %in% PLike_KOG),]
EP_KOGenes = P9K_Analysis[which(P9K_Analysis$I_GEID %in% EP_KOG),]


##################################################################################################################
############ Save lists ##########################################################################################
##################################################################################################################

setwd("Results/")
save(g,file="UnAdj_MutSub.rdt")
save(g_1,file="Adj_MutSub.rdt")
save(g2,file="WtSub.rdt")
save(g_H,file="HetSub.rdt")
load("UnAdj_MutSub.rdt")
load("WtSub.rdt")
write.csv(Sp_KOGenes, file="Sp_KOGenes.csv")
write.csv(Pl_KOGenes, file="Pl_KOGenes.csv")
write.csv(EL_KOGenes, file="EL_KOGenes.csv")
write.csv(LLZ_KOGenes, file="LLZ_KOGenes.csv")
write.csv(PLike_KOGenes, file="PLike_KOGenes.csv")
write.csv(EP_KOGenes, file="EP_KOGenes.csv")

save(Sp_KOGenes, file="Sp_KOGenes.rdt")
save(Pl_KOGenes, file="Pl_KOGenes.rdt")
save(EL_KOGenes, file="EL_KOGenes.rdt")
save(LLZ_KOGenes, file="LLZ_KOGenes.rdt")
save(PLike_KOGenes, file="PLike_KOGenes.rdt")
save(EP_KOGenes, file="EP_KOGenes.rdt")

setwd("..")
