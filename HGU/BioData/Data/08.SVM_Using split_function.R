#SVM Using split function
# 1/2 train, 1/2 test
#################################################### VarianceTest ####################################################
GetVar<-function(genes){
  VAR<-apply(genes,2,sd)
  return(VAR)
}

TopVar <- function(genes, feature){
  VAR <- GetVar(genes)
  gene_ch <- rbind(genes, VAR)
  gene_ch_VAR <- gene_ch[,rev(order(gene_ch[nrow(gene_ch),]))]
  gene_sub <-gene_ch_VAR[-nrow(gene_ch_VAR), 1:feature]
  return(gene_sub)
}

#################################################### DiffTest ####################################################
GetDiff <- function(genes, result){
  negative <- apply(genes[result==0,],2,mean)
  positive <- apply(genes[result==1,],2,mean)
  Diff <- abs(positive - negative)
  return(Diff)
}

TopDiff <- function(genes, result, feature){
  Diff <- GetDiff(genes, result)
  gene_ch <- rbind(genes, Diff)
  gene_ch_Diff <- gene_ch[,rev(order(gene_ch[nrow(gene_ch),]))]
  gene_sub <-gene_ch_Diff[-nrow(gene_ch_Diff), 1:feature]
  return(gene_sub)
}

#################################################### Coefficient of Variance ##########################################
GetCV <- function(genes){
  CV<-apply(genes, 2, sd)/apply(genes, 2, mean)
  return(CV)
}

TopCV <- function(genes, feature){
  CV <- GetCV(genes)
  gene_ch <- rbind(genes, CV)
  gene_ch_CV <- gene_ch[,rev(order(gene_ch[nrow(gene_ch),]))]
  gene_sub <-gene_ch_CV[-nrow(gene_ch_CV), 1:feature]
  return(gene_sub)
}

#################################################### VarianceTest ####################################################
GetMean<-function(genes){
  MEAN<-apply(genes,2,mean)
  return(MEAN)
}

TopMean <- function(genes, feature){
  MEAN <- GetMean(genes)
  gene_ch <- rbind(genes, MEAN)
  gene_ch_MEAN <- gene_ch[,rev(order(gene_ch[nrow(gene_ch),]))]
  gene_sub <-gene_ch_MEAN[-nrow(gene_ch_MEAN), 1:feature]
  return(gene_sub)
}
#####################

##data
CV <- read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_CV_4000.csv",header = T, sep = ",")
Mean <- read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_Mean_4000.csv",header = T, sep = ",")
Var <- read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_VAR_4000.csv",header = T, sep = ",")
Annotated_308<-read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_foundation_308.csv",header = T,sep = ",")
Annotated_2267<-read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_foundation_2267.csv",header = T, sep = ',')
