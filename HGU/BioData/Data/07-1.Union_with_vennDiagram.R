##Union count
# unique(c(A,B))
## ensemble model data.
# 1. CV
# 2. Mean
# 3. Var
# 4. annotation 308
# 5. annotation 2267

library(VennDiagram)
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


CV <- read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_CV_4000.csv",header = T, sep = ",")
Mean <- read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_Mean_4000.csv",header = T, sep = ",")
Var <- read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_VAR_4000.csv",header = T, sep = ",")
Annotated_308<-read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_foundation_308.csv",header = T,sep = ",")
Annotated_2267<-read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_foundation_2267.csv",header = T, sep = ',')

CV_union <- subset(CV,select = -c(index,result,cancer_code,patient))
Mean_union <- subset(Mean,select = -c(index,result,cancer_code,patient))
Var_union <- subset(Var,select = -c(index,result,cancer_code,patient))
Annotated_308union <- subset(Annotated_308,select = -c(index,result,cancer_code,patient))
Annotated_308colnames<-colnames(Annotated_308union)
Annotated_2267union <- subset(Annotated_2267,select = -c(index,result,cancer_code,patient))
Annotated_2267colnames<-colnames(Annotated_2267union)

#lists <- c(500,1000,1500,2000,2500,3000,3500,4000)
result <-data.frame()

#for (feature in lists){
feature <- 3000
  var_ <- TopVar(Var_union,feature)
  Mean_ <- TopMean(Mean_union,feature)
  CV_ <- TopCV(CV_union,feature)
  Mean_colnames <- colnames(Mean_)
  var_colnames <- colnames(var_)
  CV_colnames <- colnames(CV_)
  feature_sum <- sum(feature * 3, 308, 2267)
  #df_merge<-merge(Mean_colnames,var_colnames,CV_colnames,Annotated_2267colnames,Annotated_308colnames)
  Mean_colnames <- as.data.frame(Mean_colnames)
  var_colnames <- as.data.frame(var_colnames)
  CV_colnames <- as.data.frame(CV_colnames)
  
  df <- merge(Mean_colnames,var_colnames)
  
  df <- union(Mean_colnames,var_colnames)
  df <- union(df,CV_colnames)
  df <- union(df,Annotated_308colnames)
  df <- union(df,Annotated_2267colnames)
  Union <- length(df)
  df <- data.frame(feature,feature_sum,Union,stringsAsFactors = F)
  result <- rbind(result,df)
  print(feature)
  print(Union)
#}


write.csv(result,"D:/biodatalab/2018-1/Result/Union.csv",row.names = F)
