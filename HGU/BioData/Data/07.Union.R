##Union count
# unique(c(A,B))

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


CV <- read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_CV_2500.csv",header = T, sep = ",")
Mean <- read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_Mean_2500.csv",header = T, sep = ",")
Var <- read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_VAR_2500.csv",header = T, sep = ",")
Annotated<-read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_foundation_308.csv",header = T,sep = ",")

CV <- subset(CV,select = -c(index,result,cancer_code,patient))
Mean <- subset(Mean,select = -c(index,result,cancer_code,patient))
Var <- subset(Var,select = -c(index,result,cancer_code,patient))
Annotated <- subset(Annotated,select = -c(index,result,cancer_code,patient))
Annotated_colnames<-colnames(Annotated)
lists <- c(50,100,200,500,1000,2500)
result <-data.frame()
for (feature in lists){
  var_ <- TopVar(Var,feature)
  Mean_ <- TopMean(Mean,feature)
  CV_ <- TopCV(CV,feature)
  Mean_colnames <- colnames(Mean_)
  var_colnames <- colnames(var_)
  CV_colnames <- colnames(CV_)
  df <- union(Mean_colnames,var_colnames)
  df <- union(df,CV_colnames)
  df <- union(df,Annotated_colnames)
  Union <- length(df)
  df <- data.frame(feature,Union,stringsAsFactors = F)
  result <- rbind(result,df)
  print(feature)
  print(Union)
}

a<-read.csv("C:/Users/sungmin/Documents/카카오톡 받은 파일/foundation_308.csv")
a
result
length(union(a$x,Annotated_colnames))
length(intersect(a$x,Annotated_colnames))
