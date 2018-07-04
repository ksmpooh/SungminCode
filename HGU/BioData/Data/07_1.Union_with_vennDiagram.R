##Union count
# unique(c(A,B))
## ensemble model data.
# 1. CV
# 2. Mean
# 3. Var
# 4. annotation 308
# 5. annotation 2267

#library(VennDiagram)
library(limma)
library(gplots)

CV <- read.csv("D:/biodatalab/2018-1/TCGA_with_GEO/union/names_GEO_input_ensemble_CV_3000.csv",header = T, sep = ",")
Mean <- read.csv("D:/biodatalab/2018-1/TCGA_with_GEO/union/names_GEO_input_ensemble_Mean_3000.csv",header = T, sep = ",")
Var <- read.csv("D:/biodatalab/2018-1/TCGA_with_GEO/union/names_GEO_input_ensemble_VAR_3000.csv",header = T, sep = ",")
Annotated_308<-read.csv("D:/biodatalab/2018-1/TCGA_with_GEO/union/names_GEO_input_ensemble_foundation_308.csv",header = T,sep = ",")
Annotated_2267<-read.csv("D:/biodatalab/2018-1/TCGA_with_GEO/union/names_GEO_input_ensemble_foundation_2267.csv",header = T, sep = ',')

#CV_union <- subset(CV,select = -c(index,result,cancer_code,patient))
#Mean_union <- subset(Mean,select = -c(index,result,cancer_code,patient))
#Var_union <- subset(Var,select = -c(index,result,cancer_code,patient))
#Annotated_308union <- subset(Annotated_308,select = -c(index,result,cancer_code,patient))
#Annotated_308colnames<-colnames(Annotated_308union)
#Annotated_2267union <- subset(Annotated_2267,select = -c(index,result,cancer_code,patient))
#Annotated_2267colnames<-colnames(Annotated_2267union)

#lists <- c(500,1000,1500,2000,2500,3000,3500,4000)
#colnames(CV) <- "CV"
#colnames(Mean) <- "Mean"
#colnames(Var) <- "Var"
#colnames(Annotated_308) <- "Foundation_308"
#colnames(Annotated_2267) <- "Foundation_2267"

result <-union(CV$x,Mean$x)
result <-union(result,Var$x)
result <-union(result,Annotated_2267$x)
result <-union(result,Annotated_308$x)

all_gene <- Reduce(rbind,result)
all_gene <- as.data.frame(all_gene)

all_gene <- t(all_gene)
all_gene <-as.data.frame(all_gene)
colnames(all_gene)<-result

all_gene$index <- "all"


CV_ <-t(CV)
CV_ <-as.data.frame(CV_)
colnames(CV_)<- CV$x
CV_$index <- "CV"

Mean_ <-t(Mean)
Mean_ <-as.data.frame(Mean_)
colnames(Mean_) <- Mean$x
Mean_$index <-"Mean"

Annotated_2267_ <-t(Annotated_2267)
Annotated_2267_ <-as.data.frame(Annotated_2267_)
colnames(Annotated_2267_) <- Annotated_2267$x
Annotated_2267_$index <-"Foundation_2267"

Annotated_308_ <-t(Annotated_308)
Annotated_308_ <-as.data.frame(Annotated_308_)
colnames(Annotated_308_) <- Annotated_308$x
Annotated_308_$index <-"Foundation_308"

Var_ <-t(Var)
Var_ <-as.data.frame(Var_)
colnames(Var_) <- Var$x
Var_$index <-"Var"

df <-data.frame()

df <- bind_rows(all_gene,CV_)
df <- bind_rows(df,Mean_)
df <- bind_rows(df,Var_)
df <- bind_rows(df,Annotated_308_)
df <- bind_rows(df,Annotated_2267_)

rownames(df)<-df$index
df<-subset(df,select = -index)

df_ <- t(df)
df_ <-as.data.frame(df_)
#colnames(df_)
#rownames(df_)
df_ <-subset(df_,select = -all)

id <- (df_!="")
id[is.na(id)]<-FALSE 
id.df<-as.data.frame(id)

vennCounts(id.df)
venn(id.df)

