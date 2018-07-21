

CV <- read.csv("D:/biodatalab/2018-1/TCGA_with_GEO/union/CV_3000.csv")
Mean <-  read.csv("D:/biodatalab/2018-1/TCGA_with_GEO/union/Mean_3000.csv")
Var <-  read.csv("D:/biodatalab/2018-1/TCGA_with_GEO/union/VAR_3000.csv")
Annotated_2267<- read.csv("D:/biodatalab/2018-1/TCGA_with_GEO/union/foundation_2267.csv")
Annotated_308 <- read.csv("D:/biodatalab/2018-1/TCGA_with_GEO/union/foundation_308.csv")

CV_union <- subset(CV,select = -c(index,result,cancer_code,patient))
Mean_union <- subset(Mean,select = -c(index,result,cancer_code,patient))
Var_union <- subset(Var,select = -c(index,result,cancer_code,patient))
Annotated_308union <- subset(Annotated_308,select = -c(index,result,cancer_code,patient))
Annotated_308colnames<-colnames(Annotated_308union)
Annotated_2267union <- subset(Annotated_2267,select = -c(index,result,cancer_code,patient))
Annotated_2267colnames<-colnames(Annotated_2267union)

#lists <- c(500,1000,1500,2000,2500,3000,3500,4000)
feature <-3000
result <-data.frame()
Mean_colnames <- colnames(Mean_union)
var_colnames <- colnames(Var_union)
CV_colnames <- colnames(CV_union)
feature_sum <- sum(feature * 3, 308, 2267)
df <- union(Mean_colnames,var_colnames)
df <- union(df,CV_colnames)
df <- union(df,Annotated_308colnames)
df <- union(df,Annotated_2267colnames)
df_n<-df
Union <- length(df)
df <- data.frame(feature,feature_sum,Union,stringsAsFactors = F)
result <- rbind(result,df)
print(feature)
print(Union)

df_n

write.csv(df_n,"D:/biodatalab/2018-1/TCGA_with_GEO/union/union.csv",row.names = F)
