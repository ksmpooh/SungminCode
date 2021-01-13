setwd("c:/Users/user/Desktop/KCDC/imputation.tool/")
## box plot

library(ggplot2)

df <- read.table("result/time.total.txt",header = T)
#colnames(df)[5] <- ("time")
head(df)
table(df$type)
type_order <- c("impute4.1-5000","impute4.5001-10000","mergeGen","mergeGen.remove.NA","info_score","gen2vcf")

df$type <- factor(df$type,levels = type_order)
#df <- df[order(match(df$type,type_order)),]
#df$type <- factor(df$type,levels = c("impute4.1-5000","impute4.5001-10000","info_score","mergeGen","mergeGen.remove.NA","info_score","gen2vcf"))

p <- ggplot(df,aes(y=time,x=type))  + 
  geom_boxplot(outlier.colour = "black",outlier.shape = 1
               ,outlier.size = 2, notch = F)
p
?ggplot


