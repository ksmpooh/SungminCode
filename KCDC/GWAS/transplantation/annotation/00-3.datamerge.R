library(tidyverse)
library(stringr)
setwd("~/Desktop/KCDC/transplantation/02.annotation/")

df <- read.csv("KCHIP_V1.1/VEP/KCHIPV1.1_sorted_VEP.vcf_summary_dataprocessing.txt",header = T)
df1 <- read.csv("QCed.KCHIP/VEP/imputed.Kchip.marker_VEP.vcf_summary_dataprocessing.txt",header = T)
head(df)
head(df1)



out <- merge(df,df1,by="index",all = T)
head(out)


##########
setwd("/Users/ksmpooh/Desktop/KCDC/transplantation/02.annotation/imputed/VEP/pro")

#df <- read.csv("JG.KR.imputation.mergeGen.processing.chr1.chunkMerge_VEP.vcf_summary_dataprocessing.txt",header = T)
#df <- read.csv("JJG.KR.imputation.mergeGen.processing.chr1.chunkMerge_VEP.vcf_summary_dataprocessing.txt",header = T)
df <- read.csv("JG.KR.imputation.mergeGen.processing.chr1.chunkMerge_checkREFallele_VEP.vcf_summary_dataprocessing.txt",header = T)
head(df)
for (i in 2:22) {
  df1 <- read.csv(paste0("JG.KR.imputation.mergeGen.processing.chr",i,".chunkMerge_checkREFallele_VEP.vcf_summary_dataprocessing.txt"),header = T)
  df <- merge(df,df1,by="index",all=T)
  print(i)
}

out <- merge(out,df,by="index",all = T)
head(out)
out$theme <- str_split_fixed(out$index,"-",2)[,1]
out$type <- str_split_fixed(out$index,"-",2)[,2]
str_split_fixed(out$index,"-",2)[,2]


out <- out[,c(26,27,2:25)]
head(out)
write.csv(out,"/Users/ksmpooh/Desktop/KCDC/transplantation/02.annotation/RESULTs/VEP.annotation.Result_20220317_checkREFALLELE.csv",col.names = T,row.names = F,quote = F)
