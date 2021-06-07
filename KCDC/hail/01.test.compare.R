setwd("~/Desktop/KCDC/hail/")

d1 <- read.table("gaws.output.tsv",header = T)
head(d1)
d1 <- d1[,c(1,3,4,5)]
summary(d1)

d2 <- read.table("JG.KR.ESRD.association.sub_Total.b.firth.chr12.132000001_133841815.epacts",header = T)
head(d2)
d2$locus <- paste0(d2$CHROM,":",d2$BEGIN)
d2 <- d2[,c("locus","BETA","CHISQ","PVALUE")]
colnames(d1) <- colnames(d2)
head(d1)
head(d2)

library(stringr)

d1$POS <- str_split_fixed(d1$locus,":",2)[2]
head(d1)
d2$POS <- str_split_fixed(d2$locus,":",2)[2]
head(d2)
