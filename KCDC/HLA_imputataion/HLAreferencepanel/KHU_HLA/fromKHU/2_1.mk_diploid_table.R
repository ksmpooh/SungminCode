rm(list=ls())
library(dplyr);library(tidyr)
library(optparse)
library(bigreadr)
library(doParallel)

setwd( "~~~/wholesample_imputation_result" )

dosage=fread2("./dosage_imputed_snps_only.txt", header=T,sep="\t")
dim(dosage)
dosage[1:10,1:10]
colnames(dosage)[1] <- "[1]ID"
colnames(dosage) = gsub("\\[\\d+\\]", "", colnames(dosage)) %>% gsub(":DS", "", ., fixed=TRUE)
dosage = select(dosage, !c("CHROM", "REF", "ALT", "POS"))

tableC4 = dosage[grep("copy", dosage$ID),]
rownames(tableC4) = tableC4$ID
tableC4 = tableC4[-1]
tableC4=t(tableC4)
tableC4=as.data.frame(tableC4)
tableC4$ID=rownames(tableC4)
tableC4= relocate(tableC4, ID)
rownames(tableC4)=1:nrow(tableC4)

#
plot(density(tableC4$C4_copy2))
dim(tableC4)
head(tableC4)

# dosage normalize
c4sum <- apply( tableC4[,grep("C4_", colnames(tableC4) )]  , 1, sum)
c4asum <- apply( tableC4[,grep("C4A_", colnames(tableC4) )]  , 1, sum)
c4bsum <- apply( tableC4[,grep("C4B_", colnames(tableC4) )]  , 1, sum)
hervsum <- apply( tableC4[,grep("HERV_", colnames(tableC4) )]  , 1, sum)
length(c4asum);length(c4asum);length(c4bsum);length(hervsum)

tableC4[,grep("C4_", colnames(tableC4) )] <- tableC4[,grep("C4_", colnames(tableC4) )]/c4sum*2
tableC4[,grep("C4A_", colnames(tableC4) )] <- tableC4[,grep("C4A_", colnames(tableC4) )]/c4asum*2
tableC4[,grep("C4B_", colnames(tableC4) )] <- tableC4[,grep("C4B_", colnames(tableC4) )]/c4bsum*2
tableC4[,grep("HERV_", colnames(tableC4) )] <- tableC4[,grep("HERV_", colnames(tableC4) )]/hervsum*2


head(tableC4)

C4_diploid_imputed =  data.frame( ID=tableC4$ID,
            C4=apply( tableC4[ ,grep("C4_", colnames(tableC4) ) ] ,1, function(x) sum( x[1]*2 + x[2]*1 + x[3]*3 + x[4]*4 ) ) ,
            C4A=apply( tableC4[ ,grep("C4A_", colnames(tableC4) ) ] ,1, function(x) sum( x[1]*2 + x[2]*0 + x[3]*1 + x[4]*3 ) ) ,
            C4B=apply( tableC4[ ,grep("C4B_", colnames(tableC4) ) ] ,1, function(x) sum( x[1]*2 + x[2]*0 + x[3]*1 + x[4]*3 ) ) ,
            HERV=apply( tableC4[ ,grep("HERV_", colnames(tableC4) ) ] ,1, function(x) sum( x[1]*2 + x[2]*0 + x[3]*1 + x[4]*3 + x[5]*4 ) )
            )


C4_diploid_imputed[1:10,]
getwd()

write.table(C4_diploid_imputed , "./tableC4_diploid.txt", col.names = T, row.names=F, sep="\t", quote=F)

head(C4_diploid_imputed)
hist(C4_diploid_imputed$C4)
hist(C4_diploid_imputed$C4A)
hist(C4_diploid_imputed$C4B)
hist(C4_diploid_imputed$HERV)


# %