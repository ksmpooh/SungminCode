setwd("c:/Users/user/Desktop/KCDC/HLAimputation/after_filter/")

bim <- read.table("JG.QCed.MHC_rmAmbg.bim")
head(bim)
result <- read.table("hglft_genome_ab75_77d980.bed")

head(bim)
head(result)

library(stringr)



result$V2 <- str_split_fixed(result$V1,"-",2)[,2]
df <- cbind(bim,result)
head(df)

write.table(df[,c(1,2,3,8,5,6)],"JG.QCed.MHC_rmAmbg.bim",col.names = F,row.names = F,quote = F,sep = '\t')



####rs ID

ori <-read.table()