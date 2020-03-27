#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/after_filter/")
setwd("E:/hlaimputation20200324/")
bim <- read.table("test/JG.QCed.HLA_rmAmbiguous.bim.hg19")
head(bim)
result <- read.table("test/hglft_genome_23bcd_b12390.bed")

head(bim)
head(result)

library(stringr)



result$V2 <- str_split_fixed(result$V1,"-",2)[,2]
df <- cbind(bim,result)
head(df)

write.table(df[,c(1,2,3,8,5,6)],"JG.QCed.HLA_rmAmbiguous.bim",col.names = F,row.names = F,quote = F,sep = '\t')



####rs ID

ori <-read.table()