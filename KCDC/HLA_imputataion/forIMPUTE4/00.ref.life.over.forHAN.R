
setwd("C:/Users/user/Desktop/KCDC/HLAimputation/IMPUTE4/HAN.ref")
##### reference panel ref/alt position

library(BiocManager)
library(BiocVersion)



#BiocManager::install("BSgenome")
#BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')

library(BSgenome.Hsapiens.UCSC.hg19)
refallele <- getSeq(Hsapiens,"chr6",28477833,33448188)
#refallele <- getSeq(Hsapiens,"chr6",28477833,28477835)
refallele <- as.character(refallele)
refallele[1](1)
substr(refallele,start = 1,stop = 1)


library(stringr)
bim <-read.table("HAN.MHC.reference.panel.fixed.bim")
head(bim)

markers <- read.table("HAN.MHC.reference.panel.fixed.markers")
head(markers)

ref <- read.table("hglft_genome_151e2_774cc0.bed")
head(ref)

bim$hg19 <- str_split_fixed(ref$V1,'-',2)[,2]
head(bim)
tail(bim)

nsnp <- nrow(bim)
out <-data.frame(pos = NULL, ref = NULL)
out <- matrix(nrow = nsnp,ncol = 2)
out

k = 0
for (i in bim$hg19){
  k = k + 1
  j <- strtoi(i)
  r <- j - 28477832
  #print(i)
  #print(j)
  #print(j - 28477832)
  out[k,1] <- j
  out[k,2] <- substr(refallele,start = r,stop = r)
}
out<-as.data.frame(out)

colnames(out) <- c("pos","ref")
head(out)
head(bim)
df <-cbind(bim,out)
head(df)
#df <-merge(bim,out,by.x = "hg19",by.y = "pos",all.y = TRUE)
colnames(df)[1:6]<-c("chr","ID","0","hg18","a1","a2")
head(df)
write.table(df[,c(1,2,4,5,6,7,9)],"ref.allele.txt",col.names = T,row.names = F,quote = F)


