setwd("~/Desktop/KCDC/HLAimputation/1000genome/")
library(BiocManager)
library(BiocVersion)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install(version = "3.17")

#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("BSgenome")
#BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')
#BiocManager::install('SNPlocs.Hsapiens.dbSNP144.GRCh37')
library(BSgenome.Hsapiens.UCSC.hg19)
#library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
refallele <- getSeq(Hsapiens,"chr6",24913126,34855591)
#refallele <- getSeq(Hsapiens,"chr6",28477833,28477835)
refallele <- as.character(refallele)
refallele[1](1)
substr(refallele,start = 1,stop = 5)

library(stringr)
bim <-read.table("1kgp.phase3.chr6.MHC.bim")
head(bim)
tail(bim)

nsnp <- nrow(bim)
out <-data.frame(pos = NULL, ref = NULL)
out <- matrix(nrow = nsnp,ncol = 2)
out

k = 0
for (i in bim$V4){
  k = k + 1
  j <- strtoi(i)
  r <- j - 24913126 + 1
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
#df <-cbind(bim,out)

bim %>% merge(out,by.x="V4",by.y="pos") %>% 
  mutate("alt" = ifelse(ref == V5,V6,V5)) %>% 
  mutate("id" = paste0(V1,"_",V4,"_",ref,"_",alt)) %>% 
  write.table("rsID.with.ref.alt.newID.txt",col.names = T,row.names = F,quote = F,sep = "\t")




