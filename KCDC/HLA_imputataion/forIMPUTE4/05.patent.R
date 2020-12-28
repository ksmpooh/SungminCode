## 특허 관련 공통 allele 찾기

setwd('c:/Users/user/Desktop/KCDC/HLAimputation/patent/')
kba <- read.table("JG.QCed.HLA.bim")
head(kba)


library(BSgenome.Hsapiens.UCSC.hg19)
refallele <- getSeq(Hsapiens,"chr6",28477833,33448188)
#refallele <- getSeq(Hsapiens,"chr6",28477833,28477835)
refallele <- as.character(refallele)

substr(refallele,start = 1,stop = 1)


library(stringr)

nsnp <- nrow(kba)
out <-data.frame(pos = NULL, ref = NULL)
out <- matrix(nrow = nsnp,ncol = 2)
out

k = 0
for (i in kba$V4){
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
head(out)


table(out$V2 == kba$V5)
table(out$V2 == kba$V6)

colnames(out) <- c("pos","ref")
write.table(out,"KBA.ref.allele.for.HLAregion.txt",row.names = F,col.names = T,quote = F)


###############################
setwd('c:/Users/user/Desktop/KCDC/HLAimputation/patent/')

df <- read.table("KBA.ref.allele.for.HLAregion.txt",header = T)
head(df)
kba <- read.table("JG.QCed.HLA.bim")
head(kba)
colnames(kba) <- c("chr","KBA_ID","V3","pos","A1","A2")

df <- merge(df,kba,by.x = "pos",by.y = "pos")
head(df)
