
#setwd("C:/Users/user/Desktop/KCDC/HLAimputation/IMPUTE4/HAN.ref")
setwd("C:/Users/user/Desktop/KCDC/HLAimputation/IMPUTE4/PAN.ref/")
##### reference panel ref/alt position

library(BiocManager)
library(BiocVersion)


markers <- read.table("Merge_panel.markers")
colnames(markers)<-c("id","pos","a1","a2")
head(markers)
tail(markers)


#BiocManager::install("BSgenome")
#BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')

library(BSgenome.Hsapiens.UCSC.hg19)
refallele <- getSeq(Hsapiens,"chr6",24898856,34886436)
#refallele <- getSeq(Hsapiens,"chr6",28477833,28477835)
refallele <- as.character(refallele)
#refallele[1](1)
substr(refallele,start = 1,stop = 1)


library(stringr)


nsnp <- nrow(markers)
#out <-data.frame(pos = NULL, ref = NULL)
out <- matrix(nrow = nsnp,ncol = 2)
out

k = 0

for (i in markers$pos){
  k = k + 1
  j <- strtoi(i)
  r <- j - 24898855
  #print(i)
  #print(j)
  #print(j - 28477832)
  out[k,1] <- j
  out[k,2] <- substr(refallele,start = r,stop = r)
#  print(out[0:5])
#  print(substr(refallele,start = r,stop = r))
}
out<-as.data.frame(out)
out
colnames(out) <- c("pos","ref")
head(out)
head(markers)
df <-cbind(markers,out)
head(df)
df$chr <- 6
#df <-merge(bim,out,by.x = "hg19",by.y = "pos",all.y = TRUE)
#colnames(df)[1:6]<-c("chr","ID","0","hg18","a1","a2")
head(df)
write.table(df[,c(7,1,2,3,4,6)],"ref.allele.txt",col.names = T,row.names = F,quote = F)


head(markers)
#B <- df[,c(1,2,grep("HLA_B_*",colnames(df)))]

a <- markers[grep("AA",markers$id),]
a
########################make hg17 for Pan

ref <- read.table("ref.allele.txt",header = T)
head(ref)
pan <- read.table("Merge_panel.markers.hg17")
head(pan)
#df <- merge(ref,pan,by.x = "id",by.y="V1")

df <- cbind(ref,pan)

head(df)
df <- df[,c(1,2,8,4,5,3,6)]

colnames(df)<-c("chr","ID","hg17","a1","a2","hg19","ref")

colnames(ref)

write.table(df[,c("chr","ID","hg19","a1","a2","hg17","ref")],"ref.allele.txt",col.names = T,row.names = F,quote = F)
