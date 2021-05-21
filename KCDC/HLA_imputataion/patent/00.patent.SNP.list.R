setwd("c:/Users/user/Desktop/KCDC/HLAimputation/patent/")

ref <- read.table("KBA.ref.allele.for.HLAregion.txt",header = T)
bim <- read.table("JG.QCed.HLA.bim")
head(ref)
head(bim)
colnames(bim) <- c("chr","SNP_ID","0","pos","a1","a2")


out <- merge(bim,ref,by = "pos")
head(out)
out$alt = NA
for (i in 1:nrow(out)) {
  if (out[i,"ref"] == out[i,'a1']) {
    out[i,'alt'] <- out[i,'a2']
  } else{
    out[i,'alt'] <- out[i,'a1']
  }
}
out[1,'ref']
out
head(out)

out$A <- "Chromosome"
out$last <- 1
out <- out[,c(9,2,1,7,8,10)]
head(out)

ref <- read.table("../../Æ¯Çã/data/impute4/intersect.pos.txt")
head(ref)
out <- out[out$pos %in% ref$V1,]
write.table(out,"patent.snp.list.for.nexus.txt",col.names = F,row.names = F,quote = F)




###################################3
setwd("c:/Users/user/Desktop/KCDC/HLAimputation/patent/")

df <- read.table("KBA.ref.allele.for.HLAregion.txt",header = T)
han <- read.table("Han.legend",header = T)
head(df)
head(han)

out <- df[df$V1 %in% han$position,]
