setwd("~/Desktop/KCDC/HLAimputation/patent/marker_list/")

library(BSgenome.Hsapiens.UCSC.hg19)
ref <- read.csv("HLA typing_patent.rsID.csv")
gene = "A"
gene = "B"
gene = "DRB1"

df <- read.table("JG.QCed.HLA_intersect_HLA.A_pruning.bim")
df <- read.table(paste0("JG.QCed.HLA_intersect_HLA.",gene,"_pruning.bim"))

#No TRAIT CHR POS

head(df)
colnames(df) <- c("chr","KCHIP_ID","tmp","Position","A1","A2")
df$SEQ <- NA

for(i in 1:dim(df)[1]){
  tempseq <- as.character(getSeq(Hsapiens,paste0("chr",df[i,1]),df[i,4]-5,df[i,4]+5))
  df[i,]$SEQ <- tempseq
}
head(ref)
head(df)
out <- merge(df[,c("Position","SEQ")],ref,by = "Position")

out <- out[,c("rs_ID","Chromosome","Position","REF.Allele","ALT.Allele","SEQ")]
head(out)

write.table(out,paste0("HLApatent_Purning.markerList_HLA.",gene,".txt"),sep = "\t",col.names = T,row.names = F,quote = F)
