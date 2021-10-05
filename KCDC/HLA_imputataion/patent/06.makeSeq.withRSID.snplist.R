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

##########



#20211005 추가
setwd("~/Desktop/KCDC/HLAimputation/patent/marker_list/")
ref <- read.table("nexus_output.txt")
head(ref)


A <- read.table("HLApatent_Purning.markerList_HLA.A.txt",header = T)
head(A)
B_ref <- read.table("JG.QCed.HLA_intersect_HLA.B_pruning.bim")
B <- read.table("HLApatent_Purning.markerList_HLA.B.txt",header = T)
head(B)
head(B_ref)
B_add <- B_ref[!(B_ref$V4 %in% B$Position),]
head(B_add)
DRB1 <- read.table("HLApatent_Purning.markerList_HLA.DRB1.txt",header = T)
head(DRB1)


DRB1_ref <- read.table("JG.QCed.HLA_intersect_HLA.DRB1_pruning.bim")
DRB1 <- read.table("HLApatent_Purning.markerList_HLA.DRB1.txt",header = T)
head(DRB1)
head(DRB1_ref)
DRB1_add <- DRB1_ref[!(DRB1_ref$V4 %in% DRB1$Position),]
head(DRB1_add)
