setwd("/Volumes/DATA/JG/annotation/imputed/preData/")
library(BiocManager)
library(BiocVersion)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg19)
for (chr in 7:10) {
  chrom = as.character(chr)
  df <- read.table(paste0("JG.KR.imputation.mergeGen.processing.chr",chrom,".chunkMerge.vcf"))
  refallele <- getSeq(Hsapiens,paste0("chr",chrom),min(df$V2),max(df$V2))
  refallele <- as.character(refallele)
  df$ref <- 0
  head(df)
  k = 0
  tmp <- strtoi(min(df$V2))
  for (i in df$V2){
    k = k + 1
    j <- strtoi(i)
    r <- j - tmp + 1
    #print(i)
    #print(j)
    #print(j - 28477832)
    #df[k,1] <- j
    df[k,'ref'] <- substr(refallele,start = r,stop = r)
  }
  df<-as.data.frame(df)
  write.table(df,paste0("./withREF/JG.KR.imputation.mergeGen.processing.chr",chrom,".chunkMerge_withREF.txt"),sep = "\t",col.names = F,row.names = F,quote = F)
}


