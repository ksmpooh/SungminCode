##### reference panel ref/alt position
source("http://bioconductor.org/biocLite.R")
library(BiocManager)
library(BiocVersion)



BiocManager::install("BSgenome")
BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')

library(BSgenome.Hsapiens.UCSC.hg19)

refallele <- getSeq(Hsapiens,"chr6",28477833,33448188)
refallele <- as.character(refallele)


