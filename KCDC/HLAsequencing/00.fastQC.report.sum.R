library(fastqcr)
setwd("/Volumes/DATA/HLA_seq/01.fastaQC/")
setwd("~/")
qc.files <- list.files("./short/","*zip")
qc.files

qc <- qc_aggregate("./short/",qc.files)
qc <- qc_aggregate(paste0("./short/",qc.files))
qc <- qc_aggregate("/Volumes/DATA/HLA_seq/01.fastaQC/short/*zip")
qc_short <- qc_aggregate("/Volumes/DATA/HLA_seq/01.fastaQC/short/")
head(qc)
library(tidyverse)

library(writexl)
qc_stats(qc)
write_xlsx(qc,"/Users/ksmpooh/Desktop/KCDC/HLA_seq/01.fastQC/HLA_shortread.fastQC.Table.xlsx")
write_xlsx(qc_stats(qc_short),"/Users/ksmpooh/Desktop/KCDC/HLA_seq/01.fastQC/HLA_shortread.fastQC_qcstat.Table.xlsx")


qc <- qc_aggregate("/Volumes/DATA/HLA_seq/01.fastaQC/long/")
write_xlsx(qc,"/Users/ksmpooh/Desktop/KCDC/HLA_seq/01.fastQC/HLA_longread.fastQC.Table.xlsx")
write_xlsx(qc_stats(qc),"/Users/ksmpooh/Desktop/KCDC/HLA_seq/01.fastQC/HLA_longread.fastQC_qcstat.Table.xlsx")
