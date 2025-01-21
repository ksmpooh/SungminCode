## Supplementary table 3
rm(list = ls())
library(openxlsx)
library(tidyverse)

setwd("~/Dropbox/SMC_AD_WGS_working/Data/")


#cate = read_tsv(cate1) %>% dplyr::rename(sample = ...1)

TableS4 = createWorkbook()

# Table S4. Association Test for Rare Noncoding Variants with AD
# Key
openxlsx::addWorksheet(TableS4, "Table S4 Overview")
openxlsx::addWorksheet(TableS4, "Table S4A")
openxlsx::addWorksheet(TableS4, "Table S4B")
openxlsx::addWorksheet(TableS4, "Table S4C")
openxlsx::addWorksheet(TableS4, "Table S4D")
openxlsx::addWorksheet(TableS4, "Table S4E")
openxlsx::addWorksheet(TableS4, "Table S4F")

############# Table S4A Copy number variation association analysis   ################


res<-data.table::fread("logistic_regression_result")

new_column_names <- c("Gene", "CNV_type", "CytoBand", "chr", "start", "end", "length", "AnnotCNV_ID", "pvals", "BONF", "FDR", "BON_c")

# 열 이름 변경
colnames(res) <- new_column_names

res = res %>% select("chr", "start", "end", "CNV_type", "length", "CytoBand","Gene" , "pvals", "BONF", "FDR", "BON_c")
res$chr <- as.numeric(res$chr)
res$start <- as.numeric(res$start)
res$end <- as.numeric(res$end)

res <- res[order(res$chr, res$start, res$end), ]

res = res %>% filter(pvals < 0.05)

res = res %>% mutate(chr = paste0('chr',chr))

openxlsx::writeData(wb = TableS4, sheet = "Table S4A", x = res, startCol = 1, startRow = 1, rowNames = F)


# 
# openxlsx::writeData(wb = TableS4, sheet = "Table S4A", x = S4A, startCol = 1, startRow = 1, rowNames = F)

############### Table S4B Result of STR length test  ##############

S4A = read_tsv("STR/eh.v5_bed.tsv")

S4A = S4A %>% dplyr::rename(`Reference length` = ref_length) %>% select("chr", "start", "end", "motif","Reference length", "Genomic annotation", "Distance to nearest TSS" ,"Gene ENSEMBL ID"        ,"Gene symbol"  )

S4B = read_tsv("STR/eh_RCL40_visual_single_str_test_inv_noPC_240419_PC_rm.txt") %>% left_join(S4A %>% mutate(id = paste(chr,start,motif,sep = '-')))
S4B <- S4B %>% select("chr","start",    "end","motif",  "allele_count_logistic_p","allele_count_effect",     "allele_count_mean_control","allele_count_mean_ad",   "allele_count_sd","allele_count_inv_norm_logistic_p" ,"allele_count_inv_norm_effect")

S4B <- S4B %>% filter(allele_count_inv_norm_logistic_p < 0.05)

openxlsx::writeData(wb = TableS4, sheet = "Table S4B", x = S4B, startCol = 1, startRow = 1, rowNames = F)


############### Table S4D	STR threshold test ##################

S4D = read_tsv('STR/eh_RCL40_visual_single_str_thr_test_ref_240419_PC_rm.txt') %>% left_join(S4A %>% mutate(id = paste(chr,start,motif,sep = '-')))
S4D <- S4D %>% select("chr","start",    "end","motif",  "thr_1_pval","thr_1_odds"  ,"thr_5_pval","thr_5_odds","thr_10_pval","thr_10_odds"    ,"thr_20_pval","thr_20_odds")

S4D <- S4D %>%
  filter(thr_1_pval < 0.05 | thr_5_pval < 0.05 | thr_10_pval < 0.05 | thr_20_pval < 0.05)

openxlsx::writeData(wb = TableS4, sheet = "Table S4D", x = S4D, startCol = 1, startRow = 1, rowNames = F)

############### Table S4E	STR outlier test ##################
S4E = read_tsv("STR/eh_RCL40_visual_expansion_count_sample_240419_PC_rm.txt") %>% select( "sample","outlier_count" ,'s_avg_depth',"DX",'RCL40_visual') %>% dplyr::rename('Amyloid beta positivity' = 'RCL40_visual',
                                                                                                                                                                         'Clinical diagnosis' = DX)


openxlsx::writeData(wb = TableS4, sheet = "Table S4E", x = S4E, startCol = 1, startRow = 1, rowNames = F)


# ############### Table S4F	GO term ##################

GO_df <- read_tsv('STR/eh_STR_over40_goterm.tsv')

GO_df <- GO_df %>% arrange(p.adjust) %>% select(ONTOLOGY, ID,Description,  pvalue,p.adjust)

openxlsx::writeData(wb = TableS4, sheet = "Table S4F", x = GO_df, startCol = 1, startRow = 1, rowNames = F)


saveWorkbook(wb = TableS4, file = "~/Dropbox/SMC_AD_WGS_paper/Tables/Supplementary Table 4 v3.1.xlsx", overwrite=T)



