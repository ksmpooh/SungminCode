library(tidyverse)


frq <- read_table("/DATA/smkim/AR/asso_fastLMM/00.rawDATA/stats/KOTRY.AR_1441.freq.frq")
hwe <- read_table("/DATA/smkim/AR/asso_fastLMM/00.rawDATA/stats/KOTRY.AR_1441.hwe.hwe")
qc <- read_table("/DATA/smkim/AR/asso_fastLMM/00.rawDATA/stats/KOTRY.AR_1441.missing.lmiss")

#fastlmm <- read_table("fastlmm_results.keep.order.txt")
#fastlmm <- read_delim("fastlmm_results.keep.order.txt",delim="\t")
bcftools <- read_table("KOTRY.AR_2025.bcftools.query.merge",col_names = F)
colnames(bcftools) <- c("ID","CHR","POS","REF","ALT","TYPE","R2","IMPUTED","HWE")
snptest <- read_table("KOTRY.AR_2025.snptest_state.query.merge")
qctools <- read_table("KOTRY.AR_2025.qctools_state.query.merge")
qctools %>% rename(SNP = rsid) %>% left_join(frq) %>% left_join(qc) %>% left_join(hwe %>% select(-X10)) -> df
liftover <-read_table("/DATA/smkim/AR/asso_fastLMM/result_merge/liftover.check",col_names=F)
#colnames(liftover) <- c("chr","pos","SNPID")
liftover %>% filter(X3 %in% fastlmm$SNP) %>% 
  separate(X3, into = c("chrom", "pos", "ref", "alt"), sep = ":", remove = FALSE) %>%
  mutate(SNPID = paste0(chrom, ":", X2, ":", ref, ":", alt)) %>%
  rename(SNPID_hg19=X3,chr=X1, position=X2) %>% filter(chr == chrom) %>%
  select(SNPID_hg19,SNPID,chr,position) -> liftover_prop
#AR_BX.fastlmm_results.keep.order.txt  
#AR_TREAT.fastlmm_results.keep.order.txt  
#GF.fastlmm_results.keep.order.txt  run.py

#certainty <- read_table("/DATA/smkim/AR/imp_kis2/04.AR_1141_SNPQC/vcf/certainty/certainty.merge.txt")

topmed <- read_table("/DATA/smkim/AR/asso_fastLMM/result_merge/withTopmed/new/matched_output/matched_output.final.txt")

topmed %>% na.omit() %>% 
  separate(new_id, into = c("chr", "position", "A0", "A1"), sep = ":", remove = FALSE) %>%
  rename(SNPID = new_id,SNPID_hg19 = old_id) %>% 
  select(SNPID,SNPID_hg19,chr,position,A0,A1) -> topmed_proc
  
  
#fastlmm <- read_table("/DATA/smkim/AR/asso_fastLMM/asso_result/01.rawResult/AR_TREAT.fastlmm_results.keep.order.txt")
fastlmm <- read_delim("/DATA/smkim/AR/asso_fastLMM/asso_result/01.rawResult/AR_TREAT.fastlmm_results.keep.order.txt",delim="\t")
fastlmm %>% 
  select(SNP, Chr, ChrPos, PValue, SnpWeight, SnpWeightSE) %>% 
  left_join(qctools %>%
      rename(SNP = rsid) %>%
      select(SNP, alleleB_frequency, alleleA, alleleB, HW_exact_p_value, total) %>%
      rename(exp_freq_a1 = alleleB_frequency,
             #A0 = alleleA, A1 = alleleB,
             HWE = HW_exact_p_value,
             n_total = total)) %>%
  left_join(bcftools %>%
              rename(SNP = ID) %>%
              select(SNP, R2, IMPUTED,TYPE) %>% 
              mutate(IMPUTED = ifelse(IMPUTED == 1, 1, 0)) %>%
              rename(info = R2, imputed = IMPUTED)) %>%  #filter(TYPE == "SNP") %>% select(-TYPE) %>%
  left_join(snptest %>%
      select(rsid, average_maximum_posterior_call) %>%
      rename(SNP = rsid, certainty = average_maximum_posterior_call)) %>%
  rename(SNPID_hg19 = SNP,
         chr_hg19 = Chr,
         position_hg19 = ChrPos,
         pvalue = PValue,
         beta = SnpWeight,
         SE = SnpWeightSE) %>% filter(SNPID_hg19 %in% topmed_proc$SNPID_hg19) %>%
  left_join(topmed_proc) %>% na.omit() %>%
  select(SNPID, chr, position, exp_freq_a1, A0, A1, beta, SE, pvalue, HWE, n_total, imputed, info, certainty) -> out

write.table(out,"/DATA/smkim/AR/asso_fastLMM/asso_result/02.processing_result/AR_TREAT.fastLMM.Result.processing.txt",col.names = T,row.names = F,quote = F,sep = "\t")


fastlmm <- read_delim("/DATA/smkim/AR/asso_fastLMM/asso_result/01.rawResult/AR_BX.fastlmm_results.keep.order.txt",delim="\t")


fastlmm %>% 
  select(SNP, Chr, ChrPos, PValue, SnpWeight, SnpWeightSE) %>% 
  left_join(qctools %>%
              rename(SNP = rsid) %>%
              select(SNP, alleleB_frequency, alleleA, alleleB, HW_exact_p_value, total) %>%
              rename(exp_freq_a1 = alleleB_frequency,
                     #A0 = alleleA, A1 = alleleB,
                     HWE = HW_exact_p_value,
                     n_total = total)) %>%
  left_join(bcftools %>%
              rename(SNP = ID) %>%
              select(SNP, R2, IMPUTED,TYPE) %>% 
              mutate(IMPUTED = ifelse(IMPUTED == 1, 1, 0)) %>%
              rename(info = R2, imputed = IMPUTED)) %>%  #filter(TYPE == "SNP") %>% select(-TYPE) %>%
  left_join(snptest %>%
              select(rsid, average_maximum_posterior_call) %>%
              rename(SNP = rsid, certainty = average_maximum_posterior_call)) %>%
  rename(SNPID_hg19 = SNP,
         chr_hg19 = Chr,
         position_hg19 = ChrPos,
         pvalue = PValue,
         beta = SnpWeight,
         SE = SnpWeightSE) %>% filter(SNPID_hg19 %in% topmed_proc$SNPID_hg19) %>%
  left_join(topmed_proc) %>% na.omit() %>%
  select(SNPID, chr, position, exp_freq_a1, A0, A1, beta, SE, pvalue, HWE, n_total, imputed, info, certainty) -> out

write.table(out,"/DATA/smkim/AR/asso_fastLMM/asso_result/02.processing_result/AR_BX.fastLMM.Result.processing.txt",col.names = T,row.names = F,quote = F,sep = "\t")


fastlmm <- read_delim("/DATA/smkim/AR/asso_fastLMM/asso_result/01.rawResult/GF.fastlmm_results.keep.order.txt",delim="\t")


fastlmm %>% 
  select(SNP, Chr, ChrPos, PValue, SnpWeight, SnpWeightSE) %>% 
  left_join(qctools %>%
              rename(SNP = rsid) %>%
              select(SNP, alleleB_frequency, alleleA, alleleB, HW_exact_p_value, total) %>%
              rename(exp_freq_a1 = alleleB_frequency,
                     #A0 = alleleA, A1 = alleleB,
                     HWE = HW_exact_p_value,
                     n_total = total)) %>%
  left_join(bcftools %>%
              rename(SNP = ID) %>%
              select(SNP, R2, IMPUTED,TYPE) %>% 
              mutate(IMPUTED = ifelse(IMPUTED == 1, 1, 0)) %>%
              rename(info = R2, imputed = IMPUTED)) %>%  #filter(TYPE == "SNP") %>% select(-TYPE) %>%
  left_join(snptest %>%
              select(rsid, average_maximum_posterior_call) %>%
              rename(SNP = rsid, certainty = average_maximum_posterior_call)) %>%
  rename(SNPID_hg19 = SNP,
         chr_hg19 = Chr,
         position_hg19 = ChrPos,
         pvalue = PValue,
         beta = SnpWeight,
         SE = SnpWeightSE) %>% filter(SNPID_hg19 %in% topmed_proc$SNPID_hg19) %>%
  left_join(topmed_proc) %>% na.omit() %>%
  select(SNPID, chr, position, exp_freq_a1, A0, A1, beta, SE, pvalue, HWE, n_total, imputed, info, certainty) -> out

write.table(out,"/DATA/smkim/AR/asso_fastLMM/asso_result/02.processing_result/GF.fastLMM.Result.processing.txt",col.names = T,row.names = F,quote = F,sep = "\t")




#indel       n
#<chr>   <int>
#1 indel 1928656
#2 snp   7402552

'''
> colnames(fastlmm)
[1] "sid_index"       "SNP"             "Chr"             "GenDist"
[5] "ChrPos"          "PValue"          "SnpWeight"       "SnpWeightSE"
[9] "EffectSize"      "SnpFractVarExpl" "Mixing"          "Nullh2"
> colnames(snptest)
[1] "alternate_ids"                  "rsid"
[3] "chromosome"                     "position"
[5] "alleleA"                        "alleleB"
[7] "index"                          "average_maximum_posterior_call"
[9] "info"                           "cohort_1_AA"
[11] "cohort_1_AB"                    "cohort_1_BB"
[13] "cohort_1_NULL"                  "all_AA"
[15] "all_AB"                         "all_BB"
[17] "all_NULL"                       "all_total"
[19] "all_maf"                        "missing_data_proportion"
[21] "comment"
> colnames(bcftools)
[1] "ID"      "CHR"     "POS"     "REF"     "ALT"     "TYPE"    "R2"
[8] "IMPUTED" "HWE"
'''


#df %>% select(P,HW_exact_p_value,HW_lrt_p_value) %>% mutate(check = P - HW_exact_p_value) %>% summary()


## fastlmm
#sid_index SNP     Chr GenDist  ChrPos  PValue SnpWeight SnpWeightSE EffectSize

## qctools snp
#alternate_ids rsid chromosome position alleleA alleleB comment HW_exact_p_value HW_lrt_p_value alleleA_count alleleB_count alleleA_frequency alleleB_frequency minor_allele_frequency minor_allele major_allele info impute_info missing_proportion A B AA AB BB NULL total

#fastlmm SNP: SNPID 	SNP ID as rs number  
#fastlmm Chr: chr 	chromosome number. Use symbols X, XY, Y and mt for non-autosomal markers. 
#fastlmm GenDist(hg19): position 	physical position for the reference sequence (indicate build 38 in readme file) 
#qctools snp MAF: exp_freq_a1	Expected allele frequency of A1
#plink --freq A2 :A0	Reference allele
#plink --freq A1 :A1	Effect allele
#fastlmm SnpWeight:beta 	beta estimate from genotype-phenotype association, “NA” if not available
#fastlmm SnpWeightSE:SE 	standard error of beta estimate, “NA” if not available
#fastlmm PValue:pvalue	p-value of test statistic, “NA” if not available 
#plink --hwe P:HWE 	Hardy-Weinberg equilibrium p-value for each SNP 
#plink --freq A1:n_total 	total sample with phenotype and genotype for SNP 
#bcftools:imputed 	1/0 coding; 1=imputed SNP, 0=if directly typed 
#bcftools: Chr:info 	Imputation quality 
#???:certainty	average certainty of best-guess genotypes

### change rs ID

library(tidyverse)
library(readr)
library(stringr)

# 1. 경로 설정
target_dir <- "/DATA/smkim/AR/asso_fastLMM/asso_result/02.processing_result/onlySNP"
sorted_files <- list.files(path = target_dir, pattern = "sorted", full.names = TRUE)

# 2. rsID 매핑 테이블 (중복 제거)
rs_ref <- read_table("/DATA/smkim/AR/asso_fastLMM/asso_result/02.processing_result/onlySNP/rsID_info/snp_to_rsid.checked.txt") %>%
  distinct(SNPID, .keep_all = TRUE)

# 3. 파일별 처리
for (file in sorted_files) {
  df <- read_table(file)
  
  df2 <- df %>%
    left_join(rs_ref, by = "SNPID") %>%
    mutate(SNPID_new = ifelse(is.na(rsID), SNPID, rsID)) %>%
    select(SNPID = SNPID_new, everything(), -rsID, -SNPID)
  
  output_path <- str_replace(file, ".txt", ".rsid.txt")
  
  write.table(df2, file = output_path, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
}


