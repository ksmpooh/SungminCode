vcf_file <- "chr20.dose.vcf.gz"
output_file <- "KOTRY.AR.SIRPA.txt"

#!/usr/bin/env Rscript

# 1. 인자 처리
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("사용법: Rscript SIRP_A_typing.R [raw.vcf.gz] [output.txt]")
}
vcf_file <- args[1]
output_file <- args[2]

# 2. 필요 라이브러리
library(tidyverse)
library(readxl)
library(stringr)

# 3. 경로 설정
ref_path <- "/DATA/smkim/AR/ref/SIPR_ref/SIPa_haplotype_20250131.xlsx"
sample_info_path <- "/DATA/smkim/AR/ref/SIPR_ref/KBA.QC.list.AR_1441pairs.20250325.txt"
tmp_txt <- tempfile(fileext = "tmp.txt")

pair_outfile <- paste0("KRKD.",output_file)


# 4. bcftools query 실행
query_cmd <- paste(
  "bcftools query -r 20:1895750-1896101 -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO[\\t%SAMPLE=%GT]\\n'",
  shQuote(vcf_file), "> ", tmp_txt
)
system(query_cmd)

# 5. 데이터 로딩
df <- readr::read_table(tmp_txt, col_names = FALSE)

# 6. 레퍼런스 및 샘플 리스트 불러오기
ref <- readxl::read_xlsx(ref_path, sheet = 6) %>%
  select(V1,V2,V4,V5,base_change,amino_change,ID,v1,v2,v4,v8,v9) %>%
  rename(CHR = V1, POS = V2, REF = V4, ALT = V5)

sample_ref <- read.table(sample_info_path, header = TRUE)
KR <- sample_ref %>% select(KBA_ID.KD) %>% rename(ID = KBA_ID.KD)
KD <- sample_ref %>% select(KBA_ID.KR) %>% rename(ID = KBA_ID.KR)

KRKD_samplelist <- bind_rows(KR, KD) %>%
  mutate(hap1 = paste0(ID, ".hap1"), hap2 = paste0(ID, ".hap2")) %>%
  pivot_longer(cols = c(hap1, hap2), names_to = "type", values_to = "hap") %>%
  select(hap)

# 7. 유사도 계산 함수
calc_similarity <- function(df, ref) {
  df <- df %>% select(-X5) %>%
    pivot_longer(cols = 5:last_col(), names_to = "original_col", values_to = "sample_geno") %>%
    separate(sample_geno, into = c("sample", "genotype"), sep = "=") %>%
    select(-original_col) %>%
    pivot_wider(names_from = sample, values_from = genotype) %>%
    rename(CHR = X1, POS = X2, REF = X3, ALT = X4)
  
  # 샘플 컬럼 이름 처리 (_ 개수 조건 분기)
  sample_cols <- colnames(df)[5:ncol(df)]
  underscore_count <- str_count(sample_cols, "_")
  new_sample_cols <- ifelse(
    underscore_count >= 3,
    str_split_fixed(sample_cols, "_", 4)[, 4],
    str_split_fixed(sample_cols, "_", 2)[, 1]
  )
  colnames(df)[5:ncol(df)] <- new_sample_cols
  
  # hap 분리
  df2 <- df
  valid_cols <- names(df)[5:ncol(df)]
  valid_cols <- valid_cols[valid_cols != ""]
  for (col_name in valid_cols) {
    df2 <- df2 %>%
      separate(col = col_name,
               into = c(paste0(col_name, ".hap1"), paste0(col_name, ".hap2")),
               sep = "\\|", remove = TRUE)
  }
  
  check_df <- df2 %>%
    mutate(ID = paste0(CHR, ":", POS, ":", REF, ":", ALT)) %>%
    filter(ID %in% ref$ID) %>%
    left_join(ref, by = "ID")
  
  hap_cols <- check_df %>% select(contains(".hap")) %>% names()
  v_cols <- check_df %>% select(starts_with("v")) %>% names()
  
  pairs_df <- expand.grid(hap = hap_cols, v = v_cols, stringsAsFactors = FALSE)
  
  result <- pairs_df %>%
    rowwise() %>%
    mutate(similarity_pct = mean(check_df[[hap]] == check_df[[v]], na.rm = TRUE) * 100) %>%
    ungroup()
  
  return(result)
}

# 8. 실행
out <- calc_similarity(df, ref)


out %>%
  group_by(hap) %>%
  slice_max(similarity_pct, n = 1, with_ties = FALSE) %>% 
  write.table(file = paste0("similarity_pct.",output_file), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

out_hap <- out %>%
  group_by(hap) %>%
  slice_max(similarity_pct, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(hap %in% KRKD_samplelist$hap) %>% 
  mutate(ID = str_split_fixed(hap, "\\.", 2)[, 1],
         hap = str_split_fixed(hap, "\\.", 2)[, 2]) %>% select(-similarity_pct) %>%
  pivot_wider(names_from = hap, values_from = v) %>%
  mutate(new_hap1 = ifelse(hap1 < hap2, hap1, hap2),
         new_hap2 = ifelse(hap1 < hap2, hap2, hap1),
         hap = paste0(new_hap1, "/", new_hap2)) %>%
  select(ID, hap1, hap2, hap)

# 9. 저장
write.table(out_hap, file = output_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
message("??? 결과 저장 완료: ", output_file)


# ???? KR/KD 쌍 기반 페어링 결과 추가 생성
message("???? KRKD SIRPα 페어링 정보 생성 중...")


#sample_info_path

sample_ref2 <- read.table(sample_info_path, header = TRUE) %>% select(KBA_ID.KR, KBA_ID.KD) 

sirp_df <- out_hap %>% select(ID, hap) %>% rename(KBA_ID.KR = ID, KR.hap = hap)
df <- sample_ref2 %>% left_join(sirp_df, by = "KBA_ID.KR")

colnames(sirp_df) <- c("KBA_ID.KD", "KD.hap")
df <- df %>% left_join(sirp_df, by = "KBA_ID.KD")

df <- df %>%
  mutate(sirp_missmatch = ifelse(KR.hap == KD.hap, 0, 1)) %>%
  select(KBA_ID.KR, KBA_ID.KD, KR.hap, KD.hap, sirp_missmatch) %>%
  rename(SIRPa_type.KR = KR.hap, SIRPa_type.KD = KD.hap)

write.table(df, file = pair_outfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
message("??? KRKD SIRPα 페어링 파일 저장 완료: ", pair_outfile)
