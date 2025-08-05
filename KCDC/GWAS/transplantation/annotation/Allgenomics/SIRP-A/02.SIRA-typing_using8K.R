### sipA
library(tidyverse)
#bcftools query -r 20:1895750-1896101 -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO[\t%SAMPLE=%GT]\n' chr20.dose.vcf.gz > sirp.region.txt

setwd("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/AR_YS/sirp/chr_20")
df <- read_table("sirp.region.txt",col_names = F)
head(df)
ref <- readxl::read_xlsx("../../../SIPa_haplotype_20250131.xlsx",sheet=6)
head(ref)

ref %>% select(V1,V2,V4,V5,base_change,amino_change,ID,v1,v2,v4,v8,v9) %>% rename(CHR = V1,POS = V2, REF = V4, ALT=V5) -> ref
head(ref)
head(df)






df %>% select(-X5)%>% pivot_longer(cols = 5:last_col(), names_to = "original_col",values_to = "sample_geno") %>% 
  separate(col  = sample_geno,into = c("sample", "genotype"),sep  = "=") %>% 
  select(-original_col) %>%
  pivot_wider(
    names_from  = sample,    # "K001_CEL01", "K002_CEL01", ...
    values_from = genotype   # "0|0", "0|1", "1|1" ...
  ) %>% rename(CHR = X1,POS = X2,REF =X3 ,ALT = X4) -> new_df

head(new_df)
str_split_fixed(colnames(new_df)[5:ncol(new_df)],"_",4)[,4] -> new_col
colnames(new_df)[5:ncol(new_df)] <- new_col 

head(new_df)
              


df2 <- new_df  # 원본 복사

# 5번째 컬럼부터 마지막 컬럼까지 반복
for (col_name in names(new_df)[5:ncol(new_df)]) {
  # 예: col_name = "sampleA"
  df2 <- df2 %>%
    separate(
      col   = col_name,                 # 현재 열
      into  = c(paste0(col_name, ".hap1"),
                paste0(col_name, ".hap2")),# 나눌 뒤의 열 이름
      sep   = "\\|",                    # '|' 기준으로 분리
      remove = TRUE                     # 원본 열 제거(기본값: TRUE)
    )
}
head(df2)

df2 %>% mutate(ID = paste0(CHR,":",POS,":",REF,":",ALT)) %>%
  filter(ID %in% ref$ID) %>% left_join(ref) -> check_df

hap_cols <- check_df %>%
  select(contains(".hap")) %>%  # 열 이름 중 ".hap"이 들어가는 것
  names()

v_cols <- check_df %>%
  select(starts_with("v")) %>%  # 열 이름 중 "v"로 시작하는 것
  names()

pairs_df <- expand.grid(hap = hap_cols, v = v_cols, stringsAsFactors = FALSE)

# 3) 모든 쌍에 대해 유사도 계산
result <- pairs_df %>%
  rowwise() %>%
  mutate(
    similarity_pct = mean(check_df[[hap]] == check_df[[v]], na.rm = TRUE) * 100
  ) %>%
  ungroup()


head(result)
result %>% group_by(hap) %>% summarise(SIRPA_hap = max(similarity_pct))
result %>% group_by(hap) %>% slice_max(similarity_pct, n = 1, with_ties = FALSE) %>% #head()
  ungroup() %>% count(v)

result %>%   pivot_wider(
  names_from  = v, 
  values_from = similarity_pct
) %>% group_by(hap) %>%
  mutate(
  max_col = c("v1", "v2","v4","v8", "v9")[which.max(c_across(c(v1, v2, v8, v9)))],
#  max_col_max = c("v1", "v2", "v8", "v9")[which.max(c(v1, v2, v8, v9))]
) %>%
  ungroup() %>%  head()

result %>% group_by(hap) %>% slice_max(similarity_pct, n = 1, with_ties = FALSE) %>% ungroup()  %>%
  count(v,similarity_pct)
  
result %>% group_by(hap) %>% slice_max(similarity_pct, n = 1, with_ties = FALSE) %>% ungroup()  %>% count(similarity_pct,v) %>%
  #ggplot(aes(x=factor(similarity_pct),fill = v)) +
  ggplot(aes(x= factor(similarity_pct),fill = v,y=n)) +
  geom_bar(stat='identity', position = 'dodge')
  #geom_histogram()
  #geom_bar(position="identity")
  #geom_histogram(stat = 'count')
  

head(result)

  ###
calc_similarity <- function(df, ref) {
  
  # 1) X5 컬럼 제거 후, pivot_longer()로 샘플 정보 펼치기
  new_df <- df %>%
    select(-X5) %>% 
    pivot_longer(
      cols = 5:last_col(), 
      names_to = "original_col",
      values_to = "sample_geno"
    ) %>%
    separate(
      col  = sample_geno,
      into = c("sample", "genotype"),
      sep  = "="
    ) %>%
    select(-original_col) %>%
    pivot_wider(
      names_from  = sample,    # "K001_CEL01", "K002_CEL01", ...
      values_from = genotype   # "0|0", "0|1", "1|1" ...
    ) %>%
    rename(
      CHR = X1,
      POS = X2,
      REF = X3,
      ALT = X4
    )
  
  # 2) 새로 만들어진 5번째 컬럼부터의 이름에서 '_' 구분자로 split하여, 4번째 조각을 새로운 컬럼 이름으로 지정
  #new_col <- str_split_fixed(colnames(new_df)[5:ncol(new_df)], "_", 4)[,4]
  #new_col <- str_split_fixed(colnames(new_df)[5:ncol(new_df)], "_", 2)[,1]
  colnames(new_df)[5:ncol(new_df)] <- new_col
  
  # 3) new_df를 복사한 df2에서, 5번째 열부터 '|' 기준으로 hap1, hap2 두 열로 분리
  df2 <- new_df
  for (col_name in names(new_df)[5:ncol(new_df)]) {
    df2 <- df2 %>%
      separate(
        col   = col_name,
        into  = c(paste0(col_name, ".hap1"), paste0(col_name, ".hap2")),
        sep   = "\\|",
        remove = TRUE
      )
  }
  
  # 4) ID 컬럼 생성 후, ref와 ID 기준으로 left_join
  check_df <- df2 %>%
    mutate(ID = paste0(CHR, ":", POS, ":", REF, ":", ALT)) %>%
    filter(ID %in% ref$ID) %>%
    left_join(ref, by = "ID")
  
  # 5) hap 계열 열 이름, v 계열 열 이름 추출
  hap_cols <- check_df %>%
    select(contains(".hap")) %>%
    names()
  
  v_cols <- check_df %>%
    select(starts_with("v")) %>%
    names()
  
  # 6) hap-열과 v-열의 모든 쌍(expand.grid) 생성
  pairs_df <- expand.grid(hap = hap_cols, v = v_cols, stringsAsFactors = FALSE)
  
  # 7) 모든 쌍에 대해 유사도(%) 계산
  result <- pairs_df %>%
    rowwise() %>%
    mutate(
      similarity_pct = mean(check_df[[hap]] == check_df[[v]], na.rm = TRUE) * 100
    ) %>%
    ungroup()
  
  return(result)
}

head(df)
out <- calc_similarity(df, ref)
head(out)


##########
df <- read_table("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/SIRPa/ori_pair/kotry.oridatasirp.region.txt",col_names = F)
head(df[1:5,1:10])
colnames(df)
df$X5
out <- calc_similarity(df, ref)
head(out)
table(out$v)
out %>% filter(v =="v4")

sample_ref <- read.table("~/Desktop/KCDC/transplantation/allogenomic/02.KR.KD.immune.cell.co-signal_targetGene.alleleMatching01.Score_Sum.txt",header = T)
sample_ref %>% select(1) -> KD
sample_ref %>% select(2) -> KR
colnames(KD) <- c("ID")
colnames(KR) <- c("ID")
KRKD_samplelist <- rbind(KR,KD) 
KRKD_samplelist %>% mutate(hap = paste0(ID,".hap1")) %>% rbind(KRKD_samplelist %>% mutate(hap = paste0(ID,".hap2"))) -> KRKD_samplelist


o
head(out)
out %>%   pivot_wider(names_from  = v, values_from = similarity_pct) %>% group_by(hap) %>%
  mutate(max_col = c("v1", "v2","v4","v8","v9")[which.max(c_across(c(v1, v2,v4, v8, v9)))]) %>%
  ungroup() %>% filter(hap %in% KRKD_samplelist$hap) -> a
#%>% filter(max_col %in% c('v4','v8','v9')) -> a
a
out %>% group_by(hap) %>% slice_max(similarity_pct, n = 1, with_ties = FALSE) %>% ungroup()  %>% 
  filter(hap %in% KRKD_samplelist$hap) %>% 
  count(v,similarity_pct) %>% pivot_wider(names_from = similarity_pct,values_from = n)

out %>% group_by(hap) %>% slice_max(similarity_pct, n = 1, with_ties = FALSE) %>% ungroup()  %>% filter(hap %in% KRKD_samplelist$hap) %>%
  mutate(ID = str_split_fixed(hap,"\\.",2)[,1],hap = str_split_fixed(hap,"\\.",2)[,2]) %>% select(-similarity_pct) %>%
  pivot_wider(names_from = hap,values_from = v) %>% mutate(new_hap1 = ifelse(hap1 < hap2,hap1,hap2),new_hap2 = ifelse(hap1 < hap2,hap2,hap1)) %>%
  mutate(hap = paste0(new_hap1,"/",new_hap2)) -> out_hap
head(out_hap)

write.table(out_hap %>% select(1,4,5,6),"~/Desktop/KCDC/transplantation/allogenomic/2025_AR/SIRPa/ori_pair/KRKD.SIRP-a.1148pair.txt",col.names = T,row.names = F)
out_hap %>% count(hap) %>% pivot_wider(names_from = hap,values_from = n) -> a
a 
out_hap %>% count(hap)

out %>% group_by(hap) %>% slice_max(similarity_pct, n = 1, with_ties = FALSE) %>% ungroup()  %>% filter(v=="v8")

out %>% group_by(hap) %>% slice_max(similarity_pct, n = 1, with_ties = FALSE) %>% ungroup()  %>%
  count(v,similarity_pct) %>%
  ggplot(aes(x=similarity_pct,y=n,fill=v)) + 
  geom_bar(stat = "identity",position = "dodge")



out %>% group_by(hap) %>% slice_max(similarity_pct, n = 1, with_ties = FALSE) %>% ungroup()  %>% #filter(hap %in% KRKD_samplelist$hap) %>%
  mutate(ID = str_split_fixed(hap,"\\.",2)[,1],hap = str_split_fixed(hap,"\\.",2)[,2]) %>% select(-similarity_pct) %>%
  pivot_wider(names_from = hap,values_from = v) %>% mutate(new_hap1 = ifelse(hap1 < hap2,hap1,hap2),new_hap2 = ifelse(hap1 < hap2,hap2,hap1)) %>%
  mutate(hap = paste0(new_hap1,"/",new_hap2)) %>% select(ID,hap) -> a

head(a)
a %>% filter(str_detect(ID,"KR")) %>% mutate(ref = str_replace_all(ID,"_KR","")) %>% rename(ID_KR = ID, hap_KR = hap) %>% left_join(
  a %>% filter(str_detect(ID,"KD")) %>% mutate(ref = str_replace_all(ID,"_KD","")) %>% rename(ID_KD = ID, hap_KD = hap)
) %>% select(ID_KR,ID_KD,hap_KR,hap_KD) %>% write.table("../KRKD.SIRP-a.31pair_byYS.txt",col.names = T,row.names = F,quote = F,sep = "\t")

