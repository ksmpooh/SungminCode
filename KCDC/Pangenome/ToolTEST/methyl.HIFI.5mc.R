####
# 라이브러리 로드
library(tidyverse)
library(reshape2)
library(pheatmap)
setwd("/Users/ksmpooh/Desktop/KU/@research/STR/methyl")
# 데이터 불러오기
jaccard_data <- read.table("jaccard_results.tsv", header = FALSE, sep = "\t", col.names = c("sample1", "sample2", "jaccard"))
head(jaccard_data)

# 데이터 전처리: 행렬 형태로 변환
jaccard_matrix <- jaccard_data %>%
  spread(key = sample2, value = jaccard, fill = NA) %>%
  column_to_rownames(var = "sample1") %>%
  as.matrix()

head(jaccard_matrix)
# 대칭 행렬을 위한 상하 복사 (NA 대체)
jaccard_matrix[lower.tri(jaccard_matrix)] <- t(jaccard_matrix)[lower.tri(jaccard_matrix)]

# 히트맵 그리기
pheatmap(jaccard_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",show_rownames = F,show_colnames = F,
         clustering_method = "average",
         main = "Jaccard Index Heatmap")

### methyl
library(tidyverse)
setwd("/Users/ksmpooh/Desktop/KU/@research/STR/methyl/long_merge_withinterBED")
tmp <- read_table("NIH23J3016218_sorted.pb_CpG.combined.merged.bed",col_names = F)
head(tmp)
colnames(tmp) <- c("chrom","start","end","MOTIF","STR_ID","STR_region",
                   "chrom_methyl","start_methyl","end_methyl","modification_score","happlotype",
                   "coverage","modified_site_count","unmodified_site_count","average_mod_score",
                   "filename","ID")
tmp <- tmp %>% select(7:15)
flist = grep(list.files("./"),pattern = "bed", value=TRUE)
flist
df <- NULL
for (i in flist) {
  tmp <- read_table(i,col_names = FALSE)
  tmp$filename = i
  #  tmp %>% filter(X6 == "PASS") %>% select(X3,X4,X5) -> tmp
  df <- rbind(df,tmp)
}
head(df)

head(df)

ref <-read.table("~/Desktop/KCDC/pangenome/00.datacheck/KBA.Long_Revio_Nanopore_short.IDmatchinagtable.txt",header = T)
head(ref)
df$ID <- 1
for (i in ref$Revio) {
  df %>% mutate(ID = ifelse(str_detect(filename,i),i,ID)) -> df
}
head(df)
df %>% count(X1 == X7)
colnames(df) <- c("chrom","start","end","MOTIF","STR_ID","STR_region",
                  "chrom_methyl","start_methyl","end_methyl","modification_score","happlotype",
                  "coverage","modified_site_count","unmodified_site_count","average_mod_score",
                  "filename","ID")
head(df)
df %>% filter(ID != "NIH20N2594890") -> df
df %>% select(ID,STR_ID,modification_score,coverage) %>% count(ID,STR_ID) -> df_methyl_count
#df %>% select(ID,STR_ID,modification_score,coverage) %>% count(ID,STR_ID) -> df_methyl_count_byID




ggplot(df,aes(x=modification_score)) + 
  geom_histogram()


#df %>% write.table("../methyl.INFO.interSTR.txt",col.names = T,row.names = F,quote = F,sep = "\t")

#####################
df <- read_table("methyl.INFO.interSTR.txt")
df %>% select(ID,STR_ID,modification_score,coverage) %>% count(ID,STR_ID) -> df_methyl_count
df %>% select(ID,STR_ID,modification_score,coverage) %>% 
  filter(modification_score >= 30) %>% count(ID,STR_ID) -> df_methyl_count_upper30

head(df)
#head(df_methyl_count)

df_methyl_count %>% count(STR_ID) %>% filter(n==66) -> df_methyl_all_site

head(df)
df_methyl_count
head(df_methyl_count_upper30)

df_methyl_count %>% head()
  

df_methyl_count_upper30 %>%  count(STR_ID) %>% dim() #13938
df_methyl_count %>%  count(STR_ID) %>% dim() #15527

df_methyl_count_upper30 
df_methyl_count_upper30 %>%
  group_by(STR_ID) %>%
  summarise(Variance = var(n)) -> df_methyl_count_upper30_vari

df_methyl_count_upper30 %>% count(STR_ID) -> df_methyl_count_upper30_strcount


df_methyl_count_upper30 %>%
  group_by(STR_ID) %>%
  summarise(Variance = var(n)) -> df_methyl_count_upper30_vari

df_methyl_count_upper30_vari  %>% left_join(df_methyl_count_upper30_strcount) %>%
  ggplot(aes(x=n,y=Variance)) + 
  geom_point() +
  labs(title = "Distribution of methyl count variant (point : STR ID)") + 
  xlab("# of Sample")


head(df_methyl_count_upper30_vari)

df_methyl_count_upper30_vari %>%
  ggplot(aes(x="",y=Variance)) + 
  geom_violin()

df_methyl_count_upper30 %>%
  ggplot(aes(x = n)) +
  geom_histogram(binwidth = 2, fill = "orange", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of methyl sites Values",
       x = "",
       y = "Frequency")


df_methyl_count_upper30 


head(df_methyl_count_upper30_vari)

df_methyl_count_upper30_vari %>% filter(Variance >= 20) -> df_methyl_count_upper30_vari_uppor20
df_methyl_count_upper30 %>% filter(STR_ID %in% df_methyl_count_upper30_vari_uppor20$STR_ID)

df_methyl_count_upper30_vari_uppor20 %>% left_join(ref_str %>% rename(STR_ID = ID)) %>%
  mutate(ref_STR_length = end-start) %>% select(-chrom,-start,-end)
#%>% select(-chrom,-start)
head(ref)
ref_str <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header = T)
head(ref_str)
