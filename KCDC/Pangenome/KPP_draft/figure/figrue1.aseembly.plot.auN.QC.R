### auN novel

library(tidyverse)
library(googledrive)
library(googlesheets4)

drive_auth()
drive_find("type='csv")
#https://docs.google.com/spreadsheets/d/1Y8cqzCizaQ9vvAzoAEhqqGW7GaplxRRnmPsQ0dK6o30/edit?usp=sharing
url <- "https://docs.google.com/spreadsheets/d/1Y8cqzCizaQ9vvAzoAEhqqGW7GaplxRRnmPsQ0dK6o30/edit?usp=sharing"
# URL에서 file id 추출
id <- sub(".*?/d/([^/]+).*", "\\1", url)

# 인증(처음 1회: 브라우저 뜸)
drive_auth()

drive_download(as_id(id), path = "", overwrite = TRUE)


url <- "https://docs.google.com/spreadsheets/d/1Y8cqzCizaQ9vvAzoAEhqqGW7GaplxRRnmPsQ0dK6o30/edit?usp=sharing"
drive_auth()

df <- read_sheet(url,sheet = 4)  # 기본: 첫 번째 시트
df

props <- sheet_properties(url)
props$name
#13,"1.3.asm_correcteness1_inspector"
#14,"1.3.asm_correcteness1_merqury"   
#11 "1.2.assembly_contiguity" 

assembly_aun <- read_sheet(url,sheet = 11)  
inspector <- read_sheet(url,sheet = 13)  # 기본: 첫 번째 시트
mercury <- read_sheet(url,sheet = 14)  # 기본: 첫 번째 시트
head(inspector)
head(mercury)


assembly_aun %>% 
  select(Batch,Sample,Sex,Haplotype,auN)


library(dplyr)
library(tidyr)
library(ggplot2)
#assembly_aun$auN
head(plot_df)
plot_df %>% filter(Batch == "T2T") %>% count(Sex)
plot_df <- assembly_aun %>%
  select(Batch, Sample, Sex, Haplotype, auN) %>%
  mutate(auN_Mb = auN / 1e6) %>%
  select(-auN, auN = auN_Mb) %>%
  pivot_wider(
    names_from  = Haplotype,
    values_from = auN,
    names_prefix = "hap"
  ) %>%
  rename(hap1 = hap1, hap2 = hap2) %>%
  mutate(
    Batch = factor(Batch),
    Sex   = factor(Sex),
    is_T2T = Batch %in% c("T2T", "T2T_level", "T2T-level")  # 필요시 Batch 값에 맞게 조정
  )

mx <- median(plot_df$hap1, na.rm = TRUE)
my <- median(plot_df$hap2, na.rm = TRUE)

ggplot(plot_df, aes(x = hap1, y = hap2)) +
  
  ## median lines
  geom_vline(xintercept = mx, linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = my, linetype = "dashed", linewidth = 0.6) +
  ## points: other (background)
  geom_point(data = ~ filter(.x, !is_T2T),aes(shape = Sex),color = "grey40",alpha = 0.7,size = 2) +
  ## points: T2T (highlight)
  geom_point( data = ~ filter(.x, is_T2T),aes(shape = Sex),color = "red",alpha = 0.90,size = 2.6) +
  
  scale_shape_manual(
    values = c("F" = 15, "M" = 4),   # 15 = square, 4 = X
    labels = c("Female", "Male")
  ) + 
  labs(
    x = "Hap1 auN (Mb)",
    y = "Hap2 auN (Mb)",
    shape = "Sex"
  ) +
  lims(x = c(20, 150), y = c(20, 150)) +
  theme_step1() +
  theme(legend.title = element_blank(),
        legend.position = c(0, 1),
        legend.justification = c("left", "top")) +
  guides(shape = guide_legend(override.aes = list(color = "black"))) -> a



head(inspector)
dim(inspector)
tail(inspector)
head(mercury)
dim(mercury)
plot_df
inspector %>% select(Sample,QV,Haplotype,Batch,`Total length`) %>%
  group_by(Sample) %>%
  summarise(QV_weighted = weighted.mean(QV, w = `Total length`, na.rm = TRUE),
    total_length = sum(`Total length`, na.rm = TRUE),
    n_hap = sum(!is.na(QV) & !is.na(`Total length`)),
    .groups = "drop") %>% rename(inspector_QV = QV_weighted) %>% select(Sample,inspector_QV) %>%
  left_join(plot_df %>% select(Sample,Sex)) %>%
  left_join(mercury %>% filter(Haplotype == "Both") %>% select(Sample,QV,Batch) %>% rename(mercury_QV = QV)) %>%
  mutate(Batch = factor(Batch),Sex   = factor(Sex),is_T2T = Batch %in% c("T2T", "T2T_level", "T2T-level"))-> df_qv

df_qv

mx <- median(df_qv$inspector_QV, na.rm = TRUE)
my <- median(df_qv$mercury_QV, na.rm = TRUE)

ggplot(df_qv, aes(x = inspector_QV, y = mercury_QV)) +
  
  ## median lines
  geom_vline(xintercept = mx, linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = my, linetype = "dashed", linewidth = 0.6) +
  
  ## non-T2T (background)
  geom_point(
    data = ~ filter(.x, !is_T2T),
    aes(shape = Sex),
    color = "grey40",
    alpha = 0.7,
    size = 2
  ) +
  
  ## T2T (highlight)
  geom_point(
    data = ~ filter(.x, is_T2T),
    aes(shape = Sex),
    color = "red",
    alpha = 0.9,
    size = 2.6
  ) +
  
  ## Sex legend only
  scale_shape_manual(
    values = c("F" = 15, "M" = 4),   # 15 = square, 4 = X
    labels = c("Female", "Male")
  ) + 
  
  labs(
    x = "Estimated QV (variant-based)",
    y = "Estimated QV (k-mer-based)",
    shape = "Sex"
  ) +
  
  ## legend icons all black
  guides(
    shape = guide_legend(
      override.aes = list(color = "black")
    )
  ) +
  lims(x = c(49, 60), y = c(49, 60)) +
  theme_step1() +
  theme(legend.position = 'none') -> b

a+b
cowplot::plot_grid(a,b)
  

###
url <- "https://docs.google.com/spreadsheets/d/1Y8cqzCizaQ9vvAzoAEhqqGW7GaplxRRnmPsQ0dK6o30/edit?usp=sharing"
drive_auth()

props <- sheet_properties(url)
props$name
#13,"1.3.asm_correcteness1_inspector"
#14,"1.3.asm_correcteness1_merqury"   
#11 "1.2.assembly_contiguity" 

novel_seq <- read_sheet(url,sheet = "1.4.novel_sequence_coverage",skip = 1)  
novel_seq
novel_seq_chm13 <- novel_seq %>% select(1:10)
novel_seq_GRCh38 <- novel_seq %>% select(11:20)

colnames(novel_seq_chm13) <- str_split_fixed(colnames(novel_seq_chm13),"\\.",2)[,1]
colnames(novel_seq_GRCh38) <- str_split_fixed(colnames(novel_seq_GRCh38),"\\.",2)[,1] 

novel_seq_chm13$cohort <- "CHM13"
novel_seq_GRCh38$cohort <- "GRCh38"

novel_seq_GRCh38 

novel_seq_chm13 %>% rbind(novel_seq_GRCh38) %>% #head()
  mutate(fraction = ((Total_length-Unaligned_length)/Total_length)*100) %>%
  select(Assembly,fraction,Unaligned_length,cohort) %>%
  ggplot(aes(x=Unaligned_length,y=fraction,color=cohort)) + 
  geom_point()


