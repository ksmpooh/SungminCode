
#KPPD assembly plot
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/")

theme_step1 <- function(base_size = 11, base_family = "",
                        base_line_size = base_size / 22,
                        base_rect_size = base_size / 22) {
  theme(title = element_text(family = 'Arial', size = 18, color = 'black'), text = element_text(family = 'Arial', size = 16, color = 'black'),
        axis.title = element_text(family = 'Arial', size = 18, color = 'black'), axis.text = element_text(family = 'Arial', size = 16, color = 'black'), 
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = NA), axis.line = element_line(colour = "black", size = rel(1)),
        #legend.background = element_rect(color = 'black'), 
        legend.title = element_text(family = 'Arial', size = 16),
        legend.text = element_text(family = 'Arial', size = 14),
        legend.direction = "vertical", 
        legend.box = c("horizontal", "vertical"),
        legend.spacing.x = unit(0.1, 'cm'),
        plot.margin = unit(c(0.25, 1, 1, 0.5), 'cm'),
        axis.title.y = element_text(margin = margin(r = 10, unit = "pt"))) }

library(tidyverse)
library(ggbeeswarm)
library(ggplot2)
#library(plyr)
#library(dplyr)


kppd <- read_table("draft_korean_pangenome_collection/FIgure1B_Assembly_Contiguity_2/kppd.merge.fai.txt",col_names = F)
hprc <- read.csv("/Users/ksmpooh/Desktop/KCDC/pangenome/construction/03.hifiasm/merge.HPRC.phaseI.NGx.length.csv", header = F) %>% filter(!(V1 %in% c("CHM13Y.fa","human_GRCh38_no_alt_analysis_set.fasta")))
#hprc <- read.csv("/Users/ksmpooh/Desktop/KCDC/pangenome/construction/03.hifiasm/merge.HPRC.phaseI.NGx.length.csv", header = F) %>% filter(!(V1 %in% c("CHM13Y.fa")))

#chm <- read_table("draft_korean_pangenome_collection/FIgure1B_Assembly_Contiguity_2/chm13v2.0_maskedY_rCRS.fa.fai",col_names = F)     # <- 너의 첫 번째 .fai 경로
#grch <- read_table("draft_korean_pangenome_collection/FIgure1B_Assembly_Contiguity_2/GCA_000001405.15_GRCh38_no_alt_analysis_set.cleaned.fna.fai",col_names = F)     # <- 너의 두 번째 .fai 경로
#grch <- read.csv("/Users/ksmpooh/Desktop/KCDC/pangenome/construction/03.hifiasm/lnegth_grch38_plamary_scaf.csv", header = F)
grch <- read_table("Figure1/GRCh38_info.txt",col_names = F)
chm <- read_table("Figure1/chm13_info.txt",col_names = F)
genome_length = 3000000000
head(hprc)
hprc %>% filter(V1 == "HG002.1.fa") %>% tail()
head(grch)
head(chm)
head(kppd)
head(grch)
colnames(kppd) <- c("name","len","offset","line_bases","line_width")
#colnames(chm) <- c("name","len","offset","line_bases","line_width")
#colnames(grch) <- c("name","len","offset","line_bases","line_width")

colnames(chm) <- c("name","len")
colnames(grch) <- c("name","len")

head(hprc)

grch
colnames(hprc) <- c("name","len","assembly","Coverage","cum_len")
#colnames(grch) <- c("name","len","assembly","Coverage","cum_len")
hprc %>% filter(str_detect(name,"GRC"))
hprc %>% filter(str_detect(name,"CHM"))

head(chm)

hprc %>% mutate(cohort = ifelse(str_detect(name,"GRCh38"),"GRCh38","HPRC")) %>% count(cohort)

chm %>% mutate(cohort = "CHM13v2",name = "CHM13v2") %>% 
  rbind(grch %>% mutate(cohort = "GRCh38",name="GRCh38")) %>% #head()
  arrange(cohort,desc(len)) %>% 
  group_by(cohort) %>%
  dplyr::mutate(cum_len = cumsum(len)) %>% #tail
  mutate(start = cum_len - len) %>% #filter(name == "KPPD132.1")
  mutate(Coverage = start /genome_length ) %>% rbind(hprc %>% mutate(cohort = "HPRC")) %>%
  select(cohort,name,len,cum_len,start,Coverage) -> ref_contig
  
  




ref_contig
  
ref_contig
kppd %>% mutate(name = paste0(str_split_fixed(name,"#",3)[,1],".",paste0(str_split_fixed(name,"#",3)[,2]))) %>%
  arrange(name, desc(len)) %>% #head()
  group_by(name) %>%
  dplyr::mutate(cum_len = cumsum(len)) %>% #tail
  mutate(start = cum_len - len) %>% #filter(name == "KPPD132.1")
  mutate(Coverage = start /genome_length ) %>% #head()
  select(name,len,Coverage,start,cum_len) %>% mutate(cohort = "KPPD") -> kppd_contig

kppd_contig
df
kppd_contig %>% mutate(cohort = ifelse(name %in% c("KPPD132.1","KPPD131.1","KPPD130.1","KPPD129.1","KPPD132.2","KPPD131.2","KPPD130.2","KPPD129.2"),"KPPD (T2T level)","KPPD")) %>%
  rbind(ref_contig) %>% ungroup() %>% #dplyr::count(cohort)
  mutate(Coverage = Coverage*3) -> df

df %>%
  ggplot(aes(x=Coverage, y=len/1000000, color=cohort,group=name)) +
  geom_vline(xintercept = 1.5, linetype="dotted", linewidth=0.5) +
  geom_step(alpha=0.8) +
  #scale_color_manual(values = colors) +
  #scale_x_continuous(limits = c(0, 1.25), breaks = seq(0, 1.25, by = 0.25)) +
  labs(x = "Cumulative coverage (Gb)", y = "Length (Mb)") +
  theme_bw()

head(df)

table(df$cohort)

head(df)


write.table(df %>% select(name,len,cohort),"Figure1/Figure1.data.cummulativeCoverage_KPPD_HPRC_CHM13_GRCH38nsplit.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#head(df)

## v1
ggplot() +
  # Plot HPRC and KoGES first with lower alpha
  geom_step(
    data = df %>% filter(cohort %in% c("HPRC", "KPPD")),
    aes(x=Coverage, y=len/1000000, color=cohort,group=name),
    size = 0.5, alpha = 0.2
  ) +
  geom_step(
    data = df %>% filter(cohort == "KPPD (T2T level)"),
    aes(x=Coverage, y=len/1000000, color=cohort,group=name),
    size = 0.5, alpha = 0.8
  ) +
  
  # Plot GRCh38 and CHM13 with higher alpha and black color, no legend
  geom_step(
    data = df %>% filter(cohort %in% c("GRCh38", "CHM13v2")),
    aes(x = Coverage, y = len / 1000000, group = name),
    color = "black", size = 0.5, alpha = 1, show.legend = FALSE
  ) +
  geom_vline(xintercept = 1.5, linetype = "dotted", linewidth = 0.5) +
  scale_color_manual(values = c("HPRC" = "red", "KPPD" = "blue","KPPD (T2T level)" = "green")) +
  labs(x = "Cumulative coverage (Gb)", y = "Length (Mb)") +
  theme_step1() + 
  theme(
    legend.position = c(0.87, 0.8),
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  # Annotations for CHM13 and GRCh38
  annotate("text", x = 0.8, y = 200, label = "CHM13v2", color = "black", size = 5, hjust = 0) +
  annotate("text", x = 0.8, y = 100, label = "GRCh38", color = "black", size = 5, hjust = 0) -> p
p

### v2

head(df_raw)
df_raw <- read.table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/Figure1/Figure1.data.cummulativeCoverage_KPPD_HPRC_CHM13_GRCH38nsplit.txt",sep = "\t",header = T)
df_raw %>% filter(name %in% c("KPPD132.1","KPPD131.1","KPPD130.1","KPPD129.1","KPPD132.2","KPPD131.2","KPPD130.2","KPPD129.2"))
table(df_raw$cohort)

genome_length = 3000000000
df_raw %>%
  arrange(name, desc(len)) %>% #head()
  group_by(name) %>%
  dplyr::mutate(cum_len = cumsum(len)) %>% #tail
  mutate(start = cum_len - len) %>% #filter(name == "KPPD132.1")
  mutate(Coverage = start /genome_length ) %>% #head()
  select(name,len,Coverage,start,cum_len,cohort) -> df

#df_v2 %>% filter(name=="CHM13v2")
df

head(df)

ggplot() +
  # Plot HPRC and KoGES first with lower alpha
  geom_step(
    data = df %>% filter(cohort %in% c("HPRC", "KPPD")),
    aes(x=Coverage*3, y=len/1000000, color=cohort,group=name),
    size = 0.5, alpha = 0.2
  ) +
  geom_step(
    data = df %>% filter(cohort == "KPPD (T2T level)"),
    aes(x=Coverage*3, y=len/1000000, color=cohort,group=name),
    size = 0.5, alpha = 0.8
  ) +
  
  # Plot GRCh38 and CHM13 with higher alpha and black color, no legend
  geom_step(
    data = df %>% filter(cohort %in% c("GRCh38", "CHM13v2")),
    aes(x = Coverage*3, y = len / 1000000, group = name),
    color = "black", size = 0.5, alpha = 1, show.legend = FALSE
  ) +
  geom_vline(xintercept = 1.5, linetype = "dotted", linewidth = 0.5) +
  scale_color_manual(values = c("HPRC" = "red", "KPPD" = "blue","KPPD (T2T level)" = "green")) +
  labs(x = "Cumulative coverage (Gb)", y = "Length (Mb)") +
  theme_step1() + 
  theme(
    legend.position = c(0.87, 0.8),
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  # Annotations for CHM13 and GRCh38
  annotate("text", x = 0.8, y = 200, label = "CHM13v2", color = "black", size = 5, hjust = 0) +
  annotate("text", x = 0.8, y = 100, label = "GRCh38", color = "black", size = 5, hjust = 0) -> p
p


ggsave("Figure1/Figure1b.cummulativeCoverage_GRChNsplit.v2.pdf",p,width = 6,height = 5,  useDingbats = FALSE)

library(Cairo)
CairoPDF("Figure1b.cummulativeCoverage.v1.pdf", width = 6, height = 5)
print(p)
dev.off()

head(hprc)


library(dplyr)
head(df)
## 1) 샘플(name)별 N50/L50 계산 -------------------------------------------
n50_by_sample <- df %>%
  group_by(name) %>%
  summarise(
    total_bp  = sum(len, na.rm = TRUE),
    n_contigs = n(),
    N50_bp = {
      # 길이 내림차순에서 누적합이 총 길이의 50% 이상이 되는 첫 contig의 길이
      ord <- order(-len)                 # 내림차순 인덱스
      l   <- len[ord]
      cs  <- cumsum(l)
      l[which(cs >= 0.5 * sum(l))[1]]
    },
    L50   = {
      ord <- order(-len)
      l   <- len[ord]
      cs  <- cumsum(l)
      which(cs >= 0.5 * sum(l))[1]
    },
    .groups = "drop"
  )

n50_by_sample

n50_by_sample %>% filter(!(name %in% c("CHM13v2","GRCh38"))) %>%
  mutate(cohort = ifelse(str_detect(name,"KPPD"),"KPPD","HPRC")) %>% group_by(cohort) %>% 
  summarise(mean(n_contigs),mean(N50_bp),mean(L50))

n50_summary_by_cohort <- df %>%
  group_by(cohort, name) %>%
  summarise(
    N50_bp = { l <- sort(len, decreasing = TRUE); l[cumsum(l) >= 0.5*sum(l)][1] },
    .groups = "drop_last"
  ) %>%
  summarise(
    n_samples = n(),
    N50_median_bp = median(N50_bp),
    N50_IQR_bp    = IQR(N50_bp),
    N50_mean_bp   = mean(N50_bp),
    .groups = "drop"
  )

n50_summary_by_cohort

#####
# ====== 설정 ======
##


chm %>% mutate(cohort = "CHM13v2",name = "CHM13v2") %>% #dim()
  #rbind(grch %>% mutate(cohort = "GRCh38",name="GRCh38")) %>% #head()
  arrange(cohort,desc(len)) %>% 
  group_by(cohort) %>%
  dplyr::mutate(cum_len = cumsum(len)) %>% #tail
  mutate(start = cum_len - len) %>% #filter(name == "KPPD132.1")
  mutate(Coverage = start /genome_length ) %>% rbind(hprc %>% mutate(cohort = ifelse(str_detect(name,"GRCh38"),"GRCh38","HPRC"))) %>%
  rbind(grch %>% mutate(cohort = "GRCh38",name="GRCh38")) %>% #head()
  select(cohort,name,len,cum_len,start,Coverage) -> ref_contig


####



df <- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/draft_korean_pangenome_collection/FIgure1B_Assembly_Contiguity_1/Supplement_T2_Contiguity_gfastats.tab")
ref <- read_table("~/Desktop/KCDC/pangenome/00.datacheck/KPP_draft_assembly_matchingID.with.Gender.20250224.txt")
head(ref)
head(df)

df %>% left_join(ref %>% rename(Sample = KPPD_ID) %>% select(Sample,sex)) %>% 
  mutate(Batch = ifelse(Batch == "KPP","KPPD (T2T level)","KPPD")) %>% 
  select(Sample,Batch,Total_contig_length,sex) -> df

df %>% select(Sample,Batch,Total_contig_length,sex)
table(df$sex)/2




df2 <- df %>%
  mutate(total_gb = Total_contig_length / 1e9)

#write.table(df2,"Figure1/Figure1.data.TotalContigLength_bySex.txt",col.names = T,row.names = F,quote = F,sep = "\t")

# T2T level만 추출
df_t2t <- df2 %>% filter(Batch == "KPPD (T2T level)")

# ---- Plot ----


# sex → label 매핑

#ggsave("Figure1/Figure1b.TotalContigLength_bySex.pdf",p,width = 1,height = 5,  useDingbats = FALSE)


library(Cairo)
CairoPDF("Figure1/Figure1b.TotalContigLength_bySex.pdf", width = 4, height = 5)
print(p)
dev.off()


df2 <- read.table("Figure1/Figure1.data.TotalContigLength_bySex.txt",header = T,sep = "\t")


df_t2t <- df2 %>% filter(Batch == "KPPD (T2T level)")

sex_counts <- df2 %>% count(sex) %>% mutate(n = n/2)
sex_labels <- setNames(
  paste0(ifelse(sex_counts$sex == "F", "Female", "Male"),
         "\n(n=", sex_counts$n, ")"),
  sex_counts$sex
)

# T2T level 샘플별로 서로 다른 shape 부여 (최대 26개 정도까지 깔끔)
t2t_samples <- unique(df_t2t$Sample)
shape_pool  <- c(3,4,5,8)
shape_vals  <- setNames(shape_pool[seq_along(t2t_samples)], t2t_samples)

# ---- 시각화 ----
p <- ggplot(df2, aes(x = sex, y = total_gb, fill = sex)) +
  
  # 1️⃣ 더 투명한 바이올린
  geom_violin(trim = FALSE, alpha = 0.4, width = 0.9,
              color = "gray70", linewidth = 0.3) +
  # 2️⃣ 박스플롯 (IQR 강조)
  geom_boxplot(width = 0.05, outlier.shape = NA,
               color = "gray50", fill = NA, linewidth = 0.4) +
  
  # 3️⃣ 일반 샘플: 흐린 회색 점
  geom_quasirandom(
    data = df2 %>% filter(Batch != "KPPD (T2T level)"),
    aes(x = sex, y = total_gb),
    color = "grey", size = 2, alpha = 1.5, width = 0.15,shape=1
  ) +
  
  # 4️⃣ T2T 샘플: 형광 연두 포인트 (제일 강조)
  geom_point(
    data = df_t2t,
    aes(x = sex, y = total_gb, shape = Sample),
    size = 2, stroke = 1.1,
    color = "red",  # 💡 형광 연두색
    position = position_jitter(width = 0.05, height = 0)
  ) +
  
  scale_x_discrete(labels = sex_labels) +
  scale_fill_manual(values = c(F = "#F6B7A4", M = "#89A8C5")) +  # 은은한 배경색
  scale_shape_manual(values = shape_vals) +
  
  labs(
    x = NULL,
    y = "Total contig length (Gb)",
    subtitle = "* T2T level: Red"
  ) +
  
  theme_step1() +
  theme(
    legend.position = "none",
    #axis.title.y = element_text(margin = margin(r = 10), color = "black"),
    axis.text = element_text(color = "black"),
    plot.subtitle = element_text(hjust = 1, face = "italic", color = "gray30")
  )

p



###3####
head(df)

df <- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/draft_korean_pangenome_collection/FIgure1B_Assembly_Contiguity_1/Supplement_T2_Contiguity_gfastats.tab")
ref <- read_table("~/Desktop/KCDC/pangenome/00.datacheck/KPP_draft_assembly_matchingID.with.Gender.20250224.txt")
head(ref)
head(df)

df %>% left_join(ref %>% rename(Sample = KPPD_ID) %>% select(Sample,sex)) %>% 
  mutate(Batch = ifelse(Batch == "KPP","KPPD (T2T level)","KPPD")) -> df


df %>% select(Sample,Haplotype,Batch,Total_contig_length,Contigs,Contig_N50,Contig_L50,Largest_contig,Smallest_contig) -> df
head(df)
#write.table(df,"Figure1/sup_Figure1.data.BasicAssemblyQualityMetrics.txt",col.names = T,row.names = F,quote = F,sep = "\t")

library(dplyr)
library(tidyr)
library(ggplot2)

# ---- 데이터 변환 ----
df_long <- df %>%
  select(Sample, Haplotype, Batch,
         Total_contig_length, Contigs, Contig_N50,
         Contig_L50, Largest_contig, Smallest_contig) %>%
  pivot_longer(
    cols = c(Total_contig_length, Contigs, Contig_N50,
             Contig_L50, Largest_contig, Smallest_contig),
    names_to = "Metric",
    values_to = "Value"
  )

# ---- facet plot ----
ggplot(df_long, aes(x = factor(Haplotype), y = Value, fill = factor(Haplotype))) +
  geom_violin(trim = FALSE, alpha = 0.6, width = 0.8, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", alpha = 0.4) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("1" = "#E07A5F", "2" = "#3D5A80")) +
  labs(
    x = "Haplotype",
    y = NULL,
    title = "Assembly quality metrics by haplotype",
    subtitle = "Violin + boxplots for each assembly metric"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )




df_long <- df %>%
  select(Sample, Batch,
         Total_contig_length, Contigs, Contig_N50,
         Contig_L50, Largest_contig, Smallest_contig) %>%
  pivot_longer(
    cols = c(Total_contig_length, Contigs, Contig_N50,
             Contig_L50, Largest_contig, Smallest_contig),
    names_to = "Metric", values_to = "Value"
  )

# Metric 순서(원하는 순서가 있으면 조정)
df_long <- df_long %>%
  mutate(Metric = factor(
    Metric,
    levels = c("Total_contig_length","Contigs","Contig_N50",
               "Contig_L50","Largest_contig","Smallest_contig")
  ))

# 플롯: x=Metric(하나의 축), 바이올린 + 점(배치별 모양만 다르게)
ggplot(df_long, aes(x = Metric, y = Value)) +
  geom_violin(fill = "#d9dde6", color = NA, width = 0.9, alpha = 0.5, trim = FALSE) +
  geom_jitter(
    aes(shape = Batch),
    width = 0.15, height = 0,
    size = 1.6, stroke = 1.0, color = "black", alpha = 0.7
  ) +
  scale_shape_manual(values = c("KPPD" = 21, "KPPD (T2T level)" = 4)) +  # ○ vs ×
  labs(x = NULL, y = NULL, title = "Assembly quality metrics") +
  facet_wrap(~ Metric, scales = "free", ncol = 3) + 
  theme_step1() +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_blank()
    #axis.text.x = element_text(angle = 15, hjust = 1)
  ) -> p
  
CairoPDF("Figure1/F1.sup.Assemlby.Quality.metrics.pdf", width = 20, height = 10)
print(p)
dev.off()

### busco

busco <- read.table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/draft_korean_pangenome_collection/Supplement_T3.Completeness_busco.tab",header = T)
head(busco)
table(busco$Batch)
busco %>%
  mutate(Batch = ifelse(Batch == "T2T","KPPD (T2T level)",ifelse(Batch == "Non.T2T","KPPD",Batch))) %>%
  mutate(cohort = ifelse(str_detect(Batch,"KPPD"),"KPPD","Ref")) -> busco

\library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# 1) long-format + 비율(%) 계산
busco_long <- busco %>%
  select(File, Batch, cohort,
         complete_single_copy, complete_duplicated, fragmented, missing) %>%
  pivot_longer(
    c(complete_single_copy, complete_duplicated, fragmented, missing),
    names_to = "Category", values_to = "Count"
  ) %>%
  mutate(
    Fraction = Count / 11834 * 100,
    Category = factor(Category,
                      levels = c("missing","fragmented",
                                 "complete_duplicated","complete_single_copy")),
    Batch = factor(Batch,
                   levels = c("KPPD", "KPPD (T2T level)", "CHM13v2", "GRCh38"))
  )

# 2) 상단 패널 (90~100%)
p_top <- ggplot(busco_long, aes(x = File, y = Fraction, fill = Category)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  facet_grid(~ Batch, scales = "free_x", space = "free_x") +
  coord_cartesian(ylim = c(90, 100), expand = FALSE) +
  scale_y_continuous(
    breaks = seq(90, 100, by = 10),
    labels = function(x) paste0(x),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_brewer(
    palette = "Set2",
    breaks = c("missing","fragmented","complete_duplicated","complete_single_copy"),
    labels = c("Missing","Fragmented","Complete (Duplicated)","Complete (Single-copy)")
  ) +
  labs(y = " ", x = NULL) +
  theme_step1() +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x.bottom = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  guides(fill = guide_legend(reverse = TRUE))

# 3) 하단 패널 (0~10%)
p_bottom <- ggplot(busco_long, aes(x = File, y = Fraction, fill = Category)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  facet_grid(~ Batch, scales = "free_x", space = "free_x") +
  coord_cartesian(ylim = c(0, 10), expand = FALSE) +
  scale_y_continuous(
    breaks = seq(0, 10, by = 2),
    labels = function(x) paste0(x),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_brewer(
    palette = "Set2",
    breaks = c("missing","fragmented","complete_duplicated","complete_single_copy"),
    labels = c("Missing","Fragmented","Complete (Duplicated)","Complete (Single-copy)")
  ) +
  labs(y = "Percentage (%)", x = "") +
  theme_step1() +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

# 4) 결합 (상:하 = 1:3)
p <- (p_top / p_bottom) +
  plot_layout(heights = c(1, 3), guides = "collect")

p


head(busco_long)

busco_long %>% group_by(Batch,Category) %>%
  summarise(mean = mean(Fraction)) %>%
  pivot_wider(names_from = Category,values_from = mean)


### PCA
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/")

df <- read.table("draft_korean_pangenome_collection/PCA/kor132_hprc47_cpc58.chm13v2.snpgdsPCA.tsv",header = T)
df <- read.table("draft_korean_pangenome_collection/PCA/KPPD132.kor_hprc.chm13v2.snpgdsPCA.tsv",header = T)
head(df)
df %>% count(Population,Graph)

head(df)

shape_map <- c(
  "HPRCv1" = 17,  # triangle
  "KPPD"   = 19   # circle
  
)

ggplot(df, aes(x = EV1, y = EV2)) +
  geom_point(
    aes(color = Population, shape = Graph),
    size = 2.8, stroke = 0.7, alpha = 0.8
  ) +
  scale_shape_manual(values = shape_map) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = "EV1",
    y = "EV2",
    shape = "Project",
    color = "Population"
  ) +
  theme_step1() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

p



df <- read.table("draft_korean_pangenome_collection/PCA/kor132_hprc47_cpc58.chm13v2.snpgdsPCA.tsv",header = T)
df
table(df$Graph)

shape_map <- c(
  "HPRCv1" = 17,  # triangle
  "CPCv1" = 18,
  "KPPD"   = 19   # circle
  
)

ggplot(df, aes(x = EV1, y = EV2)) +
  geom_point(
    aes(color = Population, shape = Graph),
    size = 2.8, stroke = 0.7, alpha = 0.7
  ) +
  scale_shape_manual(values = shape_map) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = "EV1",
    y = "EV2",
    shape = "Project",
    color = "Population"
  ) +
  theme_step1() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )


### indspector
url <- "https://docs.google.com/spreadsheets/d/1Y8cqzCizaQ9vvAzoAEhqqGW7GaplxRRnmPsQ0dK6o30/edit?usp=sharing"
drive_auth()

props <- sheet_properties(url)
props$name
#13,"1.3.asm_correcteness1_inspector"
#14,"1.3.asm_correcteness1_merqury"   
#11 "1.2.assembly_contiguity" 

inspector <- read_sheet(url,sheet = "1.3.asm_correcteness1_inspector")  # 기본: 첫 번째 시트
inspector %>% select(File,`Base substitution`,`Small-scale collapse`,`Small-scale expansion`) %>%
  pivot_longer(2:4) %>%
  ggplot(aes(x=value,fill=name)) + 
  geom_density(alpha = 0.5) +
  labs(
    x = "Number of small-scale assembly errors",
    y = "Density",
    fill = NULL
  ) +
  theme_step1() +
  theme(
    legend.direction = "horizontal",
    legend.position = "top"   # 위에 legend
  )
