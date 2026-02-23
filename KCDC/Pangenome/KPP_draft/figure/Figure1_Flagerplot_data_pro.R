### flagger
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/flagger/compare_pipeline")

library(ggpattern)

h1 <- read.table("h1_hmm_flagger_result/prediction_summary_final.tsv")
h2 <- read.table("h2_hmm_flagger_result/prediction_summary_final.tsv")



h1_a <- read.table("h1_all_hmm_flagger_result/prediction_summary_final.tsv")
h2_a <- read.table("h2_all_hmm_flagger_result/prediction_summary_final.tsv")

head(h1)
head(h1_a)


h1 <- read.table("h1_hmm_flagger_result/transition_final.tsv")
h2 <- read.table("h2_hmm_flagger_result/transition_final.tsv")

h1_a <- read.table("h1_all_hmm_flagger_result/transition_final.tsv")
h2_a <- read.table("h2_all_hmm_flagger_result/transition_final.tsv")

head(h1)
head(h2)

dim(h1)
dim(h1_a)


h1 <- read.table("h1_hmm_flagger_result/final_flagger_prediction.bed",skip = 1)
h2 <- read.table("h2_hmm_flagger_result/final_flagger_prediction.bed",skip = 1)

h1_a <- read.table("h1_all_hmm_flagger_result/final_flagger_prediction.bed",skip = 1)
h2_a <- read.table("h2_all_hmm_flagger_result/final_flagger_prediction.bed",skip = 1)

head(h1)
head(h1_a)

dim(h1)
dim(h1_a)
head(h1)
h1 %>% count(V4) %>% mutate(prop = n / sum(n)) %>% mutate(type = "h1")
h1_a %>% count(V4) %>% mutate(prop = n / sum(n)) %>% mutate(type = "h1_a")

h2 %>% count(V4) %>% mutate(prop = n / sum(n)) %>% mutate(type = "h2")
h2_a %>% count(V4) %>% mutate(prop = n / sum(n)) %>% mutate(type = "h2_a")

##
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/flagger/final_result")
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/flagger/03.result/")
library(gridExtra)

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

library(ggplot2)
library(cowplot)  # 패치워크와 비슷하게 레이아웃 조정 가능




# 병합 대상 파일 경로 불러오기
bed_files <- list.files(
  path = ".",
  pattern = "final_flagger_prediction\\.bed$",
  recursive = TRUE,
  full.names = TRUE
)

# 파일이 제대로 인식됐는지 확인
print(bed_files)
a
head(a)
a$dir[,2]
bed_files %>% as.data.frame() %>% rename(path = ".")  %>% 
  mutate(dir = str_split_fixed(path,"/",4)[,2]) %>% #head()
  mutate(ID = str_split_fixed(dir,"_",2)[,1]) %>% count(ID)

merged_bed <- lapply(bed_files, function(f) {
  # 파일 읽기
  df <- read_tsv(f, col_names = FALSE,skip = 1)
  
  # ID 추출: 예) final_result/KPPD080_annotation/final_results/final_flagger_prediction.bed
  id <- sub("_annotation.*", "", basename(dirname(dirname(f))))
  
  # ID 열 추가
  df$ID <- id
  
  return(df)
}) %>% bind_rows()
head(merged_bed)
merged_bed %>% count((ID))

#merged_bed %>% select(ID,X1,X2,X3,X4) %>% write.table("~/Desktop/KCDC/pangenome/KPPD/flagger/final.prediction.merge.txt",col.names = F, row.names = F,quote = F,sep = "\t")


merged_bed %>% mutate(Region = X3 - X2) %>% mutate(haplotype = str_split_fixed(X1,'#',3)[,2]) %>% 
  group_by(ID,haplotype,X4) %>% summarise(sum = sum(Region)) %>% mutate(prop = prop.table(sum)) %>%
  mutate(type = ifelse(str_detect(ID,"out"),"normal","annotation_CHM13")) %>%
  mutate(hap = paste0(ID,"_",haplotype)) -> bed_prop

head(bed_prop)

head(bed_prop)

kppd <- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/draft_korean_pangenome_collection/FIgure1B_Assembly_Contiguity_2/kppd.merge.fai.txt",col_names = F)
colnames(kppd) <- c("name","len","offset","line_bases","line_width")
head(kppd)

kppd %>% mutate(hap = paste0(str_split_fixed(name,"#",3)[,1],"_",paste0(str_split_fixed(name,"#",3)[,2]))) %>% #head()
  mutate(check_100kb = ifelse(len < 100000,"pass","fail")) %>%
  count(hap,check_100kb)

kppd %>% mutate(hap = paste0(str_split_fixed(name,"#",3)[,1],"_",paste0(str_split_fixed(name,"#",3)[,2]))) %>% #head()
  group_by(hap) %>% summarise(length_fai = sum(len)) %>% 
  left_join(bed_prop %>%group_by(hap) %>% summarise(length = sum(sum)) %>% filter(!str_detect(hap,"flagger"))) %>%
  count(length_fai == length)


kppd %>% mutate(hap = paste0(str_split_fixed(name,"#",3)[,1],"_",paste0(str_split_fixed(name,"#",3)[,2]))) %>% #head()
  group_by(hap) %>% summarise(length_fai = sum(len)) %>% 
  left_join(bed_prop %>%group_by(hap) %>% summarise(length = sum(sum)) %>% 
              filter(str_detect(hap,"flagger")) %>% mutate(hap =str_replace_all(hap, "_flagger_out",""))) %>%
  count(length_fai == length)


bed_prop %>%group_by(hap) %>% summarise(length = sum(sum)) %>% 
  filter(str_detect(hap,"flagger")) %>% mutate(hap =str_replace_all(hap, "_flagger_out",""))
bed_prop %>%group_by(hap) %>% summarise(length = sum(sum)) %>% filter(!str_detect(hap,"flagger"))



bed_prop
bed_prop %>% 
  #filter(type == "annotation_CHM13") %>% 
  ggplot(aes(x = hap, y = sum,fill = X4)) +
  geom_bar(stat = "identity") +
  theme_step1() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "Sample ID", y = "Proportion", fill = "Flagger Region") + 
  scale_fill_manual(values = c(
    "Hap" = "#66c2a5",
    "Col" = "#fc8d62",
    "Dup" = "#8da0cb",
    "Err" = "#e78ac3"
  )) + coord_cartesian(ylim = c(0, 300000))

#df_all
head(bed_prop)
bed_prop %>% filter(type == "annotation_CHM13") %>% select(ID,haplotype,X4,sum) %>% 
  rename(length = sum) %>% rename(region = X4) %>% #head()
  writexl::write_xlsx("~/Desktop/KCDC/pangenome/KPPD/figure/Figure1_DATA.flagger_region.xlsx")


######
library(tidyverse)
library(patchwork)

bed_prop <- readxl::read_xlsx("~/Desktop/KCDC/pangenome/KPPD/figure/Figure1_DATA.flagger_region.xlsx")
dim(bed_prop)
bed_prop


# 0) 색/라벨 고정
fill_colors <- c("Hap"="#66c2a5","Col"="#fc8d62","Dup"="#8da0cb","Err"="#e78ac3")
fill_labels <- c("Hap"="Haploid","Col"="Collapsed","Dup"="Duplicated","Err"="Erroneous")
head(df_all)
# 1) 공통 데이터 (팩터 레벨을 4종 모두로 고정)
df_all <- bed_prop %>%
  #filter(type == "annotation_CHM13") %>%
  mutate(
    hap = paste0(ID,"_",haplotype),
    sum_mb = length/1e6,
    region = factor(region, levels = c("Hap","Col","Dup","Err"))
  )

# x축 정렬(원하면 기준 변경 가능)
ord <- df_all %>% filter(region=="Hap") %>% arrange(desc(sum_mb)) %>% pull(hap)
df_all <- df_all %>% mutate(hap = factor(hap, levels = unique(c(ord, setdiff(hap, ord)))))
head(df_all)
# 2) 상단 = Haploid만
df_top <- df_all# %>% filter(X4 == "Hap")

p_top <- ggplot(df_top, aes(hap, sum_mb, fill = region)) +
  geom_col(width = 0.9, show.legend = FALSE) +
  coord_cartesian(ylim = c(2800, 3100), expand = FALSE) +
  scale_fill_manual(values = fill_colors, breaks = names(fill_labels),
                    labels = fill_labels, drop = FALSE) +  # 범례 유지용
  labs(y = "Total length (Mb)                       ") +
  theme_step1() +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.line.x =  element_blank())
#p_top
# 3) 하단 = Haploid 제외 (Col→Dup→Err 순으로 누적)
df_bottom <- df_all %>%
  #filter(X4 != "Hap") %>%
  mutate(region = forcats::fct_relevel(region,"Hap","Col","Dup","Err"))

p_bottom <- ggplot(df_bottom, aes(hap, sum_mb, fill = region)) +
  geom_col(width = 0.9, position = "stack") +
  coord_cartesian(ylim = c(0, 400), expand = FALSE) +
  scale_fill_manual(values = fill_colors, breaks = names(fill_labels),
                    labels = fill_labels, drop = FALSE) +
  labs(y = "",x="Assemblies") +
  theme_step1() +
  theme(#axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "right")

# 4) 결합
p <- (p_top / p_bottom) + plot_layout(heights = c(2, 1), guides = "collect")
p


head(bed_prop)
bed_prop %>% group_by(ID,haplotype) %>% mutate(prop = prop.table(length)) %>%
  group_by(region) %>% 
  summarise(mean(length),mean(prop)) %>% t()






setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/")
ggsave("Figure1X.Flagger_region.v1.pdf",p,width = 10,height = 5,  useDingbats = FALSE)

library(Cairo)
CairoPDF("Figure1b.cummulativeCoverage.v1.pdf", width = 6, height = 5)
print(p)
dev.off()

ggsave(
  filename = "Figure1X.Flagger_region.v1.pdf",
  plot     = p,
  width    = 12, height = 5,
  device   = cairo_pdf,
  #useDingbats = FALSE,
  limitsize  = FALSE
)



###

library(tidyverse)
library(ggbeeswarm)

# 팔레트/라벨 고정
fill_colors <- c("Hap"="#66c2a5","Col"="#fc8d62","Dup"="#8da0cb","Err"="#e78ac3")
fill_labels <- c("Hap"="Haploid","Col"="Collapsed","Dup"="Duplicated","Err"="Erroneous")

# 데이터 준비 (region 레벨 고정, cohort 모양 지정용)
plot_df <- bed_prop %>%
  mutate(
    cohort = ifelse(ID %in% c("KPPD132","KPPD131","KPPD130","KPPD129"),
                    "KPPD (T2T level)", "KPPD"),
    region = forcats::fct_relevel(region, "Hap","Col","Dup","Err")
    # length 단위 변환 필요하면: length = length/1e6
  )

p_q <- ggplot(
  plot_df,
  aes(x = region, y = length, fill = region, color = region, shape = cohort)
) +
  geom_quasirandom(
    groupOnX    = TRUE,
    dodge.width = 0.7,
    size        = 1.9,
    alpha       = 0.8,
    stroke      = 0.25
  ) +
  # y축만 자유 척도
  facet_grid(~ region, scales = "free") +
  # 색/라벨/범례 고정 (drop=FALSE로 모든 항목 유지)
  scale_fill_manual(values = fill_colors,
                    breaks = names(fill_labels),
                    labels = fill_labels,
                    name   = "Flagger region",
                    drop   = FALSE) +
  scale_color_manual(values = fill_colors, guide = "none", drop = FALSE) +
  scale_shape_manual(values = c("KPPD" = 16, "KPPD (T2T level)" = 17),
                     name   = "Cohort") +
  labs(x = NULL, y = "Length (Mb)") +
  theme_step1() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_blank(),   # facet 스트립이 라벨 역할
    axis.ticks.x       = element_blank(),
    legend.position    = "right",
    strip.text         = element_text(face = "bold")
  ) +
  guides(
    fill  = guide_legend(order = 1,
                         override.aes = list(shape = 22, size = 5, colour = NA)),
    shape = guide_legend(order = 2)
  )

p_q
ggsave(
  filename = "Figure1X.Flagger_region_distribution.v1.pdf",
  plot     = p_q,
  width    = 12, height = 5,
  device   = cairo_pdf,
  #useDingbats = FALSE,
  limitsize  = FALSE
)


### data check
head(merged_bed)

#merged_bed %>% select(X1,X2,X3,X4) %>% write.table("~/Desktop/KCDC/pangenome/KPPD/flagger/final.prediction.merge.txt",col.names = F, row.names = F,quote = F,sep = "\t")
merged_bed  <- read_table("~/Desktop/KCDC/pangenome/KPPD/flagger/final.prediction.merge.txt",col_names = F)
head(merged_bed)



merged_bed %>% mutate(Region = X4 - X3) %>% mutate(haplotype = str_split_fixed(X2,'#',3)[,2]) %>% #head()
  mutate(contig = str_split_fixed(X2,'#',3)[,3]) %>%
  mutate(type = ifelse(str_detect(X1,"out"),"normal","annotation_CHM13")) %>% count(type)
  filter(type == "annotation_CHM13") %>%
  group_by(X2) %>%
  mutate(contig_length = sum(Region)) %>% #head()
  group_by(X2,X5) %>% summarise(region_sum = sum(Region))



merged_bed %>% mutate(Region = X4 - X3) %>% mutate(haplotype = str_split_fixed(X2,'#',3)[,2]) %>% #head()
  mutate(contig = str_split_fixed(X2,'#',3)[,3]) %>%
  mutate(type = ifelse(str_detect(X1,"out"),"normal","annotation_CHM13")) %>%
  filter(type == "annotation_CHM13") %>%
  group_by(X2,X5) %>% summarise(region_sum = sum(Region)) %>% #head()
  group_by(X2) %>% mutate(contig_length = sum(region_sum)) %>%
  mutate(prop.byContig = prop.table(region_sum)) -> contig_region_check

contig_region_check
head(contig_region_check)
contig_region_check %>% select(X2, contig_length) %>% unique() 
#%>% summarise(quantile(contig_length))
contig_region_check %>% select(X2, contig_length) %>% unique() %>% #head()
  ggplot(aes(x="",y=contig_length)) + 
  geom_violin()


contig_region_check %>%
  distinct(X2, contig_length) %>%  # 중복 제거
  ggplot(aes(x = contig_length / 1e6)) +  # Mb 단위 변환
  geom_histogram(
    bins = 100, 
    fill = "#66c2a5", 
    color = "black", 
    alpha = 0.8
  ) +
  scale_x_continuous(
    trans = "log10",
    labels = scales::comma
  ) +
  labs(
    x = "Contig length (Mb, log10 scale)",
    y = "Count",
    title = "Distribution of contig lengths"
  ) +
  theme_minimal(base_size = 14)




contig_region_check %>% filter(X5 == "Dup") #%>% head()

library(tidyverse)
contig_region_check %>% filter(X5 == "Dup") %>%  
  mutate(length_group = ntile(contig_length, 4))
# contig_length를 기준으로 4분위수 구간 만들기
plot_df <- contig_region_check %>% filter(X5 == "Dup") %>% ungroup() %>%
  mutate(
    length_group = ntile(contig_length, 4),  # 1~4분위 구간
    length_group = factor(length_group,
                          levels = 1:4,
                          labels = c("Q1 (shortest)", "Q2", "Q3", "Q4 (longest)"))
  )

plot_df
plot_df %>% count(length_group)
# 각 구간별 Dup 비율 분포 그리기
ggplot(plot_df, aes(x = length_group, y = prop.byContig, fill = length_group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 1.5) +
#  geom_jitter(width = 0.15, alpha = 0.3, size = 1, color = "black") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Contig length quantile group",
    y = "Dup proportion (per contig)",
    title = "Distribution of Dup proportion by contig length quantile"
  ) +
  theme_step1() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )


head(plot_df)
ggplot(plot_df, aes(x = length_group, y = prop.byContig, fill = length_group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 1.5) +
  #  geom_jitter(width = 0.15, alpha = 0.3, size = 1, color = "black") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Contig length quantile group",
    y = "Dup proportion (per contig)",
    title = "Distribution of Dup proportion by contig length quantile"
  ) +
  theme_step1() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )



library(tidyverse)

# ──────────────────────────────
# ① 분위 구간 계산
# ──────────────────────────────
#contig_region_check %>% filter(X5 == "Dup") %>% ungroup() %>%
quantile_info <- contig_region_check %>% filter(X5 == "Dup") %>% ungroup() %>%
  summarise(
    q1 = quantile(contig_length, 0.25),
    q2 = quantile(contig_length, 0.50),
    q3 = quantile(contig_length, 0.75),
    min_len = min(contig_length),
    max_len = max(contig_length)
  )

# 각 분위 구간 실제 범위 추출
bounds <- tibble(
  length_group = factor(1:4, labels = c("Q1 (shortest)", "Q2", "Q3", "Q4 (longest)")),
  lower = c(quantile_info$min_len, quantile_info$q1, quantile_info$q2, quantile_info$q3),
  upper = c(quantile_info$q1, quantile_info$q2, quantile_info$q3, quantile_info$max_len)
) %>%
  mutate(label = paste0(
    scales::comma(round(lower/1e6,1)), "–", 
    scales::comma(round(upper/1e6,1)), " Mb"
  ))

# ──────────────────────────────
# ② 데이터에 분위 구간 부여
# ──────────────────────────────
plot_df <- contig_region_check %>% filter(X5 == "Dup") %>% ungroup() %>%
  mutate(
    length_group = ntile(contig_length, 4),
    length_group = factor(length_group, 
                          levels = 1:4,
                          labels = c("Q1 (shortest)", "Q2", "Q3", "Q4 (longest)"))
  )

# ──────────────────────────────
# ③ boxplot + 구간 주석
# ──────────────────────────────
p <- ggplot(plot_df, aes(x = length_group, y = prop.byContig, fill = length_group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 1.5) +
  #geom_jitter(width = 0.15, alpha = 0.3, size = 1, color = "black") +
  # 각 boxplot 아래에 텍스트 추가
  geom_text(
    data = bounds,
    aes(x = length_group, y = -0.02, label = label),   # y 위치 조정 (데이터 범위에 따라 변경)
    size = 4.2, vjust = 1, color = "gray20"
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Contig length quantile group\n(Shown range = contig length in Mb)",
    y = "Dup proportion (per contig)",
    title = "Distribution of Dup proportion by contig length quantile"
  ) +
  theme_step1() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 20, l = 10)
  )

p
head(plot_df)
head(merged_bed)
head(kppd)
head()
### dup region
dup_region <-  read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/flagger/dup_Resion/annot_Final_merge.bed",col_names = F)

dup_region
head(contig_region_check)
colnames(dup_region) <- c("chr","start","end","name","score","strand","MAPQ_max","tp_list","dv_list","top_label","top_bp","op_frac","labels_pass")


contig_region_check %>% mutate(contig = paste0(str_split_fixed(X2,"#",3)[,1],"_",str_split_fixed(X2,"#",3)[,3])) %>% 
  filter(X5 == "Dup") %>% ungroup %>% select(contig,region_sum,contig_length,prop.byContig) -> contig_region_check_dup
head(contig_region_check_dup)

contig_region_check_dup %>% 
  mutate(assembly = str_split_fixed(contig,"tg",2)[,1]) %>% 
  group_by(assembly) %>%
  summarise(Dup_region_sum = sum(region_sum),contig_length_sum = sum(contig_length)) -> contig_region_check_dup_byAssembly
  


dup_region %>% select(chr,start,end,name,top_label,top_bp,op_frac) %>%
  mutate(contig = paste0(str_split_fixed(name,"_",4)[,1],"_",str_split_fixed(name,"_",4)[,2])) %>% 
  mutate(assembly_start = str_split_fixed(name,"_",4)[,3]) %>% #head()
  mutate(assembly_end = str_split_fixed(name,"_",4)[,4]) %>% #filter(is.na(assembly_end))
  mutate(assembly_region_length = as.numeric(assembly_end) - as.numeric(assembly_start))

head(contig_region_check_dup)  

dup_region %>% select(top_label,labels_pass)

dup_region %>% select(chr,start,end,name,top_label,top_bp,op_frac)

dup_region %>% select(chr,start,end,name,top_label,top_bp,op_frac) %>%
  mutate(contig = paste0(str_split_fixed(name,"_",4)[,1],"_",str_split_fixed(name,"_",4)[,2])) %>% 
  mutate(assembly = str_split_fixed(contig,"tg",2)[,1]) %>% #head
  mutate(top_label = str_split_fixed(top_label,":",2)[,1]) %>% #head
  group_by(assembly,top_label) %>% #head()
  summarise(sum_top_bp = sum(top_bp)) %>% mutate(
    cohort = ifelse(assembly %in% c("KPPD132_h1","KPPD131_h1","KPPD130_h1","KPPD129_h1","KPPD132_h2","KPPD131_h2","KPPD130_h2","KPPD129_h2"),
                    "KPPD (T2T level)", "KPPD")) %>%
  left_join(contig_region_check_dup_byAssembly) -> dup_check_portion_vs_chm13

contig_region_check_dup %>% filter(contig == "KPPD001_h1tg000243l")


dup_region %>%
  separate_rows(labels_pass, sep = ";") %>% # 1. labels_pass를 ';' 기준으로 여러 행으로 분리
  ## 2. 분리된 labels_pass 안에서 다시 ':'로 세부항목 추출
  separate(labels_pass,into = c("sub_label", "source", "bp", "frac"),
           sep = ":", fill = "right", remove = FALSE) %>% 
  # 3. bp, frac은 수치형으로 변환  
  mutate(bp = as.numeric(bp), frac = as.numeric(frac)) %>%
  select(chr,start,end,name,sub_label,bp) %>%
  mutate(contig = paste0(str_split_fixed(name,"_",4)[,1],"_",str_split_fixed(name,"_",4)[,2])) %>% 
  mutate(assembly = str_split_fixed(contig,"tg",2)[,1]) %>% #head
  group_by(assembly,sub_label) %>% #head()
  summarise(sum_bp = sum(bp)) %>% mutate(
    cohort = ifelse(assembly %in% c("KPPD132_h1","KPPD131_h1","KPPD130_h1","KPPD129_h1","KPPD132_h2","KPPD131_h2","KPPD130_h2","KPPD129_h2"),
                    "T2T level", "Normal")) %>%
  left_join(contig_region_check_dup_byAssembly) -> dup_check_portion_vs_chm13_withSuplabelSum

head(dup_check_portion_vs_chm13_withSuplabelSum)
library(tidyverse)
head()
# ---- 비율 계산 ----
dup_summary_prop <- dup_check_portion_vs_chm13 %>%
  group_by(assembly) %>%
  mutate(prop = sum_top_bp / Dup_region_sum) %>% 
  mutate(top_label = if_else(top_label == "sd", "SD", top_label))  %>%
  ungroup()
head(dup_summary_prop)
# ---- top_label 없는 그룹(noDB) 채우기 ----
all_labels <- c("STR", "VNTR", "censat", "SD")

dup_summary_filled <- dup_summary_prop %>%
  tidyr::complete(assembly, top_label = all_labels,
                  fill = list(sum_top_bp = 0, prop = 0))

# ---- 그래프 ----
ggplot(dup_summary_filled, aes(x = assembly, y = prop, fill = top_label)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c(
      "STR" = "#F2C94C",      # 노랑
      "VNTR" = "#56CCF2",     # 하늘색
      "censat" = "#BB6BD9",   # 보라색
      "SD" = "#6FCF97"),
    name = "CHM13 DB"
  ) +
  labs(
    title = "Proportion of Dup regions by annotation type",
    x = "Assembly",
    y = "Proportion of Dup regions"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "top",
    axis.text.x = element_blank()
  ) + 
  facet_grid(~cohort, scales = "free_x", space = "free_x")



##### with sup 

head(dup_check_portion_vs_chm13_withSuplabelSum)

dup_summary_prop <- dup_check_portion_vs_chm13_withSuplabelSum %>%
  group_by(assembly) %>%
  mutate(prop = sum_bp / Dup_region_sum) %>% 
  mutate(sub_label = if_else(sub_label == "sd", "SD", sub_label))  %>%
  ungroup()
head(dup_summary_prop)
# ---- top_label 없는 그룹(noDB) 채우기 ----
all_labels <- c("STR", "VNTR", "censat", "SD")

dup_summary_filled <- dup_summary_prop %>%
  tidyr::complete(assembly, sub_label = all_labels,
                  fill = list(sum_bp = 0, prop = 0))
head(dup_summary_filled)
# ---- 그래프 ----
ggplot(dup_summary_filled, aes(x = assembly, y = prop, fill = sub_label)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c(
      "STR" = "#F2C94C",      # 노랑
      "VNTR" = "#56CCF2",     # 하늘색
      "censat" = "#BB6BD9",   # 보라색
      "SD" = "#6FCF97"),
    name = "CHM13 DB"
  ) +
  labs(
    title = "Proportion of Dup regions by annotation type",
    x = "Assembly",
    y = "Proportion of Dup regions"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "top",
    axis.text.x = element_blank()
  ) + 
  facet_grid(~cohort, scales = "free_x", space = "free_x")


###
fill_vals <- c(
  "STR" = "#F2C94C",      # 노랑
  "VNTR" = "#56CCF2",     # 하늘색
  "censat" = "#BB6BD9",   # 보라색
  "SD" = "#6FCF97"
)

base_layers <- list(
  geom_bar(stat = "identity", position = "stack"),
  scale_fill_manual(values = fill_vals,
                    breaks = c("STR","VNTR","censat","SD"),
                    name = "CHM13 DB"),
  labs(
    title = "Proportion of Dup regions by annotation type",
    x = "Assembly",
    y = "Proportion of Dup regions"
  ),
  theme_minimal(base_size = 13),
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )
) 

# 2) 각 cohort별 플롯 생성 ----

# 왼쪽 패널: 범례 유지
p_kppd <- ggplot(
  dup_summary_filled %>% filter(cohort == "Normal"),
  aes(x = assembly, y = prop, fill = sub_label)
) + base_layers +
  theme(legend.position = "top")

# 오른쪽 패널: 범례 생성 금지
p_t2t <- ggplot(
  dup_summary_filled %>% filter(cohort == "T2T level"),
  aes(x = assembly, y = prop, fill = sub_label)
) + base_layers +
  guides(fill = "none") +                # ← fill 범례 자체를 생성 안 함
  theme(
    legend.position = "none",            # ← 혹시 모를 표시도 차단
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    plot.title  = element_blank()
  )

# 결합: 공통 범례 1개만 수집
p <- (p_kppd + p_t2t) +
  plot_layout(widths = c(10, 1), guides = "collect") &&
  theme(legend.position = "top")


ggsave("../../figure/Proportion_of_DupRegion_by_CHM13annotationType.pdf",p,width = 15,height = 5,  useDingbats = FALSE)




###3



fill_vals <- c("STR"="#F2C94C","VNTR"="#56CCF2","censat"="#BB6BD9","sd"="#6FCF97","SD"="#6FCF97")

base_layers <- list(
  geom_bar(stat = "identity", position = "stack"),
  scale_fill_manual(values = fill_vals,
                    breaks = c("STR","VNTR","censat","SD"),
                    name = "CHM13 DB"),
  labs(x = "Assembly", y = "Proportion of Dup regions"),
  theme_minimal(base_size = 13),
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")
)

# 1) 왼쪽 패널: Normal
p_kppd <- ggplot(
  dup_summary_filled %>% filter(cohort == "Normal"),
  aes(x = assembly, y = prop, fill = sub_label)
) + base_layers +
  ggtitle("Normal")

# 2) 오른쪽 패널: T2T level (범례 생성 금지)
p_t2t <- ggplot(
  dup_summary_filled %>% filter(cohort == "T2T level"),
  aes(x = assembly, y = prop, fill = sub_label)
) + base_layers +
  ggtitle("T2T level") +
  guides(fill = "none") +
  theme(axis.title.y = element_blank())

p_t2t
# 3) 결합: 공통 범례 1개 + 전역 제목
(p_kppd + p_t2t) +
  plot_layout(widths = c(10, 1), guides = "collect") &
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.title.position = "plot",
    # 각 패널 제목 여백(위로 살짝 붙게)
    plot.margin = margin(t = 4, r = 5, b = 2, l = 5)
  ) +
  plot_annotation(title = "Proportion of Dup regions by annotation type")



head(dup_summary_filled)
bed_prop %>% group_by(ID,haplotype) %>% mutate(prop = prop.table(length)) %>% head()
dup_summary_filled %>% mutate(assembly = paste(ID,))

dup_summary_filled %>% group_by(assembly) %>% summarise(annotation = sum(sum_bp))
head9bed_prop

bed_prop %>% mutate(assembly = paste0(ID,"_h",haplotype)) %>% pivot_wider(values_from = length,names_from = region) %>% #head()
  left_join(dup_summary_filled %>% group_by(assembly) %>% summarise(annotation = sum(sum_bp))) %>%
  mutate(Dup = Dup - annotation, Hap = Hap + annotation) %>% #head()
  group_by(ID,haplotype) %>% select(-assembly,-annotation) %>% #head()
  pivot_longer(Col:Hap,names_to = "region",values_to = "length") %>% #head()
  mutate(prop = prop.table(length)) %>% group_by(region) %>%
  summarise(mean(length),mean(prop)) %>% t()


#####
## bam depth

merge_bam_nohaploid <- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/flagger/depth_nohap/merge.bam.depth.txt",col_names = F)
colnames(merge_bam_nohaploid) <- c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")
head(merge_bam_nohaploid)
merge_bam_nohaploid %>% mutate(contig = str_split_fixed(rname,'#',3)[,3]) %>%
  mutate(ID = str_split_fixed(rname,'#',3)[,1]) -> merge_bam_nohaploid

merged_bed  <- read_table("~/Desktop/KCDC/pangenome/KPPD/flagger/final.prediction.merge.txt",col_names = F)
head(merged_bed)



merged_bed %>% #mutate(Region = X4 - X3) %>% mutate(haplotype = str_split_fixed(X2,'#',3)[,2]) %>% #head()
  mutate(contig = str_split_fixed(X2,'#',3)[,3]) %>%
  mutate(ID = str_split_fixed(X2,'#',3)[,1]) %>%
  mutate(type = ifelse(str_detect(X1,"out"),"normal","annotation_CHM13")) %>% 
  filter(type == "annotation_CHM13") -> merge_bed

head(merge_bed)
head(merge_bam_nohaploid)

merge_bed %>% rename(startpos = X3, endpos = X4) %>% 
  left_join(merge_bam_nohaploid %>% mutate(startpos = startpos - 1)) %>% na.omit() %>%
  ggplot(aes(x=meandepth,y=meanmapq)) +
  geom_point() +
  facet_grid(~X5)


merge_bed %>% rename(startpos = X3, endpos = X4) %>% 
  left_join(merge_bam_nohaploid %>% mutate(startpos = startpos - 1)) %>% na.omit() %>%
  select(X2,X5,numreads:meanmapq) %>%
  pivot_longer(numreads:meanmapq) %>%
  ggplot(aes(x=X5,y=value)) +
  geom_violin() + 
  facet_wrap(~name,scales = 'free_y')
  

head(merge_bam_nohaploid)
merge_bed %>% rename(startpos = X3, endpos = X4) %>% 
  left_join(merge_bam_nohaploid %>% mutate(startpos = startpos - 1)) %>% 
  mutate(new_flagger_v2 = ifelse(X5 == "Hap","Hap",ifelse(meanmapq >= 20|meandepth >= 20,"Hap",X5))) %>% 
  mutate(new_flagger_v1 = ifelse(X5 == "Hap","Hap",ifelse(meanmapq >= 20,"Hap",X5))) %>% 
  
  mutate(Region = endpos-startpos) %>%
  mutate(haplotype = str_split_fixed(X2,'#',3)[,2]) %>%
  mutate(hapID = paste0(ID,"_",haplotype)) -> check
head(check)

check %>% 
  group_by(hapID,X5) %>% summarise(region_sum = sum(Region)) %>% #head()
  group_by(hapID) %>% mutate(contig_length = sum(region_sum)) %>%
  mutate(prop.byContig = prop.table(region_sum)) -> ori_flagger


check %>%   group_by(hapID,new_flagger_v1) %>% summarise(region_sum = sum(Region)) %>% #head()
  group_by(hapID) %>% mutate(contig_length = sum(region_sum)) %>%
  mutate(prop.byContig = prop.table(region_sum)) -> new_flagger_v1

check %>%   group_by(hapID,new_flagger_v2) %>% summarise(region_sum = sum(Region)) %>% #head()
  group_by(hapID) %>% mutate(contig_length = sum(region_sum)) %>%
  mutate(prop.byContig = prop.table(region_sum)) -> new_flagger_v2


head(ori_flagger)
new_flagger_v1
ori_flagger$type = "1. Original"
new_flagger_v1$type = "2. MAPQ >= 20"
new_flagger_v2$type = "3. MAPQ,depth >= 20 "

ori_flagger %>% rename(flagger = X5) %>% 
  rbind(new_flagger_v1 %>% rename(flagger = new_flagger_v1)) %>% rbind(new_flagger_v2%>% rename(flagger = new_flagger_v2)) %>% #head()
  ggplot(aes(x=flagger,y=prop.byContig,fill = type)) + 
  geom_boxplot() +
  facet_wrap(~flagger,nrow = 1,scales = "free")


ori_flagger %>% rename(flagger = X5) %>% 
  rbind(new_flagger_v1 %>% rename(flagger = new_flagger_v1)) %>% rbind(new_flagger_v2%>% rename(flagger = new_flagger_v2)) %>% #head()
  group_by(type,flagger) %>%
  summarise(mean_prop = mean(prop.byContig)) %>% pivot_wider(names_from = flagger,values_from = mean_prop)

ori_flagger %>% group_by(X5) %>%
  summarise(mean(prop.byContig))

new_flagger %>% group_by(new_flagger) %>%
  summarise(mean(prop.byContig)) 


#####
##### winnowmap 2
bed_files <- list.files(
  path = ".",
  pattern = "\\.bed$",
  recursive = TRUE,
  full.names = TRUE
)

# 파일이 제대로 인식됐는지 확인
print(bed_files)
a
head(a)
a$dir[,2]
bed_files %>% as.data.frame() %>% rename(path = ".")  %>% 
  mutate(dir = str_split_fixed(path,"/",4)[,2]) %>% #head()
  mutate(ID = str_split_fixed(dir,"_",2)[,1]) %>% count(ID)

merged_bed <- lapply(bed_files, function(f) {
  # 파일 읽기
  df <- read_tsv(f, col_names = FALSE,skip = 1)
  
  # ID 추출: 예) final_result/KPPD080_annotation/final_results/final_flagger_prediction.bed
  id <- sub("_annotation.*", "", basename(dirname(dirname(f))))
  
  # ID 열 추가
  df$ID <- id
  
  return(df)
}) %>% bind_rows()
head(merged_bed)

merged_bed %>% mutate(Region = X3 - X2) %>% mutate(haplotype = str_split_fixed(X1,'#',3)[,2]) %>%
  mutate(ID = str_split_fixed(X1,'#',3)[,1]) %>%  
  group_by(ID,haplotype,X4) %>% summarise(sum = sum(Region)) %>% mutate(prop = prop.table(sum)) %>%
  mutate(type = ifelse(ID %in% c("KPPD129","KPPD130","KPPD131","KPPD132"),"T2T","nonT2T")) %>%
  mutate(hap = paste0(ID,"_",haplotype)) -> bed_prop


bed_prop  <- read_table("~/Desktop/KCDC/pangenome/KPPD/flagger/final_flagger_prediction.KPPD.8K.bed")

head(bed_prop)
bed_prop %>% #filter(type == "annotation_CHM13") %>% select(ID,haplotype,X4,sum) %>% 
  rename(length = sum) %>% rename(region = X4) %>% #head()
  #writexl::write_xlsx("~/Desktop/KCDC/pangenome/KPPD/figure/Figure1_DATA.winnowmap.flagger_region.xlsx")
  writexl::write_xlsx("~/Desktop/KCDC/pangenome/KPPD/figure/Figure1_DATA.winnowmap.flagger_region.8K.xlsx")


######
library(tidyverse)
library(patchwork)


#bed_prop <- readxl::read_xlsx("~/Desktop/KCDC/pangenome/KPPD/figure/Figure1_DATA.winnowmap.flagger_region.xlsx")
bed_prop <- readxl::read_xlsx("~/Desktop/KCDC/pangenome/KPPD/figure/Figure1_DATA.winnowmap.flagger_region.8K.xlsx")
dim(bed_prop)
bed_prop

bed_prop %>% mutate(File = paste0(ID,"_hap",haplotype)) %>%
  mutate(Haplotype = paste0("hap",haplotype)) %>%
  rename(Sample = ID) %>% select(-prop) %>%
  pivot_wider(names_from = region,values_from = length) %>% 
  mutate(Total_unreliable = Col + Dup + Err) %>% #head()
  mutate(Total_unreliable_percent = round((Total_unreliable/(Hap + Total_unreliable)*100),2)) %>%
  rename(Haploid=Hap,Collapsed = Col,Duplicated = Dup, Erroneous = Err) %>%
  mutate(Batch = ifelse(Sample %in% c("KPPD129","KPPD130","KPPD131","KPPD132"),"T2T","Non.T2T")) %>%
  select(File,Sample,Haplotype,Batch,Haploid,Duplicated,Collapsed,Erroneous) %>%
  #erroneous	falsely_duplicated	haploid	collapsed
  mutate(Batch = ifelse(Sample %in% c("KPPD129","KPPD130","KPPD131","KPPD132"),"T2T","Non.T2T")) %>% 
  writexl::write_xlsx("~/Desktop/KCDC/pangenome/KPPD/flagger/sup.data.final.xlsx")
  
  


# 0) 색/라벨 고정
fill_colors <- c("Hap"="#66c2a5","Col"="#fc8d62","Dup"="#8da0cb","Err"="#e78ac3")
fill_labels <- c("Hap"="Haploid","Col"="Collapsed","Dup"="Duplicated","Err"="Erroneous")
head(df_all)
# 1) 공통 데이터 (팩터 레벨을 4종 모두로 고정)
df_all <- bed_prop %>%
  #filter(type == "annotation_CHM13") %>%
  mutate(
    hap = paste0(ID,"_",haplotype),
    sum_mb = length/1e6,
    region = factor(region, levels = c("Hap","Col","Dup","Err"))
  )
ord_total <- df_all %>%
  group_by(hap) %>%
  summarise(total_mb = sum(sum_mb), .groups = "drop") %>%
  arrange(desc(total_mb)) %>%
  pull(hap)

ord
# x축 정렬(원하면 기준 변경 가능)
df_all
ord <- df_all %>% filter(region=="Hap") %>% arrange(desc(sum_mb)) %>% pull(hap)
#df_all <- df_all %>% mutate(hap = factor(hap, levels = unique(c(ord, setdiff(hap, ord)))))
df_all <- df_all %>%
  mutate(
    hap = factor(hap, levels = ord_total)
  )

head(df_all)
# 2) 상단 = Haploid만
df_top <- df_all# %>% filter(X4 == "Hap")

p_top <- ggplot(df_top, aes(hap, sum_mb, fill = region)) +
  geom_col(width = 0.9, show.legend = FALSE) +
  coord_cartesian(ylim = c(2800, 3100), expand = FALSE) +
  scale_fill_manual(values = fill_colors, breaks = names(fill_labels),
                    labels = fill_labels, drop = FALSE) +  # 범례 유지용
  labs(y = "Total length (Mb)                       ") +
  #facet_grid(~type,scales = "free_x",space = "free_x") + 
  theme_step1() +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.line.x =  element_blank())

p_top

#p_top
# 3) 하단 = Haploid 제외 (Col→Dup→Err 순으로 누적)
df_bottom <- df_all %>%
  #filter(X4 != "Hap") %>%
  mutate(region = forcats::fct_relevel(region,"Hap","Col","Dup","Err"))

p_bottom <- ggplot(df_bottom, aes(hap, sum_mb, fill = region)) +
  geom_col(width = 0.9, position = "stack") +
  coord_cartesian(ylim = c(0, 400), expand = FALSE) +
  scale_fill_manual(values = fill_colors, breaks = names(fill_labels),
                    labels = fill_labels, drop = FALSE) +
  labs(y = "",x="Assemblies") +
  theme_step1() +
  theme(#axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "right")


# 4) 결합
p <- (p_top / p_bottom) + plot_layout(heights = c(2, 1), guides = "collect")
p

#padding = unit(c(top, right, bottom, left), "pt")
p <- (p_top + theme(plot.margin = margin(5, 5, 7, 5))) /
  (p_bottom + theme(plot.margin = margin(7, 5, 5, 5))) +
  plot_layout(heights = c(2, 1), guides = "collect")
p




head(bed_prop)
bed_prop %>% group_by(ID,haplotype) %>% mutate(prop = prop.table(length)) %>%
  group_by(region) %>% 
  summarise(mean(length),mean(prop),median(prop)) %>% t()


bed_prop %>% group_by(ID,haplotype) %>% mutate(prop = prop.table(length)) %>% #head()
  group_by(type,region) %>% 
  summarise(mean(length),mean(prop),median(prop)) %>% t()






setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/")
ggsave("Figure1X.Flagger_region.v1.pdf",p,width = 10,height = 5,  useDingbats = FALSE)

library(Cairo)
CairoPDF("Figure1b.cummulativeCoverage.v1.pdf", width = 6, height = 5)
print(p)
dev.off()


###3 EAS

bed_files <- list.files(
  path = ".",
  pattern = "\\.bed$",
  recursive = TRUE,
  full.names = TRUE
)

# 파일이 제대로 인식됐는지 확인
print(bed_files)
a
head(a)
a$dir[,2]
bed_files %>% as.data.frame() %>% rename(path = ".")  %>% 
  mutate(dir = str_split_fixed(path,"/",4)[,2]) %>% #head()
  mutate(ID = str_split_fixed(dir,"_",2)[,1]) %>% count(ID)

merged_bed <- lapply(bed_files, function(f) {
  # 파일 읽기
  df <- read_tsv(f, col_names = FALSE,skip = 1)
  
  # ID 추출: 예) final_result/KPPD080_annotation/final_results/final_flagger_prediction.bed
  id <- sub("_annotation.*", "", basename(dirname(dirname(f))))
  
  # ID 열 추가
  df$ID <- id
  
  return(df)
}) %>% bind_rows()
head(merged_bed)


merged_bed %>% mutate(Region = X3 - X2) %>% mutate(haplotype = str_split_fixed(X1,'#',3)[,2]) %>%
  mutate(ID = str_split_fixed(X1,'#',3)[,1]) %>%  
  group_by(ID,haplotype,X4) %>% summarise(sum = sum(Region)) %>% mutate(prop = prop.table(sum)) %>%
  #mutate(type = ifelse(ID %in% c("KPPD129","KPPD130","KPPD131","KPPD132"),"T2T","nonT2T")) %>%
  mutate(hap = paste0(ID,"_",haplotype)) -> bed_prop

bed_prop <- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/flagger/final_flagger_prediction.EAS10.8K.bed")
head(bed_prop)
bed_prop %>% #filter(type == "annotation_CHM13") %>% select(ID,haplotype,X4,sum) %>% 
  rename(length = sum) %>% rename(region = X4) %>% #head()
  write.table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/EAS_10sample.winnowmap.flagger_region.txt",col.names = T,row.names = F,quote = F,sep = "\t")


######### with EAS
library(ggpattern)
eas <- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/flagger/eas.10sample.flagger.txt")
head(eas)

bed_prop <- readxl::read_xlsx("~/Desktop/KCDC/pangenome/KPPD/figure/Figure1_DATA.winnowmap.flagger_region.xlsx")
dim(bed_prop)
bed_prop

bed_prop$cohort = "KPPD"
eas$cohort = "HPRC.EAS"

bed_prop %>% select(-type) %>%
  rbind(eas) -> bed_prop

fill_colors <- c("Hap"="#66c2a5","Col"="#fc8d62","Dup"="#8da0cb","Err"="#e78ac3")
fill_labels <- c("Hap"="Haploid","Col"="Collapsed","Dup"="Duplicated","Err"="Erroneous")

## -----------------------------
## 1) 공통 데이터 준비
## -----------------------------
df_all <- bed_prop %>%
  mutate(
    hap     = paste0(ID, "_", haplotype),
    sum_mb  = length / 1e6,
    region  = factor(region, levels = c("Hap", "Col", "Dup", "Err")),
    pattern = factor(ifelse(cohort == "HPRC.EAS", "stripe", "none"),
                     levels = c("none", "stripe"))
  )

## -----------------------------
## 2) 정렬을 "강제"로 만들기:
##    KPPD(및 기타) 전체 + spacer + HPRC.EAS 전체
##    각각 total_mb 내림차순
## -----------------------------
totals <- df_all %>%
  group_by(hap, cohort) %>%
  summarise(total_mb = sum(sum_mb), .groups = "drop")

ord_kppd <- totals %>%
  filter(cohort != "HPRC.EAS") %>%
  arrange(desc(total_mb)) %>%
  pull(hap)

ord_hprc <- totals %>%
  filter(cohort == "HPRC.EAS") %>%
  arrange(desc(total_mb)) %>%
  pull(hap)

spacer_hap <- "___SPACER___"
ord_total  <- c(ord_kppd, spacer_hap, ord_hprc)

## spacer row(0값) 추가
df_spacer <- df_all %>%
  slice(1) %>%
  mutate(
    hap     = spacer_hap,
    sum_mb  = 0,
    region  = factor("Hap", levels = levels(df_all$region)),
    pattern = factor("none", levels = levels(df_all$pattern)),
    cohort  = "SPACER"
  )

df_all2 <- bind_rows(df_all, df_spacer) %>%
  mutate(hap = factor(as.character(hap), levels = ord_total))
#df_all2 %>% 
## -----------------------------
## 3) Top panel
## -----------------------------
x_kppd
x_hprc

p_top <- ggplot(df_all2, aes(x = hap, y = sum_mb, fill = region)) +
  ggpattern::geom_col_pattern(
    aes(pattern = pattern),
    width = 0.9,
    show.legend = TRUE,          # fill legend는 필요
    colour = NA,
    pattern_fill = "black",
    pattern_colour = "black",
    pattern_density = 0.01,
    pattern_spacing = 0.05
  ) +
  coord_cartesian(ylim = c(2800, 3100), expand = FALSE) +
  scale_fill_manual(
    values = fill_colors,
    breaks = names(fill_labels),
    labels = fill_labels,
    drop = FALSE
  ) +
  scale_pattern_manual(values = c(none = "none", stripe = "stripe")) +
  labs(y = "Total length (Mb)                       ") + 
  theme_step1() +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank()
  ) +
  annotate(
    "text",
    x = 264/2,
    y = 3050,
    label = "KPPD",
    vjust = -1.5,
    fontface = "italic"
  ) +
  annotate(
    "text",
    x = 273,
    y = 3050,
    label = "HPRC\n(EAS)",
    vjust = 0,
    fontface = "italic"
  ) +
  #coord_cartesian(clip = "off") +
  ## ✅ stripe 관련 legend 완전 제거 + fill legend 키에서 stripe 제거
  guides(
    pattern = "none",
    fill = guide_legend(override.aes = list(pattern = "none"))
  )
#p_top
## -----------------------------
## 4) Bottom panel
## -----------------------------
p_bottom <- ggplot(df_all2, aes(x = hap, y = sum_mb, fill = region)) +
  ggpattern::geom_col_pattern(
    aes(pattern = pattern),
    width = 0.9,
    position = "stack",
    show.legend = TRUE,
    pattern_fill = "black",
    pattern_colour = "black",
    pattern_density = 0.001,
    pattern_spacing = 0.1
  ) +
  coord_cartesian(ylim = c(0, 400), expand = FALSE) +
  scale_fill_manual(
    values = fill_colors,
    breaks = names(fill_labels),
    labels = fill_labels,
    drop = FALSE
  ) +
  scale_pattern_manual(values = c(none = "none", stripe = "stripe")) +
  labs(x = "Assemblies", y = NULL) +
  theme_step1() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  ## ✅ stripe 관련 legend 완전 제거 + fill legend 키에서 stripe 제거
  guides(
    pattern = "none",
    fill = guide_legend(override.aes = list(pattern = "none"))
  )
#p_bottom
## -----------------------------
## 5) 결합 + gap 조정
## -----------------------------
p <- (p_top + theme(plot.margin = margin(5, 5, 7, 5))) /
  (p_bottom + theme(plot.margin = margin(7, 5, 5, 5))) +
  plot_layout(heights = c(2, 1), guides = "collect")

p

### compare defualt vs 1Mb
compare_defualtvs1Mb <- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/flagger/merge.winnowmap.comparedefaultvs1Mb.txt")
head(compare_defualtvs1Mb)

compare_defualtvs1Mb %>% select(Category_Name,Err,Dup,Hap,Col,SampleID) %>%
  mutate(type = ifelse(Category_Name == "WHOLE_GENOME_DEFAULT","Default","1MB")) %>%
  mutate(unreliable_region = Err + Dup + Col) %>% 
  mutate(prop_unreliable =unreliable_region/(unreliable_region+Hap)) %>% #head()
  ggplot(aes(x=type,y=prop_unreliable,fill=type)) + 
  geom_boxplot()



compare_defualtvs1Mb %>% 
  select(Category_Name, Err, Dup, Hap, Col, SampleID) %>%
  mutate(
    type = ifelse(Category_Name == "WHOLE_GENOME_DEFAULT", "Default", "1MB"),
    unreliable_region = Err + Dup + Col,
    sum_bp = unreliable_region + Hap
  ) %>%
  mutate(
    Err_ratio  = Err / sum_bp,
    Dup_ratio  = Dup / sum_bp,
    Hap_ratio  = Hap / sum_bp,
    Col_ratio  = Col / sum_bp,
    Unrel_ratio = unreliable_region / sum_bp
  ) %>% #head()
  select(SampleID,type,sum_bp:Unrel_ratio) %>% 
  pivot_longer(sum_bp:Unrel_ratio) -> a
a %>%
  ggplot(aes(x=type,y=value,fill=type)) +
  geom_boxplot() +
  facet_wrap(~name,scales = "free_y")
  

a %>% group_by(type,name) %>% 
  summarise(mean_pct=mean(value)) %>%
  pivot_wider(names_from = name,values_from = mean_pct)
  
