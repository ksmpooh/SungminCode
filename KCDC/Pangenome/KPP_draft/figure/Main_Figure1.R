### Figure 1 Final
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

cohort_cols <- c("KPP"="#1F4E9E", "HPRC"="#C00000", "CPC"="#2E7D32")

#KPPD (T2T) level
"KPP_T2T" = "#8F6B00"


library(Cairo)
library(patchwork)
library(ggbeeswarm)

library(tidyverse)
library(ggplot2)
library(googledrive)
library(googlesheets4)


#### PCA plot
df<- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/PCA/PCA_ethnic.kppdip_hprcdip.txt")
sample_1gkp<- read_tsv("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/PCA/1kgp.ethnic.info.tsv")
head(sample_1gkp)
sample_1gkp %>% select(1,7) %>%
  mutate(
    superpop5 = case_when(
      str_detect(`Superpopulation name`, "African")    ~ "AFR",
      str_detect(`Superpopulation name`, "American")   ~ "AMR",
      str_detect(`Superpopulation name`, "East Asia")  ~ "EAS",
      str_detect(`Superpopulation name`, "European")   ~ "EUR",
      str_detect(`Superpopulation name`, "South Asia") ~ "SAS",
      TRUE                                             ~ NA_character_
    )
  ) -> sample_1gkp
colnames(sample_1gkp) <- c("ID","rawg","GROUP") 



library(dplyr)
library(ggplot2)

pop5 <- c("AFR","AMR","EAS","EUR","SAS")

df_pca_plot <- df %>%
  mutate(
    cohort = case_when(
      str_detect(IID, "KPPD") ~ "KPP",
      str_detect(IID, "HPRC") ~ "HPRC",
      TRUE                    ~ "1KGP"
    )
  ) %>%
  left_join(sample_1gkp %>% rename(IID = ID), by = "IID") %>%
  mutate(
    GROUP = ifelse(cohort %in% c("KPP","HPRC"), cohort, GROUP),
    GROUP = factor(GROUP, levels = c("HPRC", "KPP",pop5))
  )


pop5 <- c("AFR","AMR","EAS","EUR","SAS")



ggplot(df_pca_plot, aes(PC1, PC2)) +
  
  ## 배경: 5 superpop
  geom_point(
    data = ~ filter(.x, GROUP %in% pop5),
    aes(color = GROUP),
    size = 0.7,
    alpha = 0.35
  ) +
  
  ## HPRC: 삼각형 (fill 색 + darkgrey 테두리)
  geom_point(
    data = ~ filter(.x, GROUP == "HPRC"),
    aes(fill = GROUP, shape = GROUP),
    size = 3,
    alpha = 1,
    colour = "lightgrey",
    stroke = 0.3,
    show.legend = c(fill = FALSE, shape = TRUE)
  ) +
  
  ## KPPD: 네모 (fill 색 + darkgrey 테두리)
  geom_point(
    data = ~ filter(.x, GROUP == "KPP"),
    aes(fill = GROUP, shape = GROUP),
    size = 3,
    alpha = 1,
    colour = "lightgrey",
    stroke = 0.4,
    show.legend = c(fill = FALSE, shape = TRUE)
  ) +
  
  ## 색상 정의 (background + fill 공용)
  scale_fill_manual(
    values = c("HPRC" = "#C00000","KPP" = "#1F4E9E")
  ) +
  scale_color_manual(
    values = c("AFR"  = "#E69F00","AMR"  = "#009E73","EAS"  = "#56B4E9","EUR"  = "#CC79A7","SAS"  = "#F0E442"),
    breaks = pop5
  ) +
  
  
  ## shape: 테두리 가능한 도형으로 변경
  scale_shape_manual(
    values = c("HPRC" = 24,  "KPP" = 22),
    breaks = c("HPRC","KPP")
  ) +
  guides(shape = guide_legend(order = 1,override.aes = list(alpha=1,size= 4,fill= c("#C00000", "#1F4E9E"),colour = "lightgrey")),
    color = guide_legend(order = 2,override.aes = list(alpha = 1, size = 4))
  ) + 
  
  theme_step1() +
  theme(
    legend.position = c(0.98, 0.02),
    legend.justification = c("right", "bottom"),
    legend.title = element_blank()
  ) -> f1.a2

#f1.a2
f1.a2

p
#CairoPDF("Figure1/Figure1a.EthnicPlot_KPPD_HPRC_CPP.pdf", width = 7, height = 5)
#print(p)
dev.off()

f1.a2



##### ccummulatvie coverage

df_raw <- read.table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/Figure1/Figure1.data.cummulativeCoverage_KPPD_HPRC_CHM13_GRCH38nsplit.txt",sep = "\t",header = T)
df_sum <- read.table("draft_korean_pangenome_collection/FIgure1B_Assembly_Contiguity_1/Supplement_T2_Contiguity_gfastats_haplotype.tab",sep = "\t",header = T)
df_sum_kppd <- read.table("draft_korean_pangenome_collection/FIgure1B_Assembly_Contiguity_1/Supplement_T2_Contiguity_gfastats.tab",sep = "\t",header = T)

head(df_sum_kppd)
df_sum_kppd %>% summarise(mean(Contigs),mean(Total_contig_length),mean(Contig_N50),mean(Contig_L50))
df_sum_kppd %>% group_by(Batch) %>% 
  summarise(mean(Contigs),mean(Total_contig_length),mean(Contig_N50),mean(Contig_L50))

head(df_sum)
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

#ggsave("Figure1/Figure1b.cummulativeCoverage_GRChNsplit.v2.pdf",p,width = 6,height = 5,  useDingbats = FALSE)

ggplot() +
  geom_step(data = df %>% filter(cohort %in% c("HPRC", "KPPD")),
            aes(x = Coverage * 3, y = len / 1000000,color = cohort,group = name),
            size = 0.5,alpha = 0.2,linetype = "solid" ) +
    geom_step(
    data = df %>% filter(cohort == "KPPD (T2T level)"),
    aes(x = Coverage * 3,y = len / 1000000,color = cohort,group = name),
    size = 0.7,alpha = 1,linetype = "solid") +
  
  geom_step(data = df %>% filter(cohort %in% c("GRCh38", "CHM13v2")),aes(x = Coverage * 3,y = len / 1000000,group = name),
    color = "black",size = 1,alpha = 1,show.legend = FALSE) +
  
  geom_vline(xintercept = 1.5, linetype = "dotted", linewidth = 0.5) +
  scale_color_manual(
    breaks = c("HPRC","KPPD","KPPD (T2T level)"),
    values = c("HPRC"="red","KPPD"= "blue","KPPD (T2T level)" = "#8F6B00"),
    labels = c("HPRC"="HPRC","KPPD"="KPP","KPPD (T2T level)" = "KPP (T2T level)")) + 
  
  labs(x = "Cumulative coverage (Gb)",y = "Length (Mb)") +
  theme_step1() +
  theme(
    legend.position = c(0.87, 0.8),
    legend.title = element_blank()
  ) +
  guides(
    color = guide_legend(override.aes = list(alpha = 2))
  ) +
  
  ## Annotations
  annotate("text", x = 0.8, y = 200, label = "CHM13v2",
           color = "black", size = 5, hjust = 0) +
  annotate("text", x = 0.8, y = 100, label = "GRCh38",
           color = "black", size = 5, hjust = 0) -> f1.d
f1.d



############################

#write.table(df2,"Figure1/Figure1.data.TotalContigLength_bySex.txt",col.names = T,row.names = F,quote = F,sep = "\t")

df2 <- read.table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/Figure1/Figure1.data.TotalContigLength_bySex.txt",header = T,sep = "\t")
head(df2)
df2 %>% summarise(mean(total_gb))
df2 %>% group_by(sex) %>%
  summarise(mean(total_gb))

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
ggplot(df2, aes(x = sex, y = total_gb, fill = sex)) +
  
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
    color = "#1F4E9E", size = 2, alpha = 0.8, width = 0.15,shape=1
  ) +
  
  # 4️⃣ T2T 샘플: 형광 연두 포인트 (제일 강조)
  geom_point(
    data = df_t2t,
    aes(x = sex, y = total_gb, shape = Sample),
    size = 3, stroke = 1.1,
    color = "#8F6B00",  # 💡 형광 연두색
    position = position_jitter(width = 0.05, height = 0)
  ) +
  
  scale_x_discrete(labels = sex_labels) +
  scale_fill_manual(values = c(F = "grey80", M = "grey80")) +  # 은은한 배경색
  scale_shape_manual(values = shape_vals) +
  
  labs(
    x = NULL,
    y = "Total contig length (Gb)",
    subtitle = "* T2T level: Bronze"
  ) +
  
  theme_step1() +
  theme(
    legend.position = "none",
    #axis.title.y = element_text(margin = margin(r = 10), color = "black"),
    axis.text = element_text(color = "black"),
    plot.subtitle = element_text(hjust = 1, face = "italic", color = "gray30") 
  ) -> f1.b

f1.b
plot(f1.b)





#ggsave("Figure1/Figure1b.TotalContigLength_bySex.pdf",p,width = 1,height = 5,  useDingbats = FALSE)


#library(Cairo)
#CairoPDF("Figure1/Figure1b.TotalContigLength_bySex.pdf", width = 4, height = 5)
#print(p)
#dev.off()
head(df_t2t)

df2 %>% group_by(sex) %>% summarise(mean(total_gb))
df_2 %>% group_by(sex) %>% summarise(mean(total_gb))
df2 %>% group_by(Batch,sex) %>% summarise(mean(total_gb))

############
### auN


url <- "https://docs.google.com/spreadsheets/d/1Y8cqzCizaQ9vvAzoAEhqqGW7GaplxRRnmPsQ0dK6o30/edit?usp=sharing"
#drive_auth()

props <- sheet_properties(url)
props$name
#13,"1.3.asm_correcteness1_inspector"
#14,"1.3.asm_correcteness1_merqury"   
#11 "1.2.assembly_contiguity" 

assembly_aun <- read_sheet(url,sheet = 11)  
inspector <- read_sheet(url,sheet = "1.3.asm_correcteness1_inspector")  # 기본: 첫 번째 시트
mercury <- read_sheet(url,sheet = "1.3.asm_correcteness1_merqury")  # 기본: 첫 번째 시트
head(inspector)
head(mercury)


#assembly_aun$auN
head(plot_df)
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
  geom_point(data = ~ filter(.x, !is_T2T),aes(shape = Sex),color = "#1F4E9E",alpha = 0.7,size = 3) +
  ## points: T2T (highlight)
  geom_point( data = ~ filter(.x, is_T2T),aes(shape = Sex),color = "#8F6B00",alpha = 1,size = 3) +
  
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

a

plot_df
inspector
inspector %>% select(Sample,QV,Haplotype,Batch,`Total length`) %>%
  group_by(Sample) %>%  #head()
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
    color = "#1F4E9E",
    alpha = 0.7,
    size = 3
  ) +
  
  ## T2T (highlight)
  geom_point(
    data = ~ filter(.x, is_T2T),
    aes(shape = Sex),
    color = "#8F6B00",
    alpha = 1,
    size = 3
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

b
#a+b
cowplot::plot_grid(a,b) -> f1.c
f1.c




###
inspector %>% select(File,`Base substitution`,`Small-scale collapse`,`Small-scale expansion`) %>% #head()
  rename("Small Collapse" = `Small-scale collapse`,"Small Expansion" = `Small-scale expansion`,`Base Substitution`=`Base substitution`) %>%
  pivot_longer(2:4) %>%
  ggplot(aes(x=value,fill=name)) + 
  geom_density(alpha = 0.5) +
  scale_y_continuous(
    breaks = seq(0, 0.02, by = 0.004),   # 필요 범위에 맞게 상한 조절
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Number of small-scale assembly errors",
    y = "Density",
    fill = NULL
  ) +
  theme_step1() +
  theme(
    #legend.direction = "horizontal",
    legend.position = c(0.98, 0.98),        # 오른쪽 위 (plot 내부)
    legend.justification = c("right", "top")
    
    #legend.text = element_text(size=1),
#    legend.position = "top"   # 위에 legend
  ) -> f1.e

f1.e



######## Flgger 

# url https://docs.google.com/spreadsheets/d/13ng9s0i83w5V6ckM_tNUyn22WsHl-qNunYpfRrOxJ0Q/edit?usp=sharing

#bed_prop
eas <- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/EAS_10sample.winnowmap.flagger_region.txt")
head(eas)


bed_prop <- readxl::read_xlsx("~/Desktop/KCDC/pangenome/KPPD/figure/Figure1_DATA.winnowmap.flagger_region.xlsx")
dim(bed_prop)
bed_prop

bed_prop$cohort = "KPP"
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
ord_total  <- c(ord_hprc,spacer_hap, ord_kppd)

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

cohort_cols <- c("KPP" = "#1F4E9E", "HPRC.EAS" = "#C00000")

df_all2 <- bind_rows(df_all, df_spacer) %>%
  mutate(hap = factor(as.character(hap), levels = ord_total)) %>%
  mutate(outline = ifelse(cohort %in% names(cohort_cols), cohort, NA_character_))

p_top <- ggplot(df_all2, aes(x = hap, y = sum_mb, fill = region)) +
  geom_col(
    aes(color = outline),   # ✅ 테두리만 cohort 색
    width = 0.9,
    linewidth = 0.35        # 테두리 두께(원하면 0.5~0.8로)
  ) +
  coord_cartesian(ylim = c(2850, 3100), expand = FALSE) +
  scale_y_continuous(breaks = seq(2900, 3100, by = 100)) + 
  scale_fill_manual(
    values = fill_colors,
    breaks = names(fill_labels),
    labels = fill_labels,
    drop = FALSE
  ) +
  scale_color_manual(
    values = cohort_cols,
    na.value = NA           # ✅ SPACER는 테두리 없음
  ) +
  labs(y = "Total length (Mb)                       ") +
  theme_step1() +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    legend.position = "right"
  ) +
  # ✅ cohort 테두리 색 legend는 숨기고(fill legend만 남김)
  guides(color = "none") + 
  #annotate(
   # "text",
  #  x = 145,
 #   y = 3050,
#    label = "KPPD",
   # vjust = -1.5,
  #  fontface = "italic"
  #) +
  annotate(
    "text",
    x = 10,
    y = 3050,
    label = "HPRC\n(EAS)",
    vjust = 0,
    fontface = "italic"
  )

p_top

p_bottom <- ggplot(df_all2, aes(x = hap, y = sum_mb, fill = region)) +
  geom_col(
    aes(color = outline),
    width = 0.9,
    position = "stack",
    linewidth = 0.25
  ) +
  coord_cartesian(ylim = c(0, 500), expand = FALSE) +
  scale_y_continuous(
    breaks = c(0,200,400)) + 
  scale_fill_manual(
    values = fill_colors,
    breaks = names(fill_labels),
    labels = fill_labels,
    drop = FALSE
  ) +
  scale_color_manual(
    values = cohort_cols,
    na.value = NA
  ) +
  labs(x = "Assemblies", y = NULL) +
  theme_step1() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  guides(color = "none")

f1.h <- (p_top + theme(plot.margin = margin(5, 5, 7, 5))) /
  (p_bottom + theme(plot.margin = margin(7, 5, 5, 5))) +
  plot_layout(heights = c(2, 1), guides = "collect")

f1.h




head(bed_prop)
bed_prop %>% group_by(ID,haplotype) %>% mutate(prop = prop.table(length)) %>%
  group_by(region) %>% 
  summarise(mean(length),mean(prop),median(prop)) %>% t()


bed_prop %>% group_by(ID,haplotype) %>% mutate(prop = prop.table(length)) %>% #head()
  group_by(type,region) %>% 
  summarise(mean(length),mean(prop),median(prop)) %>% t()

head(eas)
eas %>% group_by(region) %>%
  summarise(median(prop))

eas %>% select(hap,region,prop) %>%
  pivot_wider(names_from = region,values_from = prop)



eas %>% select(hap,region,prop) %>% #head()
  ggplot(aes(x=region,y=prop,fill=region)) + 
  geom_boxplot() + 
  theme_step1() + 
  facet_wrap(~region,scales = "free",nrow=1)
  

bed_prop %>% count(cohort)


bed_prop %>% mutate(cohort = ifelse(ID %in% c("KPPD132","KPPD131","KPPD130","KKPD129"),"KPP_T2T",cohort)) %>% #head()
  group_by(cohort,region) %>%
  summarise(median= median(prop)) %>%
  pivot_wider(names_from = region,values_from = median)


bed_prop %>% mutate(cohort = ifelse(ID %in% c("KPPD132","KPPD131","KPPD130","KKPD129"),"KPP_T2T",cohort)) %>% #head()
  select(hap,region,cohort,prop) %>%
  ggplot(aes(x=cohort,y=prop,fill=cohort)) + 
  geom_boxplot() + 
  theme_step1() +
  facet_wrap(~region,nrow = 2,scales = "free")



## liftoff gene annotation


url <- "https://docs.google.com/spreadsheets/d/13ng9s0i83w5V6ckM_tNUyn22WsHl-qNunYpfRrOxJ0Q/edit?usp=sharing"
props <- sheet_properties(url)
props$name
#S6,

liftoff <- read_sheet(url,sheet = "S7",skip = 1)  

liftoff
head(df2)

liftoff %>% left_join(df2 %>% select(Sample,sex) %>% unique()) %>% #head()
  pivot_longer(Protein_coding_genes:Noncoding_trans,names_to = "coding_type",values_to = "annotation_pct") %>% #head()
  mutate(coding_type = factor(coding_type, levels = c("Protein_coding_genes", "Noncoding_genes","Protein_coding_trans","Noncoding_trans"))) %>%
  ggplot(aes(x = coding_type, y = annotation_pct,color = coding_type)) +
  geom_quasirandom(shape = 1, alpha = 1, dodge.width = 0.8) +
  #labs(y="Percentage\nannotated (%)") +
  labs(y="Annotation rate (%)") +
  scale_x_discrete(labels = c("Protein_coding_genes" = "Protein-coding\ngenes", 
                              "Noncoding_genes" = "Noncoding\ngenes",
                              "Protein_coding_trans" = "Protein-coding\ntranscripts", 
                              "Noncoding_trans" = "Noncoding\ntranscripts")) +
  scale_y_continuous(
    breaks = seq(90, 100, by = 5),   # 필요 범위에 맞게 상한 조절
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_step1() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(), 
        legend.position = 'none', 
        strip.text.x = element_text(
          size = 14,        # ??Ʈ ũ??
          face = "bold",    # ??Ʈ ???? (bold, italic ??)
          family = "Arial")) -> f1.g
f1.g


liftoff %>% left_join(df2 %>% select(Sample,sex) %>% unique()) %>% #head()
  pivot_longer(Protein_coding_genes:Noncoding_trans,names_to = "coding_type",values_to = "annotation_pct") %>% #head()
  mutate(coding_type = factor(coding_type, levels = c("Protein_coding_genes", "Noncoding_genes","Protein_coding_trans","Noncoding_trans"))) %>% 
  #group_by(Batch) %>% summarise(max(annotation_pct),min(annotation_pct),median(annotation_pct),mean(annotation_pct))
  #group_by(Batch,coding_type) %>% summarise(pct=mean(annotation_pct)) %>% pivot_wider(names_from = coding_type,values_from = pct)
  group_by(coding_type) %>% summarise(pct=mean(annotation_pct)) %>% pivot_wider(names_from = coding_type,values_from = pct)
  

#### repeat masker
### repeatmasker pattern
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/repeat_masker")
df <- read_table("merge.repeat_size.txt")
hprc <- read_table("merge.repeat_size.hprc.txt")
head(df)
head(hprc)

df$cohort <- "KPP"
hprc$cohort <- "HPRC"


#####
df %>%
  rbind(hprc) %>%
  mutate(
    cohort = factor(cohort,levels = c("HPRC","KPP")),
    Category = factor(
      Category,
      levels = c(
        "LINE","SINE","LTR","Satellite",
        "DNA_transposon","Simple_repeat",
        "Low_complexity","Other"
      )
    )
  ) %>%
  ggplot(aes(x = cohort, y = bp, color = cohort, fill = cohort)) +
  
  geom_violin(
    width = 0.9,
    alpha = 0.5,
    linewidth = 0
  ) +
  
  geom_quasirandom(
    width = 0.25,
    alpha = 0.5,
    size = 1
  ) +
  
  facet_wrap(
    ~ Category,
    scales = "free",
    nrow = 2
  ) +
  
  scale_y_continuous(
    labels = function(x) x / 1e6
  ) +
  scale_color_manual(
    values = c("HPRC" = "#C00000", "KPP" = "#1F4E9E")
  ) +
  scale_fill_manual(
    values = c("HPRC" = "#C00000", "KPP" = "#1F4E9E")
  ) +
  
  labs(y = "Length (Mb)") +
  guides(color = guide_legend(override.aes = list(alpha = 1)),fill = guide_legend(override.aes = list(alpha = 1))) +
  theme_step1() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold",size = 12),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
 ) -> f1.f



f1.f
library(png)
library(Cairo)
library(grid)
library(cowplot)

img <- readPNG("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/Figure1/f1.a.png")

f1.a1 <- ggplot() +
  annotation_custom(
    rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"))
  ) +
  theme_void() + 
  theme(plot.margin = unit(c(0.25, 1, 1, 0.5), 'cm'))

f1.a1
plot_grid(f1.a1,f1.a2,f1.b,nrow = 1,rel_widths = c(1,1,1),labels = c("A","","B"),label_size = 25) -> f1.row1
#f1.row1
plot_grid(f1.c,f1.d,nrow = 1,rel_widths = c(2,1),labels = c("C","D"),label_size = 25) ->  f1.row2

f1.g+theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 30)) -> f1.g_new

plot_grid(f1.e,f1.g_new,nrow = 2,rel_heights = c(1.2,1),labels = c("E","F"),label_size = 25) -> f1.row3.left

#plot_grid(f1.e,f1.g+theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 20)),nrow =  1,rel_widths = c(1,1.5),labels = c("E","F"),label_size = 20)

plot_grid(f1.row3.left,f1.f,nrow = 1,rel_widths = c(1,1.2),labels = c("","G"),label_size = 25) -> f1.row3
plot_grid(f1.h,nrow = 1,labels = c("H"),label_size = 25) -> f1.row4
f1.row4
#f1.h


plot_grid(f1.row1,f1.row2,f1.row3,f1.row4,nrow = 4,rel_heights = c(1.1,1,1.2,1)) -> f1
f1


ggsave(
  filename = "/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/Figure1/test/Figure1.v1.pdf",
  plot = f1,
  device = cairo_pdf,
  width = 18,      # inch
  height =22,     # inch
  dpi = 300
)

ggsave(
  filename = "/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/Figure1/test/Figure1.v1.png",
  plot = f1,
  device = "png",
  width = 18,      # inch
  height = 22,     # inch
  dpi = 300
)

####################


cohort_cols <- c("KPPD"="#1F4E9E", "HPRC"="#C00000", "CPC"="#2E7D32")

KPPD (T2T) level
"KPPD_T2T" = "#8F6B00"


ggsave(
  filename = "/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/Figure1/Figure1.a2.PCA.png",
  plot = f1.a2,
  device = "png",
  width = 8,      # inch
  height = 6,     # inch
  dpi = 300
)

ggsave(
  filename = "/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/Figure1/Figure1.b.contig_length.png",
  plot = f1.b,
  device = "png",
  width = 6,      # inch
  height = 6,     # inch
  dpi = 300
)

cowplot::plot_grid(a,b) -> f1.c
ggsave(
  filename = "/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/Figure1/Figure1.c.auN.png",
  plot = a,
  device = "png",
  width = 7,      # inch
  height = 6,     # inch
  dpi = 300
)

ggsave(
  filename = "/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/Figure1/Figure1.c.QV.png",
  plot = b,
  device = "png",
  width = 7,      # inch
  height = 6,     # inch
  dpi = 300
)



ggsave(
  filename = "/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/Figure1/Figure1.d.Cummulative.png",
  plot = f1.d,
  device = "png",
  width = 8,      # inch
  height = 6,     # inch
  dpi = 300
)

ggsave(
  filename = "/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/Figure1/Figure1.e.inpector.png",
  plot = f1.e,
  device = "png",
  width = 10,      # inch
  height = 4,     # inch
  dpi = 300
) 


ggsave(
  filename = "/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/Figure1/Figure1.f.repeatmaskers.png",
  plot = f1.f,
  device = "png",
  width = 16,      # inch
  height = 10,     # inch
  dpi = 300
)




ggsave(
  filename = "/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/Figure1/Figure1.g.liftoff.png",
  plot = f1.g,
  device = "png",
  width = 10,      # inch
  height = 4,     # inch
  dpi = 300
)


ggsave(
  filename = "/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/Figure1/Figure1.h.flagger.png",
  plot = f1.h,
  device = "png",
  width = 18,      # inch
  height = 4,     # inch
  dpi = 300
)

setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/test/")


DPI   <- 300
W_STD <- 7
H_STD <- 5
OUT   <- "/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/Figure1/panels"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# --- (중요) 너의 theme_step1가 base_size를 받게 되어 있어야 함 ---
# theme_step1 <- function(base_size = 11, base_family = "Arial", ...) { ... }

save_panel_scaled <- function(p, fname, w, h,
                              base_size_std = 16,  # 7x5에서 원하는 글씨 크기
                              std_w = W_STD, std_h = H_STD,
                              family = "Arial",
                              dpi = DPI) {

  # 캔버스 면적 기준 스케일 팩터 (가장 안정적)
  s <- sqrt((w * h) / (std_w * std_h))

  # theme_step1를 base_size로 재적용
  p2 <- p + theme_step1(base_size = base_size_std * s, base_family = family)

  ggsave(
    filename = file.path(OUT, fname),
    plot     = p2,
    device   = "png",
    width    = w, height = h, units = "in",
    dpi      = dpi,
    bg       = "white"
  )
}
21
save_panel_scaled(f1.a2, "Figure1.a.png", w = W_STD,      h = H_STD)
save_panel_scaled(f1.b,  "Figure1.b.png", w = W_STD,      h = H_STD)
save_panel_scaled(f1.c,  "Figure1.c.png", w = W_STD*2,    h = H_STD)
save_panel_scaled(f1.d,  "Figure1.d.png", w = W_STD,      h = H_STD)
save_panel_scaled(f1.e,  "Figure1.e.png", w = W_STD*1.2,  h = H_STD*0.8)
save_panel_scaled(f1.f,  "Figure1.f.png", w = W_STD*1.2,  h = H_STD*1.8)
save_panel_scaled(f1.g,  "Figure1.g.png", w = W_STD*1.5,  h = H_STD*0.8)
save_panel_scaled(f1.h,  "Figure1.h.png", w = W_STD*2.5,  h = H_STD*0.8)

