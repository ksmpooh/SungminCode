### 2026 KOGO CNV

library(tidyverse)
setwd("/ADATA/pangenome/liftoff/kogo_2026")

#cpc <- read_table("/CDATA/pangenome/liftoff/merge/CPC/CPC.all.merge.with_flagger.gene_protein_coding.selected_withhgnc.txt")

hprc_gene <- read_table("/CDATA/pangenome/liftoff/hprc/gencode_GRCh38.v38_partial/merge/merge.HPRC.liftoff_afterprep.gene.txt")
hprc_trans <- read_table("/CDATA/pangenome/liftoff/hprc/gencode_GRCh38.v38_partial/merge/merge.HPRC.liftoff_afterprep.trans.txt")

hprc_gene %>% 
  mutate(ID = ifelse(str_detect(contig,"chr"),"CHM13",(paste0(str_split_fixed(contig,"#",3)[,1],"_h",str_split_fixed(contig,"#",3)[,2])))) %>%
  mutate(coding_type = ifelse(gene_type=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  select(ID,coding_type,gene_name) %>% unique() %>% 
  count(ID,coding_type) %>% write.table("HPRC_genecode.v.38.gene.protein_coding_count.txt",col.names = T,row.names = F,quote = F,sep = "\t")


hprc_trans %>% select(contig,transcript_type,transcript_name) %>%
  mutate(ID = ifelse(str_detect(contig,"chr"),"CHM13",(paste0(str_split_fixed(contig,"#",3)[,1],"_h",str_split_fixed(contig,"#",3)[,2])))) %>%
  #mutate(coding_type = ifelse(gene_type=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  mutate(coding_type = ifelse(transcript_type =="protein_coding","Protein_coding_trans","Noncoding_trans")) %>%
  select(ID,coding_type,transcript_name) %>% unique() %>% 
  count(ID,coding_type) %>% write.table("HPRC_genecode.v.38.trans.protein_coding_count.txt",col.names = T,row.names = F,quote = F,sep = "\t")


df_gene <- read_table("/ADATA/pangenome/liftoff/02.annotation/merge.genecode.v38.prep.txt_forR") %>% mutate(ID = str_split_fixed(contig,"tg",2)[,1])
df_trans <- read_table("/ADATA/pangenome/liftoff/02.annotation/merge.genecode_transcript.v38.prep.txt_forR") %>% mutate(ID = str_split_fixed(contig,"tg",2)[,1])

df_trans %>% mutate(coding_type = ifelse(transcript_type=="protein_coding","Protein_coding_trans","Noncoding_trans")) %>% 
  select(ID,coding_type,transcript_name) %>% unique() %>% #dim()
  count(ID,coding_type) %>% write.table("KPPD.transcript.annotation.byhap.txt",col.names = T,row.names = F,quote = F,sep = "\t") 


df_gene %>% mutate(coding_type = ifelse(gene_type=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>% 
  select(ID,coding_type,gene_name) %>% unique() %>%
  count(ID,coding_type) %>% write.table("KPPD.gene.annotation.byhap.txt",col.names = T,row.names = F,quote = F,sep = "\t") 


### CNV
hprc_gene %>% filter(gene_type == "protein_coding") %>% 
  mutate(ID = ifelse(str_detect(contig,"chr"),"CHM13",(paste0(str_split_fixed(contig,"#",3)[,1],"_h",str_split_fixed(contig,"#",3)[,2])))) %>%
  filter(coverage >= 0.9) %>%
  filter(sequence_ID>=0.9) %>%
  filter(level < 3) %>%
  filter(valid_ORFs > 0) -> hprc_gene_proteincoding_QC


hprc_gene_proteincoding_QC %>%
  count(ID,gene_name) %>% filter(n != 1) %>% 
  count(ID) %>% write.table("HPRC.liftoff.CNV.bysample.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#write.table("/ADATA/pangenome/liftoff/99.summary/gene.annotation.byhap.txt",col.names = T,row.names = F,quote = F,sep = "\t") 
hprc_gene_proteincoding_QC %>%
  count(ID,gene_name) %>% filter(n != 1) %>% 
  count(gene_name) %>% write.table("HPRC.liftoff.CNV.byGene.txt",col.names = T,row.names = F,quote = F,sep = "\t")

kppd <- read_table("/ADATA/pangenome/liftoff/02.annotation/merge.genecode.v38.prep.txt") %>% mutate(ID = str_split_fixed(contig,"tg",2)[,1])


kppd %>% filter(gene_type == "protein_coding") %>% 
  filter(coverage >= 0.9) %>%
  filter(sequence_ID>=0.9) %>%
  filter(level < 3) %>%
  filter(valid_ORFs > 0) -> kppd_gene_proteincoding_QC


kppd_gene_proteincoding_QC %>%
  count(ID,gene_name) %>% filter(n != 1) %>% 
  count(ID) %>% write.table("KPPD.liftoff.CNV.bysample.txt",col.names = T,row.names = F,quote = F,sep = "\t")

#write.table("/ADATA/pangenome/liftoff/99.summary/gene.annotation.byhap.txt",col.names = T,row.names = F,quote = F,sep = "\t") 

kppd_gene_proteincoding_QC %>%
  count(ID,gene_name) %>% filter(n != 1) %>% 
  count(gene_name) %>% write.table("KPPD.liftoff.CNV.byGene.txt",col.names = T,row.names = F,quote = F,sep = "\t")

cpc <- read_table("/CDATA/pangenome/liftoff/merge/CPC/CPC.all.merge.with_flagger.gene_protein_coding.select.txt")
cpc %>%
  mutate(ID = ifelse(str_detect(contig,"chr"),"CHM13",(paste0(str_split_fixed(contig,"#",3)[,1],"_h",str_split_fixed(contig,"#",3)[,2])))) %>%
  filter(coverage >= 0.9) %>%
  filter(sequence_ID>=0.9) %>%
  filter(level < 3) %>%
  filter(valid_ORFs > 0) -> cpc_gene_proteincoding_QC


cpc_gene_proteincoding_QC %>%
  count(ID,gene_name) %>% filter(n != 1) %>% 
  count(ID) %>% write.table("CPC.liftoff.CNV.bysample.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#write.table("/ADATA/pangenome/liftoff/99.summary/gene.annotation.byhap.txt",col.names = T,row.names = F,quote = F,sep = "\t") 
cpc_gene_proteincoding_QC %>%
  count(ID,gene_name) %>% filter(n != 1) %>% 
  count(gene_name) %>% write.table("CPC.liftoff.CNV.byGene.txt",col.names = T,row.names = F,quote = F,sep = "\t")



#############################################
####  plot

library(tidyverse)
library(ggbeeswarm)

theme_step1 <- function(base_size = 11, base_family = "",
                        base_line_size = base_size / 22,
                        base_rect_size = base_size / 22) {
  theme(title = element_text(family = 'Arial', size = 18, color = 'black'), text = element_text(family = 'Arial', size = 16, color = 'black'),
        axis.title = element_text(family = 'Arial', size = 18, color = 'black'), axis.text = element_text(family = 'Arial', size = 16, color = 'black'), 
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = NA), axis.line = element_line(colour = "black", size = rel(1)),
        legend.background = element_rect(color = 'black'), legend.title = element_text(family = 'Arial', size = 16),
        legend.text = element_text(family = 'Arial', size = 14),
        legend.direction = "vertical", 
        legend.box = c("horizontal", "vertical"),
        legend.spacing.x = unit(0.1, 'cm'),
        plot.margin = unit(c(0.25, 1, 1, 0.5), 'cm'),
        axis.title.y = element_text(margin = margin(r = 10, unit = "pt"))) }


setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/99.kogo/2026.winter.kogo")

hprc_sup_sex <- readxl::read_xlsx("~/Desktop/KCDC/paper/pangenome/Draft human pangenome reference41586_2023_5896_MOESM4_ESM.xlsx",sheet = 7,skip = 1)
head(hprc_sup_sex)
hprc_sup_sex %>% select(sample,Sex) %>% unique() -> hprc_sup_sex

kppd_sex <- read_table("~/Desktop/KCDC/pangenome/00.datacheck/KPP_draft_assembly_matchingID.with.Gender.20250224.txt")
head(kppd_sex)
kppd_sex %>% mutate(Sex = ifelse(sex =="M","male","female")) %>% mutate(sample = KPPD_ID) %>% select(sample,Sex) -> kppd_sex

kppd_sex
hprc_sup_sex

kppd_gene <- read_table("KPPD.gene.annotation.byhap.txt")
kppd_trans <- read_table("KPPD.transcript.annotation.byhap.txt")

hprc_gene <- read_table("HPRC_genecode.v.38.gene.protein_coding_count.txt")
hprc_trans <- read_table("HPRC_genecode.v.38.trans.protein_coding_count.txt")


head(kppd_gene)
head(kppd_trans)

head(hprc_gene)
head(kppd_trans)
table(hprc_gene$ID) %>% dim()/2

kppd_gene %>% rbind(kppd_trans) %>% mutate(cohort = "KPPD") %>%
  rbind(hprc_gene %>% rbind(hprc_trans) %>% mutate(cohort = "HPRC")) %>%
  mutate(
    total = case_when(
      coding_type == "Protein_coding_genes"  ~ 19959,
      coding_type == "Noncoding_genes"       ~ 39509,
      coding_type == "Protein_coding_trans"  ~ 86728,
      coding_type == "Noncoding_trans"       ~ 150192,
      TRUE ~ NA_real_
    ),
    annotation_pct = (n / total) * 100
  ) %>% mutate(sample = str_split_fixed(ID,"_",2)[,1]) %>% left_join(kppd_sex %>% rbind(hprc_sup_sex)) %>%
  mutate(coding_type = factor(coding_type, levels = c("Protein_coding_genes", "Noncoding_genes","Protein_coding_trans","Noncoding_trans"))) -> liftoff_QC_forannotation

liftoff_QC_forannotation %>% 
  group_by(cohort,coding_type) %>%
  summarise(mean(annotation_pct),sd(annotation_pct)) 


liftoff_QC_forannotation %>% filter(cohort == "KPPD") %>%
  filter(coding_type == "Protein_coding_genes") %>%
  group_by(cohort,Sex,coding_type) %>%
  summarise(mean(annotation_pct)) 


liftoff_QC_forannotation %>%
  mutate(cohort = factor(cohort, levels = c("KPPD", "HPRC"))) %>%
  ggplot(aes(x = coding_type,y = annotation_pct,group = cohort,shape = Sex,fill  = cohort)) +
  geom_quasirandom(
    alpha = 0.8,dodge.width = 0.6,size = 2,colour = "black"   # ✅ 테두리만 검정
  ) +
  labs(y = "Percentage annotated (%)") +
  scale_x_discrete(
    labels = c(
      "Protein_coding_genes" = "Protein-coding\ngenes",
      "Noncoding_genes"      = "Noncoding\ngenes",
      "Protein_coding_trans" = "Protein-coding\ntranscripts",
      "Noncoding_trans"      = "Noncoding\ntranscripts"
    )
  ) +
  scale_shape_manual(
    values = c(
      "female" = 24,  # triangle
      "male"   = 22   # square
    ),
    breaks = c("male", "female"),
    labels = c("Male", "Female")
  ) +
  ## ✅ cohort 색을 fill로 강제 (KPPD 먼저)
  scale_fill_manual(
    values = c("KPPD" = "#1F4E9E", "HPRC" = "#C00000"),
    breaks = c("KPPD", "HPRC")
  ) +
  theme_step1() +
  guides(
    ## cohort legend(색): alpha 영향 없이 진하게
    fill = guide_legend(override.aes = list(alpha = 1, shape = 21, colour = "black", size = 3)),
    ## Sex legend(모양): alpha 영향 없이 진하게
    shape = guide_legend(override.aes = list(alpha = 1, fill = "grey80", colour = "black", size = 3))
  ) +
  theme(
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.background     = element_blank(),
    legend.box.background = element_blank(),
    legend.key            = element_blank()
  )



######## CNV 
kppd_cnv_sample <- read_table("KPPD.liftoff.CNV.bysample.txt")
hprc_cnv_sample <- read_table("HPRC.liftoff.CNV.bysample.txt")
cpc_cnv_smaple <- read_table("CPC.liftoff.CNV.bysample.txt")

kppd_cnv_sample$cohort <- "KPPD"
hprc_cnv_sample$cohort <- "HPRC"
cpc_cnv_smaple$cohort <- "CPC"


kppd_cnv_sample %>% rbind(hprc_cnv_sample) %>% rbind(cpc_cnv_smaple) %>% 
  group_by(cohort) %>% summarise(mean(n))

kppd_cnv_sample %>% rbind(hprc_cnv_sample) %>% rbind(cpc_cnv_smaple) %>%
  mutate(cohort = factor(cohort, levels = c("KPPD", "HPRC", "CPC"))) %>%
  ggplot(aes(x = cohort, y = n, color = cohort)) +
  geom_boxplot(outlier.shape = NA) +   # outlier 점 숨김
  geom_quasirandom(
    shape = 23,alpha = 0.8,dodge.width = 0.8,size = 2) +
  scale_color_manual(
    values = c("KPPD" = "#1F4E9E","HPRC" = "#C00000","CPC"  = "#2E7D32")) +
  labs(y = "Duplicated genes for assemblies",x="\n") +
  theme_step1() +
  theme(
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    strip.text.x = element_text(
      size = 14,
      face = "bold",
      family = "Arial")) -> p1

p1
kppd_cnv <- read_table("KPPD.liftoff.CNV.byGene.txt")
hprc_cnv <- read_table("HPRC.liftoff.CNV.byGene.txt")
cpc_cnv <- read_table("CPC.liftoff.CNV.byGene.txt")
head(cpc_cnv)
#kppd_cnv %>% rename(KPPD = n)  -> kppd_cnv
#hprc_cnv %>% rename(HPRC = n) -> hprc_cnv
#cpc_cnv %>% rename(CPC = n) -> cpc_cnv


gene_sets <- list(
  KPPD = kppd_cnv$gene_name,
  HPRC = hprc_cnv$gene_name,
  CPC  = cpc_cnv$gene_name
)


library(VennDiagram)
library(grid)


names(gene_sets)

p2 <- venn.diagram(
  x = gene_sets,
  filename = NULL,              # 화면 출력
  
  category.names = names(gene_sets),
  
  ## 채움색
  fill  = c("blue", "red", "green"),
  alpha = 0.2,
  
  ## 선
  lwd = 1,
  
  ## 🔹 교집합 숫자
  cex = 1.2,                     # 숫자 크기
  fontfamily = "Arial",          # 숫자 폰트
  fontface = "plain",            # "plain", "bold", "italic"
  
  ## 🔹 집합 이름 (KPPD, HPRC, CPC 등)
  cat.cex = 1.2,
  cat.fontfamily = "Arial",
  cat.fontface = "plain",
  
  ## 위치
  cat.pos  = c(-10, -1, 1),
  cat.dist = c(0.05, 0.05, -0.45),
  
  margin = 0.1
)

grid.newpage()
grid.draw(p2)



##
kppd_cnv %>%
  arrange(desc(n)) %>%
  mutate(rank = row_number(),
         prop = n / 264,
         prop_bin = case_when(
           n == 1                      ~ "singleton\n(n=1)",
           prop < 0.01                 ~ "0~0.01",
           prop < 0.05                 ~ "0.01~0.05",
           prop < 0.10                 ~ "0.05~0.10",
           prop < 0.20                 ~ "0.10~0.20",
           prop < 0.30                 ~ "0.20~0.30",
           prop < 0.40                 ~ "0.30~0.40",
           prop < 0.50                 ~ "0.40~0.50",
           prop <= 1                   ~ "0.50~1.00",
           TRUE                        ~ NA_character_),
         prop_bin = factor(
           prop_bin,levels = c("singleton\n(n=1)",
                               "0~0.01","0.01~0.05","0.05~0.10","0.10~0.20",
                               "0.20~0.30","0.30~0.40","0.40~0.50","0.50~1.00"))) %>% 
  left_join(cpc_cnv %>% rename(cpc=n)) %>% left_join(hprc_cnv %>% rename(hprc=n)) %>%
  mutate( category = case_when(
    !is.na(n) &  is.na(cpc) &  is.na(hprc) ~ "KPPD",
    !is.na(n) & !is.na(cpc) &  is.na(hprc) ~ "KPPD ∩ CPC",
    !is.na(n) &  is.na(cpc) & !is.na(hprc) ~ "KPPD ∩ HPRC",
    !is.na(n) & !is.na(cpc) & !is.na(hprc) ~ "ALL",
    TRUE ~ NA_character_),
    category = factor(category, levels = c("ALL", "KPPD ∩ HPRC", "KPPD ∩ CPC","KPPD"))) %>% #head()
  count(prop_bin,category) %>% #head()
  ggplot(aes(x = prop_bin, y = n,fill=category)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene copy number",
       y = "Frequency") +
  scale_fill_discrete(limits = c("KPPD", "KPPD ∩ CPC", "KPPD ∩ HPRC", "ALL")) + 
  scale_color_manual(
    values = c("KPPD" = "#1F4E9E","KPPD ∩ HPRC" = "#C00000","KPPD ∩ CPC"  = "#2E7D32","ALL" = "#8338EC")) +
  theme_step1() + 
  theme(legend.title = element_blank(),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = "white",color = NA,alpha = 0.7)) -> p3

p3
kppd_cnv %>%
  arrange(desc(n)) %>%
  mutate(rank = row_number(),
         prop = n / 264,
         prop_bin = case_when(
           n == 1                      ~ "singleton (n=1)",
           prop < 0.01                 ~ "0~0.01",
           prop < 0.05                 ~ "0.01~0.05",
           prop <= 1                   ~ "0.05~1.00",
           TRUE                        ~ NA_character_),
         prop_bin = factor(
           prop_bin,levels = c("singleton (n=1)",
                               "0~0.01","0.01~0.05","0.05~1.00"))) %>% 
  left_join(cpc_cnv %>% rename(cpc=n)) %>% left_join(hprc_cnv %>% rename(hprc=n)) %>%
  mutate( category = case_when(
    !is.na(n) &  is.na(cpc) &  is.na(hprc) ~ "KPPD",
    !is.na(n) & !is.na(cpc) &  is.na(hprc) ~ "KPPD ∩ CPC",
    !is.na(n) &  is.na(cpc) & !is.na(hprc) ~ "KPPD ∩ HPRC",
    !is.na(n) & !is.na(cpc) & !is.na(hprc) ~ "ALL",
    TRUE ~ NA_character_),
    category = factor(category, levels = c("ALL", "KPPD ∩ HPRC", "KPPD ∩ CPC","KPPD"))) %>% #head()
  count(prop_bin,category) -> p3_data
p3_data %>%
  ggplot(aes(x = prop_bin, y = n,fill=category)) +
  geom_bar(stat = "identity") +
  labs(x = "Frequency of KPPD CNV",
       y = "Number of duplicated genes") +
  scale_fill_discrete(limits = c("KPPD", "KPPD ∩ CPC", "KPPD ∩ HPRC", "ALL")) + 
  scale_fill_manual(
    values = c("KPPD" = "#1F4E9E","KPPD ∩ HPRC" = "#C00000","KPPD ∩ CPC"  = "#2E7D32","ALL" = "#8338EC")) +
  theme_step1() + 
  theme(legend.title = element_blank(),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = "white",color = NA,alpha = 0.7)) -> p3
p3


library(cowplot)

p <- plot_grid(
  p1, p2, p3,
  nrow = 1,
  rel_widths = c(1, 1, 1.4),
  labels = c("A", "B", "C"),
  label_size = 32,
  label_fontface = "bold"
)
p


p3_data 
831/1001

merge_freq %>% filter(is.na(HPRC),is.na(CPC)) %>% select(gene_name) %>%
  write.table("KPPD.unique.gene.1001.txt",col.names = F,row.names = F,quote = F)
####### last prequency
kppd_cnv %>% mutate(KPPD = n/264) %>% 
  left_join(hprc_cnv %>% mutate(HPRC = n/94) %>% select(-n)) %>%
  left_join(cpc_cnv %>% mutate(CPC = n/116) %>% select(-n)) -> merge_freq
library(ggplot2)
library(smplot2)
merge_freq

merge_freq %>% select(gene_name,KPPD,HPRC) %>%
  na.omit() %>% #dim() #662
  ggplot(aes(x=HPRC,y=KPPD)) + 
  geom_point(shape = 21, fill = "#C00000", size = 2,alpha = 0.8) +
  #geom_smooth(method = "lm", se = FALSE, color = "0f993d", linetype = "solid", size = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray40", size = 0.8) +  # ??? y=x ??��??
  sm_statCorr(
    color = "#C00000", corr_method = "spearman",
    linetype = "dashed",
    text_size = 8
  ) +
  annotate("text",x = 1, y = 0,label = "KPPD ∩ HPRC : 662 CNVs",hjust = 1, vjust = 0,size = 5,color = "black" ) +
  ylim(0,1) + xlim(0,1) + 
  theme_step1() -> p1



merge_freq %>% select(gene_name,KPPD,CPC) %>%
  na.omit() %>% #dim() #738
  ggplot(aes(x=CPC,y=KPPD)) + 
  geom_point(shape = 21, fill = "#2E7D32", size = 2,alpha = 0.8) +
  #geom_smooth(method = "lm", se = FALSE, color = "0f993d", linetype = "solid", size = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray40", size = 0.8) +  # ??? y=x ??��??
  sm_statCorr(
    color = "#2E7D32", corr_method = "spearman",
    linetype = "dashed",
    text_size = 8
  ) +
  ylim(0,1) + xlim(0,1) + 
  theme_step1() + 
  annotate("text",x = 1, y = 0,label = "KPPD ∩ CPC : 738 CNVs",hjust = 1, vjust = 0,size = 5,color = "black" ) +
  theme(axis.title.y=element_blank())-> p2

plot_grid(p1, p2, rel_widths = c(1,1),nrow = 1)
  
cpc_cnv_smaple
max(cpc_cnv$n)

merge_freq %>% arrange(-KPPD) %>%
  mutate(KPPD_CPC = KPPD-CPC) %>%
  mutate(KPPD_HPRC = KPPD-HPRC) %>% #head()
  filter(KPPD > 0.1) %>% filter(HPRC < 0.1) %>% filter(KPPD_HPRC > 0.05) %>%
  arrange(-KPPD_HPRC) %>%  #head()
  slice_head(n = 5) %>% select(gene_name,KPPD,HPRC,KPPD_HPRC) -> KPPD_HPRC


merge_freq %>% arrange(-KPPD) %>%
  mutate(KPPD_CPC = KPPD-CPC) %>%
  mutate(KPPD_HPRC = KPPD-HPRC) %>%
  filter(KPPD > 0.1) %>% filter(CPC < 0.1) %>% filter(KPPD_CPC > 0.05) %>%
  arrange(-KPPD_CPC) %>% slice_head(n = 5) %>% select(gene_name,KPPD,CPC,KPPD_CPC) -> KPPD_CPC


KPPD_HPRC %>% #select(gene_name,KPPD,HPRC) %>% 
#  na.omit() %>% 
  pivot_longer(
    cols = c(KPPD, HPRC),
    names_to = "cohort",
    values_to = "freq"
  ) %>% 
  mutate(
    freq_signed = ifelse(cohort == "HPRC", -freq, freq)
  ) %>% 
  group_by(gene_name) %>%
  mutate(order_key = max(freq_signed)) %>%
  ungroup() %>%
  mutate(gene_name = reorder(gene_name, order_key)) %>%  
  ggplot(aes(x = freq_signed, y = gene_name, fill = cohort)) +
  ## bar
  geom_col(width = 0.7) +
  ## 기준선
  geom_vline(xintercept = 0, linewidth = 0.5) +
  ## 값 텍스트
  geom_text(
    aes(
      label = sprintf("%.3f", freq),
      hjust = ifelse(cohort == "HPRC", 1.05, -0.05)
    ),size = 5) +
  
  ## x축 (절댓값 표기)
  scale_x_continuous(
    labels = function(x) abs(x),
    expand = expansion(mult = c(0.10, 0.10))
  ) +
  xlim(-1,1) + 
  ## 색상 고정
  scale_fill_manual(values = c("KPPD" = "#1F4E9E","HPRC" = "#C00000")) +
  labs(x = "CNV Frequency", y = NULL, fill = NULL) +
  theme_step1() +
  theme(
    legend.position = c(0.98, 0.1),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", colour = "white"),
    ## box 제거 핵심
    panel.border = element_blank(),
    panel.background = element_blank(),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank())


KPPD_CPC %>%
  select(gene_name, KPPD, CPC) %>%
  pivot_longer(cols = c(KPPD, CPC),names_to = "cohort",values_to = "freq") %>%
  mutate(freq_signed = ifelse(cohort == "CPC", -freq, freq)) %>%
  group_by(gene_name) %>%
  mutate(order_key = max(freq_signed)) %>%
  ungroup() %>%
  mutate(gene_name = reorder(gene_name, order_key)) %>%
  ggplot(aes(x = freq_signed, y = gene_name, fill = cohort)) +
  ## bar
  geom_col(width = 0.7) +
  ## 기준선
  geom_vline(xintercept = 0, linewidth = 0.5) +
  ## 값 텍스트
  geom_text(aes(label = sprintf("%.3f", freq),hjust = ifelse(cohort == "CPC", 1.05, -0.05)),size = 5) +
  ## x축 (절댓값 표기)
  scale_x_continuous(
    labels = function(x) abs(x),
    expand = expansion(mult = c(0.10, 0.10))) +
  xlim(-1,1) +
  ## 색상 고정
  scale_fill_manual(values = c("KPPD" = "#1F4E9E","CPC" = "#2E7D32")) +
  labs(x = "CNV Frequency", y = NULL, fill = NULL) +
  theme_step1() +
  theme(
    legend.position = c(0.98, 0.1),
    legend.justification = c("right", "bottom"),
    legend.box.background = element_blank(),
    ## box 제거 핵심
    panel.border = element_blank(),
    panel.background = element_blank(),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank())




KPPD_HPRC %>% #select(gene_name,KPPD,HPRC) %>% 
#  na.omit() %>% 
  pivot_longer(
    cols = c(KPPD, HPRC),
    names_to = "cohort",
    values_to = "freq"
  ) %>% 
  mutate(
    freq_signed = ifelse(cohort == "HPRC", -freq, freq)
  ) %>% 
  group_by(gene_name) %>%
  mutate(order_key = max(freq_signed)) %>%
  ungroup() %>%
  mutate(gene_name = reorder(gene_name, order_key)) %>%  
  ggplot(aes(x = freq_signed, y = gene_name, fill = cohort)) +
  ## bar
  geom_col(width = 0.7) +
  ## 기준선
  geom_vline(xintercept = 0, linewidth = 0.5) +
  ## 값 텍스트
  geom_text(
    aes(
      label = sprintf("%.3f", freq),
      hjust = ifelse(cohort == "HPRC", 1.05, -0.05)
    ),size = 5) +
  
  ## x축 (절댓값 표기)
  scale_x_continuous(
    labels = function(x) abs(x),
    expand = expansion(mult = c(0.10, 0.10))
  ) +
  xlim(-1,1) + 
  ## 색상 고정
  scale_fill_manual(values = c("KPPD" = "#1F4E9E","HPRC" = "#C00000")) +
  labs(x = "CNV Frequency", y = NULL, fill = NULL) +
  theme_step1() +
  theme(
    legend.position = c(0.98, 0.1),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", colour = "white"),
    ## box 제거 핵심
    panel.border = element_blank(),
    panel.background = element_blank(),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank())


KPPD_CPC %>%
  select(gene_name, KPPD, CPC) %>%
  pivot_longer(cols = c(KPPD, CPC),names_to = "cohort",values_to = "freq") %>%
  mutate(freq_signed = ifelse(cohort == "CPC", -freq, freq)) %>%
  group_by(gene_name) %>%
  mutate(order_key = max(freq_signed)) %>%
  ungroup() %>%
  mutate(gene_name = reorder(gene_name, order_key)) %>%
  ggplot(aes(x = freq_signed, y = gene_name, fill = cohort)) +
  ## bar
  geom_col(width = 0.7) +
  ## 기준선
  geom_vline(xintercept = 0, linewidth = 0.5) +
  ## 값 텍스트
  geom_text(aes(label = sprintf("%.3f", freq),hjust = ifelse(cohort == "CPC", 1.05, -0.05)),size = 5) +
  ## x축 (절댓값 표기)
  scale_x_continuous(
    labels = function(x) abs(x),
    expand = expansion(mult = c(0.10, 0.10))) +
  xlim(-1,1) +
  ## 색상 고정
  scale_fill_manual(values = c("KPPD" = "#1F4E9E","CPC" = "#2E7D32")) +
  labs(x = "CNV Frequency", y = NULL, fill = NULL) +
  theme_step1() +
  theme(
    legend.position = c(0.98, 0.1),
    legend.justification = c("right", "bottom"),
    legend.box.background = element_blank(),
    ## box 제거 핵심
    panel.border = element_blank(),
    panel.background = element_blank(),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank())




library(tidyverse)
library(cowplot)

# 1) 공통 팔레트/순서 (legend에 3개가 항상 나오게)
cohort_cols   <- c("KPPD"="#1F4E9E", "HPRC"="#C00000", "CPC"="#2E7D32")
cohort_breaks <- c("KPPD","HPRC","CPC")

# 2) p1 / p2: legend 제거 + scale_fill_manual을 3개로 통일(drop=FALSE)
p1 <- KPPD_HPRC %>%
  pivot_longer(cols = c(KPPD, HPRC), names_to="cohort", values_to="freq") %>%
  mutate(freq_signed = ifelse(cohort=="HPRC", -freq, freq)) %>%
  group_by(gene_name) %>% mutate(order_key = max(freq_signed)) %>% ungroup() %>%
  mutate(gene_name = reorder(gene_name, order_key),
         cohort = factor(cohort, levels = cohort_breaks)) %>%
  ggplot(aes(x=freq_signed, y=gene_name, fill=cohort)) +
  geom_col(width=0.7) +
  geom_vline(xintercept=0, linewidth=0.5) +
  geom_text(aes(label=sprintf("%.3f", freq),
                hjust=ifelse(cohort=="HPRC", 1.05, -0.05)), size=5) +
  scale_x_continuous(labels=function(x) abs(x),
                     expand=expansion(mult=c(0.10,0.10))) +
  xlim(-1,1) +
  scale_fill_manual(values=cohort_cols, breaks=cohort_breaks, drop=FALSE) +
  labs(x="CNV Frequency", y=NULL, fill=NULL) +
  theme_step1() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

p2 <- KPPD_CPC %>%
  select(gene_name, KPPD, CPC) %>%
  pivot_longer(cols = c(KPPD, CPC), names_to="cohort", values_to="freq") %>%
  mutate(freq_signed = ifelse(cohort=="CPC", -freq, freq)) %>%
  group_by(gene_name) %>% mutate(order_key = max(freq_signed)) %>% ungroup() %>%
  mutate(gene_name = reorder(gene_name, order_key),
         cohort = factor(cohort, levels = cohort_breaks)) %>%
  ggplot(aes(x=freq_signed, y=gene_name, fill=cohort)) +
  geom_col(width=0.7) +
  geom_vline(xintercept=0, linewidth=0.5) +
  geom_text(aes(label=sprintf("%.3f", freq),
                hjust=ifelse(cohort=="CPC", 1.05, -0.05)), size=5) +
  scale_x_continuous(labels=function(x) abs(x),
                     expand=expansion(mult=c(0.10,0.10))) +
  xlim(-1,1) +
  scale_fill_manual(values=cohort_cols, breaks=cohort_breaks, drop=FALSE) +
  labs(x="CNV Frequency", y=NULL, fill=NULL) +
  theme_step1() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# 3) legend만 추출: 3개 cohort를 강제로 포함시키는 dummy 데이터 사용
legend_dummy <- tibble(
  x = 1:3, y = 1:3,
  cohort = factor(cohort_breaks, levels = cohort_breaks)
)

legend_plot <- ggplot(legend_dummy, aes(x, y, fill = cohort)) +
  geom_col() +
  scale_fill_manual(values=cohort_cols, breaks=cohort_breaks, drop=FALSE) +
  theme_void() +
  theme_step1() + 
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left"
  )

legend_plot
leg <- cowplot::get_legend(legend_plot)
leg

# 4) cowplot으로 결합: 위(p1) / 아래(p2) / 맨 아래 legend(왼쪽)
main <- plot_grid(p1, p2, nrow = 1, align = "v", axis = "lr")

plot_grid(p1, p2, leg, rel_widths = c(1, 1,0.1),nrow = 1, align = "v", axis = "lr")
plot_grid(main, leg, rel_widths = c(1,0.1),nrow = 1)

final <- plot_grid(
  main,
  leg,
  nrow = 1,
  rel_heights = c(1, 0.12)  # legend 높이 조절
)

final





p <- plot_grid(
  p1, p2, #p3,
  nrow = 1,
  rel_widths = c(1, 1),
  labels = c("A", "B"),
  label_size = 32,
  label_fontface = "bold"
)
p
