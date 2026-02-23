## liftoff final data pro
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


setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/liftoff/with_flagger/99.summary/")

gene_anno <- read_table("gene.annotation.byhap.txt")
trans_anno <- read_table("transcript.annotation.byhap.txt")
head(gene_anno)
head(trans_anno)
#gene coding = 19959, non= 19983
#trasn 86728 150192
trans_anno %>% mutate(pct = ifelse(coding_type == "Noncoding_trans",n/150192,n/86728)) %>% count(pct > 1)
gene_anno %>% mutate(pct = ifelse(coding_type == "Noncoding_genes",n/39483,n/19959)) %>% count(pct > 1)


trans_anno %>% mutate(pct = ifelse(coding_type == "Noncoding_trans",n/150192,n/86728)) %>% #head()
  rbind(gene_anno %>% mutate(pct = ifelse(coding_type == "Noncoding_genes",n/39483,n/19959))) %>% #head()
  mutate(pct = round(pct*100,2)) %>% select(-n) %>%
  pivot_wider(names_from = coding_type,values_from = pct) %>% 
  mutate(File = str_replace_all(ID,"h","hap")) %>% 
  mutate(Sample = str_split_fixed(File,"_",2)[,1]) %>%
  mutate(Haplotype = str_split_fixed(File,"_",2)[,2]) %>%
  mutate(Batch = ifelse(Sample %in% c("KPPD129","KPPD130","KPPD131","KPPD132"),"T2T","Non.T2T")) %>% 
  select(File,Sample,Haplotype,Batch,Protein_coding_genes,Noncoding_genes,Protein_coding_trans,Noncoding_trans) %>%
  writexl::write_xlsx("liftoff_annotation.gencodev38.xlsx")
  



multi_exon <- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/liftoff/genocode.v38.multiexon_genes.tsv") %>% filter(max_exon_count > 1)
head(multi_exon)
multi_exon %>% filter(gene_name == "SPATA13")
multi_exon %>%
  group_by(gene_name) %>%
  filter(
    if (any(gene_type == "protein_coding")) {
      gene_type == "protein_coding"
    } else {
      TRUE
    }
  ) %>% ungroup() -> multi_exon

kppd_bysample <- read_table("liftoff.CNV.bysample.txt")
kppd_bysample %>% 
  mutate(Sample = str_split_fixed(ID,"_",2)[,1]) %>%
  mutate(Haplotype = str_split_fixed(ID,"_",2)[,2]) %>% 
  mutate(Haplotype = str_replace(Haplotype,"h","hap")) %>% 
  rename(Number_of_Extra_copy_number_gene = n) %>%
  select(Sample,Haplotype,Number_of_Extra_copy_number_gene) %>%
  writexl::write_xlsx("liftoff.CNV.bysample.forTable.xlsx")

head(kppd_bysample)
kppd_bysample %>% summarise(mean = mean(n),median = median(n),max  = max(n), min  = min(n))
#1  49.8     49    83    30
kppd_bysample %>%
  mutate(sampleID = str_split_fixed(ID, "_", 2)[,1]) %>%
  mutate(cohort = ifelse(
    sampleID %in% c("KPPD132","KPPD131","KPPD130","KPPD129"),
    "T2T", "nonT2T"
  )) %>% 
  ggplot(aes(x="",y=n)) +
  #geom_violin(fill = "steelblue",alpha = 0.35,color = NA) +
  geom_boxplot() +
  geom_quasirandom(shape = 23, alpha = 0.8, dodge.width = 0.8,size=2,fill = "steelblue") +
  labs(y="Duplicated genes for genome") + 
  #geom_hline(yintercept = 54,color = "blue",linetype = "dashed",linewidth = 1) +
  geom_hline(yintercept = 36,color = "red",linetype = "dashed",linewidth = 1) +
  geom_hline(yintercept = 53, color = "green4",linetype = "dashed",linewidth = 1) +
  #annotate("text", x = 0.4, y = 50, label = "KPPD:49",color = "blue", hjust = 0, size = 3) +
  #annotate("text", x = 1.3, y = 35, label = "HPRC:36",color = "red", hjust = 0, size = 3) +
  annotate("text", x = 1.3, y = 35, label = "HPRC",color = "red", hjust = 0, size = 4) +
  annotate("text", x = 1.4, y = 52, label = "CPC",color = "green4", hjust = 0, size = 4) + 
  #facet_wrap(~ref,scales = "free_x") + 
  theme_step1() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(), 
        legend.position = 'none', 
        strip.text.x = element_text(
          size = 14,        # ??Ʈ ũ??
          face = "bold",    # ??Ʈ ???? (bold, italic ??)
          family = "Arial"))
# 36 33

kppd_bysample_plot <- kppd_bysample %>%
  mutate(sampleID = str_split_fixed(ID, "_", 2)[,1]) %>%
  mutate(cohort = ifelse(
    sampleID %in% c("KPPD132","KPPD131","KPPD130","KPPD129"),
    "T2T", "nonT2T"
  )) %>%
  mutate(cohort = factor(cohort, levels = c("nonT2T","T2T"))) %>%  # 핵심
  arrange(cohort, n) %>%
  mutate(ID = factor(ID, levels = unique(ID)))

kppd_bysample_plot %>% 
  summarise(mean = mean(n),median = median(n),max  = max(n), min  = min(n))
kppd_bysample_plot %>% group_by(cohort) %>%
  summarise(mean = mean(n),median = median(n),max  = max(n), min  = min(n))

kppd_bysample_plot
kppd_bysample_plot %>% group_by(cohort) %>%
  summarise(mean = mean(n),median = median(n),max  = max(n), min  = min(n))

ggplot(kppd_bysample_plot, aes(x = ID, y = n, fill = cohort)) +
  geom_bar(stat = "identity") +
  facet_grid(
    ~ cohort,
    scales = "free_x",
    space  = "free_x",
    switch = "x"          # ← 핵심: strip을 아래쪽으로
  ) +
  labs(y = "Duplicated genes for Assembly") +
  theme_step1() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    
    strip.background = element_blank(),
    strip.placement  = "outside",  # ← 패널 바깥으로
    strip.text.x     = element_text(size = 12, face = "bold")
  )


kppd_byGene <- read_table("liftoff.CNV.byGene.txt")
head(kppd_byGene)
dim(kppd_byGene)

##kppd_byGene %>% arrange(-n) %>% rename(Gene = gene_name,Frequency = n) %>%  writexl::write_xlsx("liftoff.CNV.byGene.forTable.xlsx")

summary(kppd_byGene$n)

kppd_byGene %>% 
  count(n) 
  

kppd_byGene %>% 
  count(n) %>%
  ggplot(aes(x=n,y=nn)) +
  geom_bar(stat = "identity")


kppd_byGene %>%
  mutate(n_bin = ifelse(n >= 15, "15+", as.character(n))) %>%
  count(n_bin) %>%
  mutate(
    n_bin = factor(
      n_bin,
      levels = c(as.character(1:14), "15+")
    )
  ) %>% 
  ggplot(aes(x = n_bin, y = n,fill="#F6B7A4",alpha = 1.5)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene copy number",
       y = "Frequency") +
  theme_step1() + 
  theme(legend.position = "none")

kppd_byGene %>% arrange(-n) %>% 
  mutate(rank = row_number()) %>% 
  mutate(prop = n/264) %>% head()
  filter(rank %in% 10:1) %>%
  mutate(gene_name = factor(gene_name, levels = gene_name)) %>%
  ggplot(aes(x=prop,y=gene_name)) + 
  geom_bar(stat = "identity")

kppd_byGene %>%
    arrange(desc(n)) %>%
    mutate(rank = row_number(),
      prop = n / 264,
      prop_bin = case_when(
        n == 1                      ~ "singleton (n=1)",
        prop < 0.01                 ~ "0–0.01",
        prop < 0.05                 ~ "0.01–0.05",
        prop < 0.10                 ~ "0.05–0.10",
        prop < 0.20                 ~ "0.10–0.20",
        prop < 0.30                 ~ "0.20–0.30",
        prop < 0.40                 ~ "0.30–0.40",
        prop < 0.50                 ~ "0.40–0.50",
        prop <= 1                   ~ "0.50–1.00",
        TRUE                        ~ NA_character_
      ),
      prop_bin = factor(
        prop_bin,levels = c("singleton (n=1)",
          "0–0.01","0.01–0.05","0.05–0.10","0.10–0.20",
          "0.20–0.30","0.30–0.40","0.40–0.50","0.50–1.00"))) %>% 
    count(prop_bin) %>% 
    ggplot(aes(x=prop_bin,y=n)) +
    geom_bar(stat = "identity")


head(kppd_byGene)

hprc <-  readxl::read_xlsx("~/Desktop/KCDC/paper/pangenome/Draft human pangenome reference41586_2023_5896_MOESM4_ESM.xlsx",sheet = 9,skip = 1)
head(hprc)
hprc %>% pivot_longer(2:ncol(hprc)) %>% na.omit()  %>% #dim()
  filter(!name %in% c("GRCh38","CHM13_combined")) %>% #head()
  filter(!(value %in% c(0,1))) %>% unique() %>% count(name) %>% #head()
  summarise(mean = mean(n),max  = max(n), min  = min(n))

hprc %>% pivot_longer(2:ncol(hprc)) %>% na.omit()  %>% #dim()
  filter(!name %in% c("GRCh38","CHM13_combined")) %>% #head()
  filter(value != 0) %>% unique() %>% count(name) %>% #head()
  summarise(mean = mean(n),max  = max(n), min  = min(n))

hprc %>% 
  pivot_longer(2:ncol(hprc)) %>% na.omit()  %>% #dim()
  filter(!name %in% c("GRCh38","CHM13_combined")) %>% #head()
  filter((gene %in% mutli_exon_count$gene_name)) %>%
  filter(value > 0) %>% unique() %>% count(name) %>% #head()
  summarise(mean = mean(n),max  = max(n), min  = min(n))

  #filter((gene %in% mutli_exon_count$gene_name)) 


hprc %>% pivot_longer(2:ncol(hprc)) %>% na.omit()  %>% #head()
  filter(!name %in% c("GRCh38","CHM13_combined")) %>% #head()
  mutate(sample = str_split_fixed(name,"_",2)[,1]) %>%
  mutate(hap = str_split_fixed(name,"_",2)[,2]) %>% mutate(hap = paste0("hap",hap)) %>%
  select(-name) %>%
  pivot_wider(names_from = hap,values_from = value) %>% #head()
  filter(hap1 > 0, hap2 > 0) %>% count(sample) %>%
  summarise(mean = mean(n),max  = max(n), min  = min(n))

hprc  
hprc %>% pivot_longer(2:ncol(hprc)) %>% na.omit()  %>% #dim()
  filter(!name %in% c("GRCh38","CHM13_combined")) %>% #head()
  filter(value != 0) %>% count(gene) %>% filter(n == 1) %>% dim()/1115


788/1115
  

arab <- readxl::read_xlsx("~/Desktop/KCDC/paper/pangenome/arab_suptable41467_2025_61645_MOESM2_ESM.xlsx",skip = 1,sheet = 20)
head(arab)
dim(arab)
arab %>% pivot_longer(2:ncol(arab)) %>% na.omit()  %>% #head()
  filter(value != 0) %>% #unique() %>% 
  count(name) %>% #head()
  summarise(mean = mean(n),max  = max(n), min  = min(n))
dim(arab)

cpc <- readxl::read_xlsx("~/Desktop/KCDC/paper/pangenome/CPC_cnv.xlsx",skip = 1)
head(cpc)
tail(cpc)
col(cpc)

cpc %>% select(Gene...1,`CPC\nCount (n=116)...2`) %>% rename(gene = Gene...1,count = `CPC\nCount (n=116)...2`) %>%
  rbind(cpc %>% select(Gene...7,`CPC\nCount (n=116)...8`) %>% rename(gene = Gene...7,count = `CPC\nCount (n=116)...8`)) %>% na.omit() -> cpc
head(cpc)

#cpc %>% write_tsv("../../db/cpc.gene.list.txt")
cpc <- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/liftoff/db/cpc.gene.list.txt")
cpc

hprc %>% mutate(Zero_count = "")
dim(hprc)
hprc %>%  na.omit()  %>% dim()
head(hprc)
ncol(hprc)/2
91 - 3

hprc %>% count(GRCh38)
hprc %>% filter(GRCh38 != 0) %>% pivot_longer(2:ncol(hprc)) %>% #na.omit()  %>% dim()
  filter(!name %in% c("CHM13_combined")) %>% #head()
  filter(value != 0) %>%
  group_by(gene) %>%
  mutate(
    GRCh38_ref = value[name == "GRCh38"][1]
  ) %>% ungroup() %>% filter(GRCh38_ref != value) %>% 
  mutate(check = value - GRCh38_ref) %>% 
  filter(check > 0)

hprc %>% pivot_longer(2:ncol(hprc)) %>% #na.omit()  %>% dim()
  filter(!name %in% c("CHM13_combined")) %>% count()
  filter(value != 0) %>% count(gene) %>% arrange(-n)


hprc %>% pivot_longer(2:ncol(hprc)) %>% #na.omit()  %>% dim()
  filter(!name %in% c("GRCh38","CHM13_combined")) %>% 
  filter(value != 0) %>% count(gene) %>% arrange(-n)





hprc %>% pivot_longer(2:ncol(hprc)) %>% #na.omit()  %>% dim()
  filter(!name %in% c("GRCh38","CHM13_combined")) %>% 
  filter(value != 0) %>%
  count(gene) -> hprc
head(hprc)
tail(hprc)
head(kppd_byGene)
head(hprc)
head(cpc)


kppd_byGene %>% rename(KPPD = "gene_name")
hprc %>% rename(HPRC = "gene")
cpc %>% rename(CPC = "gene")


library(ggVennDiagram)

gene_sets <- list(
  KPPD = kppd_byGene$gene_name,
  HPRC = hprc$gene,
  CPC  = cpc$gene
)


#kppd_byGene %>% filter(!(gene_name %in% hprc$gene)) %>% filter(!(gene_name %in% cpc$gene)) %>% select(gene_name) %>% write.table("liftoff.CNV.byGene.KPPDunique.txt",col.names = F,row.names = F,quote = F)
kppd_byGene %>% filter(!(gene_name %in% hprc$gene)) %>% filter(!(gene_name %in% cpc$gene)) %>% filter(n>1) %>% #dim()
  select(gene_name) %>% write.table("liftoff.CNV.byGene.KPPDunique.notsingleton.txt",col.names = F,row.names = F,quote = F)

library(VennDiagram)
library(grid)


names(gene_sets)
venn.plot <- venn.diagram(
  x = gene_sets,
  filename = NULL,              # 화면에 그릴 때는 NULL
  category.names = names(gene_sets),
  
  # 채움색(연하게)
  fill  = c("blue", "red", "green"),
  alpha = 0.2,
  
  # 원/라벨/숫자 스타일
  lwd = 1,
  cex = 1.2,                     # 영역 숫자 크기
  cat.cex = 1.2,                 # 집합 이름 크기
  cat.pos = c(-10, -1, 1),      # 집합 이름 위치(필요시 조정)
  cat.dist = c(0.05, 0.05, -0.45),
  
  margin = 0.1
)

grid.newpage()
grid.draw(venn.plot)


dim(hprc)
dim(cpc)
dim(kppd_byGene)
dim(arab)
head(hprc)
head(cpc)
hprc %>% count(n) %>% mutate(prop = prop.table(nn))
cpc %>% count(count) %>% mutate(prop = prop.table(n))
kppd_byGene %>% count(n) %>% mutate(prop = prop.table(nn))


library(ggVennDiagram)



kppd_byGene %>%
  arrange(desc(n)) %>%
  mutate(rank = row_number(),
         prop = n / 264,
         prop_bin = case_when(
           n == 1                      ~ "singleton (n=1)",
           prop < 0.01                 ~ "0–0.01",
           prop < 0.05                 ~ "0.01–0.05",
           prop < 0.10                 ~ "0.05–0.10",
           prop < 0.20                 ~ "0.10–0.20",
           prop < 0.30                 ~ "0.20–0.30",
           prop < 0.40                 ~ "0.30–0.40",
           prop < 0.50                 ~ "0.40–0.50",
           prop <= 1                   ~ "0.50–1.00",
           TRUE                        ~ NA_character_),
         prop_bin = factor(
           prop_bin,levels = c("singleton (n=1)",
                               "0–0.01","0.01–0.05","0.05–0.10","0.10–0.20",
                               "0.20–0.30","0.30–0.40","0.40–0.50","0.50–1.00"))) %>% 
  left_join(cpc %>% rename(gene_name = gene,cpc=count)) %>% left_join(hprc %>% rename(gene_name = gene,hprc=n)) %>%
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
  theme_step1() + 
  theme(legend.title = element_blank(),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = "white",color = NA,alpha = 0.7))


kppd_byGene %>%
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
  left_join(cpc %>% rename(gene_name = gene,cpc=count)) %>% left_join(hprc %>% rename(gene_name = gene,hprc=n)) %>%
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
  geom_bar(stat = "identity", width = 0.7) +
  labs(x = "Frequency of KPPD CNV",
       y = "Duplicated genes") +
  scale_fill_discrete(limits = c("KPPD", "KPPD ∩ CPC", "KPPD ∩ HPRC", "ALL")) + 
  scale_fill_manual(
    values = c("KPPD" = "#1F4E9E","KPPD ∩ HPRC" = "#C00000","KPPD ∩ CPC"  = "#2E7D32","ALL" = "#8338EC")) +
  theme_step1() + 
  theme(legend.title = element_blank(),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = "white",color = NA,alpha = 0.7)) 



####
head(kppd_byGene)
head(hprc)
head(cpc)

kppd_byGene %>% mutate(KPPD = n/264) %>% 
  left_join(hprc %>% mutate(HPRC = n/88) %>% select(-n) %>% rename(gene_name = gene)) %>% 
  left_join(cpc %>% mutate(CPC = count/116) %>% select(-count) %>% rename(gene_name = gene)) -> merge_freq

merge_freq

merge_freq %>% select(gene_name,KPPD,HPRC) %>%
  na.omit() %>% #dim() #255
  ggplot(aes(x=HPRC,y=KPPD)) + 
  geom_point(shape = 21, fill = "#C00000", size = 2,alpha = 0.8) +
  #geom_smooth(method = "lm", se = FALSE, color = "0f993d", linetype = "solid", size = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray40", size = 0.8) +  # ??? y=x ??��??
  sm_statCorr(
    color = "#C00000", corr_method = "spearman",
    linetype = "dashed",
    text_size = 5
  ) +
  annotate("text",x = 1, y = 0,label = "KPPD ∩ HPRC : 255 CNVs",hjust = 1, vjust = 0,size = 5,color = "black" ) +
  ylim(0,1) + xlim(0,1) + 
  theme_step1() -> p1



merge_freq %>% select(gene_name,KPPD,CPC) %>%
  na.omit() %>% #dim() #501
  ggplot(aes(x=CPC,y=KPPD)) + 
  geom_point(shape = 21, fill = "#2E7D32", size = 2,alpha = 0.8) +
  #geom_smooth(method = "lm", se = FALSE, color = "0f993d", linetype = "solid", size = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray40", size = 0.8) +  # ??? y=x ??��??
  sm_statCorr(
    color = "#2E7D32", corr_method = "spearman",
    linetype = "dashed",
    text_size = 5
  ) +
  ylim(0,1) + xlim(0,1) + 
  theme_step1() + 
  annotate("text",x = 1, y = 0,label = "KPPD ∩ CPC : 501 CNVs",hjust = 1, vjust = 0,size = 5,color = "black" ) +
  theme(axis.title.y=element_blank())-> p2

plot_grid(p1, p2, rel_widths = c(1,1),nrow = 1)





### multi exon

kppd_byGene_multi <- read_table("liftoff.CNV_multiexon.byGene.txt")
head(kppd_byGene_multi)
dim(kppd_byGene_multi)

multi_exon
hprc %>% left_join(multi_exon %>% rename(gene = gene_name)) %>% count(gene_type)
hprc %>% left_join(multi_exon %>% rename(gene = gene_name)) %>% filter(gene_type == "protein_coding") %>%
  select(gene,n) -> hprc_qc
hprc %>% left_join(multi_exon %>% rename(gene = gene_name)) %>% filter(gene_type != "protein_coding") %>% #head()
  select(gene,gene_type)

kppd_byGene %>% filter(gene_name %in% multi_exon$gene_name) -> kppd_byGene_multi_exon_qc
cpc %>% filter(gene %in% multi_exon$gene_name) -> cpc_multi_exon_qc

dim(kppd_byGene_multi_exon_qc)
hprc_qc %>% count(n) %>% mutate(prop = prop.table(nn))
cpc_multi_exon_qc %>% count(count) %>% mutate(prop = prop.table(n))
kppd_byGene_multi_exon_qc %>% count(n) %>% mutate(prop = prop.table(nn))


gene_sets <- list(
  KPPD = kppd_byGene_multi_exon_qc$gene_name,
  HPRC = hprc_qc$gene,
  CPC  = cpc_multi_exon_qc$gene
)


library(VennDiagram)
library(grid)


names(gene_sets)
venn.plot <- venn.diagram(
  x = gene_sets,
  filename = NULL,              # 화면에 그릴 때는 NULL
  category.names = names(gene_sets),
  
  # 채움색(연하게)
  fill  = c("blue", "red", "green"),
  alpha = 0.2,
  
  # 원/라벨/숫자 스타일
  lwd = 1,
  cex = 1.2,                     # 영역 숫자 크기
  cat.cex = 1.2,                 # 집합 이름 크기
  cat.pos = c(-10, -1, 1),      # 집합 이름 위치(필요시 조정)
  cat.dist = c(0.05, 0.05, -0.45),
  
  margin = 0.1
)

grid.newpage()
grid.draw(venn.plot)

head(kppd_byGene_multi_exon_qc)

hprc_qc

head(cpc)
hprc
hprc_qc
kppd_byGene
kppd_byGene %>%
  arrange(desc(n)) %>%
  mutate(rank = row_number(),
         prop = n / 264,
         prop_bin = case_when(
           n == 1                      ~ "singleton (n=1)",
           prop < 0.01                 ~ "0–0.01",
           prop < 0.05                 ~ "0.01–0.05",
           prop < 0.10                 ~ "0.05–0.10",
           prop < 0.20                 ~ "0.10–0.20",
           prop < 0.30                 ~ "0.20–0.30",
           prop < 0.40                 ~ "0.30–0.40",
           prop < 0.50                 ~ "0.40–0.50",
           prop <= 1                   ~ "0.50–1.00",
           TRUE                        ~ NA_character_),
         prop_bin = factor(
           prop_bin,levels = c("singleton (n=1)",
                               "0–0.01","0.01–0.05","0.05–0.10","0.10–0.20",
                               "0.20–0.30","0.30–0.40","0.40–0.50","0.50–1.00"))) %>% 
  left_join(cpc %>% rename(gene_name = gene,cpc=count)) %>% left_join(hprc %>% rename(gene_name = gene,hprc=n)) %>%
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
  theme_step1() + 
  theme(legend.title = element_blank(),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = "white",color = NA,alpha = 0.7))


kppd_byGene_multi_exon_qc %>% left_join(cpc_multi_exon_qc %>% rename(gene_name = gene,cpc = count)) %>%
  left_join(hprc_qc %>% rename(gene_name = gene,hprc = n)) %>%
  mutate( category = case_when(
      !is.na(n) &  is.na(cpc) &  is.na(hprc) ~ "KPPD",
      !is.na(n) & !is.na(cpc) &  is.na(hprc) ~ "KPPD ∩ CPC",
      !is.na(n) &  is.na(cpc) & !is.na(hprc) ~ "KPPD ∩ HPRC",
      !is.na(n) & !is.na(cpc) & !is.na(hprc) ~ "ALL",
      TRUE ~ NA_character_),
      category = factor(category, levels = c("ALL", "KPPD ∩ HPRC", "KPPD ∩ CPC","KPPD"))) %>%
  mutate(n_bin = ifelse(n >= 15, "15+", as.character(n))) %>%
  count(n_bin,category) %>% 
  mutate(
    n_bin = factor(
      n_bin,
      levels = c(as.character(1:14), "15+")
    )
  ) %>% 
  ggplot(aes(x = n_bin, y = n,fill=category)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene copy number",
       y = "Frequency") +
  scale_fill_discrete(limits = c("KPPD", "KPPD ∩ CPC", "KPPD ∩ HPRC", "ALL")) + 
  theme_step1() + 
  theme(legend.title = element_blank(),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = "white",color = NA,alpha = 0.7))



kppd_byGene_multi_exon_qc %>% left_join(cpc_multi_exon_qc %>% rename(gene_name = gene,cpc = count)) %>%
  left_join(hprc_qc %>% rename(gene_name = gene,hprc = n))

cpc_multi_exon_qc %>% rename(gene_name = gene) %>%
  mutate(CPC = count/116) %>% select(gene_name,CPC) -> cpc_multi_exon_qc_freq

hprc_qc %>% rename(gene_name = gene) %>% arrange(-n)

hprc_qc %>% rename(gene_name = gene) %>% #head()
  mutate(HPRC = n/88) %>% select(gene_name,HPRC) -> hprc_qc_freq
hprc_qc_freq

kppd_byGene_multi_exon_qc %>% mutate(KPPD = n/264) %>% 
  select(gene_name,KPPD) ->  kppd_byGene_multi_exon_qc_freq
  
library(smplot2)
kppd_byGene_multi_exon_qc_freq %>% left_join(hprc_qc_freq) %>% #head()
  na.omit() %>% dim() #255
  ggplot(aes(x=KPPD,y=HPRC)) + 
  geom_point(shape = 21, fill = "#0f993d", size = 2,alpha = 0.8) +
  #geom_smooth(method = "lm", se = FALSE, color = "0f993d", linetype = "solid", size = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray40", size = 0.8) +  # ??? y=x ??��??
  sm_statCorr(
    color = "#0f993d", corr_method = "spearman",
    linetype = "dashed",
    text_size = 10
  ) +
  ylim(0,1) + xlim(0,1) + 
  theme_step1()


kppd_byGene_multi_exon_qc_freq %>% left_join(hprc_qc_freq) %>% 
  na.omit() %>% filter(HPRC != 1) %>% #dim() #198
  ggplot(aes(x=KPPD,y=HPRC)) +
  geom_point(shape = 21, fill = "#0f993d", size = 2,alpha = 0.8) +
  #geom_smooth(method = "lm", se = FALSE, color = "0f993d", linetype = "solid", size = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray40", size = 0.8) +  # ??? y=x ??��??
  sm_statCorr(
    color = "#0f993d", corr_method = "spearman",
    linetype = "dashed",
    text_size = 10
  ) +
  ylim(0,1) + xlim(0,1) + 
  theme_step1()


kppd_byGene_multi_exon_qc_freq %>% left_join(cpc_multi_exon_qc_freq) %>% #head()
  na.omit() %>% #dim #400
  ggplot(aes(x=KPPD,y=CPC)) +
  geom_point(shape = 21, fill = "#0f993d", size = 2,alpha = 0.8) +
  #geom_smooth(method = "lm", se = FALSE, color = "0f993d", linetype = "solid", size = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray40", size = 0.8) +  # ??? y=x ??��??
  sm_statCorr(
    color = "#0f993d", corr_method = "spearman",
    linetype = "dashed",
    text_size = 10
  ) +
  theme_step1()

cpc %>% filter(gene == "USP17L11")

kppd_byGene_multi_exon_qc_freq %>% left_join(cpc_multi_exon_qc_freq) %>% 
  left_join(hprc_qc_freq) -> cnv_allmerge

cnv_allmerge %>% na.omit() %>%
  ggplot(aes(x=KPPD,y=CPC)) +
  geom_point(shape = 21, fill = "#0f993d", size = 2,alpha = 0.8) +
  #geom_smooth(method = "lm", se = FALSE, color = "0f993d", linetype = "solid", size = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray40", size = 0.8) +  # ??? y=x ??��??
  sm_statCorr(
    color = "#0f993d", corr_method = "spearman",
    linetype = "dashed",
    text_size = 10
  ) +
  ylim(0,0.3) + xlim(0,0.3) + 
  theme_step1()

cnv_allmerge %>% na.omit() %>%
  ggplot(aes(x=KPPD,y=HPRC)) +
  geom_point(shape = 21, fill = "#0f993d", size = 2,alpha = 0.8) +
  #geom_smooth(method = "lm", se = FALSE, color = "0f993d", linetype = "solid", size = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray40", size = 0.8) +  # ??? y=x ??��??
  sm_statCorr(
    color = "#0f993d", corr_method = "spearman",
    linetype = "dashed",
    text_size = 10
  ) +
  ylim(0,0.3) + xlim(0,0.3) + 
  theme_step1()

cpc %>% mutate(CPC = count/116) %>%  rename(gene_name = gene) %>% select(gene_name,CPC)
hprc %>% rename(gene_name = gene) %>% mutate(HPRC = n/88)  %>% select(gene_name,HPRC)
  
kppd_byGene %>% mutate(KPPD = n/264) %>% select(-n) %>%
  left_join(cpc %>% mutate(CPC = count/116) %>%  rename(gene_name = gene) %>% select(gene_name,CPC)) %>% 
  left_join(hprc %>% rename(gene_name = gene) %>% mutate(HPRC = n/88)  %>% select(gene_name,HPRC)) -> cnv_allmerge_freq

cnv_allmerge_freq%>% filter(gene_name == "PNMA6A")

cnv_allmerge_freq %>% arrange(-KPPD) %>%
  mutate(KPPD_CPC = KPPD-CPC) %>%
  mutate(KPPD_HPRC = KPPD-HPRC) %>%
  filter(KPPD > 0.1) %>% arrange(-KPPD_HPRC) #%>% slice_head(n = 5) %>% select(gene_name,KPPD,HPRC,KPPD_HPRC) -> KPPD_HPRC

cnv_allmerge_freq %>% arrange(-KPPD) %>%
  mutate(KPPD_CPC = KPPD-CPC) %>%
  mutate(KPPD_HPRC = KPPD-HPRC) %>%
  filter(KPPD > 0.1) %>% arrange(-KPPD_CPC) %>% slice_head(n = 10) %>% select(gene_name,KPPD,CPC,KPPD_CPC) -> KPPD_CPC

KPPD_HPRC

cnv_allmerge_freq %>% arrange(-KPPD) %>%
  mutate(KPPD_CPC = KPPD-CPC) %>%
  mutate(KPPD_HPRC = KPPD-HPRC) %>%
  filter(KPPD > 0.1) %>% arrange(-KPPD_CPC) #%>% slice_head(n = 5) %>% select(gene_name,KPPD,HPRC,KPPD_HPRC) -> KPPD_HPRC


cnv_allmerge_freq %>% filter(gene_name == "AMY2A")

cnv_allmerge_freq %>% filter(gene_name == "OR2A1")



library(dplyr)
library(tidyr)
library(ggplot2)

KPPD_HPRC %>%
  select(gene_name, KPPD, HPRC) %>%
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

