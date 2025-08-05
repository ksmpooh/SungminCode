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

#HPRC_GRCh38_gene_merge.txt
#HPRC_GRCh38_transcript_merge.txt
#HPRC_T2T_gene_merge.txt
#HPRC_T2T_transcripts_merge.txt
#KPPD_GRCh38_gene_merge.txt
#KPPD_GRCh38_transcript_merge.txt
#KPPD_T2T_gene_merge.txt
#KPPD_T2T_transcripts_merge.txt
##205:/ADATA/pangenome/annotation/merge


KPPD_GRCh38_gene <- read_table("KPPD_GRCh38_gene_merge.txt",col_names = F)
KPPD_T2T_gene <- read_table("KPPD_T2T_gene_merge.txt",col_names = F)

KPPD_GRCh38_transcripts <- read_table("KPPD_GRCh38_transcript_merge.txt",col_names = F)
KPPD_T2T_transcripts <- read_table("KPPD_T2T_transcripts_merge.txt",col_names = F)

CPC_GRCh38_gene <- read_table("CPC_GRCh38_gene_merge.cov0.9.txt",col_names = F)

KPPD_GRCh38_gene %>% mutate(ID = str_split_fixed(X1,"tg",2)[,1]) %>% 
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  select(ID,coding_type,X2) %>% unique() %>% 
  group_by(ID) %>% count(ID,coding_type) %>% write.table("preprocessing/KPPD_GRCh38_gene.protein_coding_count.txt",col.names = T,row.names = F,quote = F,sep = "\t")

KPPD_T2T_gene %>% mutate(ID = str_split_fixed(X1,"tg",2)[,1]) %>% 
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  select(ID,coding_type,X2) %>% unique() %>% 
  group_by(ID) %>% count(ID,coding_type) %>% write.table("preprocessing/KPPD_T2T_gene.protein_coding_count.txt",col.names = T,row.names = F,quote = F,sep = "\t")


KPPD_GRCh38_transcripts %>% mutate(ID = str_split_fixed(X1,"tg",2)[,1]) %>% 
  mutate(X2 = str_split_fixed(X2,"_",2)[,1]) %>%
  mutate(X4 = str_split_fixed(X4,"_",2)[,1]) %>%
  mutate(coding_type = ifelse(str_detect(X3,"protein_coding"),"Protein_coding_trans","Noncoding_trans")) %>%
  select(ID,coding_type,X2,X4) %>% unique() %>% 
  group_by(ID) %>% count(ID,coding_type) %>% write.table("preprocessing/KPPD_GRCh38_transcripts.protein_coding_count.txt",col.names = T,row.names = F,quote = F,sep = "\t")


KPPD_T2T_transcripts %>% mutate(ID = str_split_fixed(X1,"tg",2)[,1]) %>% 
  mutate(X2 = str_split_fixed(X2,"_",2)[,1]) %>%
  mutate(X4 = str_split_fixed(X4,"_",2)[,1]) %>%
  mutate(coding_type = ifelse(str_detect(X3,"protein_coding"),"Protein_coding_trans","Noncoding_trans")) %>%
  select(ID,coding_type,X2,X4) %>% unique() %>% 
  group_by(ID) %>% count(ID,coding_type) %>% write.table("preprocessing/KPPD_T2T_transcripts.protein_coding_count.txt",col.names = T,row.names = F,quote = F,sep = "\t")

HPRC_GRCh38_gene <- read_table("HPRC_GRCh38_gene_merge.txt",col_names = F)
HPRC_T2T_gene <- read_table("HPRC_T2T_gene_merge.txt",col_names = F)

HPRC_GRCh38_transcripts <- read_table("HPRC_GRCh38_transcript_merge.txt",col_names = F)
HPRC_T2T_transcripts <- read_table("HPRC_T2T_transcripts_merge.txt",col_names = F)



HPRC_GRCh38_gene %>% #filter(!str_detect(X1,"chr")) %>% 
  mutate(ID = ifelse(str_detect(X1,"chr"),"CHM13",(paste0(str_split_fixed(X1,"#",3)[,1],"_h",str_split_fixed(X1,"#",3)[,2])))) %>%
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  select(ID,coding_type,X2) %>% unique() %>% 
  group_by(ID) %>% count(ID,coding_type) %>% write.table("preprocessing/HPRC_GRCh38_gene.protein_coding_count.txt",col.names = T,row.names = F,quote = F,sep = "\t")

HPRC_T2T_gene %>% 
  mutate(ID = ifelse(str_detect(X1,"chr"),"CHM13",(paste0(str_split_fixed(X1,"#",3)[,1],"_h",str_split_fixed(X1,"#",3)[,2])))) %>%
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  select(ID,coding_type,X2) %>% unique() %>% 
  group_by(ID) %>% count(ID,coding_type) %>% write.table("preprocessing/HPRC_T2T_gene.protein_coding_count.txt",col.names = T,row.names = F,quote = F,sep = "\t")


HPRC_GRCh38_transcripts %>%
  mutate(X2 = str_split_fixed(X2,"_",2)[,1]) %>%
  mutate(X4 = str_split_fixed(X4,"_",2)[,1]) %>%
  mutate(ID = ifelse(str_detect(X1,"chr"),"CHM13",(paste0(str_split_fixed(X1,"#",3)[,1],"_h",str_split_fixed(X1,"#",3)[,2])))) %>%
  mutate(coding_type = ifelse(str_detect(X3,"protein_coding"),"Protein_coding_trans","Noncoding_trans")) %>%
  select(ID,coding_type,X2,X4) %>% unique() %>% 
  group_by(ID) %>% count(ID,coding_type) %>% write.table("preprocessing/HPRC_GRCh38_transcripts.protein_coding_count.txt",col.names = T,row.names = F,quote = F,sep = "\t")


HPRC_T2T_transcripts %>% 
  mutate(X2 = str_split_fixed(X2,"_",2)[,1]) %>%
  mutate(X4 = str_split_fixed(X4,"_",2)[,1]) %>%
  mutate(ID = ifelse(str_detect(X1,"chr"),"CHM13",(paste0(str_split_fixed(X1,"#",3)[,1],"_h",str_split_fixed(X1,"#",3)[,2])))) %>%
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_trans","Noncoding_trans")) %>%
  select(ID,coding_type,X2,X4) %>% unique() %>% 
  group_by(ID) %>% count(ID,coding_type) %>% write.table("preprocessing/HPRC_T2T_transcripts.protein_coding_count.txt",col.names = T,row.names = F,quote = F,sep = "\t")


#### server copy number by ID
KPPD_GRCh38_gene %>% mutate(ID = str_split_fixed(X1,"tg",2)[,1]) %>% 
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  filter(X4 != 0) %>% 
  select(ID,coding_type,X2) %>% unique() %>% 
  group_by(ID) %>% count(ID,coding_type) %>% write.table("cnv_count/KPPD_GRCh38_gene.copy_number_count.txt",col.names = T,row.names = F,quote = F,sep = "\t")

KPPD_T2T_gene %>% mutate(ID = str_split_fixed(X1,"tg",2)[,1]) %>% 
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  filter(X4 != 0) %>% 
  select(ID,coding_type,X2) %>% unique() %>% 
  group_by(ID) %>% count(ID,coding_type) %>% write.table("cnv_count/KPPD_T2T_gene.copy_number_count.txt",col.names = T,row.names = F,quote = F,sep = "\t")

HPRC_GRCh38_gene %>% #filter(!str_detect(X1,"chr")) %>% 
  mutate(ID = ifelse(str_detect(X1,"chr"),"CHM13",(paste0(str_split_fixed(X1,"#",3)[,1],"_h",str_split_fixed(X1,"#",3)[,2])))) %>%
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  filter(X4 != 0) %>% 
  select(ID,coding_type,X2) %>% unique() %>% 
  group_by(ID) %>% count(ID,coding_type) %>% write.table("cnv_count/HPRC_GRCh38_gene.copy_number_count.txt",col.names = T,row.names = F,quote = F,sep = "\t")

HPRC_T2T_gene %>% 
  mutate(ID = ifelse(str_detect(X1,"chr"),"CHM13",(paste0(str_split_fixed(X1,"#",3)[,1],"_h",str_split_fixed(X1,"#",3)[,2])))) %>%
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  filter(X4 != 0) %>% 
  select(ID,coding_type,X2) %>% unique() %>% 
  group_by(ID) %>% count(ID,coding_type) %>% write.table("cnv_count/HPRC_T2T_gene.copy_number_count.txt",col.names = T,row.names = F,quote = F,sep = "\t")

CPC_GRCh38_gene %>% #filter(!str_detect(X1,"chr")) %>% 
  mutate(ID = ifelse(str_detect(X1,"chr"),"CHM13",(paste0(str_split_fixed(X1,"#",3)[,1],"_h",str_split_fixed(X1,"#",3)[,2])))) %>%
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  filter(X4 != 0) %>% 
  select(ID,coding_type,X2) %>% unique() %>% 
  group_by(ID) %>% count(ID,coding_type) %>% write.table("cnv_count/CPC_GRCh38_gene.copy_number_count.txt",col.names = T,row.names = F,quote = F,sep = "\t")

a %>% filter(coding_type == "Protein_coding_genes") %>% filter(X4 != 0) %>% select(X2,X4) %>% unique() %>% dim
#2684    2

a %>% filter(coding_type == "Protein_coding_genes") %>% filter(X4 != 0) %>% select(X2) %>% unique() %>% dim

#1395
a %>% filter(coding_type == "Protein_coding_genes") %>% filter(X4 == 1) %>% select(X2) %>% unique() %>% dim()
#1395
a %>% filter(coding_type == "Protein_coding_genes") %>% filter(X4 != 0,X7 != 0) %>% select(X2) %>% unique() %>% dim
#1324

a %>% filter(coding_type == "Protein_coding_genes") %>% filter(X4 == 1) %>% select(X2,X6) %>% unique() %>% count(X6)
#<dbl> <int>
#1   227
#2  1154
#3    14
#227 + 1154 + 14 = 1395
#1395 - 14  = 1381

a %>% filter(coding_type == "Protein_coding_genes") %>% filter(X4 == 1,X7 > 0) %>% select(X2) %>% unique() %>% dim

#592 + 76 +396 + 331 = 1395 #
#10 + 113 +165 + 1079 = 1367 # CPC paper


cpc_sup_excel <- readxl::read_xlsx("~/Desktop/KCDC/paper/pangenome/CPC_cnv.xlsx",skip = 1) %>% mutate(across(everything(), ~ str_replace_all(., " ", ""))) %>%
  mutate(across(everything(), ~ str_replace_all(., "\n", ""))) 
head(cpc_sup_excel)
dim(cpc_sup_excel)

#1,367
cpc_sup_excel %>% select(1:6) -> a
cpc_sup_excel %>% select(7:12) -> b
colnames(a) <- c("Gene","CPC_count","HPRC.EAS.count","HPRC,nEAS.count","prop","Tajima")
colnames(b) <- c("Gene","CPC_count","HPRC.EAS.count","HPRC,nEAS.count","prop","Tajima")
a %>% rbind(b) -> a
head(a)
dim(a)
a %>% filter(is.na(CPC_count))
a %>% na.omit() %>% dim()
head(df)
df %>% filter(cohort=="CPC") %>% filter(!(X2 %in% a$Gene)) -> b
df %>% filter(cohort=="CPC") -> c
b$X2
head(a)
a %>% filter(!(Gene %in% c$X2))

#### server copy number with level by ID

KPPD_GRCh38_gene %>% mutate(ID = str_split_fixed(X1,"tg",2)[,1]) %>% 
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  filter(X4 != 0) %>% filter(X6 != 3) %>%
  select(ID,coding_type,X2) %>% unique() %>% 
  group_by(ID) %>% count(ID,coding_type) %>% #rename("level" = X6) %>%
  write.table("cnv_count/KPPD_GRCh38_gene.copy_number_count_withoutlevel3.txt",col.names = T,row.names = F,quote = F,sep = "\t")


HPRC_GRCh38_gene %>% #filter(!str_detect(X1,"chr")) %>% 
  mutate(ID = ifelse(str_detect(X1,"chr"),"CHM13",(paste0(str_split_fixed(X1,"#",3)[,1],"_h",str_split_fixed(X1,"#",3)[,2])))) %>%
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  filter(X4 != 0) %>% filter(X6 != 3) %>%
  select(ID,coding_type,X2) %>% unique() %>% 
  group_by(ID) %>% count(ID,coding_type) %>% #rename("level" = X6) %>%
  write.table("cnv_count/HPRC_GRCh38_gene.copy_number_count_withoutlevel3.txt",col.names = T,row.names = F,quote = F,sep = "\t")

CPC_GRCh38_gene %>% #filter(!str_detect(X1,"chr")) %>% 
  mutate(ID = ifelse(str_detect(X1,"chr"),"CHM13",(paste0(str_split_fixed(X1,"#",3)[,1],"_h",str_split_fixed(X1,"#",3)[,2])))) %>%
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  filter(X4 != 0) %>% filter(X6 != 3) %>%
  select(ID,coding_type,X2) %>% unique() %>% 
  group_by(ID) %>% count(ID,coding_type) %>% #rename("level" = X6) %>%
  write.table("cnv_count/CPC_GRCh38_gene.copy_number_count_withoutlevel3.txt",col.names = T,row.names = F,quote = F,sep = "\t")



### copy number gene freq
KPPD_GRCh38_gene %>% mutate(ID = str_split_fixed(X1,"tg",2)[,1]) %>% 
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  filter(coding_type == "Protein_coding_genes") %>%
  filter(X4 != 0) %>% 
  select(ID,coding_type,X2) %>% unique() %>% 
  count(X2) %>% write.table("cnv_count/KPPD_GRCh38_gene.copy_number_gene_freq.txt",col.names = T,row.names = F,quote = F,sep = "\t")

KPPD_T2T_gene %>% mutate(ID = str_split_fixed(X1,"tg",2)[,1]) %>% 
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  filter(coding_type == "Protein_coding_genes") %>%
  filter(X4 != 0) %>% 
  select(ID,coding_type,X2) %>% unique() %>% 
  count(X2) %>% write.table("cnv_count/KPPD_T2T_gene.copy_number_gene_freq.txt",col.names = T,row.names = F,quote = F,sep = "\t")

HPRC_GRCh38_gene %>% #filter(!str_detect(X1,"chr")) %>% 
  mutate(ID = ifelse(str_detect(X1,"chr"),"CHM13",(paste0(str_split_fixed(X1,"#",3)[,1],"_h",str_split_fixed(X1,"#",3)[,2])))) %>%
  filter(ID != "CHM13") %>%
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  filter(coding_type == "Protein_coding_genes") %>%
  filter(X4 != 0) %>% 
  select(ID,coding_type,X2) %>% unique() %>% 
  count(X2) %>% write.table("cnv_count/HPRC_GRCh38_gene.copy_number_gene_freq.txt",col.names = T,row.names = F,quote = F,sep = "\t")

HPRC_T2T_gene %>% 
  mutate(ID = ifelse(str_detect(X1,"chr"),"CHM13",(paste0(str_split_fixed(X1,"#",3)[,1],"_h",str_split_fixed(X1,"#",3)[,2])))) %>%
  filter(ID != "CHM13") %>%
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  filter(coding_type == "Protein_coding_genes") %>%
  filter(X4 != 0) %>% 
  select(ID,coding_type,X2) %>% unique() %>% 
  count(X2) %>% write.table("cnv_count/HPRC_T2T_gene.copy_number_gene_freq.txt",col.names = T,row.names = F,quote = F,sep = "\t")


CPC_GRCh38_gene %>% #filter(!str_detect(X1,"chr")) %>% 
  mutate(ID = ifelse(str_detect(X1,"chr"),"CHM13",(paste0(str_split_fixed(X1,"#",3)[,1],"_h",str_split_fixed(X1,"#",3)[,2])))) %>% #head()
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  filter(coding_type == "Protein_coding_genes") %>%
  filter(X4 != 0) %>% 
  select(ID,coding_type,X2) %>% unique() %>% 
  count(X2) %>% write.table("cnv_count/CPC_GRCh38_gene.copy_number_gene_freq.txt",col.names = T,row.names = F,quote = F,sep = "\t")


# with level
KPPD_GRCh38_gene %>% mutate(ID = str_split_fixed(X1,"tg",2)[,1]) %>% 
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  filter(coding_type == "Protein_coding_genes") %>%
  filter(X4 != 0) %>% 
  select(ID,coding_type,X2,X6) %>% unique() %>% 
  count(X6,X2) %>% write.table("cnv_count/KPPD_GRCh38_gene.copy_number_gene_freq_withlevel.txt",col.names = T,row.names = F,quote = F,sep = "\t")

HPRC_GRCh38_gene %>% #filter(!str_detect(X1,"chr")) %>% 
  mutate(ID = ifelse(str_detect(X1,"chr"),"CHM13",(paste0(str_split_fixed(X1,"#",3)[,1],"_h",str_split_fixed(X1,"#",3)[,2])))) %>%
  filter(ID != "CHM13") %>%
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  filter(coding_type == "Protein_coding_genes") %>%
  filter(X4 != 0) %>% 
  select(ID,coding_type,X2,X6) %>% unique() %>% 
  count(X6,X2) %>% write.table("cnv_count/HPRC_GRCh38_gene.copy_number_gene_freq_withlevel.txt",col.names = T,row.names = F,quote = F,sep = "\t")


CPC_GRCh38_gene %>% #filter(!str_detect(X1,"chr")) %>% 
  mutate(ID = ifelse(str_detect(X1,"chr"),"CHM13",(paste0(str_split_fixed(X1,"#",3)[,1],"_h",str_split_fixed(X1,"#",3)[,2])))) %>%
  filter(ID != "CHM13") %>%
  mutate(coding_type = ifelse(X3=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>%
  filter(coding_type == "Protein_coding_genes") %>%
  filter(X4 != 0) %>% 
  select(ID,coding_type,X2,X6) %>% unique() %>% 
  count(X6,X2) %>% write.table("cnv_count/CPC_GRCh38_gene.copy_number_gene_freq_withlevel.txt",col.names = T,row.names = F,quote = F,sep = "\t")



#######
#######



setwd("~/Desktop/KCDC/pangenome/KPPD/liftoff/")

#genome@genome205:/ADATA/pangenome/db/chm13$ head chm13.draft_v2.0.gene_annotation.gff3_transcript
#protein_coding 156386
#other 78517
#genome@genome205:/ADATA/pangenome/db/chm13$ cat chm13.draft_v2.0.gene_annotation.gff3_gene
#protein_coding 20067
#other 44146
flist <- list.files("count/",pattern = "HPRC")
flist
df <- NULL
for (i in flist) {
  tmp <- read_table(paste0("count/",i))
  tmp$filename <- i
  df <- rbind(df,tmp)
  
}
hprc <- df
hprc %>% count(filename)
head(df)
df %>% filter(!str_detect(ID,"chr")) %>% head()

flist <- list.files("count/",pattern = "KPPD")
flist
df <- NULL
for (i in flist) {
  tmp <- read_table(paste0("count/",i))
  tmp$filename <- i
  df <- rbind(df,tmp)
  
}
kppd <- df
df <- NULL
head(kppd)
head(hprc)
hprc %>% filter(coding_type == "Protein_coding_genes")
hprc %>% filter(!str_detect(ID,"chr")) %>% count(coding_type)
hprc %>% filter(!str_detect(ID,"chr")) %>% #count(coding_type)
#  mutate(ID = paste0(str_split_fixed(ID,"#",3)[,1],"_h",str_split_fixed(ID,"#",3)[,2])) %>% #head()
  mutate(ref = str_split_fixed(filename,"\\.",3)[,1]) %>% #count(ref)
  group_by(ref,ID,coding_type) %>% 
  summarise(n=sum(n)) %>% 
  ungroup() -> hprc_pro

head(hprc_pro)
hprc_pro %>% count(ref,coding_type)
  
head(kppd)
kppd %>% mutate(ref = str_split_fixed(filename,"\\.",3)[,1]) %>% select(-filename) %>%
  rbind(hprc_pro) %>% 
  mutate(cohort = str_split_fixed(ref,"_",3)[,1],type = str_split_fixed(ref,"_",3)[,3],ref = str_split_fixed(ref,"_",3)[,2]) -> df
head(hprc_pro)
table(hprc_pro$ref)
hprc_pro %>% filter(ref == "HPRC_GRCh38_transcripts",coding_type != "Protein_coding_trans")
head(df)
df %>% count(ref,cohort,type)

df %>% #filter(type == "gene") %>%
  mutate(annotation_pct = ifelse(ref=="T2T",ifelse(coding_type == "Protein_coding_genes", n/19866, n/39626),
         ifelse(coding_type == "Protein_coding_genes", n/20089, n/57606))) %>% #head()
  mutate(annotation_pct = ifelse(type == "gene",annotation_pct,ifelse(ref=="T2T",ifelse(coding_type == "Protein_coding_trans", n/154726, n/75345),
                                 ifelse(coding_type == "Protein_coding_trans", n/89861, n/298083)))) %>% #head()
  select(ID, coding_type, ref,cohort,annotation_pct) %>% #count(ref)
  mutate(coding_type = factor(coding_type, levels = c("Protein_coding_genes", "Noncoding_genes"))) %>% #head
  ggplot(aes(x = coding_type, y = annotation_pct, fill = cohort)) +
    geom_quasirandom(shape = 23, alpha = 0.7, dodge.width = 0.8) +
    labs(y="Percentage annotated (%)") + 
    facet_wrap(~ref) + 
    scale_x_discrete(labels = c("Protein_coding_genes" = "Protein-coding\ngenes", 
                                "Noncoding_genes" = "Noncoding\ngenes")) +
    theme_step1() +
    theme(axis.title.x = element_blank(),
          legend.title = element_blank(),
#          legend.position = 'none', 
          strip.text.x = element_text(
            size = 14,        # 폰트 크기
            face = "bold",    # 폰트 굵기 (bold, italic 등)
            family = "Arial"))
  


df %>% filter(type == "gene") %>%
  mutate(annotation_pct = ifelse(ref=="T2T",ifelse(coding_type == "Protein_coding_genes", n/19866, n/39626),
                                 ifelse(coding_type == "Protein_coding_genes", n/20089, n/57606))) %>% 
  filter(ref == "T2T",cohort == "HPRC",annotation_pct > 0.96) -> c


hprc_sup_sex <- readxl::read_xlsx("~/Desktop/KCDC/paper/pangenome/Draft human pangenome reference41586_2023_5896_MOESM4_ESM.xlsx",sheet = 7,skip = 1)
head(hprc_sup_sex)
hprc_sup_sex %>% select(sample,Sex) %>% unique() -> hprc_sup_sex

kddp_sex <- read_table("~/Desktop/KCDC/pangenome/00.datacheck/KPP_draft_assembly_matchingID.with.Gender.20250224.txt")
head(kddp_sex)
kddp_sex %>% mutate(Sex = ifelse(sex =="M","male","female")) %>% mutate(sample = KPPD_ID) %>% select(sample,Sex) -> kddp_sex

kddp_sex %>% rbind(hprc_sup_sex) -> sex_info
head(sex_info)
 
hprc_sup_sex %>% mutate(sample = str_split_fixed(ID,"_",2)[,1]) %>% left_join(hprc_sup %>% select(sample,Sex)) %>% count(Sex)

hprc_sup <- readxl::read_xlsx("~/Desktop/KCDC/paper/pangenome/Draft human pangenome reference41586_2023_5896_MOESM4_ESM.xlsx",sheet = 8,skip = 1) %>% 
  filter(!(`...1` %in% c("median","min","max"))) %>% filter(!is.na(`...1`))
head(hprc_sup)
dim(hprc_sup)
tail(hprc_sup)
hprc_sup %>% select(1:5) %>% mutate(ID = str_replace_all(`...1`,"\\.","_")) -> hprc_sup
hprc_sup
colnames(hprc_sup) <- c("sample","Protein_coding_genes","Protein_coding_trans","Noncoding_genes","Noncoding_trans","ID")
head(hprc_sup)
hprc_sup %>% select(-sample) %>% pivot_longer(1:4,names_to = "coding_type",values_to = "annotation_pct") %>% 
  mutate(annotation_pct = ifelse(annotation_pct < 1,annotation_pct*100,annotation_pct)) -> hprc_sup

head(hprc_sup)
  
head(hprc_sup)
hprc_sup %>% count(ID) %>% count(n)
hprc_sup %>% filter(str_detect(ID,"GRC"))
hprc_sup %>% select(ID) %>% unique() -> a

head(df)
df %>% filter(ref == "T2T",cohort == "KPPD") %>% #head()
  mutate(cohort = ifelse(ID %in% c("KPPD129_h1","KPPD129_h2","KPPD130_h1","KPPD130_h2","KPPD131_h1","KPPD131_h2","KPPD132_h1","KPPD132_h2"),"KPP",cohort)) %>%
  #mutate(annotation_pct = ifelse(coding_type == "Protein_coding_genes", n/19866, n/39626)) %>% 
  mutate(annotation_pct = ifelse(type == "gene",ifelse(coding_type == "Protein_coding_genes", n/19866, n/39626),
                                 ifelse(coding_type == "Protein_coding_trans", n/154726, n/75345))) %>% 
  mutate(annotation_pct = annotation_pct*100) %>% select(ID,coding_type,cohort,annotation_pct) %>% 
  rbind(hprc_sup %>% mutate(cohort = "HPRC")) %>% mutate(sample = str_split_fixed(ID,"_",2)[,1],hap = str_split_fixed(ID,"_",2)[,2]) %>% left_join(sex_info) %>% #head()
  mutate(cohort = ifelse(ID == "CHM13","CHM13",cohort),hap = ifelse(ID == "CHM13","CHM13",hap),Sex = ifelse(ID == "CHM13","CHM13",Sex)) %>% #count(cohort)
  mutate(g = ifelse(ID == "CHM13","HPRC",cohort)) %>%
  #group_by(cohort,coding_type) %>% count(cohort,coding_type)#summarise(mean(annotation_pct))
  mutate(coding_type = factor(coding_type, levels = c("Protein_coding_genes", "Noncoding_genes","Protein_coding_trans","Noncoding_trans"))) %>% #head
  ggplot(aes(x = coding_type, y = annotation_pct, group = g,shape= Sex,color = cohort)) +
  #geom_quasirandom(aes(shape= Sex,fill = cohort),alpha = 0.8, dodge.width = 0.8) +
  geom_quasirandom(alpha = 0.7, dodge.width = 0.8,size=2) +
  labs(y="Percentage annotated (%)") + 
  scale_x_discrete(labels = c("Protein_coding_genes" = "Protein-coding\ngenes", 
                              "Noncoding_genes" = "Noncoding\ngenes",
                              "Protein_coding_trans" = "Protein-coding\ntranscripts", 
                              "Noncoding_trans" = "Noncoding\ntranscripts")) +
  theme_step1() +
  theme(axis.title.x = element_blank(),
        #legend.position = 'none', 
        legend.title = element_blank(),
        strip.text.x = element_text(
          size = 14,        # 폰트 크기
          face = "bold",    # 폰트 굵기 (bold, italic 등)
          family = "Arial"))



####
table(df$coding_type)
head(df)
df %>% #filter(type == "gene") %>%
  mutate(annotation_pct = ifelse(ref=="T2T",ifelse(coding_type == "Protein_coding_genes", n/19866, n/39626),
                                 ifelse(coding_type == "Protein_coding_genes", n/20089, n/57606))) %>% #head()
  mutate(annotation_pct = ifelse(type == "gene",annotation_pct,ifelse(ref=="T2T",ifelse(coding_type == "Protein_coding_trans", n/154726, n/75345),
                                                                      ifelse(coding_type == "Protein_coding_trans", n/116490, n/271454)))) %>% #head()
  mutate(cohort = ifelse(ID == "CHM13","CHM13",cohort)) %>% 
  mutate(g = ifelse(cohort == "KPPD","KPPD","HPRC")) %>% #head()
  select(ID, coding_type, ref,cohort,annotation_pct,g) %>% #head()
  mutate(coding_type = factor(coding_type, levels = c("Protein_coding_genes", "Noncoding_genes","Protein_coding_trans","Noncoding_trans"))) %>%#-> a
  ggplot(aes(x = coding_type, y = annotation_pct, fill = cohort,group = g)) +
  geom_quasirandom(shape = 23, alpha = 0.7, dodge.width = 0.8) +
  labs(y="Percentage annotated (%)") + 
  facet_wrap(~ref,ncol = 1) + 
  scale_x_discrete(labels = c("Protein_coding_genes" = "Protein-coding\ngenes", 
                              "Noncoding_genes" = "Noncoding\ngenes",
                              "Protein_coding_trans" = "Protein-coding\ntranscripts", 
                              "Noncoding_trans" = "Noncoding\ntranscripts")) +
  theme_step1() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(), 
        #          legend.position = 'none', 
        strip.text.x = element_text(
          size = 14,        # 폰트 크기
          face = "bold",    # 폰트 굵기 (bold, italic 등)
          family = "Arial"))


head(a)
a %>% group_by(ref,cohort,coding_type) %>% summarise(annotation_pct = mean(annotation_pct)) %>% #$head()
  filter(cohort != "CHM13") %>% #head()
  pivot_wider(names_from = cohort,values_from = annotation_pct) %>%
  p

head(df)
df %>% filter(type == "transcripts") %>% filter(ref =="T2T")

###### cnv count by ID
setwd("~/Desktop/KCDC/pangenome/KPPD/liftoff/cnv_count/")

flist <- list.files(pattern = "copy_number_count.txt")
flist
df <- NULL
for (i in flist) {
  tmp <- read_table(i)
  tmp$filename <- i
  df <- rbind(df,tmp)
  
}
head(df)
df %>% filter(max(n))
df %>% filter(coding_type == "Protein_coding_genes") %>% group_by(filename) %>% summarise(max(n))
df %>% filter(n==217) #HG02572_h1 #HG03516_h1

df %>% 
  mutate(cohort = case_when(
    str_detect(ID, "CHM13") ~ "CHM13",
    str_detect(filename, "HPRC") ~ "HPRC",
    str_detect(filename, "KPPD") ~ "KPPD",
    str_detect(filename, "CPC") ~ "CPC",
    TRUE ~ "other"
  )) %>% mutate(ref = ifelse(str_detect(filename,"GRCh38"),"GRCh38","T2T")) %>% #filter(ref !="T2T") %>%
  filter(coding_type == "Protein_coding_genes") -> df_cnv_count


df_cnv_count %>% filter(ref !="T2T") %>%
  group_by(ref) %>%  arrange(n) %>% mutate(ID_ordered = factor(ID, levels = unique(ID))) %>%  # n 오름차순으로 level 지정
  ungroup() %>%
  ggplot(aes(x = ID_ordered, y = n, fill = cohort)) +
  #ggplot(aes(x=fct_reorder(ID, n),y=n,fill=cohort)) + 
  geom_bar(stat = "identity",position = "stack") + 
  facet_wrap(~ref,ncol = 1)

df_cnv_count %>% 
  filter(cohort !="CHM13") %>% #-> cnv_all_level
  mutate(cohort = factor(cohort, levels = c("HPRC", "CPC", "KPPD"))) %>% #-> a
  ggplot(aes(x=cohort,y=n,fill=cohort)) +
  geom_violin() + 
  geom_quasirandom(shape = 23, alpha = 0.7, dodge.width = 0.8) +
  labs(y="Duplicated genes for genome") + 
  facet_wrap(~ref,scales = "free_x") + 
  theme_step1() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(), 
        legend.position = 'none', 
        strip.text.x = element_text(
          size = 14,        # 폰트 크기
          face = "bold",    # 폰트 굵기 (bold, italic 등)
          family = "Arial"))

head(df_cnv_count)
head(a)

df_cnv_count %>% filter(cohort !="CHM13") %>% group_by(ref,cohort) %>% 
  summarise(mean(n),max(n),min(n)) -> a
a %>% filter()

df %>% #mutate(cohort = ifelse(ID == "CHM13",ID,ifelse(str_detect(filename,"HPRC"),"HPRC","KPPD"))) %>%
  mutate(cohort = case_when(
    str_detect(ID, "CHM13") ~ "CHM13",
    str_detect(filename, "HPRC") ~ "HPRC",
    str_detect(filename, "KPPD") ~ "KPPD",
    str_detect(filename, "CPC") ~ "CPC",
    TRUE ~ "other"
  )) %>%
  mutate(ref = ifelse(str_detect(filename,"GRCh38"),"GRCh38","T2T")) %>% 
  filter(coding_type == "Protein_coding_genes") %>% filter(cohort !="CHM13") -> cnv_all_level
  

###### cnv count by ID (without level 3)
setwd("~/Desktop/KCDC/pangenome/KPPD/liftoff/cnv_count/")

flist <- list.files(pattern = "without")
flist
df <- NULL
for (i in flist) {
  tmp <- read_table(i)
  tmp$filename <- i
  df <- rbind(df,tmp)
  
}
head(df)
df %>% filter(max(n))
df %>% filter(coding_type == "Protein_coding_genes") %>% group_by(filename) %>% summarise(max(n))
#df %>% filter(n==217) #HG02572_h1 #HG03516_h1


df %>% 
  mutate(cohort = case_when(
    str_detect(ID, "CHM13") ~ "CHM13",
    str_detect(filename, "HPRC") ~ "HPRC",
    str_detect(filename, "KPPD") ~ "KPPD",
    str_detect(filename, "CPC") ~ "CPC",
    TRUE ~ "other"
  )) %>% mutate(ref = ifelse(str_detect(filename,"GRCh38"),"GRCh38","T2T")) %>% #filter(ref !="T2T") %>%
  filter(coding_type == "Protein_coding_genes") -> df_cnv_count_withoutlevel3

df_cnv_count_withoutlevel3 %>% filter(ref !="T2T") %>%
  #filter(n=)
  group_by(ref) %>%  arrange(n) %>% mutate(ID_ordered = factor(ID, levels = unique(ID))) %>%  # n 오름차순으로 level 지정
  ungroup() %>%
  ggplot(aes(x = ID_ordered, y = n, fill = cohort)) +
  #ggplot(aes(x=fct_reorder(ID, n),y=n,fill=cohort)) + 
  geom_bar(stat = "identity",position = "stack") + 
  facet_wrap(~ref,ncol = 1)

df_cnv_count_withoutlevel3 %>% filter(ref !="T2T") %>% 
  filter(cohort !="CHM13") %>%
  mutate(cohort = factor(cohort, levels = c("HPRC", "CPC", "KPPD"))) %>%
  ggplot(aes(x=cohort,y=n,fill=cohort)) +
  geom_violin() + 
  geom_quasirandom(shape = 23, alpha = 0.7, dodge.width = 0.8) +
  labs(y="Duplicated genes for genome") + 
  facet_wrap(~ref) + 
  theme_step1() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(), 
        legend.position = 'none', 
        strip.text.x = element_text(
          size = 14,        # 폰트 크기
          face = "bold",    # 폰트 굵기 (bold, italic 등)
          family = "Arial"))




###############
setwd("~/Desktop/KCDC/pangenome/KPPD/liftoff/freq/")

flist <- list.files(pattern = "withlevel.txt")
flist
df <- NULL
for (i in flist) {
  tmp <- read_table(i)
  tmp$filename <- i
  df <- rbind(df,tmp)
  
}
head(df)
head(df)
df %>% #mutate(cohort = ifelse(str_detect(filename,"HPRC"),"HPRC","KPPD")) %>%
  mutate(cohort = case_when(
#    str_detect(ID, "CHM13") ~ "CHM13",
    str_detect(filename, "HPRC") ~ "HPRC",
    str_detect(filename, "KPPD") ~ "KPPD",
    str_detect(filename, "CPC") ~ "CPC",
    TRUE ~ "other")) %>%
  mutate(ref = ifelse(str_detect(filename,"GRCh38"),"GRCh38","T2T")) %>%
  filter(ref == "GRCh38") -> df

head(df)
df %>% count(cohort)
#1 CPC     1395
#2 HPRC    2217
#3 KPPD    2110

df %>% select(cohort,X2) %>% unique() %>% count(cohort)
#1 CPC     1395
#2 HPRC    2217
#3 KPPD    2110

#592 + 76 +396 + 331 = 1395 #
#10 + 113 +165 + 1079 = 1367 # CPC paper

head(df) 
df %>% filter(n != 3) %>%select(cohort,X2) %>% unique() %>% count(cohort) 
dim(df)
head(df)
df %>% 
  count(X2,cohort) %>% #dim #5722
  pivot_wider(names_from = cohort,values_from = n) %>% #head()
  mutate(HPRC = ifelse(HPRC == 1,X2,NA)) %>%
  mutate(CPC = ifelse(CPC == 1,X2,NA)) %>%
  mutate(KPPD = ifelse(KPPD == 1,X2,NA)) -> df_cnv_venn

head(df_cnv_venn)
df_cnv_venn %>% select(HPRC) %>% na.omit() %>% dim #2217
df_cnv_venn %>% select(KPPD) %>% na.omit() %>% dim #2110
df_cnv_venn %>% select(CPC) %>% na.omit() %>% dim #1395

ggVennDiagram::ggVennDiagram(list(HPRC = df_cnv_venn$HPRC,KPPD = df_cnv_venn$KPPD))
#1481 + 737
library(ggVennDiagram)

library(ggvenn)

ggVennDiagram(
  list(HPRC = df_cnv_venn$HPRC, KPPD = df_cnv_venn$KPPD),
  label_alpha = 0
) +
    scale_fill_gradient(low = "#d95f02", high = "#1b9e77") +  
  theme(legend.position = "none") +
  coord_flip()

ggvenn(
  list(HPRC = df_cnv_venn$HPRC, KPPD = df_cnv_venn$KPPD),
#  fill_color = c("HPRC" = "#1b9e77", "KPPD" = "#d95f02"),
  stroke_size = 1,
  text_size = 5,
  show_percentage = FALSE
)

library(VennDiagram)

# 1. 집합 데이터 준비
venn_data <- list(
  HPRC = df_cnv_venn$HPRC %>% na.omit(),
  CPC = df_cnv_venn$CPC %>% na.omit(),
  KPPD = df_cnv_venn$KPPD %>% na.omit()
)

# 2. 이미지로 바로 저장하지 않고, 그리기만 하기 위해 grid 사용
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = names(venn_data),
  fill = c("lightblue", "red","lightgreen"),
  lwd = 1.5,
  fontfamily = "sans",
  cex = 2,               # ?????? 숫자 표시 제거
  cat.cex = 3,         # 범주(label) 크기
  #cat.pos = c(-20, 20),  # label 위치 조정
  cat.fontface = "bold",
  scaled = TRUE,         # ??? 면적 비율 반영
  filename = NULL        # ?????? 화면에 바로 출력 (파일 저장 X)
)

# 3. 화면에 그리기
grid.newpage()
grid.draw(venn.plot)

###
df %>% mutate(cohort = ifelse(str_detect(filename,"HPRC"),"HPRC","KPPD")) %>%
  mutate(ref = ifelse(str_detect(filename,"GRCh38"),"GRCh38","T2T")) %>% filter(X6==1)
  filter(ref == "GRCh38") %>% count(X2,cohort) %>% #count(n)
  pivot_wider(names_from = cohort,values_from = n) %>% #head()
  mutate(HPRC = ifelse(HPRC == 1,X2,NA)) %>% 
  mutate(KPPD = ifelse(KPPD == 1,X2,NA)) -> df_cnv_venn
head(df_cnv_venn)
df_cnv_venn %>% select(HPRC) %>% na.omit() %>% dim #2217
df_cnv_venn %>% select(KPPD) %>% na.omit() %>% dim #2110
df_cnv_venn %>% select(CPC) %>% na.omit() %>% dim #1395

df

###
head(df)
table(df$cohort)
library(smplot2)
df %>% filter(cohort == "CPC")
df %>% #mutate(cohort = ifelse(str_detect(filename,"HPRC"),"HPRC","KPPD")) %>%
  mutate(ref = ifelse(str_detect(filename,"GRCh38"),"GRCh38","T2T")) %>% #head()
  filter(ref == "GRCh38") %>% #count(cohort)
  #mutate(freq = ifelse(cohort == "KPPD",n/(132*2),n/(47*2))) %>% rename("gene" = X2) %>%
  mutate(freq = ifelse(cohort == "KPPD",n/(132*2),ifelse(cohort == "HPRC",n/(47*2),n/(58*2)))) %>% rename("gene" = X2) %>% 
  select(gene,freq,cohort) -> df_cnv_freq


df_cnv_freq %>% filter(cohort != "CPC") %>%
  pivot_wider(names_from = cohort,values_from = freq) %>% na.omit() %>% #dim() # 736
  ggplot(aes(x=KPPD,y=HPRC)) +
  geom_point(shape = 21, fill = "#0f993d", size = 2,alpha = 0.8) +
  #geom_smooth(method = "lm", se = FALSE, color = "0f993d", linetype = "solid", size = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray40", size = 0.8) +  # ??? y=x 보조선
  sm_statCorr(
    color = "#0f993d", corr_method = "spearman",
    linetype = "dashed",
    text_size = 10
  ) +
  theme_step1()
  #labs(x="KPPD copy number gene frequency")
  
df_cnv_freq %>% filter(cohort != "HPRC") %>%
  pivot_wider(names_from = cohort,values_from = freq) %>% na.omit() %>% #dim() # 736
  ggplot(aes(x=KPPD,y=CPC)) +
  geom_point(shape = 21, fill = "#0f993d", size = 2,alpha = 0.8) +
  #geom_smooth(method = "lm", se = FALSE, color = "0f993d", linetype = "solid", size = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray40", size = 0.8) +  # ??? y=x 보조선
  sm_statCorr(color = "#0f993d", corr_method = "spearman",linetype = "dashed",text_size = 10) +
  theme_step1()

df_cnv_freq %>% filter(cohort != "HPRC") %>% pivot_wider(names_from = cohort,values_from = freq) %>% na.omit() %>% 
  dim() # 727
df_cnv_freq %>% filter(cohort != "CPC") %>% pivot_wider(names_from = cohort,values_from = freq) %>% na.omit() %>% 
  dim() # 736

  

df %>% mutate(cohort = ifelse(str_detect(filename,"HPRC"),"HPRC","KPPD")) %>%
  mutate(ref = ifelse(str_detect(filename,"GRCh38"),"GRCh38","T2T")) %>% #head()
  filter(ref == "GRCh38") %>% 
  mutate(freq = ifelse(cohort == "KPPD",n/(132*2),n/(47*2))) %>% rename("gene" = X2) %>%
  filter(str_detect(gene,"ENSG"))
  #filter(freq > 0.9,cohort == "KPPD")
  
  
setwd("~/Desktop/KCDC/pangenome/KPPD/liftoff/freq/")
kppd <- read_table("KPPD_GRCh38_gene.copy_number_gene_freq_withlevel.txt")
head(kppd)
colnames(kppd) <- c("level","gene","count")
kppd$cohort <- "KPPD"
hprc <- read_table("HPRC_GRCh38_gene.copy_number_gene_freq_withlevel.txt")
head(hprc)
colnames(hprc) <- c("level","gene","count")

hprc$cohort <- "HPRC"
hprc %>% filter(is.na(level))
hprc %>% filter(gene == "ADGRB2")
hprc %>% rbind(kppd) %>% 
  count(cohort,level)


###########################################################################################
head(df)
df <- read_table("KPPD.liftoff_genecode_chm13.protein_coding_count.txt")
df_09 <- read_table("KPPD.liftoff_genecode_chm13.protein_coding_count.coverge0.9.txt")
df_tr <- read_table("KPPD.liftoff_genecode_chm13.protein_coding_count_transcript.txt")
head(df)
head(df_09)
head(df_tr)

df %>% 
  left_join(df_09 %>% rename(n90 = n)) %>% 
  mutate(annotation_pct = ifelse(coding_type == "Protein_coding_genes", n/19866, n/39626)) %>% 
  mutate(annotation_pct_n90 = ifelse(coding_type == "Protein_coding_genes", n90/19866, n90/39626)) %>% 
  select(ID, coding_type, annotation_pct_n90, annotation_pct) %>% 
  pivot_longer(annotation_pct:annotation_pct_n90) %>% 
  # coding_type을 팩터로 변환하여 순서를 지정합니다.
  mutate(coding_type = factor(coding_type, levels = c("Protein_coding_genes", "Noncoding_genes"))) %>% 
  ggplot(aes(x = coding_type, y = value, fill = coding_type)) +
  geom_quasirandom(shape = 21, alpha = 0.7, dodge.width = 0.8) +
  labs(y="Percentage annotated") + 
  facet_grid(~name,labeller = as_labeller(
    c("annotation_pct" = "coverage (Overall)",
      "annotation_pct_n90" = "coverage (>=0.9)"))) +
  scale_x_discrete(labels = c("Protein_coding_genes" = "Protein-coding\ngenes", 
                              "Noncoding_genes" = "Noncoding\ngenes")) +
  theme_step1() +
  theme(axis.title.x = element_blank(),
        legend.position = 'none', 
        strip.text.x = element_text(
          size = 14,        # 폰트 크기
          face = "bold",    # 폰트 굵기 (bold, italic 등)
          family = "Arial"))

head(df_tr)
df_tr %>% 
  mutate(pct = ifelse(coding_type == "Protein_coding_trans", n/154726, n/75345)) %>%  #head()
  ggplot(aes(x = coding_type, y = pct, fill = coding_type)) +
  geom_quasirandom(shape = 21, alpha = 0.7, dodge.width = 0.8) +
  labs(y="Percentage annotated") + 
  scale_x_discrete(labels = c("Protein_coding_trans" = "Protein-coding\ntranscripts", 
                              "Noncoding_trans" = "Noncoding\ntranscripts")) +
  theme_step1() +
  theme(axis.title.x = element_blank(),
        legend.position = 'none')


head(df_tr)
df %>% mutate(pct = ifelse(coding_type == "Protein_coding_genes", n/19866, n/39626)) %>%  #head()
  rbind(df_tr %>% mutate(pct = ifelse(coding_type == "Protein_coding_trans", n/154726, n/75345))) %>% #count(coding_type)
  mutate(coding_type = factor(coding_type, levels = c("Protein_coding_genes", 
                                                      "Noncoding_genes", 
                                                      "Protein_coding_trans", 
                                                      "Noncoding_trans"))) %>%  
  ggplot(aes(x = coding_type, y = pct, fill = coding_type)) +
  geom_quasirandom(shape = 21, alpha = 0.7, dodge.width = 0.8) +
  labs(y="Percentage annotated") + 
  scale_x_discrete(labels = c("Protein_coding_genes" = "Protein-coding\ngenes", 
                              "Noncoding_genes" = "Noncoding\ngenes",
                              "Protein_coding_trans" = "Protein-coding\ntranscripts", 
                              "Noncoding_trans" = "Noncoding\ntranscripts")) +
  theme_step1() +
  theme(axis.title.x = element_blank(),
        legend.position = 'none')

  

### copy number
cp_ref <- readxl::read_xlsx("~/Desktop/KCDC/paper/pangenome/Draft human pangenome reference41586_2023_5896_MOESM4_ESM.xlsx",sheet = 9,skip = 1)


cp_num <- read_table("KPPD.liftoff_genecode_chm13.change_copy_number_INFO.txt")
head(cp_num)
table(cp_num$coding_type)
hist(cp_num$copy_num)
cp_num %>% filter(coding_type == "Protein_coding_genes") %>%
  ggplot(aes(x=copy_num)) + 
  geom_histogram()

hist(cp_num$copy_num)

head(cp_num)
cp_num %>% filter(coding_type == "Protein_coding_genes") %>% 
  select(ID,gene) %>% 
  unique() %>% count(gene) %>%
  ggplot(aes(x=n)) + 
  geom_histogram()

head(cp_num)
cp_num %>% filter(gene %in% cp_ref$gene) %>% count(gene)  

cp_num %>% filter(coding_type == "Protein_coding_genes") %>% 
  select(ID,gene) %>% unique() %>% count(gene) %>%
  ggplot(aes(x=n)) +
  geom_histogram() + 
  theme_classic() + 
  labs(y = "Freqeuncy", x = "gene copy number")

head(cp_num)


df_summary <- df_all %>%
  group_by(assembly) %>%
  summarise(duplicated_genes = sum(copies > 1))

cp_num %>% filter(coding_type == "Protein_coding_genes") %>% 
  group_by(ID) %>%
  summarise(duplicated_genes = sum(copy_num>0)) %>%  #head()
  mutate(sampleID = str_split_fixed(ID,"_",2)[,1]) %>% #head()
  mutate(asssembly_level = ifelse(sampleID %in% c("KPPD129","KPPD130","KPPD131","KPPD132"),"T2T-level","HIFI")) %>%
  mutate(asssembly_level = factor(asssembly_level, levels = rev(unique(asssembly_level)))) %>%  # ???? 순서 반전
  ggplot(aes(y = reorder(ID, duplicated_genes), x = duplicated_genes,fill=asssembly_level)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(y = "Assembly", x = "Duplicated genes per genome", title = "Duplicated Genes Across Assemblies") +
  theme_minimal() + 
  theme(axis.text.x = element_blank())
#  theme(legend.position = 'none')

  
cp_num %>% filter(coding_type == "Protein_coding_genes") %>% 
  group_by(ID) %>%
  summarise(duplicated_genes = sum(copy_num > 1)) %>% mutate(sampleID = str_split_fixed(ID,"_",2)[,1]) %>%
  mutate(asssembly_level = ifelse(sampleID %in% c("KPPD129","KPPD130","KPPD131","KPPD132"),"T2T-level","HIFI")) %>% 
  ggplot(aes(x = reorder(assembly, duplicated_genes), y = duplicated_genes, fill = assembly)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Assembly", y = "Duplicated genes per genome", title = "Duplicated Genes Across Assemblies") +
  theme_minimal()


cp_num %>% filter(coding_type == "Protein_coding_genes") %>% 
  group_by(gene) %>%
  summarise(total_copies = sum(copy_num)) %>%
  arrange(desc(total_copies)) %>%
  head(25) -> top25_cp_num

ggplot(top25_cp_num, aes(x = reorder(gene, total_copies), y = total_copies, fill = gene)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Gene family", y = "Total additional copies", title = "Top Gene Families with Duplications") +
  theme_minimal() + 
  theme(legend.position = 'none')
