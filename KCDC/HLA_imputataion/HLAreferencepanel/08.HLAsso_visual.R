library(tidyverse)
library(gridExtra)
library(ggpubr)
library(ggthemes)
library(reshape2)

setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/asso/")

df <-read.table("130K.HLAimpt.AR2_DR2_AF.txt")
ref_freq <- readxl::read_xlsx("supplementary_table_ddac016.xlsx",sheet = 2,skip = 2) %>% select(1,3)
head(ref_freq)
colnames(ref_freq) <- c("HLAtype","Kim 2022")
ref_freq$`Kim 2022` <- ref_freq$`Kim 2022`/100
head(df)
colnames(df)<-c("ID","AR2","DR2","AF","AC","AN")
df %>% mutate(HLAtype = str_split_fixed(ID,"HLA_",2)[,2]) %>% 
  mutate(Gene = str_split_fixed(HLAtype,"\\*",2)[,1]) %>% #head()
  filter(grepl(":",HLAtype)) %>% #head()
  select(Gene,HLAtype,AC,AR2,DR2) %>% 
  group_by(Gene) %>%
  mutate(Frequency = prop.table(AC)) -> df_freq

head(df_freq)
head(ref_freq)
table(df_freq$Gene)
table(ref_freq$Gene)
df_freq %>% count(Gene) %>% rename("KBA_13K"=n) -> a
ref_freq %>% mutate(Gene = str_split_fixed(HLAtype,"\\*",2)[,1]) %>% count(Gene) %>% rename("Kim"=n) -> b
df_freq %>% select(Gene,HLAtype)->a
ref_freq %>% mutate(Gene = str_split_fixed(HLAtype,"\\*",2)[,1]) %>% select(Gene,HLAtype) -> b

a %>% inner_join(b) %>% count(Gene) -> t

a %>% inner_join(b) -> t
ref_freq
df_freq %>% mutate(KBA_130K = Frequency) %>%
  select(Gene,HLAtype,KBA_130K) %>% inner_join(ref_freq) %>% #head()
  ggscatter(.,x='KBA_130K',y="Kim 2022",color='Gene',
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            xlab = "KBA 130K",
            ylab = "Kim 2022",
            xlim = c(0,0.5),
            ylim = c(0,0.5),
            #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
            cor.coeff.args = list(method = "pearson", label.sep = "\n")) + 
  theme(legend.position = "right")


df_freq %>% mutate(KBA_130K = Frequency) %>%
  select(Gene,HLAtype,KBA_130K) %>% inner_join(ref_freq) %>% #head()
  ggscatter(.,x='KBA_130K',y="Kim 2022",color='Gene',
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            xlab = "KBA 130K",
            ylab = "Kim 2022",
            xlim = c(0,0.5),
            ylim = c(0,0.5),
            #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
            cor.coeff.args = list(method = "pearson", label.sep = "\n")) + 
  theme(legend.position = "right")
  #facet_grid(~Gene)


HDL <- read.table("KCHIP_130K_HLA_TAPAS_asso_HDL_z.assoc.linear",header = T)
HDL %>% filter(grepl("HLA",SNP)) %>% filter(grepl(":",SNP)) %>% mutate(SNP = str_split_fixed(SNP,"HLA_",2)[,2]) %>% 
   mutate(Trait = "HDL") %>% select(Trait,SNP,BETA,P) -> HDL

TC <- read.table("KCHIP_130K_HLA_TAPAS_asso_TC_z.assoc.linear",header = T)
TC %>% filter(grepl("HLA",SNP)) %>% filter(grepl(":",SNP)) %>% mutate(SNP = str_split_fixed(SNP,"HLA_",2)[,2]) %>% 
  mutate(Trait = "TC") %>% select(Trait,SNP,BETA,P) -> TC

TG <- read.table("KCHIP_130K_HLA_TAPAS_asso_TG_logz.assoc.linear",header = T)
TG %>% filter(grepl("HLA",SNP)) %>% filter(grepl(":",SNP)) %>% mutate(SNP = str_split_fixed(SNP,"HLA_",2)[,2]) %>% 
  mutate(Trait = "TG") %>% select(Trait,SNP,BETA,P) -> TG


head(HDL)
colnames(HDL_ref)
HDL_ref <- read_tsv("phenotypes/HDL_cholesterol/step_01.result.tsv")
HDL_ref %>% filter(grepl("HLA",marker_name_pub)) %>% select(phenotype_name_publication,marker_name_pub,coef,P) %>% 
  mutate(marker_name_pub = str_split_fixed(marker_name_pub,"HLA-",2)[,2]) %>% mutate(phenotype_name_publication = "HDL")-> HDL_ref
head(HDL_ref)
#colnames(HDL_ref) <- c("Trait","HLAtype","Z_HDL_ref","P_HDL_ref")
colnames(HDL_ref) <- c("Trait","HLAtype","coef_ref","P_ref")
TC_ref <- read_tsv("phenotypes/Total_cholesterol/step_01.result.tsv")
TC_ref %>% filter(grepl("HLA",marker_name_pub)) %>% select(phenotype_name_publication,marker_name_pub,coef,P) %>%
  mutate(phenotype_name_publication = "TC") %>%
  mutate(marker_name_pub = str_split_fixed(marker_name_pub,"HLA-",2)[,2]) -> TC_ref
colnames(TC_ref) <- c("Trait","HLAtype","coef_ref","P_ref")
TG_ref <- read_tsv("phenotypes/Triglyceride/step_01.result.tsv")
TG_ref %>% filter(grepl("HLA",marker_name_pub)) %>% select(phenotype_name_publication,marker_name_pub,coef,P) %>%
  mutate(phenotype_name_publication = "TG") %>%
  mutate(marker_name_pub = str_split_fixed(marker_name_pub,"HLA-",2)[,2]) -> TG_ref
colnames(TG_ref) <- c("Trait","HLAtype","coef_ref","P_ref")


head(HDL)
head(HDL_ref)
HDL_ref %>% rbind(TC_ref) %>% rbind(TG_ref) -> ref
head(ref)

HDL %>% rbind(TC) %>% rbind(TG) %>% rename("HLAtype" = SNP) %>% inner_join(ref) %>% #head()
  ggscatter(.,x='P',y="P_ref",color='Trait',
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            xlab = "KBA 130K (P)",
            ylab = "Kim 2022 (P)",
            #xlim = c(0,0.5),
            #ylim = c(0,0.5),
            #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
            cor.coeff.args = list(method = "pearson", label.sep = "\n")) + 
  theme(legend.position = "right") + 
  facet_grid(~Trait)

HDL %>% rbind(TC) %>% rbind(TG) %>% rename("HLAtype" = SNP) %>% inner_join(ref) %>% #head()
  ggscatter(.,x='BETA',y="coef_ref",color='Trait',
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            xlab = "KBA 130K (BETA)",
            ylab = "Kim 2022 (coef)",
            #xlim = c(0,0.5),
            #ylim = c(0,0.5),
            #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
            cor.coeff.args = list(method = "pearson", label.sep = "\n")) + 
  theme(legend.position = "right") + 
  facet_grid(~Trait)


HDL %>% rbind(TC) %>% rbind(TG) %>% rename("HLAtype" = SNP) %>% inner_join(ref) -> out
colnames(out) <- c("Trait","HLAtype","KBA130K_BETA","KBA130K_P-value","Kim2022_BETA","Kim2022_P-avlue")  
writexl::write_xlsx(out,"~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/HLAassociation_Result.xlsx")

#### HLA freqeuncy vs Park

setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/")
park <-readxl::read_xlsx("Result/HLAtype.freq_park2016.xlsx",sheet = 3)
park %>% mutate(Frequency = Frequency/100)%>% 
  rename('Park 2016'=Frequency) %>% mutate(Gene = str_split_fixed(value,"\\*",2)[,1]) %>% 
  select(Gene,value,`Park 2016`) -> park


df <- readxl::read_xlsx("Result/KMHC.HLAtype.freq.xlsx")
head(df)
df %>% rename('KMHC'=prop) %>% select(Gene,value,KMHC) -> df
head(park)

df %>% inner_join(park) %>% 
  ggscatter(.,x='KMHC',y="Park 2016",color='Gene',
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            xlab = "KMHC",
            ylab = "Park 2016",
#            xlim = c(0,0.5),
#            ylim = c(0,0.5),
            #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
            cor.coeff.args = list(method = "pearson", label.sep = "\n")) + 
  theme(legend.position = "right")
#  facet_grid(~Gene)

head(df)
df %>% inner_join(park) %>% count(Gene) %>% summarise(sum(n))
head(park)



park %>% filter(Gene =="DQB1")


