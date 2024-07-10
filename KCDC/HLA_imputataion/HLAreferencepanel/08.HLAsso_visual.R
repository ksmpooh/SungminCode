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
  select(Gene,HLAtype,KBA_130K) %>% inner_join(ref_freq) %>% head()
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


######## after KHU hema

## freq
head(ref_freq)

df <- read.table("after_khu_hema/vcf/130K.HLAimp.minimac4.new.txt")
head(df)
colnames(df) <- c("ID","MAF","R2","AC","AN")


df %>% mutate(HLAtype = str_split_fixed(ID,"HLA_",2)[,2]) %>% 
  mutate(Gene = str_split_fixed(HLAtype,"\\*",2)[,1]) %>% #head()
  filter(grepl(":",HLAtype)) %>% #head()
  select(Gene,HLAtype,AC,R2) %>% 
  group_by(Gene) %>%
  mutate(Frequency = prop.table(AC)) -> df_freq

head(df_freq)
head(ref)

df_freq %>% select(Gene,HLAtype) %>% inner_join(ref) %>%  filter(freq != 'rare') %>% dim()
df_freq %>% select(Gene,HLAtype) %>% inner_join(ref) %>%  filter(freq != 'rare') %>% count(Gene)
df_freq %>% select(Gene,HLAtype) %>% inner_join(ref) %>%  filter(freq != 'rare') %>% #count(Gene)
  filter(HLAtype %in% ref_freq$HLAtype) %>% #dim()
  count(Gene) 


head(df_freq) 
head(ref_freq)
dim(df_freq)
df_freq %>% select(Gene,HLAtype) %>% inner_join(ref) %>%  filter(freq != 'rare') %>% filter(!(HLAtype %in% ref_freq$HLAtype)) -> only_kmhc

  
df_freq %>% inner_join(ref) %>%  filter(freq != 'rare') %>% #head()
  mutate(KBA_130K = Frequency) %>%
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
  
### pheno number of sample


setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/asso/")
#WBC RBC PLAT MCV MCHC MCH HCT HB
phenos <- c("White_blood_cell_count",'Red_blood_cell_count',"Platelet_count",'Mean_corpuscular_volume',"Mean_corpuscular_hemoglobin_concentration","Mean_corpuscular_hemoglobin","Hematocrit","Hemoglobin")


tmp <- read_tsv(paste0("phenotypes/",pheno,"/step_01.result.tsv"))
head(tmp)
df <- NULL
for (pheno in phenos) {
  tmp <- read_tsv(paste0("phenotypes/",pheno,"/step_01.result.tsv"))
  tmp %>% filter(grepl("HLA",marker_name_pub)) %>% select(phenotype_name_publication,`samples(case/control)`) %>% 
    #mutate(marker_name_pub = str_split_fixed(marker_name_pub,"HLA-",2)[,2]) %>% 
    mutate(phenotype_name_publication = pheno)-> tmp
  head(tmp)
  colnames(tmp) <- c("Trait","nsample_ref")
  df <- rbind(df,tmp)
}

ref_kim <- df
head(ref_kim)
ref_kim %>% unique() -> ref_kim

flist = grep(list.files("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/asso/after_khu_hema/minimac4_imp/"),pattern = "assoc", value=TRUE)
#tmp <- read_table(paste0("./after_khu_hema/minimac4_imp/",flist[1]))
head(tmp)
df <- NULL
for (i in flist) {
  tmp <- read_table(paste0("./after_khu_hema/minimac4_imp/",i))
  tmp %>% #filter(grepl("HLA",SNP)) %>% filter(grepl(":",SNP)) %>% mutate(SNP = str_split_fixed(SNP,"HLA_",2)[,2]) %>% 
    mutate(Trait = str_remove(i,"KCHIP_130K.KMHC_v1_ref.eagle_phasing.minimac4_imputation.HLA_TAPAS_asso.")) %>%
    mutate(Trait = str_remove(Trait,".assoc.linear")) %>% select(Trait,NMISS) %>% unique()-> tmp
  df <- rbind(df,tmp)
}
head(df)
table(df$Trait)
table(ref_kim$Trait)

df %>% mutate(Trait = recode(Trait,"WBC" = "White_blood_cell_count",
                             "RBC" = "Red_blood_cell_count",
                             "PLAT" = "Platelet_count",
                             "MCV" = "Mean_corpuscular_volume",
                             "MCHC" = "Mean_corpuscular_hemoglobin_concentration",
                             "MCH" = "Mean_corpuscular_hemoglobin",
                             "HCT" = "Hematocrit",
                             "HB" = "Hemoglobin")) %>%
  rename("nsample" = NMISS)-> df
head(df)
head(ref_kim)

df %>% inner_join(ref_kim) -> a


####asso

head(only_kmhc)
setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/asso/")
#WBC RBC PLAT MCV MCHC MCH HCT HB
phenos <- c("White_blood_cell_count",'Red_blood_cell_count',"Platelet_count",'Mean_corpuscular_volume',"Mean_corpuscular_hemoglobin_concentration","Mean_corpuscular_hemoglobin","Hematocrit","Hemoglobin")


tmp <- read_tsv(paste0("phenotypes/",pheno,"/step_01.result.tsv"))
head(tmp)
df <- NULL
for (pheno in phenos) {
  tmp <- read_tsv(paste0("phenotypes/",pheno,"/step_01.result.tsv"))
  tmp %>% filter(grepl("HLA",marker_name_pub)) %>% select(phenotype_name_publication,marker_name_pub,coef,P) %>% 
    mutate(marker_name_pub = str_split_fixed(marker_name_pub,"HLA-",2)[,2]) %>% mutate(phenotype_name_publication = pheno)-> tmp
  head(tmp)
  colnames(tmp) <- c("Trait","HLAtype","coef_ref","P_ref")
  df <- rbind(df,tmp)
}

ref_kim <- df


flist = grep(list.files("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/asso/after_khu_hema/minimac4_imp/"),pattern = "assoc", value=TRUE)

df <- NULL
for (i in flist) {
  tmp <- read_table(paste0("./after_khu_hema/minimac4_imp/",i))
  tmp %>% filter(grepl("HLA",SNP)) %>% filter(grepl(":",SNP)) %>% mutate(SNP = str_split_fixed(SNP,"HLA_",2)[,2]) %>% 
    mutate(Trait = str_remove(i,"KCHIP_130K.KMHC_v1_ref.eagle_phasing.minimac4_imputation.HLA_TAPAS_asso.")) %>%
    mutate(Trait = str_remove(Trait,".assoc.linear")) %>% select(Trait,SNP,BETA,P) -> tmp
  df <- rbind(df,tmp)
}
head(df)
table(df$Trait)
table(ref_kim$Trait)

df %>% mutate(Trait = recode(Trait,"WBC" = "White_blood_cell_count",
                             "RBC" = "Red_blood_cell_count",
                             "PLAT" = "Platelet_count",
                             "MCV" = "Mean_corpuscular_volume",
                             "MCHC" = "Mean_corpuscular_hemoglobin_concentration",
                             "MCH" = "Mean_corpuscular_hemoglobin",
                             "HCT" = "Hematocrit",
                             "HB" = "Hemoglobin")) -> df

ref <- read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLAtype_check/KMHC_HLAtype_frequency.txt",header = T)
ref_after_imp <- read.table("after_khu_hema/vcf/130K.HLAimp.minimac4.new.txt")
head(ref_after_imp)
ref_after_imp %>% filter(grepl("HLA",V1)) %>% mutate(V1 = str_remove(V1,"HLA_")) %>% select(V1,V3) %>%
  rename("HLAtype" = V1)


head(ref)
ref %>% select(value,freq) %>% rename("HLAtype" =value) -> ref


head(df)
head(ref_kim)
df %>% rename("HLAtype" = SNP) %>% inner_join(ref_kim) -> out
#colnames(out) <- c("Trait","HLAtype","KBA130K_BETA","KBA130K_P-value","Kim2022_BETA","Kim2022_P-avlue")  
out <- out %>% left_join(ref)

head(out)
table(out$freq)
out %>% filter(freq != 'rare') %>% #head()
  filter(P <=0.05) %>% #head()
  mutate(P = -log10(P),P_ref = -log10(P_ref)) %>%
  ggscatter(.,x='P',y="P_ref",color='freq',
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            xlab = "KBA 130K (-logP)",
            ylab = "Kim 2022 (-logP)",
            #xlim = c(0,0.5),
            #ylim = c(0,0.5),
            #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
            cor.coeff.args = list(method = "pearson", label.sep = "\n")) + 
  theme(legend.position = "right") + 
  #facet_grid(~Trait,)
  facet_wrap(~Trait,ncol = 4,scales = 'free')


out %>% filter(freq != 'rare') %>%
  ggscatter(.,x='BETA',y="coef_ref",color='freq',
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
  #facet_grid(~Trait,)
  facet_wrap(~Trait,ncol = 4,scales = 'free')

head(out)
# only 4 type
head(df)
df %>% rename("HLAtype" = SNP) %>% mutate(only = ifelse(HLAtype %in% only_kmhc$HLAtype,HLAtype,"Intersect")) %>% #arrange(desc(only)) %>% #head()
  #mutate(size = ifelse(only == "Unique",1,0.5)) %>%
  ggplot(aes(x=BETA,y=-log10(P),color = only)) +
  geom_point()+
  #scale_color_manual(values=c("grey","blue","red","green","purple")) +
  scale_color_manual(values=c("blue","red","green","purple",rgb(0.6,0.6,0.6,0.2))) + 
  scale_size_continuous(range = c(1,3,3,3,3)) + 
  theme(legend.position = "bottom",
        legend.title = element_blank()) + 
  facet_wrap(~Trait,ncol = 4,scales = 'free')
  

