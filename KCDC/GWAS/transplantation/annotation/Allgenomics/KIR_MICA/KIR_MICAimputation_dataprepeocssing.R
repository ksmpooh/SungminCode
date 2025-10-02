### KIR imputation data preprocessing
library(tidyverse)
setwd("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/KIR/KIRIMP_UK")
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


kir_allele <- read_csv("alleles_kir.csv")
head(kir_allele)
dim(kir_allele)
kir_allele %>% count(locus)
dim(kir_allele)

snps_allele <- read_csv("alleles_snp.csv")
head(snps_allele)
colnames(snps_allele)
cor_val <- cor(snps_allele$allele1_frequency_reference,
               snps_allele$allele1_frequency_input,
               method = "pearson")



library(smplot2)

snps_allele %>% ggplot(aes(x=allele1_frequency_reference,y=allele1_frequency_input)) +
  geom_point() + 
#  geom_smooth(method = "lm", se = FALSE, color = "red") +   # 상관 직선
    labs(title = "SNPs frequency comparison\n (common SNPs: 170)",x="KOTRY allele frequency(n=3,817)",y="Reference allele frequency: KIR*IMP") + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray40", size = 0.8) +  # ??? y=x ??????????????????
  sm_statCorr(
    color = "#0f993d", corr_method = "spearman",
    linetype = "dashed",
    text_size = 5
  ) +  theme_step1()

kir_imp_result <- read.csv("imputations.csv")
head(kir_imp_result)
kir_imp_result %>% count(locus)

kir_imp_result %>% #filter(!(locus %in% c("AvsB","KIRhaplotype"))) %>% 
  ggplot(aes(x=locus,y=posteriorProbability,fill=locus)) + 
  geom_violin() + 
  theme_step1() + 
  theme(legend.position = "none",
        axis.text = element_text(family = 'Arial', size = 13, color = 'black',angle = 90))

kir_imp_result %>% filter(locus == "KIRhaplotype") %>% count(imputedType)
kir_imp_result %>% filter(locus == "AvsB") %>% count(imputedType)
kir_imp_result %>% count(locus)
kir_imp_result %>% filter(locus == "KIR2DS2") %>% count(imputedType)
kir_imp_result %>% filter(!(locus %in% c("KIRhaplotype","AvsB"))) %>% count(imputedType)

kir_imp_result %>% head()
kir_imp_result %>% mutate(haplotype = paste0("hap.",str_split_fixed(haplotypeID,"\\.",3)[,3])) %>% rename(ID = ID_1) %>% 
  select(ID,locus,haplotype,imputedType,posteriorProbability) -> kir_imp_result


ref <- read_table("~/Desktop/KCDC/transplantation/00.sampleInfo/JG.IDupdate.NIHtobCODE.txt")
prof_yang <- read.table("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/AR_YS/pair_33/KRKD.YS.33pair.SIRPa.list.txt",header = T)
head(prof_yang)
prof_yang %>% select(KBA_ID.KR) %>% rename(ID = KBA_ID.KR) %>% rbind(prof_yang %>% select(KBA_ID.KD) %>% rename(ID = KBA_ID.KD)) -> prof_yang_sapmle66



head(ref)
kir_imp_result %>% filter(ID %in% prof_yang_sapmle66$ID) %>% write.table("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/AR_YS/pair_33/KRKD.YS.33pair.KIRIMP.posteriorProbability.txt",col.names = T,row.names = F,quote = F,sep = "\t")

kir_imp_result %>% select(-posteriorProbability) %>%  
  filter(ID %in% prof_yang_sapmle66$ID) %>% mutate(type = str_split_fixed(ID,"_",3)[,3]) %>% #count(type)
  #mutate(ID = str_split_fixed(ID,"_K",2)[,1]) %>% #head()
  pivot_wider(names_from = locus:haplotype,values_from = imputedType) %>% 
  arrange(ID) %>% writexl::write_xlsx("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/AR_YS/pair_33/KRKD.YS.33pair.KIRIMP.reformat.xlsx")
  #pivot_wider(names_from = locus:haplotype,values_from = imputedType) %>% head()
  #filter(grepl("kotry_2884",ID)) %>% select(ID)

### MICA imptation

setwd("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/MICA/aftermichiganIMP/")
mica_after_PMRA <- read.table("results_PMRA/MICA_imputed.tsv",header = T)
mica_after_V <- read.table("results_V/MICA_imputed.tsv",header = T)
mica_after_VI <- read.table("results_VI/MICA_imputed.tsv",header = T)
mica_after_VII <- read.table("results_PMRA/MICA_imputed.tsv",header = T)


setwd("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/MICA/beforeIMP/")
mica_PMRA <- read.table("results_PMRAmodel/MICA_imputed.tsv",header = T)
mica_V <- read.table("results_V/MICA_imputed.tsv",header = T)
mica_VI <- read.table("results_VI/MICA_imputed.tsv",header = T)
mica_VII <- read.table("results_1000G_FIN_VII/MICA_imputed.tsv",header = T)


head(mica_after_PMRA)

mica_after_PMRA$batch = "MICA_PMRA_afterIMP"
mica_after_V$batch = "MICA_V_afterIMP"
mica_after_VI$batch = "MICA_VI_afterIMP"
mica_after_VII$batch = "MICA_VII_afterIMP"

mica_PMRA$batch = "MICA_PMRA"
mica_V$batch = "MICA_V"
mica_VI$batch = "MICA_VI"
mica_VII$batch = "MICA_VII"
head(mica_after_PMRA)

mica_after_PMRA %>% rbind(mica_after_V) %>% rbind(mica_after_VI) %>% rbind(mica_after_VII) %>%
  rbind(mica_PMRA) %>% rbind(mica_V) %>% rbind(mica_VI) %>% rbind(mica_VII) %>% #count(batch)
  #filter(grepl("IMP",batch))
  mutate(QC = ifelse(grepl("IMP",batch),"Imputation","Genotype")) %>% #head()
  mutate(model = str_split_fixed(batch,"_",3)[,2]) %>% #count(model)
  ggplot(aes(x=model,y=prob,fill=QC)) + 
  geom_violin() + 
  theme_step1() + 
  theme(#legend.position = "none",
        axis.text = element_text(family = 'Arial', size = 13, color = 'black'))

  


setwd("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/MICA/")
mica <- read.table("aftermichiganIMP/results_PMRA/MICA_imputed.tsv",header = T) 
mica  %>% rename(ID = sample.id) %>% rename(MICA.hap1 = allele1,MICA.hap2 = allele2) %>%   filter(ID %in% prof_yang_sapmle66$ID) %>%  
  arrange(ID) %>% writexl::write_xlsx("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/AR_YS/pair_33/KRKD.YS.33pair.MICA_IMP.rawResult.xlsx")
micb <- read.table("aftermichiganIMP/results_PMRA/MICB_imputed.tsv",header = T)
micb %>% rename(ID = sample.id) %>% rename(MICB.hap1 = allele1,MICB.hap2 = allele2) %>%
  filter(ID %in% prof_yang_sapmle66$ID) %>%  
  arrange(ID) %>% writexl::write_xlsx("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/AR_YS/pair_33/KRKD.YS.33pair.MICB_IMP.rawResult.xlsx")

mica %>% rename(ID = sample.id) %>% rename(MICA.hap1 = allele1,MICA.hap2 = allele2) %>% select(ID,MICA.hap1,MICA.hap2) -> mica
micb %>% rename(ID = sample.id) %>% rename(MICB.hap1 = allele1,MICB.hap2 = allele2) %>% select(ID,MICB.hap1,MICB.hap2) -> micb
  
mica %>% left_join(micb) %>% #head()
  filter(ID %in% prof_yang_sapmle66$ID) %>% 
  mutate(type = str_split_fixed(ID,"_",3)[,3]) %>%
  mutate(ID = str_split_fixed(ID,"_K",2)[,1]) %>% #head()
  arrange(ID) %>% 
  select(c("ID", "type"), everything()) %>% writexl::write_xlsx("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/AR_YS/pair_33/KRKD.YS.33pair.MICAB_IMP.prep.xlsx")


####KOTRY all sample QC
#sampleID QC
ys_cel_info <- read.table("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/sample_info/ys.sampleID.final_20250529.txt")
head(ys_cel_info)
tail(ys_cel_info)
ys_cel_info %>% mutate(KBA_ID = ifelse(grepl("NIH",V1),str_replace(str_split_fixed(V1,"_",6)[,6],".CEL",""),V3)) %>% #head()
  select(KBA_ID,V3) %>% rename(ID = V3) -> ys_id_info
  






library(tidyverse)
setwd("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/KIR/KIRIMP_UK")
ref <- read_table("~/Desktop/KCDC/transplantation/00.sampleInfo/JG.IDupdate.NIHtobCODE.txt")
final_id <- read.table("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/SIRPa/new_pair/KOTRY.AR.SIRPA.txt",header = T)
final_pair <- read.table("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/SIRPa/new_pair/KRKD.KOTRY.AR.SIRPA.txt",header = T) %>% select(KBA_ID.KR,KBA_ID.KD)

head(final_pair)
dim(final_pair)
pheno <- readxl::read_xlsx("~/Desktop/KCDC/transplantation/allogenomic/phenotype/DATASET_20250315_Final.xlsx") 
head(pheno)
head(ref)
head(final_id)
final_id %>% rename(KBA_ID = ID) %>% left_join(ref) %>% select(KBA_ID,bCODE) ->final_id
pheno %>% filter(bCODE_R %in% final_id$bCODE) ->  pheno
dim(pheno)

kir_imp_result <- read.csv("imputations.csv")

head(kir_imp_result)
kir_imp_result %>% filter(!grepl("NIH",ID_1))
head(ys_id_info)
kir_imp_result %>% filter(ID_1 == 2)


final_id %>% filter(!(KBA_ID %in% check$KBA_ID))

#KBA_ID            bCODE  IID->2
#1 NIH19KT6374 41373386DNA01102

kir_imp_result %>% left_join(ys_id_info %>% rename(ID_1 = ID)) %>% #head()
  mutate(KBA_ID = ifelse(is.na(KBA_ID),ID_1,KBA_ID)) %>% mutate(KBA_ID = ifelse(KBA_ID == 2,"NIH19KT6374",KBA_ID)) %>% #head()
  filter(KBA_ID %in% final_id$KBA_ID) %>% #count(KBA_ID) %>% dim()
  left_join(ref) %>% mutate(haplotype = paste0("hap.",str_split_fixed(haplotypeID,"\\.",3)[,3])) %>%
  select(bCODE,locus,haplotype,imputedType,posteriorProbability) ->  kir_imp_result_1441pair_raw
  
  
kir_imp_result_1441pair_raw %>% write.table("KOTRY.2882sample.KIRimp.rawData.txt",col.names = T,row.names = F,quote = F,sep = "\t")
kir_imp_result_1441pair_raw %>% writexl::write_xlsx("KOTRY.2882sample.KIRimp.rawData.xlsx")
head(final_pair)
dim(final_pair)
dim(ref)
head(ref)
#final_pair
kir_imp_result_1441pair_raw %>% head()
  
head(pheno)
head(final_pair)
head(ref)
head(pheno)
ref %>% filter(bCODE == "02910712DNA01102") #NIH19KT0023
final_pair %>% filter(KBA_ID.KR == "NIH19KT0023")

kir_imp_result_1441pair_raw %>% select(-posteriorProbability) %>% #head()
  pivot_wider(names_from = locus:haplotype,values_from = imputedType) %>% write.table("KOTRY.2882sample.KIRimp.prepocessing.txt",col.names = T,row.names = F,quote = F,sep = "\t")

kir_imp_result_1441pair_raw %>% select(-posteriorProbability) %>% #head()
  pivot_wider(names_from = locus:haplotype,values_from = imputedType) %>% writexl::write_xlsx("KOTRY.2882sample.KIRimp.prepocessing.xlsx")
  

#### MICA
setwd("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/MICA/")
mica <- read.table("aftermichiganIMP/results_PMRA/MICA_imputed.tsv",header = T)
micb <- read.table("aftermichiganIMP/results_PMRA/MICB_imputed.tsv",header = T)
head(mica)

mica %>% mutate(gene = "MICA") %>% left_join(micb %>% mutate(gene = "MICB")) %>% 
  left_join(ys_id_info %>% rename(sample.id = ID)) %>% #head()
  mutate(KBA_ID = ifelse(is.na(KBA_ID),sample.id,KBA_ID)) %>% mutate(KBA_ID = ifelse(KBA_ID == 2,"NIH19KT6374",KBA_ID)) %>% #head()
  filter(KBA_ID %in% final_id$KBA_ID) %>% #count(KBA_ID) %>% #dim()
  left_join(ref) %>% 
  select(bCODE,gene,prob,matching,allele1,allele2) %>% writexl::write_xlsx("KOTRY.2882sample.MICAimp.prepocessing.xlsx")
