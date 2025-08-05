library(tidyverse)


df <- read.table("~/Downloads/KBAv1_CEL.txt",header = T)
ref <- read.table("~/Desktop/KCDC/KCHIP_open/reQC/KBAv1_KoGES_QCed_list.txt",header = T)
head(ref)
dim(ref)
head(df)
table(df$USE)
df %>% filter(USE == "USE(7K)")
df %>% filter(USE %in% c("USE(7K)","USE")) %>% count(COHORT,YEAR)
df %>% count(COHORT,YEAR)


dim(ref)

xmtable(df$DATATYPE)
df %>% filter(USE %in% c("USE(7K)","USE")) %>% filter(ID %in% ref$ID) %>% count(COHORT)
120 + 15406 + 7017 + 80
15406 + 7017 + 80
df %>% filter(USE %in% c("USE(7K)","USE")) %>% filter(ID %in% ref$ID) %>% count(COHORT,YEAR)

v2 <- read.table("~/Desktop/KCDC/DATA_open/KBA/KBAv2.0AB_42K_CEL_ID.txt")
ref_2 <- readxl::read_xlsx("~/Desktop/KCDC/DATA_open/KBA/KBAv1.v2_47K_KC°³¹ß°ú_240614.xlsx")
head(ref_2)
dim(ref_2)
head(v2)
dim(v2)
v2 %>% filter(grepl("HC",V4)) %>% dim()
v2 %>% filter(grepl("DL",V4)) %>% dim()
v2 %>% filter(!grepl("DL",V4) & !grepl("HC",V4)) %>% head()


## 
setwd("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/AR_ori/SIRPa/")

kba <- read.table("KBAv2_annotation_v1_202410.SIRPA.txt")
kgp <- read.table("KGP3KRG.SIRP.txt")
kis_8k <- read.table("8Kref.SIRP.txt")

head(kba)
head(kgp)
head(kis_8k)
colnames(kba) <- c("KBAversion", "VariantID", "HG38", "WGS8K_MAF", "gnomAD_EAS_MAF", 
  "FST", "rmDUP", "HG19", "GENE", "FUNC", "dbNSFP", "CLN", 
  "CLINVAR", "CLINVAR_DISEASE", "LOVD", "CA130", "OMIM", "ACMG", 
  "ADME", "PUBDB", "HLA", "GWASGRID", "dbSNP_HG38", "dbSNP_HG19")


head(kba)
head(str_split_fixed)
kba %>% mutate(hg19_pos = str_split_fixed(HG19,":",4)[,2]) %>% select(KBAversion,VariantID,HG19,hg19_pos) %>% head()

head(kgp)
kgp %>% mutate(ref = str_split_fixed(V1,":",4)[,3],alt = str_split_fixed(V1,":",4)[,4]) %>% filter(ref != V3)
kgp %>% mutate(ref = str_split_fixed(V1,":",4)[,3],alt = str_split_fixed(V1,":",4)[,4]) %>% count(ref != V3)
kgp %>% mutate(ref = str_split_fixed(V1,":",4)[,3],alt = str_split_fixed(V1,":",4)[,4]) %>% count(alt != V4)

kgp %>% mutate(kgp_ID = paste0("20:",str_split_fixed(V1,":",2)[,2])) %>% head()

kgp %>% mutate(kgp_ID = paste0("20:",V2,":",V3,":",V4)) %>% filter(str_length(V3) != 1 | str_length(V4) != 1)

kgp %>% mutate(kgp_ID = paste0("20:",V2,":",V3,":",V4)) %>% head()
head(kis_8k)
kis_8k %>% mutate(kis_8k_ID = paste0(V1,":",V3,":",V4)) %>% head()

kba %>% mutate(hg19_pos = str_split_fixed(HG19,":",4)[,2]) %>% select(KBAversion,VariantID,HG19,hg19_pos) -> kba_1
kgp %>% mutate(kgp_ID = paste0("20:",V2,":",V3,":",V4)) -> kgp_1
kis_8k %>% mutate(kis_8k_ID = paste0(V1,":",V3,":",V4)) -> kis_8k_1

head(kba_1)
head(kgp_1)
head(kis_8k_1)


kba_1 %>% mutate(ID = HG19) %>% merge(kgp_1 %>% mutate(ID = kgp_ID,kgp = "kgp"),all=T,by="ID") %>% #head()
  merge(kis_8k_1 %>% mutate(ID = kis_8k_ID,kis_8k = "kis_8k"),all=T,by="ID") %>% #head()
  select(ID,KBAversion,kgp,kis_8k) %>% write.table("SIRPA.variant.list.merge.txt",col.names = T,row.names = F,quote = F,sep = "\t")

#ANN=ALLELE|EFFECT|IMPACT|GENE|GENE_ID|FEATURE_TYPE|FEATURE_ID|TRANSCRIPT_BIOTYPE|...|PROTEIN_CHANGE|...
anno <- read.table("anno/output.vcf")
head(anno)
dim(anno)

anno %>% mutate(V8 = str_split_fixed(V8,"=",2)[,2]) %>%
  separate_rows(V8,sep = ",") %>% #head()
  filter(grepl("NM_080792",V8)) %>% 
  mutate(variant_type = str_split_fixed(V8,"\\|",20)[,2]) %>%
  mutate(anno = str_split_fixed(V8,"\\|",20)[,3]) -> a
head(a)
a %>% filter(grepl("Ser",V8))
a %>% filter(grepl("intron",anno))
a %>% count(variant_type)
a %>% filter(variant_type %in% c("missense_variant","missense_variant&splice_region_variant","synonymous_variant")) %>% head()
a %>% mutate(base_change = str_split_fixed(V8,"\\|",20)[,10]) %>%
  mutate(amino_change = str_split_fixed(V8,"\\|",20)[,11]) %>% filter(amino_change != "") -> a

'''
MEPAG PAPGR LGPLL CLLLA ASCAW 25 SGVAG EEELQ VIQPD KSVLV AAGET
ATLRC TATSL IPVGP IQWFR GAGPG RELIY NQKEG HFPRV TTVSD LTKRNN M
DFSIRIGNITPADAGTYYCVKFRKGSPDDVEFKSGAGTELSVRAKPSAPVVSGPAARA
TPQHTVSFTCESHGFSPRDITLKWFKNGNELSDFQTNVDPVGESVSYSIHSTAKVVLT
REDVHSQVICEVAHVTLQGDPLRGTANLSETIRVPPTLEVTQQPVRAENQVNVTCQVR
KFYPQRLQLTWLENGNVSRTETASTVTENKDGTYNWMSWLLVNVSAHRDDVKLTCQVE
HDGQPAVSKSHDLKVSAHPKEQGSNTAAENTGSNERNIYIVVGVVCTLLVALLMAALY
LVRIRQKKAQGSTSSTRLHEPEKNAREITQDTNDITYADLNLPKGKKPAPQAAEPNNH
TEYASIQTSPQPASEDTLTYADLDMVHLNRTPKQPAPKPEPSFSEYASVQVPRK
'''
#27
# 31 -> 1
sirp <- "MEPAGPAPGRLGPLLCLLLAASCAWSGVAGEEELQVIQPDKSVLVAAGETATLRCTATSLIPVGPIQWFRGAGPGRELIYNQKEGHFPRVTTVSDLTKRNNMDFSIRIGNITPADAGTYYCVKFRKGSPDDVEFKSGAGTELSVRAKPSAPVVSGPAARATPQHTVSFTCESHGFSPRDITLKWFKNGNELSDFQTNVDPVGESVSYSIHSTAKVVLTREDVHSQVICEVAHVTLQGDPLRGTANLSETIRVPPTLEVTQQPVRAENQVNVTCQVRKFYPQRLQLTWLENGNVSRTETASTVTENKDGTYNWMSWLLVNVSAHRDDVKLTCQVEHDGQPAVSKSHDLKVSAHPKEQGSNTAAENTGSNERNIYIVVGVVCTLLVALLMAALYLVRIRQKKAQGSTSSTRLHEPEKNAREITQDTNDITYADLNLPKGKKPAPQAAEPNNHTEYASIQTSPQPASEDTLTYADLDMVHLNRTPKQPAPKPEPSFSEYASVQVPRK"

substring(sirp, 27, 35)

substring(sirp, 31, 35) #E #3 #E or G
a %>% filter(grepl("31",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) 
# Glu31fs -> frameshift
substring(sirp, 42, 42) #S  #12 #S
a %>% filter(grepl("42",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) %>% unique()

substring(sirp, 44, 45) #LV  #14 #L or S
a %>% filter(grepl("44",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) %>% unique()
#c.131T>C    p.Leu44Ser
#c.132G>A    p.Leu44Leu

substring(sirp, 50, 57) #TATLRCTA T:20 A T:22 L R:24 CT A:27
a %>% filter(grepl("50",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) %>% unique()
# T or S
#c.148A>T    p.Thr50Ser 
#c.150A>G    p.Thr50Thr  
a %>% filter(grepl("52",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) %>% unique()
# T or I
#c.155C>T p.Thr52Ile

a %>% filter(grepl("54",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) %>% unique()
#R H L
#c.161G>A    p.Arg54His #R H

a %>% filter(grepl("57",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) %>% unique()
# A or V
#c.170C>T    p.Ala57Val

substring(sirp, 75, 75) #G:45
a %>% filter(grepl("75",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) %>% unique()
# G or A
#c.224G>C    p.Gly75Ala  


substring(sirp, 95, 99) #DLTKR D: L T: K R:
a %>% filter(grepl("95",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) %>% unique()
#c.284_285i¡¦ p.Asp95fs 
#c.285C>G    p.Asp95Glu 

a %>% filter(grepl("96",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) %>% unique()

#c.287_288d¡¦ p.Leu96fs
#c.286C>T    p.Leu96Ph
#c.287T>C    p.Leu96Pro


substring(sirp, 100,100)
# N -> E
a %>% filter(grepl("100",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) %>% unique()
# c.298A>G    p.Asn100Asp
# c.300_301d¡¦ p.Asn100fs
# c.298A>G    p.Asn100Asp
# c.300C>A    p.Asn100Lys


substring(sirp, 105,109) # SIRIG
# 105 107 109
# S-P R-S G-S
a %>% filter(grepl("105",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) %>% unique()
a %>% filter(grepl("107",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) %>% unique()
#c.319C>A    p.Arg107Ser 
a %>% filter(grepl("109",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) %>% unique()
#c.325G>A    p.Gly109Ser
substring(sirp, 129, 130) #PD -> .P
a %>% filter(grepl("129",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) %>% unique()
a %>% filter(grepl("130",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) %>% unique()

#c.387_389delCGA p.Asp130del 

substring(sirp, 132, 132) #V -> T
a %>% filter(grepl("132",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) %>% unique()
#c.394G>A    p.Val132Met
#c.395T>C    p.Val132Ala

substring(sirp, 146, 146) #A -> G
a %>% filter(grepl("146",amino_change)) %>% select(V1,V2,V4,V5,variant_type,anno,base_change,amino_change) %>% unique()
#c.436G>A    p.Ala146Thr 

####
variant_list <- read_table("SIRPA.variant.list.merge.txt")
aa <- c(
  "p.Glu31fs", "p.Leu44Ser", "p.Leu44Leu",
  "p.Thr50Ser", "p.Thr50Thr",
  "p.Thr52Ile", "p.Arg54His",
  "p.Ala57Val", "p.Gly75Ala",
  "p.Asp95fs", "p.Asp95Glu",
  "p.Leu96fs", "p.Leu96Phe","p.Leu96Pro",
  "p.Asn100Asp", "p.Asn100Asp", "p.Asn100fs", 
  "p.Asn100Asp", "p.Asn100Lys","p.Arg107Ser", "p.Gly109Ser",
  "p.Asp130del", "p.Val132Met", "p.Val132Ala","p.Ala146Thr"
)

print(mutation_vector)
head(a)
head(variant_list)
a %>% filter(amino_change %in% aa) %>% unique() %>% #head()-> cc
  mutate(ID = paste0(V1,":",V2,":",V4,":",V5)) %>% left_join(variant_list) %>% unique() %>% #dim()
  writexl::write_xlsx("SIRPA.variant.list.merge.withAnnoInfo.xlsx")

getwd()
