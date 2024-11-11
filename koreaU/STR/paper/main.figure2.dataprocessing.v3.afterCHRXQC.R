##### figure 2  re analysis chrX with gender
#### circus plot + aano
library(tidyverse)
setwd("/Users/ksmpooh/Desktop/KU/@research/STR/figure/figure2_withchrX/")

final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR

final_ref %>% select(MOTIFS,ID) %>% rename(STR_ID = ID) %>% 
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> final_ref_pro

head(final_ref_pro)
head(final_ref)
final_ref %>% filter(chrom == "chrX") -> final_ref_chrX
head(final_ref_chrX)

sex_info <- read_table("~/Desktop/KCDC/pangenome/00.datacheck/Revio.WGS.sex.info.txt") %>% select(Revio,sex) %>% rename(ID = Revio)
head(sex_info)
sex_info %>% filter(sex == "F") -> sex_info_female
sex_info %>% count(sex)


eh_trgt_merge_simple_pass_intersect_forConcordance <- read_table("~/Desktop/KU/@research/STR/figure/figure2/eh_trgt_merge_simple_pass_intersect_forConcordance.txt")
eh_trgt_simpleSTR_prep_chrX.txt <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/figure/chrX/eh_trgt_simpleSTR_prep_chrX.txt")
head(eh_trgt_simpleSTR_prep_chrX.txt)
#eh_trgt_simpleSTR_prep_chrX.txt %>% count(ID)
eh_trgt_merge_simple_pass_intersect_forConcordance %>% count(EH_STR1 <= EH_STR2)
#eh_trgt_simpleSTR_prep_chrX.txt %>% count(EH_STR1 <= EH_STR2)
head(eh_trgt_merge_simple_pass_intersect_forConcordance_withoutchrX)
eh_trgt_merge_simple_pass_intersect_forConcordance_withoutchrX %>% count(EH_STR1 <= EH_STR2)
#eh_trgt_merge_simple_pass_intersect_forConcordance_withoutchrX %>% count(TRGT_STR1 <= TRGT_STR2)

head(eh_trgt_merge_simple_pass_intersect_forConcordance)
head(eh_trgt_simpleSTR_prep_chrX.txt)
eh_trgt_merge_simple_pass_intersect_forConcordance %>% filter(!(STR_ID %in% final_ref_chrX$ID)) -> eh_trgt_merge_simple_pass_intersect_forConcordance_withoutchrX
eh_trgt_merge_simple_pass_intersect_forConcordance %>% filter((STR_ID %in% final_ref_chrX$ID)) %>% filter(ID %in% sex_info_female$ID) -> eh_trgt_merge_simple_pass_intersect_forConcordance_onlychrX_female
head(eh_trgt_merge_simple_pass_intersect_forConcordance_onlychrX_female)
head(eh_trgt_merge_simple_pass_intersect_forConcordance_withoutchrX)
#eh_trgt_merge_simple_pass_intersect_forConcordance_onlychrX_female %>% count(TRGT_STR1 <=  TRGT_STR2)
eh_trgt_merge_simple_pass_intersect_forConcordance_onlychrX_female %>% rbind(eh_trgt_merge_simple_pass_intersect_forConcordance_withoutchrX) %>%
  select(ID,STR_ID,TRGT_STR1,EH_STR1,TRGT_AM1,TRGT_AP1) %>% 
  rename(TRGT_STR = TRGT_STR1,EH_STR= EH_STR1,TRGT_AM = TRGT_AM1,TRGT_AP = TRGT_AP1) %>% mutate(allele = 1) -> df1


eh_trgt_merge_simple_pass_intersect_forConcordance_onlychrX_female %>% rbind(eh_trgt_merge_simple_pass_intersect_forConcordance_withoutchrX) %>%
  select(ID,STR_ID,TRGT_STR2,EH_STR2,TRGT_AM2,TRGT_AP2) %>% 
  rename(TRGT_STR = TRGT_STR2,EH_STR= EH_STR2,TRGT_AM = TRGT_AM2,TRGT_AP = TRGT_AP2) %>% mutate(allele = 2) -> df2

eh_trgt_merge_simple_pass_intersect_forConcordance_onlychrX_female %>% rbind(eh_trgt_merge_simple_pass_intersect_forConcordance_withoutchrX) %>% 
  filter(str_detect(STR_ID,"chrX")) %>% count(ID) %>% dim()

head(df1)
head(df2)
head(eh_trgt_simpleSTR_prep_chrX.txt)
eh_trgt_simpleSTR_prep_chrX.txt %>% count(ID) %>% count(n)
df1 %>% rbind(df2) %>% filter(str_detect(STR_ID,"chrX")) %>% count(ID) %>% dim()


eh_trgt_simpleSTR_prep_chrX.txt %>% rename(TRGT_AP = AP,TRGT_AM = AM) %>% #head()
  select(ID,STR_ID,TRGT_STR,EH_STR,TRGT_AM,TRGT_AP) %>% mutate(allele = 1) %>% dim() #241410

eh_trgt_simpleSTR_prep_chrX.txt %>% rename(TRGT_AP = AP,TRGT_AM = AM) %>% #head()
  select(ID,STR_ID,TRGT_STR,EH_STR,TRGT_AM,TRGT_AP) %>% mutate(allele = 1) %>% unique() %>% dim() #241410

eh_trgt_simpleSTR_prep_chrX.txt[duplicated(eh_trgt_simpleSTR_prep_chrX.txt),]


head(eh_trgt_simpleSTR_prep_chrX.txt)
eh_trgt_simpleSTR_prep_chrX.txt %>% rename(TRGT_AP = AP,TRGT_AM = AM) %>% #head()
  select(ID,STR_ID,TRGT_STR,EH_STR,TRGT_AM,TRGT_AP) %>% mutate(allele = 1) %>% rbind(df1) %>% rbind(df2) %>% 
  write.table("eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")

####
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC <- read_table("eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC.txt")
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
#eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% 
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID != "RFC1") %>%
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% group_by(ID) %>% 
  summarise(concordance_rate = mean(check)) %>% write.table("f2.concordance_byID_simpleSTR.afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")



eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID != "RFC1") %>%
  select(ID:EH_STR) %>% 
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% group_by(STR_ID) %>% 
  summarise(concordance_rate = mean(check)) -> concordacne_bySTR
head(concordacne_bySTR)
write.table(concordacne_bySTR,"f2.concordance_bySTR_simpleSTR.afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")

head(concordacne_bySTR)
head(concordacne_bySTR)
#concordacne_bySTR %>% select(STR_ID) %>% unique() %>% dim()
concordacne_bySTR %>%
  mutate(concordance_rate_range = case_when(
  concordance_rate == 0 ~ "0",
  concordance_rate > 0 & concordance_rate < 0.1 ~ "(0,0.1)",
  concordance_rate >= 0.1 & concordance_rate < 0.2 ~ "[0.1,0.2)",
  concordance_rate >= 0.2 & concordance_rate < 0.3 ~ "[0.2,0.3)",
  concordance_rate >= 0.3 & concordance_rate < 0.4 ~ "[0.3,0.4)",
  concordance_rate >= 0.4 & concordance_rate < 0.5 ~ "[0.4,0.5)",
  concordance_rate >= 0.5 & concordance_rate < 0.6 ~ "[0.5,0.6)",
  concordance_rate >= 0.6 & concordance_rate < 0.7 ~ "[0.6,0.7)",
  concordance_rate >= 0.7 & concordance_rate < 0.8 ~ "[0.7,0.8)",
  concordance_rate >= 0.8 & concordance_rate < 0.9 ~ "[0.8,0.9)",
  concordance_rate >= 0.9 & concordance_rate < 1 ~ "[0.9,1)",
  concordance_rate == 1 ~ "1")) %>% count(concordance_rate_range) %>% #write.table("f2.concordance.range.bySTR_simpleSTR.afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")
  mutate(facet = ifelse(concordance_rate_range %in% c("[0.9,1)","1"),2,1)) %>% #head()
  ggplot(aes(x=factor(concordance_rate_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                                       "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                                       "[0.8,0.9)", "[0.9,1)","1")),y=n, fill = n)) + 
  geom_bar(stat = 'identity') + 
  theme_bw() + 
  scale_y_continuous(labels = scales::comma) + 
  labs(y="# of STRs",x="Concordance") + 
  geom_text(aes(label = scales::comma(n), y = n), # ???????????? ???????????? ???????????????? ???????????????? ???????????? y?????? n???????????? ?????? ???????? ????????????
            size = 5, family = "Arial", vjust = 0) +
  facet_row(~ facet, scales = "free", space = "free") + 
  theme(strip.text = element_blank(),
        legend.position = "none")


#eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC <- read_table("eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC.txt")
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID != "RFC1") %>%
  select(ID:EH_STR) %>%
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% group_by(STR_ID) %>% 
  summarise(concordance_rate = mean(check)) -> concordacne_bySTR

write.table(concordacne_bySTR,"f2.concordance_bySTR_simpleSTR.afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")

### length
dim(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
head(final_ref_chrX)
head(final_ref_pro)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% select(ID:EH_STR) %>% filter(STR_ID != "RFC1") %>%
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% left_join(final_ref_pro) %>%
  mutate(STR_length = TRGT_STR * RU.length) %>%
  mutate(STR_length = case_when(
    STR_length < 50 ~ "[0~50)",
    STR_length < 100 ~ "[50~100)",
    STR_length < 150 ~ "[100~150)",
    TRUE ~ '[150~Inf)')) %>% count(STR_length) -> a

head(a)
  
a %>% write.table("f2.STR.allele.count.bySTR_length.txt",col.names = T,row.names = F,quote = F,sep = "\t")

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% select(ID:EH_STR) %>% filter(STR_ID != "RFC1") %>%
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% left_join(final_ref_pro) %>%
  mutate(STR_length = TRGT_STR * RU.length) %>%
  mutate(STR_length = case_when(
    STR_length < 50 ~ "[0~50)",
    STR_length < 100 ~ "[50~100)",
    STR_length < 150 ~ "[100~150)",
    TRUE ~ '[150~Inf)')) %>% count(STR_length,check) -> b
b %>% group_by(STR_length) %>%
  mutate(prop = prop.table(n)) %>% write.table("f2.STR.allele.count.match.bySTR_length.txt",col.names = T,row.names = F,quote = F,sep = "\t")



eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% select(ID:EH_STR) %>% filter(STR_ID != "RFC1") %>%
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% left_join(final_ref_pro) %>%
  mutate(STR_length = TRGT_STR * RU.length) %>%
  mutate(STR_length = case_when(
    STR_length < 50 ~ "[0~50)",
    STR_length < 100 ~ "[50~100)",
    STR_length < 150 ~ "[100~150)",
    TRUE ~ '[150~Inf)')) %>% group_by(ID,STR_length) %>% 
  summarise(concordance_rate = mean(check)) -> c

c %>% write.table("f2.sample.concordance.rate.bySTR_length.txt",col.names = T,row.names = F,quote = F,sep = "\t")


eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% select(ID:EH_STR) %>% filter(STR_ID != "RFC1") %>%
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% left_join(final_ref_pro) %>%
  mutate(STR_length = TRGT_STR * RU.length) %>%
  mutate(STR_length = case_when(
    STR_length < 50 ~ "[0~50)",
    STR_length < 100 ~ "[50~100)",
    STR_length < 150 ~ "[100~150)",
    TRUE ~ '[150~Inf)')) %>%
  mutate(GC = case_when(
  GC < 0.25 ~ "[0~0.25)",
  GC < 0.5 ~ "[0.25~0.5)",
  GC < 0.75 ~ "[0.5~0.75)",
  TRUE ~ '[0.75~1]')) %>% #count(GC)
  group_by(ID,STR_length,GC) %>% 
  summarise(concordance_rate = mean(check)) -> d
d %>% write.table("f2.sample.concordance.rate.bySTR_length_GC.txt",col.names = T,row.names = F,quote = F,sep = "\t")
head(d)

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% select(ID:EH_STR) %>% filter(STR_ID != "RFC1") %>%
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% left_join(final_ref_pro) %>%
  mutate(STR_length = TRGT_STR * RU.length) %>%
  mutate(STR_length = case_when(
    STR_length < 50 ~ "[0~50)",
    STR_length < 100 ~ "[50~100)",
    STR_length < 150 ~ "[100~150)",
    TRUE ~ '[150~Inf)')) %>%
  mutate(GC = case_when(
    GC < 0.25 ~ "[0~0.25)",
    GC < 0.5 ~ "[0.25~0.5)",
    GC < 0.75 ~ "[0.5~0.75)",
    TRUE ~ '[0.75~1]')) %>% #count(GC)
  group_by(STR_length,GC) %>% 
  count(check) -> e
head(e)
e %>% write.table("f2.match.count.bySTR_length_GC.txt",col.names = T,row.names = F,quote = F,sep = "\t")


### pandepth GC
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC <- read_table("eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC.txt")
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
head(common_STR)
head(eh_notinter_pandepth)
eh_pandepth_common <- read_table("~/Desktop/KU/@research/STR/eh/pandepth.merge.from.eh.oribam.txt") %>% filter(STR_ID %in% common_STR$STR_DB) %>% select(ID,STR_ID,type,`GC(%)`,MeanDepth)
trgt_pandepth_common <- read_table("~/Desktop/KU/@research/STR/trgt/pandepth/pandepth.merge.from.trgt.oribam.txt") %>% filter(STR_ID %in% common_STR$STR_DB) %>% select(ID,STR_ID,type,`GC(%)`,MeanDepth)
head(eh_pandepth_common)
head(trgt_pandepth_common)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% left_join(trgt_pandepth_common) %>% dim() # 40257880
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% left_join(trgt_pandepth_common) %>% na.omit() %>% dim()

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% left_join(trgt_pandepth_common) -> eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
#head(final_ref_pro)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% 
  left_join(final_ref_pro) %>%
  mutate(STR_length = TRGT_STR * RU.length) %>%
  mutate(STR_length = case_when(
    STR_length < 50 ~ "[0~50)",
    STR_length < 100 ~ "[50~100)",
    STR_length < 150 ~ "[100~150)",
    TRUE ~ '[150~Inf)')) %>%
  mutate(GC = case_when(
  #  is.na(`GC(%)`) ~ "NA",
    `GC(%)` < 0.25 ~ "[0~0.25)",
    `GC(%)` < 0.5 ~ "[0.25~0.5)",
    `GC(%)` < 0.75 ~ "[0.5~0.75)",
    TRUE ~ '[0.75~1]')) %>% #count(GC)
  group_by(ID,STR_length,GC) %>% 
  summarise(concordance_rate = mean(check)) -> d

head(d)

d %>%
  mutate(GC = case_when(
    GC == "[0~0.25)" ~ "0~25",
    GC == "[0.25~0.5)" ~ "25~50",
    GC == "[0.5~0.75)" ~ "50~75",
    TRUE ~ '75~100')) %>% #count(GC)
  ggplot(aes(x=factor(GC,levels=c("0~25","25~50","50~75","75~100")),
             y=concordance_rate,fill=factor(GC,levels=c("0~25","25~50","50~75","75~100")))) + 
  geom_violin() + 
  labs(x="GC content (%)",y="Concordance",title = "STR Length") + 
  theme(legend.position = 'none') + 
  facet_grid(~factor(STR_length,levels=c("[0~50)","[50~100)","[100~150)","[150~Inf)"))) +
  theme_bw() +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +  # y=0.5 ���� �߰�
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 18,family = 'Arial'),
        plot.margin = unit(c(10, 5.5, 5.5, 20), "pt"),
        plot.title = element_text(family = 'Arial', size = 16, color = 'black',hjust = 0.5))


####
### circus plot + aano
#df <- read_table("eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC.txt")
table(duplicated(df))

ref <- read.table("~/Desktop/KU/@research/STR/db/annovar/output_bed_ref")
anno <- read_table("~/Desktop/KU/@research/STR/db/annovar/Raw.anno.merge.processing.onlyneed.txt")
#concordance_rate <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/figure/figure2/f2.concordance_bySTR_simpleSTR.txt")
concordance_rate <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/figure/figure2_withchrX/f2.concordance_bySTR_simpleSTR.afterchrXQC.txt")
head(concordance_rate)
concordance_rate %>% filter(STR_ID == "RFC1")
head(concordance_rate)

final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR

final_ref %>% select(MOTIFS,ID) %>% rename(STR_ID = ID) %>% 
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> final_ref_pro


head(ref)
ref %>% unique()-> ref
colnames(ref) <- c("chrom","start","end","STR_ID")
head(a)
head(concordance_rate)
concordance_rate %>% left_join(anno) %>% left_join(ref) -> a
head(a)

cyto <- read_table("~/Desktop/KU/@research/STR/db/annovar/hg38_cytoBand.txt",col_names = F)
head(cyto)
colnames(cyto) <- c("chrom","start","end","cytoband","anno")
cyto %>% mutate(cytoband = paste0(str_replace_all(chrom,"chr",""),cytoband)) -> cyto
head(cyto)
cbind(cyto$start, cyto$end)

head(anno)
head(ref)
head(final_ref)
head(cyto)



data<-a[1:100000,] %>% filter(concordance_rate != 1)
data<-a

#head(data,20)

# Circos �÷� ����
##### Good

library(circlize)


concordance_rate %>% left_join(anno) %>% left_join(ref) %>% left_join(final_ref_pro) -> b
b %>% select(concordance_rate,chrom,cytoband,main_cyto_type) %>%mutate(concordance_range = ifelse(concordance_rate != 1,0,1)) %>% #left_join(cyto)
  group_by(concordance_range,chrom,cytoband,main_cyto_type) %>%  #head()
  summarise(concordance_rate = mean(concordance_rate)) %>% left_join(cyto) %>% 
  ungroup()-> b_cyto_bychr
head(b_cyto_bychr)

#fill = c("grey", "#F8766D", "#619CFF"),
#### �̰ɷ� ��
circos.clear()
circos.initializeWithIdeogram(chromosome.index = c(paste0("chr", c(1:22)),"chrX"),species = "hg38")
bed_list = list(b_cyto_bychr %>% filter(concordance_rate == 1) %>% select(chrom,start,end,concordance_rate) %>% arrange(chrom, start) %>% as.data.frame(),
                b_cyto_bychr %>% filter(concordance_rate != 1) %>% select(chrom,start,end,concordance_rate) %>% arrange(chrom, start) %>% as.data.frame())
head(bed_list)
circos.genomicTrack(bed_list, 
                    panel.fun = function(region, concordance_rate, ...) {
                      i = getI(...)
                      circos.genomicLines(region, concordance_rate, col = i, ...)
                    })

circos.genomicDensity(b %>% select(chrom,start,end,concordance_rate) %>% filter(concordance_rate == 1),col = c("#619CFF"), track.height = 0.1)
circos.genomicDensity(b %>% select(chrom,start,end,concordance_rate) %>% filter(concordance_rate != 1),col = c("#F8766D"), track.height = 0.1)


b_cyto_bychr 
####
png("circos_plot.png", width = 2000, height = 2000, res = 300)

circos.clear()
circos.initializeWithIdeogram(chromosome.index = c(paste0("chr", c(1:22)),"chrX"),species = "hg38")

bed_list = list(b_cyto_bychr %>% select(chrom,start,end,concordance_rate) %>% arrange(chrom, start) %>% as.data.frame())
head(bed_list)
circos.genomicTrack(bed_list, 
                    panel.fun = function(region, concordance_rate, ...) {
                      i = getI(...)
                      circos.genomicLines(region, concordance_rate, col = i, ...)
                    })
circos.genomicDensity(b %>% select(chrom,start,end,concordance_rate) %>% filter(concordance_rate != 1),col = c("#F8766D"), track.height = 0.1)
circos.genomicDensity(b %>% select(chrom,start,end,concordance_rate) %>% filter(concordance_rate == 1),col = c("#619CFF"), track.height = 0.1)
dev.off()
# cytoband

# cytoband
# ����: Cent
# Ǫ����: stalk

###################
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% left_join(final_ref_pro) #%>% 
head(anno)
head(cyto)
#### asso in my com
concordance_rate %>% left_join(anno) %>% left_join(ref) %>% left_join(final_ref_pro) -> a
concordance_rate %>% filter(STR_ID == "RFC1")
str(a)
head(a)


a %>% mutate(arm = ifelse(str_detect(cytoband,"p"),"p","q")) %>% 
  mutate(type = ifelse(type == "upstream","promoter",type)) %>%
  mutate(type = ifelse(str_detect(type,"RNA"),"ncRNA",type)) -> a

colnames(a)
a %>% select(STR_ID,concordance_rate,type,main_cyto_type,distance_bycen,chrom,GC) %>% str()
a %>% mutate(arm = as.factor(arm)) -> b
head(a)
asso_result <- NULL
linear_model <- lm(concordance_rate ~ type, data = a)
summary(linear_model)$`Adjusted R-squared`
summary(linear_model)$coefficients %>% as.data.frame() -> asso_result
summary(linear_model)$coefficients[2,]
asso_result$type <- row.names(asso_result)
head(asso_result)
colnames(asso_result) <- c("Estimate","Std.Error","t.value","P","type")


linear_model <- lm(concordance_rate ~ main_cyto_type, data = a)
linear_model <- lm(concordance_rate ~ cyto_type, data = a)
linear_model <- lm(concordance_rate ~ distance_bycen, data = a)
linear_model <- lm(concordance_rate ~ chrom, data = a)
linear_model <- lm(concordance_rate ~ GC, data = a)
linear_model <- lm(concordance_rate ~ arm, data = b)
summary(linear_model)


linear_model <- lm(concordance_rate ~ type, data = a)
linear_model <- lm(concordance_rate ~ type + main_cyto_type + distance_bycen + chrom + GC, data = a)
linear_model <- lm(concordance_rate ~ type + main_cyto_type + distance_bycen + chrom + GC, data = a)
linear_model <- lm(concordance_rate ~ type + main_cyto_type + distance_bycen + chrom + GC, data = a)







library(randomForest)

# ���� ������Ʈ �� ����
random_forest_model <- randomForest(concordance_rate ~ type + main_cyto_type + distance_bycen + chrom + GC, data = a, importance = TRUE)

# ���� �߿䵵 ���
importance(random_forest_model)
#saveRDS(random_forest_model, file = "random_forest_model.rds")
# ���� �߿䵵 �ð�ȭ
varImpPlot(random_forest_model)

importance_df <- as.data.frame(importance(random_forest_model))
importance_df$Variables <- rownames(importance_df)
head(importance_df)
row.names(importance_df) <- NULL
head(importance_df)



loaded_model <- readRDS("/Users/ksmpooh/Desktop/KU/@research/STR/figure/figure2_withchrX/random_forest_model.rds")
predicted_values <- predict(loaded_model, newdata = a)

# �������� �������� ����Ͽ� R�� ���
# R�� = 1 - (SSE / SST)
SSE <- sum((a$concordance_rate - predicted_values)^2)  # SSE: ���� ������
SST <- sum((a$concordance_rate - mean(a$concordance_rate))^2)  # SST: �� ������
R2 <- 1 - (SSE / SST)

# R�� �� ���
print(paste("R�� ��: ", R2))

linear_model <- lm(concordance_rate ~ type + main_cyto_type + distance_bycen + chrom + GC, data = a)

# �� ���
summary(linear_model)

varImpPlot(loaded_model)

library(ggplot2)
library(ggplot)

importance_df <- as.data.frame(importance(loaded_model))
importance_df$Variables <- rownames(importance_df)
head(importance_df)
row.names(importance_df) <- NULL
str(importance_df)
head(importance_df)
head(importance_df)
colnames(importance_df)[1] <- c("IncMSE")
#importance_df$IncNodePurity
# ggplot�� ����� �߿䵵 �ð�ȭ

importance_df <- data.frame(
  incMSE = c(45.15906, 15.34429, 20.75179, 17.49869, 84.38882),
  incNodePurity = c(17.51795, 11.29153, 73.22416, 28.94505, 164.09260),
  variables = c("type", "main_cyto_type", "distance_bycen", "chrom", "GC")
)

importance_df

importance_df %>%
  ggplot(aes(x = variables, y = incMSE)) +
  geom_point() + 
  theme_step1()


importance_df %>% #head()
  ggplot(aes(x = factor(Variables), y = `%IncMSE`)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # ���η� ������
  theme_step1()



head(a)
cor_data <- data.frame(
  concordance_rate = a$concordance_rate,
  distance_bycen = a$distance_bycen,
  RU_length = a$RU.length,
  GC = a$GC
)
correlation_matrix <- cor(cor_data, use = "complete.obs")
print(correlation_matrix)


ggplot(a, aes(x=type, y=concordance_rate, fill=type)) + 
  geom_boxplot() +
  labs(title="Type�� concordance_rate ��", x="Type", y="Concordance Rate") +
  theme_minimal()

# K-means clustering ���� (������ ����)
kmeans_result <- kmeans(scale(cor_data), centers = 3)
a$cluster <- as.factor(kmeans_result$cluster)

# Ŭ�����ͺ� �ð�ȭ
ggplot(a, aes(x=RU.length, y=concordance_rate, color=cluster)) + 
  geom_point() +
  labs(title="K-means Clustering", x="RU Length", y="Concordance Rate") +
  theme_minimal()

colnames(a)
library(rpart)
library(rpart.plot)
decision_tree <- rpart(concordance_rate ~ type + cytoband + cyto_type+GC+distance_bycen+main_cyto_type+chrom, data = a)
rpart.plot(decision_tree)
head(decision_tree)

#### setver asso
setwd("/BDATA/smkim/STR/asso")
library(tidyverse)
df <- read_table("eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC.txt")
cyto <- read_table("hg38_cytoBand.txt",col_names = F)

ref <- read.table("output_bed_ref")
anno <- read_table("Raw.anno.merge.processing.onlyneed.txt")
concordance_rate <- read_table("f2.concordance_bySTR_simpleSTR.afterchrXQC.txt")

final_ref <- read.table("Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR

final_ref %>% select(MOTIFS,ID) %>% rename(STR_ID = ID) %>% 
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> final_ref_pro


head(ref)
ref %>% unique()-> ref
colnames(ref) <- c("chrom","start","end","STR_ID")

head(anno)
dim(anno)
anno %>% unique() %>% dim()
anno %>% filter(STR_ID %in% common_STR$STR_DB) %>%
  filter(!str_detect(STR_ID,"chrX")) %>% dim()

head(df)
head(ref)
head(anno)
head(final_ref_pro)
final_ref_pro %>% left_join(ref) %>% left_join(anno)
df %>% left_join(anno) %>% left_join(final_ref_pro) %>% left_join(ref) -> merge_all
concordance_rate %>% left_join(anno) %>% left_join(final_ref_pro) %>% left_join(ref) -> merge_all_meanSTR
merge_all %>% mutate(STR_length = TRGT)


merge_all_meanSTR %>%
  mutate(arm = ifelse(str_detect(cytoband,"p"),"p","q")) %>% 
  mutate(type = ifelse(type == "upstream","promoter",type)) %>%
  mutate(type = ifelse(str_detect(type,"RNA"),"ncRNA",type)) -> merge_all_meanSTR



df %>% left_join(anno) %>% left_join(final_ref_pro) %>% left_join(ref) -> merge_all
head(merge_all)
dim(df) # 40257880
df %>% unique() %>% dim() #40257880
head(merge_all)
merge_all %>%  mutate(STR_diff_logi = ifelse(TRGT_STR == EH_STR,1,0)) %>%
  mutate(STR_length = TRGT_STR * RU.length) %>% 
  mutate(arm = ifelse(str_detect(cytoband,"p"),"p","q")) %>% 
  mutate(type = ifelse(type == "upstream","promoter",type)) %>%
  mutate(type = ifelse(str_detect(type,"RNA"),"ncRNA",type)) -> merge_all_1


str_anno <- read_table("~/Desktop/KU/@research/STR/figure/sup.figure/STR.concordance.INFO.withanno.simpleSTR.txt")

glm(STR_diff_logi ~ GC,data=merge_all_1,family = binomial) -> concordance_model

lm(concordance_rate~type+main_cyto_type+cyto_type+distance_bycen+GC+chrom,data = str_anno) -> a
glm(STR_diff_logi~type+main_cyto_type+cyto_type+distance_bycen+GC+chrom+arm,data = merge_all_1,family = binomial) -> concordance_model
lm(concordance_rate~arm,data = merge_all_1) -> concordance_model_arm_lm
head(merge_all_meanSTR)
lm(concordance_rate~type+main_cyto_type+cyto_type+distance_bycen+GC+chrom+arm,data = merge_all_meanSTR) -> a
lm(concordance_rate~arm,data = merge_all_meanSTR) -> a
lm(concordance_rate~chrom,data = merge_all_meanSTR) -> a
glm(concordance_range~chrom,data = merge_all_meanSTR) -> a

table(merge_all_meanSTR$concordance_range)

lm(concordance_rate~arm,data = merge_all_1) -> concordance_model_arm_lm_rate

head(merge_all_1)
summary(concordance_model_arm_lm)

summary(concordance_model)$coefficients %>% as.data.frame() -> asso_result
summary(concordance_model)$coefficients[2,]
asso_result$type <- row.names(asso_result)
head(asso_result)
colnames(asso_result) <- c("Estimate","Std.Error","t.value","P","type")
asso_result$type
asso_result$`Pr(>|t|)`
write.table(asso_result,"Asso.annovar.annotation.txt",col.names = T,row.names = F,quote = F,sep = "\t")



#######3 STR length diff
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/figure2_withchrX/eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC.txt")
simple_concordance <- read.table("~/Desktop/KU/@research/STR/figure/figure2_withchrX/f2.concordance_bySTR_simpleSTR.afterchrXQC.txt",header = T)

head(final_ref_pro)
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% 
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% group_by(STR_ID) %>%
  mean



eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% #
  group_by(STR_ID) %>%
  summarise(mean_TRGT_STR = mean(TRGT_STR),max_TRGT_STR = max(TRGT_STR),min_TRGT_STR = min(TRGT_STR),sd_TRGT_STR =sd(TRGT_STR),
            mean_EH_STR = mean(EH_STR),max_EH_STR = max(EH_STR),min_EH_STR = min(EH_STR),sd_EH_STR =sd(EH_STR)) -> eh_trgt_merge_simple_pass_intersect_STRcountINFO


eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% mutate(STR_count_diff = TRGT_STR-EH_STR) %>%
  group_by(STR_ID) %>%
  summarise(mean_TRGT_STR = mean(TRGT_STR),mean_EH_STR = mean(EH_STR),sd_STR_count_diff =sd(STR_count_diff)) -> eh_trgt_merge_simple_pass_intersect_STRcount_withdiffSD

          
eh_trgt_merge_simple_pass_intersect_STRcountINFO
eh_trgt_merge_simple_pass_intersect_STRlengthINFO

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% left_join(final_ref_pro %>% select(STR_ID,RU.length)) %>%
  mutate(TRGT_STR = TRGT_STR*RU.length,EH_STR=EH_STR*RU.length) %>%
  group_by(STR_ID) %>%
  summarise(mean_TRGT_STR_legnth = mean(TRGT_STR),max_TRGT_STR_legnth = max(TRGT_STR),min_TRGT_STR_legnth = min(TRGT_STR),sd_TRGT_STR_legnth =sd(TRGT_STR),
            mean_EH_STR_legnth = mean(EH_STR),max_EH_STR_legnth = max(EH_STR),min_EH_STR_legnth = min(EH_STR),sd_EH_STR_legnth =sd(EH_STR)) -> eh_trgt_merge_simple_pass_intersect_STRlengthINFO

head(eh_trgt_merge_simple_pass_intersect_STRcountINFO)
head(eh_trgt_merge_simple_pass_intersect_STRlengthINFO)
head(simple_concordance)

eh_trgt_merge_simple_pass_intersect_STRcountINFO
eh_trgt_merge_simple_pass_intersect_STRlengthINFO
#write.table(eh_trgt_merge_simple_pass_intersect_STRcountINFO,"~/Desktop/KU/@research/STR/figure/figure2_withchrX/eh_trgt_merge_simple_pass_intersect_meanSTRcountINFO_afterQCchrX.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#write.table(eh_trgt_merge_simple_pass_intersect_STRlengthINFO,"~/Desktop/KU/@research/STR/figure/figure2_withchrX/eh_trgt_merge_simple_pass_intersect_meanSTRlengthINFO_afterQCchrX.txt",col.names = T,row.names = F,quote = F,sep = "\t")


eh_trgt_merge_simple_pass_intersect_STRlengthINFO %>% #head(1000) %>% #head()
  pivot_longer(mean_TRGT_STR_legnth:sd_EH_STR_legnth) %>% mutate(value = log(value)) %>%
  pivot_wider(names_from = name,values_from = value) %>%
  ggplot(aes(x = mean_TRGT_STR_legnth, y = mean_EH_STR_legnth)) +
  geom_point() +  # ?????? ??????????????
  #geom_hex(bins=15) + 
  #geom_errorbar(aes(ymin = min_EH_STR_legnth, ymax = max_EH_STR_legnth), width = 0,alpha = 0.5) +  # y ???????????? ??????????????????
  #geom_errorbarh(aes(xmin = min_TRGT_STR_legnth, xmax = max_TRGT_STR_legnth), height = 0,alpha = 0.5) +  # x ???????????? ??????????????????
  labs(x = "Mean of STR length by LRS", y = "Mean of STR length by SRS") +
#  theme_bw() + 
  coord_fixed(ratio = 1) 

max(final_ref_pro$RU.length)
#final_ref_pro %>% filter(RU.length > 20)


eh_trgt_merge_simple_pass_intersect_STRcountINFO %>% left_join(final_ref_pro) %>% 
  mutate(RU.length = case_when(
    RU.length < 6 ~"[2,5]",
    RU.length < 11 ~"[6,10]",
    TRUE ~ "[11~24]"
  )) %>%
  ggplot(aes(x = mean_TRGT_STR, y = mean_EH_STR,color=RU.length)) +
  geom_point(alpha=0.3) +  # ?????? ??????????????
  #geom_hex(bins=15) + 
  #geom_errorbar(aes(ymin = min_EH_STR, ymax = max_EH_STR), width = 0,alpha = 0.5) +  # y ???????????? ??????????????????
  #geom_errorbarh(aes(xmin = min_TRGT_STR, xmax = max_TRGT_STR), height = 0,alpha = 0.5) +  # x ???????????? ??????????????????
  labs(x = "Mean of STR Repeat Count by LRS", y = "Mean of STR Repeat Count by SRS") +
  #  theme_bw() + 
  coord_fixed(ratio = 1) 

eh_trgt_merge_simple_pass_intersect_STRcount_withdiffSD

eh_trgt_merge_simple_pass_intersect_STRcount_withdiffSD %>%
  pivot_longer(mean_TRGT_STR:sd_STR_count_diff) %>% mutate(value = log(value)) %>%
  pivot_wider(names_from = name,values_from = value) %>%
  ggplot(aes(x = mean_TRGT_STR, y = mean_EH_STR,color=sd_STR_count_diff)) +
  geom_point() + 
  coord_fixed(ratio = 1) 

  