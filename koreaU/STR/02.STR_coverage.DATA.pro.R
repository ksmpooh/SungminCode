
## coverage
library(tidyverse)
library(data.table)

ru <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.pass.IDregion.txt")
ref <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/db/str-analysis/str_analysis/variant_catalogs/catalog.GRCh38.with_adjacent_repeats.TRGT.bed",col_names = F)
final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
final_ref %>% select(MOTIFS,ID) %>% rename(STR_ID = ID) %>% 
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> final_ref_pro

head(ref)
head(final_ref)
dim(final_ref)
final_ref %>% count(ID  == STR_region)
final_ref %>% rename(STR_ID = ID)
head(venn_rawdata)

venn_rawdata %>% na.omit() %>% select(STR_DB) %>% rename("STR_ID" = STR_DB) %>% left_join(final_ref %>% rename(STR_ID = ID) %>% select(STR_ID,STR_region)) -> ref_strgion


setwd("/Users/ksmpooh/Desktop/KU/@research/STR/eh/coverage")
flist = grep(list.files("./"),pattern = "txt", value=TRUE)
flist
df <- NULL
for (i in flist) {
  tmp <- read_table(i,col_names = FALSE)
  tmp$ID <- str_replace(i,".EH.coverage_LC.txt","")
  df <- rbind(df,tmp)
}
#df <- rbind(df,tmp)
head(df)
colnames(df) <- c("STR_ID","MOTIFS","LC","ID")

df -> df_guo
head(df_guo)

setwd("/Users/ksmpooh/Desktop/KU/@research/STR/eh/coverage_pathgenicSTRanalysisDB/")
flist = grep(list.files("./"),pattern = "txt", value=TRUE)
flist
df <- NULL
for (i in flist) {
  tmp <- read_table(i,col_names = FALSE)
  tmp$ID <- str_replace(i,".EH.pathgenicSTRanalysisDB.coverage_LC.txt","")
  df <- rbind(df,tmp)
}
#df <- rbind(df,tmp)
head(df)
colnames(df) <- c("STR_ID","MOTIFS","LC","ID")

colnames(df_patho) <- c("STR_ID","MOTIFS","LC","ID")

df -> df_patho
head(df_patho)

##############
setwd("/Users/ksmpooh/Desktop/KU/@research/STR/trgt/coverage_guo/")
flist = grep(list.files("./"),pattern = "merge", value=TRUE)
flist
df <- NULL
for (i in flist) {
  tmp <- read_table(i,col_names = FALSE)
  tmp$ID <- str_replace(str_replace(i,".pbmm2_hg38_withunmapped_trgt_genotype_gangstr.sorted.bam.coverage.merge",""),"_sorted_trgt_genotype_gangstr.sorted.bam.coverage.merge","")
  df <- rbind(df,tmp)
}
colnames(df) = c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq","ID")
df %>% select(ID) %>% unique() %>% tail()
df %>% mutate(ID = ifelse(str_detect(ID,".pbmm2"),str_replace(ID,".pbmm2_hg38_withunmapped_trgt_genotype_gangstr.spanning.sorted.bam.coverage.merge",""),
                          str_replace(ID,"_sorted_trgt_genotype_gangstr.spanning.sorted.bam.coverage.merge",""))) -> trgt_guo


head(ref_strgion)
head(trgt_guo)
trgt_guo %>% mutate(STR_region = paste0(rname,"_",startpos,"_",endpos)) %>% select(STR_region) %>% unique() %>% count(STR_region %in% ref_strgion$STR_region)

ref_strgion %>% filter(str_detect(STR_ID,"chr")) -> ref_strgion_guo
trgt_guo %>% mutate(STR_region = paste0(rname,"_",startpos,"_",endpos)) %>% filter(STR_region %in% ref_strgion_guo$STR_region) %>% 
  left_join(ref_strgion_guo) %>% select(-rname,-startpos,-endpos,-STR_region) -> trgt_guo

trgt_guo %>% mutate(ID = str_replace_all(ID,"\\.merge","")) -> trgt_guo

setwd("/Users/ksmpooh/Desktop/KU/@research/STR/trgt/coverage_patho/")
flist = grep(list.files("./"),pattern = "merge", value=TRUE)
flist
df <- NULL
for (i in flist) {
  tmp <- read_table(i,col_names = FALSE)
  tmp$ID <- str_replace(str_replace(i,".pbmm2_hg38_withunmapped_trgt_genotype_stranalysisDB.spanning.sorted.bam.coverage.merge",""),"_sorted_trgt_genotype_stranalysisDB.spanning.sorted.bam.coverage.merge","")
  tmp$ID <- str_replace(tmp$ID,"\\.merge","")
  df <- rbind(df,tmp)
}
colnames(df) = c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq","ID")
#df %>% select(ID) %>% filter(str_detect(ID,"_"))
head(df)
trgt_path <- df

ref_strgion %>% filter(!str_detect(STR_ID,"chr")) -> ref_strgion_path

trgt_path %>% mutate(STR_region = paste0(rname,"_",startpos,"_",endpos)) %>% #filter(STR_region %in% ref_strgion_path$STR_region) %>% 
  left_join(ref_strgion_path) %>% select(-rname,-startpos,-endpos,-STR_region) -> trgt_path

head(trgt_path)
trgt_path %>% select(STR_ID) %>% unique() 
##############
head(df_guo)
head(df_patho)
head(trgt_guo)
head(trgt_path)


ref <- read.table("~/Desktop/KCDC/pangenome/00.datacheck/KBA.Long_Revio_Nanopore_short.IDmatchinagtable.txt",header = T)

ref %>% select(Revio,Illumina) -> ref
head(ref)
colnames(ref) <- c("ID","EH_ID")


filter(ID != "NIH23J3904558") 
eh_pass_common_gt %>% left_join(ref) #%>% #filter(ID %in% trgt_pass_common_gt$ID) %>% count(ID)
  
df_guo %>% rename(EH_ID = ID) %>%
  left_join(ref) %>% filter(ID != "NIH23J3904558")  %>% select(-EH_ID) -> df_guo

df_patho %>% rename(EH_ID = ID) %>%
  left_join(ref) %>% filter(ID != "NIH23J3904558")  %>% select(-EH_ID) -> df_patho

df_guo %>% filter(STR_ID %in% ref_strgion_guo$STR_ID) -> df_guo
#df_patho %>% filter(STR_ID %in% ref_strgion_path$STR_ID) -> df_patho

head(df_patho)
head(df_guo)
eh_complex <- read.table("/Users/ksmpooh/Desktop/KU/@research/STR/figure/extra_info/complex_STR.v2.txt",header = T)
head(eh_complex)
#eh_complex %>% select(STR_ID) %>% unique() -> eh_complex

head(eh_complex)
head(df_guo)
df_guo %>% filter(!str_detect(STR_ID,"chr"))
df_guo %>% rbind(df_patho) %>% filter(STR_ID %in% eh_complex$STR_ID) -> eh_complex_LC
df_guo %>% rbind(df_patho) %>% filter(!(STR_ID %in% eh_complex$STR_ID)) -> eh_simple_LC

head(eh_complex_LC);dim(eh_complex_LC)
eh_complex_LC %>% select(STR_ID) %>% unique() %>% dim()

write.table(eh_complex_LC,"~/Desktop/KU/@research/STR/eh/eh_complexSTR_LC.txt",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(eh_simple_LC,"~/Desktop/KU/@research/STR/eh/eh_simpleSTR_LC.txt",col.names = T,row.names = F,quote = F,sep = "\t")

########
head(trgt_guo);dim(trgt_guo)/66
head(trgt_path);dim(trgt_path)/66
head(eh_complex)
trgt_guo %>% rbind(trgt_path) %>% filter(ID != "NIH23J3904558") %>% filter(!(STR_ID %in% eh_complex$STR_ID)) -> trgt_simple_coverage
trgt_path %>% filter(ID != "NIH23J3904558") %>% filter((STR_ID %in% eh_complex$STR_ID)) -> trgt_complex_coverage

trgt_simple_coverage %>% select(STR_ID) %>% unique()
trgt_complex_coverage %>% select(STR_ID) %>% unique() 
write.table(trgt_complex_coverage,"~/Desktop/KU/@research/STR/trgt/trgt_complexSTR_coverage.txt",col.names = T,row.names = F,quote = F,sep = "\t")

write.table(trgt_simple_coverage,"~/Desktop/KU/@research/STR/trgt/trgt_simpleSTR_coverage.txt",col.names = T,row.names = F,quote = F,sep = "\t")


#### short-read samtools coverage

setwd("/Users/ksmpooh/Desktop/KU/@research/STR/eh/samtools_coverage/coverage/")
flist = grep(list.files("./"),pattern = "merge", value=TRUE)
flist
df <- NULL
for (i in flist) {
  tmp <- read_table(i,col_names = FALSE)
  tmp$ID <- str_replace(i,".bam.coverage.merge","")
  df <- rbind(df,tmp)
}
colnames(df) = c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq","ID")
srs <- df
head(srs)

setwd("/Users/ksmpooh/Desktop/KU/@research/STR/eh/samtools_coverage/coverage_path/")
flist = grep(list.files("./"),pattern = "bam", value=TRUE)
flist
df <- NULL
#tmp <- read_table("NIH20N2215198.bam.ABCD3")
head(tmp)
for (i in flist) {
  tmp <- read_table(i)
  tmp %>% mutate(ID = str_split_fixed(i,"\\.",3)[,1]) %>% mutate(STR_ID = str_split_fixed(i,"\\.",3)[,3]) -> tmp
  df <- rbind(df,tmp)
}
colnames(df) = c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq","ID","STR_ID")
#df %>% select(ID) %>% filter(str_detect(ID,"_"))
head(df)
srs_path <- df
head(srs_path)
head(srs)
head(final_ref)
srs %>% mutate(STR_region = paste0(rname,"_",startpos,"_",endpos)) %>% 
  left_join(final_ref %>% select(-chrom,-start,-end,-MOTIFS) %>% rename(STR_ID = ID)) %>% 
  filter(str_detect(STR_ID,"chr")) -> srs

head(srs)
head(srs_path)
trgt_path %>% mutate(STR_region = paste0(rname,"_",startpos,"_",endpos)) %>% #filter(STR_region %in% ref_strgion_path$STR_region) %>% 
  left_join(ref_strgion_path) %>% select(-rname,-startpos,-endpos,-STR_region) -> trgt_path


ref <- read.table("~/Desktop/KCDC/pangenome/00.datacheck/KBA.Long_Revio_Nanopore_short.IDmatchinagtable.txt",header = T)

ref %>% select(Revio,Illumina) -> ref
head(ref)
colnames(ref) <- c("ID","EH_ID")


filter(ID != "NIH23J3904558") 
eh_pass_common_gt %>% left_join(ref) #%>% #filter(ID %in% trgt_pass_common_gt$ID) %>% count(ID)

srs %>% rename(EH_ID = ID) %>%
  left_join(ref) %>% filter(ID != "NIH23J3904558")  %>% select(-EH_ID) -> srs

srs_path %>% rename(EH_ID = ID) %>%
  left_join(ref) %>% filter(ID != "NIH23J3904558")  %>% select(-EH_ID) -> srs_path


head(srs)
head(srs_path)
srs %>% select(-STR_region) %>% rbind(srs_path) -> srs_merge

head(srs_merge)
final_ref %>% filter(str_detect(MOTIFS,",")) -> complex_str
complex_str
srs_merge %>% filter(!str_detect(STR_ID,"chr")) #%>% write.table("~/Desktop/KU/@research/STR/eh/eh_complexSTR_comp.txt")

srs_merge %>% filter(!str_detect(STR_ID,"chr")) %>% filter(STR_ID %in% complex_str$ID) %>% write.table("~/Desktop/KU/@research/STR/eh/eh_complexSTR_oribam_coverage.txt",col.names = T,row.names = F,quote = F,sep = "\t")
srs_merge %>% filter(!(STR_ID %in% complex_str$ID)) %>% select(-rname,-startpos,-endpos) %>%
  write.table("~/Desktop/KU/@research/STR/eh/eh_simpleSTR_oribam_coverage.txt",col.names = T,row.names = F,quote = F,sep = "\t")


#####
eh_simple_coverage <- read_table("~/Desktop/KU/@research/STR/eh/eh_simpleSTR_oribam_coverage.txt") #%>% filter(STR_ID %in% common_STR$STR_DB)
trgt_simple_coverage <- read_table("~/Desktop/KU/@research/STR/trgt/trgt_simpleSTR_coverage.txt") #%>% select(ID,STR_ID,meandepth)

head(eh_simple_coverage)
head(trgt_simple_coverage)

##### pandepth GC 
setwd("/Users/ksmpooh/Desktop/KU/@research/STR/trgt/pandepth/normal")
setwd("/BDATA/smkim/STR/trgt/pandepth/normal")

header = c("Chr","Length","CoveredSite","TotalDepth","Coverage(%)","MeanDepth")
## sh
flist = grep(list.files("./"),pattern = "stat.gz", value=TRUE)
flist
head(flist)
df <- NULL
for (i in flist) {
  tmp <- read_table(i)
  colnames(tmp)<-header
  tmp$filename <- i
  tmp$platform <- "trgt"
  df <- rbind(df,tmp)
}
head(df)
colnames(df) <- c("Chr","Start","End","GeneID","Length","CoveredSite","TotalDepth","GC(%)","Coverage(%)","MeanDepth","filename","platfrom")

trgt_a <- df

setwd("/Users/ksmpooh/Desktop/KU/@research/STR/trgt/pandepth/patho/")
setwd("/BDATA/smkim/STR/trgt/pandepth/patho/")
flist = grep(list.files("./"),pattern = "stat.gz", value=TRUE)
flist
head(flist)
df <- NULL
for (i in flist) {
  tmp <- read_table(i)
#  colnames(tmp)<-header
  tmp$filename <- i
  tmp$platform <- "trgt"
  df <- rbind(df,tmp)
}
head(df)
colnames(df) <- c("Chr","Start","End","GeneID","Length","CoveredSite","TotalDepth","GC(%)","Coverage(%)","MeanDepth","filename","platfrom")

trgt_b <- df

head(trgt_a)
head(trgt_b)

trgt_a %>% filter(str_detect(GeneID,"chr")) %>% rbind(trgt_b) %>% select(-Length,-CoveredSite,-TotalDepth) %>%
  rename(STR_ID = GeneID) %>% 
  mutate(type = ifelse(str_detect(filename,"normal"),"normal","pathogenic")) %>% #filter(!(type == "normal" & !str_detect(STR_ID,"chr"))) %>%
  mutate(ID = ifelse(str_detect(filename,"sorted"),str_split_fixed(filename,"_",2)[,1],str_split_fixed(filename,"\\.",2)[,1])) %>%
  filter(!str_detect(Chr,"#")) %>% select(-Chr,-Start,-End,-filename) -> trgt

head(trgt)

trgt %>% select(-filename) %>% write.table("~/Desktop/KU/@research/STR/trgt/pandepth/pandepth.merge.from.trgt.oribam.txt",col.names = T,row.names = F,quote = F,sep = "\t")


setwd("/Users/ksmpooh/Desktop/KU/@research/STR/eh/pandepth")
flist = grep(list.files("./"),pattern = "stat.gz", value=TRUE)
flist
head(flist)
df <- NULL
for (i in flist) {
  tmp <- read_table(i)
  #  colnames(tmp)<-header
  #tmp$filename <- i
  tmp$platform <- "eh"
  tmp$ID <- str_split_fixed(i,"\\.",2)[,1]
  tmp$type <- ifelse(str_detect(i,"normal"),"normal","pathogenic")
  df <- rbind(df,tmp)
}
head(df)
colnames(df) <- c("Chr","Start","End","STR_ID","Length","CoveredSite","TotalDepth","GC(%)","Coverage(%)","MeanDepth","platfrom","ID","type")

df %>% select(-Length,-CoveredSite,-TotalDepth) %>% filter(!(type == "normal" & !str_detect(STR_ID,"chr"))) %>%
  write.table("~/Desktop/KU/@research/STR/eh/pandepth.merge.from.eh.oribam.txt",col.names = T,row.names = F,quote = F,sep = "\t")


eh_pandepth <- read_table("~/Desktop/KU/@research/STR/eh/pandepth.merge.from.eh.oribam.txt")
trgt_pandepth<- read_table("~/Desktop/KU/@research/STR/trgt/pandepth/pandepth.merge.from.trgt.oribam.txt")

tail(trgt_pandepth)

trgt_pandepth %>% filter(!str_detect(Chr,"#")) %>% select(-Chr,-Start,-End) %>% mutate(ID = str_replace_all(ID,"pandepth","")) %>% 
  write.table("~/Desktop/KU/@research/STR/trgt/pandepth/pandepth.merge.from.trgt.oribam.txt",col.names = T,row.names = F,quote = F,sep = "\t")

eh_pandepth  %>% filter(!str_detect(Chr,"#")) %>% select(-Chr,-Start,-End) %>% 
  write.table("~/Desktop/KU/@research/STR/eh/pandepth.merge.from.eh.oribam.txt",col.names = T,row.names = F,quote = F,sep = "\t")


#### ori bam sammotols coverage
setwd("/BDATA/smkim/STR/pandepth/coverage_oribam")
flist = grep(list.files("./"),pattern = "coverage", value=TRUE)
flist
head(flist)
df <- NULL
header <- c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")
for (i in flist) {
  tmp <- read_table(i,col_names = F)
  colnames(tmp)<-header
  tmp$filename <- i
  tmp$platform <- "trgt"
  df <- rbind(df,tmp)
}
head(df)

df %>% mutate(STR_region = paste0(rname,"_",startpos,"_",endpos)) %>%
  mutate(type = ifelse(str_detect(filename,"normal"),"normal","pathogenic")) %>% 
  mutate(ID = ifelse(str_detect(filename,"sorted"),str_split_fixed(filename,"_",2)[,1],str_split_fixed(filename,"\\.",2)[,1])) -> trgt


trgt %>% select(-filename,-rname,-startpos,-endpos) %>% write.table("trgt_normal_patho_oribam_coverage.txt",col.names = T,row.names = F,quote = F,sep = "\t")



## mac
trgt <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/trgt/oribam_cover/trgt_normal_patho_oribam_coverage.txt")

head(ref_strgion)
ref_strgion_guo %>% unique() %>% dim()
ref_strgion %>% filter(str_detect(STR_ID,"chr")) -> ref_strgion_guo
ref_strgion %>% filter(!str_detect(STR_ID,"chr")) -> ref_strgion_path

head(trgt)
head(ref_strgion_guo)

table(trgt$type)
trgt %>% filter(type == "pathogenic" & STR_region %in% ref_strgion_path$STR_region) %>% left_join(ref_strgion_path) %>% 
  rbind(trgt %>% filter(type == "normal" & STR_region %in% ref_strgion_guo$STR_region) %>% left_join(ref_strgion_guo)) -> trgt

head(trgt)

#write.table(trgt,"/Users/ksmpooh/Desktop/KU/@research/STR/trgt/oribam_cover/trgt_normal_patho_oribam_coverage.txt",col.names = T,row.names = F,quote = F,sep = "\t")

##### sample covergae
setwd("~/Desktop/KCDC/pangenome/bam.stats/pandepth/")

#flist <- list()

flist = grep(list.files("./"),pattern = ".bam.pandepth.chr.stat.gz", value=TRUE)
flist[!str_detect(flist,'pbmm2') & !str_detect(flist,'sort')] -> flist
head(flist)
df <- NULL
for (i in flist) {
  tmp <- read_table(i)
  #colnames(tmp)<-header
  #tmp$file <- i
  tmp$ID <- str_replace(i,".bam.pandepth.chr.stat.gz","")
  df <- rbind(df,tmp)
}
head(df) #NIH23J3904558 'NIH20N2594890'

df %>%na.omit() %>% filter(!str_detect(`#Chr`,"RegionLength")) %>% filter(ID != "NIH20N2594890") %>% 
  filter(`#Chr` %in% c(paste0('chr',1:22),'chrX')) %>%
  group_by(ID) %>%
  summarise(check = sum(TotalDepth)/sum(Length))  %>% mutate(platfrom = "SRS") -> srs


'NIH20N2594890'


setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/bam.stats/pandepth/")

#flist <- list()

flist = grep(list.files("./"),pattern = ".bam.pandepth.chr.stat.gz", value=TRUE)
flist[str_detect(flist,'NIH23F')] -> flist
flist
head(flist)
df <- NULL
for (i in flist) {
  tmp <- read_table(i)
  #colnames(tmp)<-header
  tmp$file <- i
  #tmp$ID <- i
  df <- rbind(df,tmp)
}
head(df) #NIH23J3904558 'NIH20N2594890'

df %>%na.omit() %>% filter(!str_detect(`#Chr`,"RegionLength")) %>% #filter(ID != "NIH20N2594890") %>% 
  filter(`#Chr` %in% c(paste0('chr',1:22),'chrX')) %>%
  group_by(file) %>%
  summarise(check = sum(TotalDepth)/sum(Length))  %>% mutate(platfrom = "LRS") -> df1


head(df1)


setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/bam.stats/2023.pro.ref_panel.Revio/")
header <- c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")

#flist <- list()

flist = grep(list.files("./"),pattern = ".bam.coverage", value=TRUE)
length(flist)
head(flist)
df <- NULL
for (i in flist) {
  tmp <- read_table(i)
  colnames(tmp)<-header
  tmp$file <- i
  #tmp$ID <- str_replace(i,".bam.coverage","")
  df <- rbind(df,tmp)
}
df2 <- df
head(df2)
df2 %>% filter(rname %in% c(paste0('chr',1:22),'chrX')) %>%
  group_by(file) %>% filter(!str_detect(file,"NIH23J3904558")) %>% #head()
  summarise(check = sum(meandepth * endpos)/sum(endpos)) %>% mutate(platfrom = "LRS") %>% rbind(df1) -> lrs


head(lrs)
head(srs)
lrs %>% rename(ID = file)%>% rbind(srs) %>% #write.table("~/Desktop/KU/@research/STR/figure/sup.figure/mapping.depth.txt",col.names = T,row.names = F,quote = F,sep = "\t")
  ggplot(aes(x=platfrom,y=check,fill=platfrom))+ 
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black",alpha=0.8) + 
  labs(y="Mean of mapping depth (X)") +
  theme_bw() +
  theme(legend.position = 'none')

  
  

### patho
setwd("/Users/ksmpooh/Desktop/KU/@research/STR/eh/samtools_coverage/coverage_path")
header <- c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")

#flist <- list()

flist = grep(list.files("./"),pattern = "bam", value=TRUE,)
length(flist)
head(flist)
df <- NULL
for (i in flist) {
  tmp <- read_table(i)
  colnames(tmp)<-header
  tmp$ID <- str_split_fixed(i,"\\.",3)[1]
  tmp$STR_ID <- str_split_fixed(i,"\\.",3)[3]
  #tmp$
  #tmp$ID <- str_replace(i,".bam.coverage","")
  df <- rbind(df,tmp)
}
eh <- df

trgt <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/trgt/oribam_cover/trgt_normal_patho_oribam_coverage.txt")
head(trgt)

eh$platform <- "SRS"
trgt$platform <- "LRS"
head(eh)
head(trgt)
head(ref)
dim(eh)
trgt %>% filter(type == "pathogenic") %>% select(-STR_region) -> trgt


head(eh)
head(trgt)

df_guo %>% rename(EH_ID = ID) %>%
  left_join(ref) %>% filter(ID != "NIH23J3904558")  %>% select(-EH_ID) -> df_guo

head(eh)
eh %>%  filter(STR_ID %in% trgt$STR_ID) %>% select(meandepth:platform) %>% 
  rename(EH_ID = ID) %>%
  left_join(ref) %>% select(-EH_ID) %>% #head()
  rbind(trgt %>% select(-type)) %>%  filter(ID != "NIH23J3904558")  %>%
  write.table("~/Desktop/KU/@research/STR/figure/qc/trgt_eh_patho_oribam.txt",col.names = T,row.names = F,quote = F,sep = '\t')


##### bam mapQ bam
head(trgt_eh_patho_oribam)
trgt_eh_patho_oribam<- read_table("~/Desktop/KU/@research/STR/figure/qc/trgt_eh_patho_oribam.txt")
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/figure2_withchrX/eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC.txt") %>%
  select(-TRGT_AM)
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
head(final_ref_pro)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% left_join(final_ref_pro %>% select(STR_ID,RU.length)) %>%
  mutate(TRGT_STR_length = TRGT_STR*RU.length,EH_STR_length = EH_STR*RU.length) %>%
  filter(TRGT_STR_length >= 100 & EH_STR_length >=100 ) -> eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_strlength100


eh_simple_coverage <- read_table("~/Desktop/KU/@research/STR/eh/eh_simpleSTR_oribam_coverage.txt") %>% select(ID,STR_ID,meandepth,meanbaseq,meanmapq)
head(eh_simple_coverage)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_strlength100 %>% filter(str_detect(STR_ID,"chr")) %>%
  left_join(eh_simple_coverage) -> simpleSTR.strlength100
rm(eh_simple_coverage)
head(simpleSTR.strlength100)
trgt_coverage <- read_table("~/Desktop/KU/@research/STR/trgt/oribam_cover/trgt_normal_patho_oribam_coverage.txt") %>% select(ID,STR_ID,meandepth,meanbaseq,meanmapq)

simpleSTR.strlength100 %>% rename(eh_meandepth = meandepth,eh_meanbaseq = meanbaseq,eh_meanmapq = meanmapq) %>%
  left_join(trgt_coverage %>% rename(trgt_meandepth = meandepth,trgt_meanbaseq = meanbaseq,trgt_meanmapq = meanmapq)) -> simpleSTR.strlength100

head(simpleSTR.strlength100)  
write.table(simpleSTR.strlength100,"~/Desktop/KU/@research/STR/figure/figure4/normalSTR.STRlength100.with.bamstats.txt",col.names = T,row.names = F,quote = F,sep = "\t")


trgt_eh_patho_oribam<- read_table("~/Desktop/KU/@research/STR/figure/qc/trgt_eh_patho_oribam.txt")
head(trgt_eh_patho_oribam)

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/figure2_withchrX/eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC.txt")
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% select(-TRGT_AM) -> eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC

trgt_cov <- read_table("~/Desktop/KU/@research/STR/trgt/oribam_cover/trgt_normal_patho_oribam_coverage.txt") %>% select(-STR_region,-type,-platform)
  
head(trgt_cov)

trgt_cov %>% select(STR_ID) %>% filter(!str_detect(STR_ID,"chr")) %>% unique()
