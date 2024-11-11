#### 2024 KOGO
library(tidyverse)
library(readr)
library(cowplot)
### stat 
## bam stats
setwd("~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/")
ref <- read.table("~/Desktop/KCDC/pangenome/00.datacheck/KBA.Long_Revio_Nanopore_short.IDmatchinagtable.txt",header = T) %>% na.omit()
header = c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")
## sh
setwd("~/Desktop/KCDC/pangenome/bam.stats/wgs/")
flist <- list()

flist = grep(list.files("./"),pattern = "coverage", value=TRUE)
flist
head(flist)
df <- NULL
for (i in flist) {
  tmp <- read_table(i)
  colnames(tmp)<-header
  tmp$ID <- str_replace(i,".bam.coverage","")
  df <- rbind(df,tmp)
}
#df$batch <- "NovaSeq6000"
df$batch <- "Illumina"
df %>% filter(ID %in% ref$Illumina)  -> novaseq


########


setwd("~/Desktop/KCDC/pangenome/bam.stats/2023.pro.ref_panel.nanopore/")

flist <- system("ls | grep coverage",intern = T)
df <- NULL
for (i in flist) {
  tmp <- read.table(paste0(i))
  #sample_id <- str_replace(unlist(strsplit(i, split = "\\.", perl=T))[4],"_mapped","")
  colnames(tmp)<-header
  #tmp$ID <- sample_id
  tmp$ID <- str_replace(i,"_sorted.bam.coverage","")
  df <- rbind(df,tmp)
}
#df$batch <- "Nanopore"
df$batch <- "ONT"
head(df)

df -> nanopore


setwd("~/Desktop/KCDC/pangenome/bam.stats/2023.pro.KCHIP.Revio/bam.withunmapped/")

flist <- system("ls | grep coverage",intern = T)
revio <- NULL
df <- NULL
for (i in flist) {
  tmp <- read.table(paste0(i))
  #sample_id <- str_replace(unlist(strsplit(i, split = "\\.", perl=T))[4],"_mapped","")
  colnames(tmp)<-header
  #tmp$ID <- sample_id
  tmp$ID <- str_replace(str_replace(i,".pbmm2_hg38.bam.coverage",""),".merge","")
  df <- rbind(df,tmp)
}
#df$batch <- "Revio"
df$batch <- "PacBio"
revio <- rbind(revio,df)
head(revio)
#table(revio$ID)


setwd("~/Desktop/KCDC/pangenome/bam.stats/2023.pro.ref_panel.Revio/")
flist <- system("ls | grep coverage",intern = T)
flist
df <- NULL
for (i in flist) {
  tmp <- read.table(paste0(i))
  colnames(tmp)<-header
  tmp$ID <- str_replace(i,"_sorted.bam.coverage","")
  df <- rbind(df,tmp)
}
#df$batch <- "Revio"
df$batch <- "PacBio"
revio <- rbind(revio,df)
head(revio)

revio %>% filter(ID %in% ref$Revio)  -> revio
head(revio)

revio %>% filter(rname %in% "chrY") %>% mutate(sex=ifelse(coverage>10,"M","F")) %>% select(ID,sex) %>%
  write.table("~/Desktop/KCDC/pangenome/00.datacheck/revio.sex.info.txt",col.names = T,row.names = F,quote = F,sep = "\t")

#######

head(revio)
head(nanopore)
head(novaseq)


df_ori <- rbind(revio,nanopore) %>% rbind(novaseq)

df_ori %>% filter(!(rname %in% c("chrX","chrY","chrM"))) %>% group_by(batch,ID) %>%
  summarise(coverage=mean(coverage),meandepth=mean(meandepth),meanbaseq=mean(meanbaseq),meanmapq=mean(meanmapq)) -> df_1

df_ori %>% filter(rname != "chrM") %>% group_by(batch,ID) %>%
  summarise(coverage=mean(coverage),meandepth=mean(meandepth),meanbaseq=mean(meanbaseq),meanmapq=mean(meanmapq)) %>% 
  summarise(numreads = sum(numreads),covbases=sum(covbases)) -> df_2

df_ori %>% filter(rname != "chrM") %>% group_by(batch) %>%
  summarise(numreads = sum(numreads),coverage=mean(coverage),meandepth=mean(meandepth),meanbaseq=mean(meanbaseq),meanmapq=mean(meanmapq)) -> df_stats
df_ori %>% filter(rname != "chrM") %>% group_by(ID,batch) %>% summarise(numreads = sum(numreads)) %>% group_by(batch) %>% summarise(numreads=mean(numreads))


df_ori %>% filter(rname != "chrM") %>% group_by(ID,batch) %>% summarise(numreads = sum(numreads),covbases=sum(covbases)) %>% 
  group_by(batch) %>% 
  summarise(numreads=mean(numreads),covbases=mean(covbases))


df_ori %>% filter(!(rname %in% c("chrX","chrY","chrM"))) %>% group_by(batch,ID) %>%
  summarise(coverage=mean(coverage),meandepth=mean(meandepth),meanbaseq=mean(meanbaseq),meanmapq=mean(meanmapq)) -> df_auto_1

head(df_auto_1)

df_ori %>% filter(rname != "chrM") %>% group_by(batch) %>%
  summarise(numreads = sum(numreads),coverage=mean(coverage),meandepth=mean(meandepth),meanbaseq=mean(meanbaseq),meanmapq=mean(meanmapq)) -> df_stats
colnames해서, 수행

df_ori %>% filter(!(rname %in% c("chrX","chrY","chrM"))) %>% group_by(batch,ID) %>% 
  summarise(numreads = sum(numreads),covbases = sum(covbases),endpos = sum(endpos)) %>%
  mutate(coverage = covbases/endpos) %>% head()#-> df_auto_1


head(df_ori)
df_ori %>% filter(!(rname %in% c("chrX","chrY","chrM"))) %>% #group_by(batch,ID) %>%  head()
  #mutate(depth_weight = meandepth * endpos) %>%
  group_by(batch,ID) %>%
  summarise(check = sum(meandepth * endpos)/sum(endpos)) 
  
  #group_by(batch,ID) %>%
  #summarise(check = sum(meandepth * endpos))

df_ori %>% filter(!(rname %in% c("chrX","chrY","chrM"))) %>% #group_by(batch,ID) %>%  head()
  #mutate(depth_weight = meandepth * endpos) %>%
  group_by(batch,ID) %>%
  summarise(check = sum(meandepth * endpos)/sum(endpos))  %>% group_by(batch) %>% summarise(max(check),min(check)) -> a


#a  

#### pandepth
setwd("~/Desktop/KCDC/pangenome/bam.stats/pandepth/")
#flist <- list()

flist = grep(list.files("./"),pattern = "gz", value=TRUE)
flist
head(flist)
df <- NULL
for (i in flist) {
  tmp <- read_table(i)
  #colnames(tmp)<-header
  tmp$file <- i
  #tmp$ID <- str_replace(i,".bam.coverage","")
  df <- rbind(df,tmp)
}
head(df)
head(tmp)
tail(tmp)
head(ref)
ref %>% select(-KBAv2) %>% pivot_longer(cols = 1:3) -> ref_platform

head(ref_platform)
colnames(ref_platform) <- c("KBAv1","batch","ID")

df %>% mutate(ID = str_split_fixed(file,'_',2)[,1]) %>% mutate(ID = str_split_fixed(ID,'\\.',2)[,1]) %>% #head()
  mutate(mapQ = str_replace(str_split_fixed(file,"pandepth",2)[,2],".chr.stat.gz","")) %>%  
  mutate(mapQ = ifelse(str_detect(mapQ,"_"),str_replace(mapQ,'_',""),"q0")) %>% 
  left_join(ref_platform) %>% na.omit() %>% filter(!str_detect(`#Chr`,"RegionLength")) -> df_pandepth


1700/4/17
head(df_pandepth)
df_pandepth %>% filter(!(`#Chr` %in% c("chrX","chrY","chrM"))) %>% #group_by(batch,ID) %>%  head()
  #mutate(depth_weight = meandepth * endpos) %>%
  #group_by(batch,ID) %>%
  mutate(calculate_depth = TotalDepth/Length) %>% #head()
  mutate(calculate_cover = as.numeric(CoveredSite)/Length) -> b 
  #summarise(check = sum(meandepth * endpos)/sum(endpos))


head(df_ori)

df_ori %>% filter(!(rname %in% c("chrX","chrY","chrM"))) %>% #group_by(batch,ID) %>%  head()
  #mutate(depth_weight = meandepth * endpos) %>%
  group_by(batch,ID) %>% head()
  summarise(check = sum(meandepth * endpos)/sum(endpos)) -> a

head(a)
tail(a)
table(a$batch)
head(b)
colnames(b)[1] <- "rname"

b %>% group_by(batch,ID) %>% filter(mapQ == 'q0') %>%
  summarise(sum_TotalDepth = sum(TotalDepth),sum_Length = sum(Length)) %>%
  mutate(check2 = sum_TotalDepth/sum_Length) %>% 
  mutate(batch = ifelse(batch == "Illumina","NovaSeq6000","batch")) %>%
  left_join(a)

## 계산 값이 동일함
###
############## pandepth  by Maq
df_pandepth %>% filter(!(`#Chr` %in% c("chrX","chrY","chrM"))) %>% #group_by(batch,ID) %>%  head()
  mutate(batch = ifelse(batch == "Illumina","NovaSeq6000",ifelse(batch == "Revio","Revio","PromethION"))) %>%
  group_by(batch,mapQ) %>% #head()
  summarise(Length = sum(Length),TotalDepth=sum(TotalDepth),CoveredSite=sum(as.numeric(CoveredSite))) %>%
  mutate(`mean Depth(X)` = TotalDepth/Length) %>% #head()
  mutate(`Coverage(%)` = CoveredSite/Length*100) %>% #head()
  select(batch,mapQ,`Coverage(%)`,`mean Depth(X)`) -> a

table(df_pandepth$batch)
df_pandepth %>% filter(!(`#Chr` %in% c("chrX","chrY","chrM"))) %>% #group_by(batch,ID) %>%  head()
  #mutate(calculate_depth = TotalDepth/Length) %>% #head()
  #mutate(calculate_cover = as.numeric(CoveredSite)/Length) %>% #head()
  mutate(batch = ifelse(batch == "Illumina","NovaSeq6000",ifelse(batch == "Revio","Revio","PromethION"))) %>%
  group_by(batch,mapQ) %>% #head()
  summarise(Length = sum(Length),TotalDepth=sum(TotalDepth),CoveredSite=sum(as.numeric(CoveredSite))) %>%
  mutate(`mean Depth(X)` = TotalDepth/Length) %>% #head()
  mutate(`Coverage(%)` = CoveredSite/Length*100) %>% #head()
  select(batch,mapQ,`Coverage(%)`,`mean Depth(X)`) %>% #head()
  pivot_longer(3:4) %>% #head()
  mutate(batch = factor(batch, levels=platform_factor)) %>%
  ggplot(aes(x=mapQ,y=value,color=batch,group=batch)) +
  geom_line() + 
  theme(axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x =  element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom") + 
  facet_wrap(~name,scales = 'free')
  

  

df_pandepth %>% filter(!(`#Chr` %in% c("chrX","chrY","chrM"))) %>% #group_by(batch,ID) %>%  head()
  mutate(batch = ifelse(batch == "Illumina","NovaSeq6000",ifelse(batch == "Revio","Revio","PromethION"))) %>%
  group_by(batch,mapQ) %>% #head()
  #summarise(Length = sum(Length),TotalDepth=sum(TotalDepth),CoveredSite=sum(as.numeric(CoveredSite))) %>%
  mutate(`mean Depth(X)` = TotalDepth/Length) %>% #head()
  mutate(`Coverage(%)` = as.numeric(CoveredSite)/Length*100) %>% #head()
  select(batch,ID,mapQ,`Coverage(%)`,`mean Depth(X)`) %>% #head()
  pivot_longer(4:5) %>% #head()
  mutate(batch = factor(batch, levels=platform_factor)) %>% #head()
  ggplot(aes(x=mapQ,y=value,fill=batch)) + 
  geom_boxplot() + 
  theme(axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x =  element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom") + 
  facet_wrap(~name,scales = 'free')
  
  




#platform_factor <- c("NovaSeq6000","Revio","PromethION")
platform_factor <- c("Illumina","PacBio","ONT")

df_ori %>% filter(!(rname %in% c("chrX","chrY","chrM"))) %>% #group_by(batch,ID) %>%  head()
  mutate(batch = ifelse(batch == "Nanopore","PromethION",batch)) %>% #head()
  mutate(batch = factor(batch, levels=platform_factor)) %>%
  group_by(batch,ID) %>% #head()
  summarise(`mean Depth(X)` = sum(meandepth * endpos)/sum(endpos),
            `Coverage(%)` = sum(covbases)/sum(endpos)*100,
            `mean BaseQ` = sum(meanbaseq * endpos)/sum(endpos),
            `mean MapQ` = sum(meanmapq * endpos)/sum(endpos)) -> a

df_ori %>% filter(!(rname %in% c("chrX","chrY","chrM"))) %>% #group_by(batch,ID) %>%  head()
  mutate(batch = ifelse(batch == "Nanopore","PromethION",batch)) %>% #head()
  mutate(batch = factor(batch, levels=platform_factor)) %>%
  group_by(batch,ID) %>% #head()
  summarise(numreads = sum(numreads),
            covbases = sum(covbases)) -> a2


head(a)
head(a)
table(a$batch)
a %>% pivot_longer(3:6) %>% #head()
  ggplot(aes(x=batch,y=value,fill=batch)) + 
  geom_boxplot() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x =  element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom") + 
  facet_wrap(~name,scales = "free",nrow = 1) 

#write.table(df_ori,"~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/bam.stats.merge.txt",col.names = T,row.names = F,quote = F,sep = "\t")




a %>% #pivot_longer(3:6) %>% #head()
  #mutate(batch = ifelse(batch == "Nanopore","PromethION",batch)) %>% #head()
  mutate(batch = factor(batch, levels=platform_factor)) -> a1

a1 %>% left_join(a2) %>%
  pivot_longer(3:8) %>%
  group_by(batch,name) %>% 
  summarise(mean = mean(value)) %>% pivot_wider(names_from = name,values_from = mean) -> a2

head(df_pandepth)
df_pandepth %>% filter(!(`#Chr` %in% c("chrX","chrY","chrM"))) %>% #group_by(batch,ID) %>%  head()
  filter(mapQ == 'q0') %>%  mutate(`GC(%)` = as.numeric(`GC(%)`)) %>%
  group_by(batch,ID) %>% summarise(check = mean(`GC(%)`))
  #summarise(`GC(%)` = sum(`GC(%)` * Length)/sum(Length)) %>% head()#->b
            
  #mutate(`GC (%)` = `GC(%)`/Length) %>% head()
  #mutate(calculate_cover = as.numeric(CoveredSite)/Length) -> b 
#summarise(check = sum(meandepth * endpos)/sum(endpos))

#####
  

head(df_pandepth)
df_pandepth %>%

head(df)
head(df_ori)
#################################
head(df_all)
df %>% pivot_longer(3:7) %>%
  ggplot(aes(y=value,x=batch,fill=batch)) + 
  geom_boxplot() +
  labs(title = "Sample-level Mapping Quality metrics by batch") +
  theme(plot.title = element_text(size=15,hjust = 0.5)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size=10)) +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(size=15)) + 
  theme(legend.position = "bottom",
        legend.text = element_text(size=15)) +
  facet_wrap(~name,scales = "free_y",nrow = 1)

chrom_factor <- paste0("chr",seq(1:22))

df_ori %>% filter(!(rname %in% c("chrX","chrY","chrM"))) %>% #group_by(batch,ID) %>%  head()
  mutate(rname = factor(rname,levels=chrom_factor)) %>%
  ggplot(aes(y=meandepth,x=batch,fill=batch)) + 
  geom_boxplot() + 
  labs(title = "Sample-level Mapping Depth by batch") +
  geom_hline(yintercept = 30,linetype="dotted",color='red') +
  theme(plot.title = element_text(size=15,hjust = 0.5)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size=10)) +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(size=10)) + 
  theme(legend.position = "bottom",
        legend.text = element_text(size=15)) +
  facet_wrap(~rname,nrow = 4)



####### VCF stats

setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/99.kogo/2024.kogo/vcf.pass.stats/")
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/99.kogo/2024.kogo/vcf.pass.stats/rmMXY.stats")
df <- read.table("filter.info.txt",header = T)
ref <- read.table("~/Desktop/KCDC/pangenome/00.datacheck/KBA.Long_Revio_Nanopore_short.IDmatchinagtable.txt",header = T) %>% na.omit()
head(ref)
#colnames(ref)[1:4]<-c("Revio","PromethION","NovaSeq6000","KBAv1")
colnames(ref)[1:4]<-c("PacBio","ONT","Illumina","KBAv1")
ref %>% select(-KBAv2) %>% pivot_longer(cols = 1:3,names_to = "platform",values_to = "ID") ->ref

head(df)
df %>% mutate(ID = str_split_fixed(file,"\\.",2)[,1]) %>% filter(ID == "Pangenome") -> a
df %>% mutate(ID = str_split_fixed(file,"\\.",2)[,1]) %>% filter(ID != "Pangenome") %>% 
  mutate(ID = str_replace(ID,".sorted","")) %>%
#  mutate(platform = ifelse(str_detect(file,"pbmm2"),"Revio",ifelse(str_detect(file,"bwa"),"NovaSeq6000","PromethION"))) -> df
  mutate(platform = ifelse(str_detect(file,"pbmm2"),"PacBio",ifelse(str_detect(file,"bwa"),"Illumina","ONT"))) -> df

head(df)
platform_factor <- c("Illumina","PacBio","ONT")
#mutate(platform = factor(platform, levels=platform_factor)) %>%
  


library(ggpattern)
library(ggplot2)

head(df)
head(ref)

df %>% left_join(ref) %>% select(KBAv1) %>% unique() %>% mutate(platform = "z",SNPs=0,indels=0) -> df_add

head(df)
df %>% left_join(ref) %>%  select(KBAv1,platform,records,SNPs,indels,others,multiallelic_sites,multiallelic_SNP_sites,ts_tv) %>% #head
  pivot_longer(3:9) %>%
  group_by(platform,name) %>%
  summarise(mean=mean(value)) %>%
  pivot_wider(names_from = name,values_from = mean) %>%
  select(platform,records,SNPs,indels,others,multiallelic_sites,multiallelic_SNP_sites,ts_tv) -> a

df %>% left_join(ref) %>%  select(KBAv1,platform,SNPs,indels) %>% #head
  #rbind(df_add) %>%
  mutate(platform = factor(platform, levels=platform_factor)) %>%
  pivot_longer(cols = SNPs:indels, names_to = "type", values_to = "value") %>% group_by(platform,type) %>%
  summarise(value = mean(value)) -> vcf_mean_df
head(vcf_mean_df)
vcf_mean_df$theme <- c("Mean")
vcf_mean_df$KBAv1 <- 'ALL'
head(vcf_mean_df)

library(ggforce)

df %>% left_join(ref) %>%  select(KBAv1,platform,SNPs,indels) %>% #head
  rbind(df_add) %>%
  mutate(theme = "Count by Sample") %>% 
  mutate(platform = factor(platform, levels=platform_factor)) %>%
  pivot_longer(cols = SNPs:indels, names_to = "type", values_to = "value") %>% #head()
  rbind(vcf_mean_df) %>%
  ggplot(aes(x = interaction(platform,KBAv1), y = value, fill = platform, pattern = type)) +
  geom_bar_pattern(stat = "identity", 
                   position = "stack", 
                   pattern_density = 0.01, 
                   #pattern_fill = "white",
                   pattern_fill = "black",
                   pattern_spacing = 0.01,width = 0.7) +
  scale_pattern_manual(values = c(SNPs = "none", indels = "stripe")) +
  labs(x = "Sample", y = "Count", title = "Comparison of SNPs and Indels by platform") +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  guides(pattern = guide_legend(override.aes = list(fill = "grey")), fill = guide_legend(override.aes = list(pattern = "none"))) + 
  facet_row(~theme,space = "free")

#c(platform_factor,"z")
df %>% left_join(ref) %>%  select(KBAv1,platform,SNPs,indels) %>% #head
  #rbind(df_add) %>%
  mutate(platform = factor(platform, levels=c(platform_factor,"z"))) %>%
  pivot_longer(cols = SNPs:indels, names_to = "type", values_to = "value") %>% #head()
  ggplot(aes(x = platform, y = value, fill = platform, pattern = type)) +
  geom_bar_pattern(stat = "identity", 
                   position = "stack", 
                   pattern_density = 0.01, 
                   #pattern_fill = "white",
                   pattern_fill = "black",
                   pattern_spacing = 0.01,width = 0.7) +
  scale_pattern_manual(values = c(SNPs = "none", indels = "stripe")) +
  labs(x = "Sample", y = "Count", title = "Comparison of SNPs and Indels by platform") +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  guides(pattern = guide_legend(override.aes = list(fill = "grey")), fill = guide_legend(override.aes = list(pattern = "none"))) + 
  facet_grid(~KBAv1)


df %>% left_join(ref) %>%  select(KBAv1,platform,SNPs,indels) %>% #head
  #mutate(platform = factor(platform, levels=c(platform_factor,"z"))) %>%
  mutate(platform = factor(platform, levels=platform_factor)) %>%
  pivot_longer(cols = SNPs:indels, names_to = "type", values_to = "value") %>% #head()
  ggplot(aes(x = platform, y = value, fill = platform)) +
  geom_barplot()
  



df %>% left_join(ref) %>%  select(KBAv1,platform,SNPs,indels) %>% #head
  rbind(df_add) %>%
  mutate(platform = factor(platform, levels=c(platform_factor,"z"))) %>%
  pivot_longer(cols = SNPs:indels, names_to = "type", values_to = "value") %>% #head()
  ggplot(aes(x = interaction(platform,KBAv1), y = value, fill = platform, pattern = type)) +
  geom_bar_pattern(stat = "identity", 
                   position = "stack", 
                   pattern_density = 0.01, 
                   #pattern_fill = "white",
                   pattern_fill = "black",
                   pattern_spacing = 0.01,width = 0.7) +
  scale_pattern_manual(values = c(SNPs = "none", indels = "stripe")) +
  labs(x = "Sample", y = "Count", title = "Comparison of SNPs and Indels by platform") +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none") + 
  guides(pattern = guide_legend(override.aes = list(fill = "grey")), fill = guide_legend(override.aes = list(pattern = "none"))) -> p1

p1

vcf_mean_df %>% #head()
  ggplot(aes(x = interaction(platform,KBAv1), y = value, fill = platform, pattern = type)) +
  geom_bar_pattern(stat = "identity", 
                   position = "stack", 
                   pattern_density = 0.01, 
                   #pattern_fill = "white",
                   pattern_fill = "black",
                   pattern_spacing = 0.01,width = 0.7) +
  scale_pattern_manual(values = c(SNPs = "none", indels = "stripe")) +
  labs(x='Platform',title = "Mean of Count by Platform") +
  theme(#axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") -> p11
p11

p2
library(ggplot2)
df %>% left_join(ref) %>%  select(KBAv1,platform,SNPs,indels) %>% #head
  pivot_longer(cols = SNPs:indels, names_to = "type", values_to = "value") %>% #head()
  mutate(platform = factor(platform, levels=platform_factor)) %>% #head()
  ggplot(aes(x = interaction(platform,KBAv1), y = value, fill = platform)) +
  geom_bar(stat = "identity") + 
  theme(legend.position = "bottom",
        legend.title = element_blank()) -> p2
p2 <- get_plot_component(p2, "guide-box", return_all = TRUE)
#p2 <- get_legend(p2)
p2
p2[[3]]


df %>% left_join(ref) %>%  select(KBAv1,platform,SNPs,indels) %>% #head
  mutate(INDELs = indels) %>% select(-indels) %>%
  pivot_longer(cols = INDELs:SNPs, names_to = "type", values_to = "value") %>% #head()
  ggplot(aes(x = interaction(platform,KBAv1), y = value, pattern = type)) +
  geom_bar_pattern(stat = "identity", 
                   position = "stack", 
                   pattern_density = 0.01, 
                   #pattern_fill = "white",
                   pattern_spacing = 0.01,width = 0.7) +
  scale_pattern_manual(values = c(SNPs = "none", INDELs = "stripe")) +
  labs(x = "Sample", y = "Count", title = "Comparison of SNPs and Indels by platform") +
  theme(axis.text.x = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  guides(pattern = guide_legend(override.aes = list(fill = "grey"),reverse=T)) -> p3

p3 <- get_plot_component(p3, "guide-box", return_all = TRUE)
p3


plot_grid(p2[[3]],p3[[3]],nrow = 1,rel_widths = c(3,2)) -> p4


plot_grid(p1,p4,ncol=1,rel_heights = c(5,1))



plot_grid(p1,p5,nrow = 1,rel_widths = c(10,1)) -> p15
p15

plot_grid(p15,p4,ncol=1,rel_heights = c(5,1))

red = '#D55E00'

green = '#009E73'



df %>% left_join(ref) %>%  select(KBAv1,platform,SNPs,indels) %>% #head
  mutate(platform = factor(platform, levels=c(platform_factor))) %>%
  pivot_longer(cols = SNPs:indels, names_to = "type", values_to = "value") %>% #head()
  mutate(type = factor(type, levels=c("SNPs","indels"))) %>%
  ggplot(aes(x = platform, y = value, fill = platform)) + 
  geom_boxplot() + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank()) + 
  facet_wrap(~type,scales = 'free')
  


############VC

setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/99.kogo/2024.kogo/concordance/stats")

#Pangenome.Revio_kchip.same_17sample.pbmm2_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlySNP_withtag.vcf.gz.AF_stat.txt.query
#Pangenome.Revio_kchip.same_17sample.pbmm2_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlyindels_withtag.vcf.gz.AF_stat.txt.query
#Pangenome.nanopore.same_17sample.minimap2_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlySNP_withtag.vcf.gz.AF_stat.txt.query
#Pangenome.nanopore.same_17sample.minimap2_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlyindels_withtag.vcf.gz.AF_stat.txt.query
#Pangenome.wgs.same_17sample.bwa_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlySNP_withtag.vcf.gz.AF_stat.txt.query
#Pangenome.wgs.same_17sample.bwa_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlyindels_withtag.vcf.gz.AF_stat.txt.query

revio_snp <- read_table("Pangenome.Revio_kchip.same_17sample.pbmm2_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlySNP_withtag.vcf.gz.AF_stat.txt.query",col_names = F)
revio_indel <- read_table("Pangenome.Revio_kchip.same_17sample.pbmm2_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlyindels_withtag.vcf.gz.AF_stat.txt.query",col_names = F)

wgs_snp <- read_table("Pangenome.wgs.same_17sample.bwa_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlySNP_withtag.vcf.gz.AF_stat.txt.query",col_names = F)
wgs_indel <- read_table("Pangenome.wgs.same_17sample.bwa_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlyindels_withtag.vcf.gz.AF_stat.txt.query",col_names = F)


nanopore_snp <- read_table("Pangenome.nanopore.same_17sample.minimap2_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlySNP_withtag.vcf.gz.AF_stat.txt.query",col_names = F)
nanopore_indel <- read_table("Pangenome.nanopore.same_17sample.minimap2_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlyindels_withtag.vcf.gz.AF_stat.txt.query",col_names = F)

head(wgs_indel)
head(revio_snp)
head(revio_indel)
head(wgs_snp)
head(wgs_indel)
head(nanopore_snp)
head(nanopore_indel)

colnames(revio_snp) <- c("ID","QUAL","AC")
colnames(revio_indel) <- c("ID","QUAL","AC")
colnames(wgs_snp) <- c("ID","QUAL","AC")
colnames(wgs_indel) <- c("ID","QUAL","AC")
colnames(nanopore_snp) <- c("ID","QUAL","AC")
colnames(nanopore_indel) <- c("ID","QUAL","AC")

library(ggplot2)
library(ggVennDiagram)
library(VennDiagram)
#platform_factor <- c("NovaSeq6000","Revio","PromethION")
platform_factor <- c("Illumina","PacBio","ONT")
x=list(A = wgs_snp$ID, B = revio_snp$ID, C = nanopore_snp$ID)
x=list(Illumina = wgs_snp$ID, PacBio = revio_snp$ID, ONT = nanopore_snp$ID)


ggVennDiagram(x, label = "count",set_size = TRUE) + 
  scale_fill_gradient(low = "lightblue", high ="white") + 
  theme(legend.position = 'none')


vd <- VennDiagram::venn.diagram(x, filename = NULL,)
grid::grid.draw(vd)
draw.triple.venn()

wgs_snp_only<-wgs_snp %>% filter(!(ID %in% revio_snp$ID)) %>% filter(!(ID %in% nanopore_snp$ID)) %>% dim()
revio_snp_only<-revio_snp %>% filter(!(ID %in% wgs_snp$ID)) %>% filter(!(ID %in% nanopore_snp$ID)) %>% dim()
nanopore_snp_only<-nanopore_snp %>% filter(!(ID %in% wgs_snp$ID)) %>% filter(!(ID %in% revio_snp$ID)) %>% dim()

wgs_revio_snp <- wgs_snp %>% filter(ID %in% revio_snp$ID) %>% filter(!(ID %in% nanopore_snp$ID)) %>% dim()
wgs_nanopore_snp <- wgs_snp %>% filter(ID %in% nanopore_snp$ID) %>% filter(!(ID %in% revio_snp$ID)) %>% dim()
revio_nanopore_snp <- revio_snp %>% filter(ID %in% nanopore_snp$ID) %>% filter(!(ID %in% wgs_snp$ID)) %>% dim()

all_snp <- revio_snp %>% filter(ID %in% nanopore_snp$ID) %>% filter(ID %in% wgs_snp$ID) %>% dim()



x=c()
library(venneuler)
v <-venneuler::venneuler(c(A=wgs_snp_only[1],B=revio_snp_only[1],C=nanopore_snp_only[1],"A&B"=wgs_revio_snp[1],"A&C"=wgs_nanopore_snp[1],"C&B"=revio_nanopore_snp[1],"A&B&C"=all_snp[1]))
plot(v)
#v$labels <- c("")  
#platform_factor <- c("NovaSeq6000","Revio","PromethION")
#x=list(A = wgs_indel$ID, B = revio_snp$ID, C = nanopore_snp$ID)
#x=list(NovaSeq6000 = wgs_indel$ID, Revio = revio_indel$ID, PromethION = nanopore_indel$ID)

x=list(Illumina = wgs_snp$ID, PacBio = revio_snp$ID, ONT = nanopore_snp$ID)

ggVennDiagram(x, label = "count",set_size = 0) + 
  scale_fill_gradient(low = "white", high ="lightblue") + 
  #scale_fill_distiller(palette = "RdBu") + 
  theme(legend.position = 'none') -> p1



x=list(Illumina = wgs_indel$ID, PacBio = revio_indel$ID, ONT = nanopore_indel$ID)

ggVennDiagram(x, label = "count",set_size = 0) + 
  scale_fill_gradient(low = "white", high ="lightblue") + 
  #scale_fill_distiller(palette = "RdBu") + 
  theme(legend.position = 'none') -> p2



cowplot::plot_grid(p1,p2,labels = c("A","B"))







head(wgs_snp)
head(revio_snp)
head(nanopore_snp)

head(revio_indel)
head(wgs_indel)
head(nanopore_indel)

revio_snp %>% select(ID) %>% rbind(nanopore_snp %>% select(ID)) %>% unique() %>% count(ID %in% wgs_snp$ID) #1116195
revio_indel %>% select(ID) %>% rbind(nanopore_indel %>% select(ID)) %>% unique() %>% count(ID %in% wgs_indel$ID) #342591

revio_snp %>% select(ID) %>% rbind(nanopore_snp %>% select(ID)) %>% unique() %>% filter(!(ID %in% wgs_snp$ID)) -> a1
revio_indel %>% select(ID) %>% rbind(nanopore_indel %>% select(ID)) %>% unique() %>% filter(!(ID %in% wgs_indel$ID)) -> a2
a1$type <- 'snp'
a2$type <- 'indel' 

a1%>% rbind(a2) %>% count(type)
a1%>% rbind(a2) %>% dim() #1458786
a1%>% rbind(a2) -> a

head(a)
dim(a)
revio_snp %>% filter(ID %in% nanopore_snp$ID) %>% filter(ID %in% a$ID) %>% rename(Revio_AC = AC) %>% select(-QUAL)-> revio_snp_uniqinterNanopore
nanopore_snp %>% filter(ID %in% revio_snp$ID) %>% filter(ID %in% a$ID) %>% rename(nanopore_AC = AC)%>% select(-QUAL) -> Nanopore_snp_uniqinterrevio


revio_indel %>% filter(ID %in% nanopore_indel$ID) %>% filter(ID %in% a$ID) %>% rename(Revio_AC = AC) %>% select(-QUAL)-> revio_indel_uniqinterNanopore
nanopore_indel %>% filter(ID %in% revio_indel$ID) %>% filter(ID %in% a$ID) %>% rename(nanopore_AC = AC)%>% select(-QUAL) -> Nanopore_indel_uniqinterrevio


dim(revio_snp_uniqinterNanopore) #423840
dim(revio_indel_uniqinterNanopore) #60687

#423840 + 60687 : 484527

revio_snp_uniqinterNanopore %>% left_join(Nanopore_snp_uniqinterrevio) %>% #count(Revio_AC == nanopore_AC)
  mutate(revio_maf = Revio_AC/34,nanopore_maf = nanopore_AC/34) %>% 
  count(revio_maf > 0.05)


revio_snp_uniqinterNanopore %>% left_join(Nanopore_snp_uniqinterrevio) %>% #count(Revio_AC == nanopore_AC)
  rbind(revio_indel_uniqinterNanopore %>% left_join(Nanopore_indel_uniqinterrevio)) %>% #dim()
  mutate(revio_maf = Revio_AC/34,nanopore_maf = nanopore_AC/34) %>% 
  count(revio_maf > 0.1)
#  count(nanopore_maf > 0.05)
 
328131/484527

revio_indel_uniqinterNanopore %>% left_join(Nanopore_indel_uniqinterrevio) %>% #dim()
  mutate(revio_maf = Revio_AC/34,nanopore_maf = nanopore_AC/34) %>% 
  count(revio_maf > 0.05)
  
#: revio frequency 기준, revio+ nanopore 동시에 있는것
## SNP
359417/423840 #0.848

## indel
48781/60687 #0.804
## totl : snp +indel
408198/484527 #0.842


revio_snp %>% select(ID) %>% rbind(nanopore_snp %>% select(ID)) %>% dup %>% count(ID %in% wgs_snp$ID) #567329
revio_snp %>% count(ID %in% wgs_snp$ID)

revio_indel_uniqinterNanopore %>% left_join(Nanopore_indel_uniqinterrevio) %>% 
  mutate(ref = str_split_fixed(ID,"_",4)[,3],alt = str_split_fixed(ID,"_",4)[,4]) %>%
  mutate(ref_length = str_length(ref),alt_length = str_length(alt)) %>%
  mutate(indel_length = ifelse(ref_length > alt_length,ref_length,alt_length)) %>% 
  #summary(indel_length)
  summarise(mean(indel_length),min(indel_length),max(indel_length))


revio_indel_uniqinterNanopore %>% left_join(Nanopore_indel_uniqinterrevio) %>%  #dim()
  mutate(ref = str_split_fixed(ID,"_",4)[,3],alt = str_split_fixed(ID,"_",4)[,4]) %>%
  mutate(ref_length = str_length(ref),alt_length = str_length(alt)) %>%
  mutate(indel_length = ifelse(ref_length > alt_length,ref_length,alt_length)) %>% count(indel_length %in% c(2:4))

# 2: 23203
23203/60687 # 0.38
# 2~3
33567/60687 # 0.55

revio_indel_uniqinterNanopore %>% left_join(Nanopore_indel_uniqinterrevio) %>% 
  mutate(ref = str_split_fixed(ID,"_",4)[,3],alt = str_split_fixed(ID,"_",4)[,4]) %>%
  mutate(ref_length = str_length(ref),alt_length = str_length(alt)) %>%
  mutate(indel_length = ifelse(ref_length > alt_length,ref_length,alt_length)) %>% 
  ggplot(aes(y=indel_length)) +
  geom_boxplot()

head(a)
head(a2)
head(revio_indel)
wgs_indel
revio_indel %>% rbind(nanopore_indel) %>% select(ID) %>% unique() %>% filter(!(ID %in% wgs_indel$ID)) %>% #dim()
  mutate(ref = str_split_fixed(ID,"_",4)[,3],alt = str_split_fixed(ID,"_",4)[,4]) %>%
  mutate(ref_length = str_length(ref),alt_length = str_length(alt)) %>%
  mutate(indel_length = ifelse(ref_length > alt_length,ref_length,alt_length)) -> indel_length


summary(indel_length)
dim(indel_length)
head(indel_length)
indel_length %>% count(indel_length == 2) #145088
145088/342591 # 42.3%

indel_length %>% count(indel_length %in% c(2,3)) #198954
198954/342591 # 58.1 %  

indel_length %>% count(indel_length >= 50)  #20024

20024/342591 # 0.05%

##
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/99.kogo/2024.kogo/concordance/concordance/")




setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/99.kogo/2024.kogo/concordance/concordance/snp/")
flist <- list()

flist = grep(list.files("./"),pattern = "by_sample", value=TRUE)
flist
head(flist)
df <- NULL

#revio <- read.table("concordance_revio_nanopore.by_sample.txt",header = T)
#wgs <- read.table("concordance_wgs_nanopore.by_sample.txt",header = T)
#nanopore <- read.table("concordance_wgs_nanopore.by_sample.txt",header = T)

#revio$batch <- "Revio"
#wgs$batch <- "NovaSeq6000"
#nanopore$batch <- "PromethION"


df <- NULL
flist

for (i in flist) {
  data <- read.table(i,header = T)
  data <- cbind(data, TP=(data$REF.REF + data$ALT_1.ALT_1 + data$ALT_2.ALT_2))
  data <- cbind(data, FN=(data$ALT_1.REF + data$ALT_2.REF + data$ALT_2.ALT_1))
  data <- cbind(data, FP=(data$REF.ALT_1 + data$REF.ALT_2 + data$ALT_1.ALT_2))
  
  data$Sensitivity <- data$TP/(data$TP+data$FN)
  data$Precision <- data$TP/(data$TP+data$FP)
  data$Accuracy <- data$TP/(data$TP+data$FP+data$FN)
  
  a <- data %>% #inner_join(maf1) %>% 
    select(sample,Sensitivity,Precision,Accuracy) %>% 
    pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% 
    mutate(ID = str_replace(i,".by_sample.txt",""))
  df <- rbind(df,a)
    #mutate(title = str_remove(str_remove(i,".txt"),"QulityMetrix_"))
  df <- rbind(df,a)
}
head(df)
head(df$variant)
df$variant <- "SNP"
df_snp <- df



setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/99.kogo/2024.kogo/concordance/concordance/indel/")
flist <- list()

flist = grep(list.files("./"),pattern = "by_sample", value=TRUE)
flist
head(flist)
df <- NULL


for (i in flist) {
  data <- read.table(i,header = T)
  data <- cbind(data, TP=(data$REF.REF + data$ALT_1.ALT_1 + data$ALT_2.ALT_2))
  data <- cbind(data, FN=(data$ALT_1.REF + data$ALT_2.REF + data$ALT_2.ALT_1))
  data <- cbind(data, FP=(data$REF.ALT_1 + data$REF.ALT_2 + data$ALT_1.ALT_2))
  
  data$Sensitivity <- data$TP/(data$TP+data$FN)
  data$Precision <- data$TP/(data$TP+data$FP)
  data$Accuracy <- data$TP/(data$TP+data$FP+data$FN)
  
  a <- data %>% #inner_join(maf1) %>% 
    select(sample,Sensitivity,Precision,Accuracy) %>% 
    pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% 
    mutate(ID = str_replace(i,".by_sample.indel.txt",""))
  df <- rbind(df,a)
  #mutate(title = str_remove(str_remove(i,".txt"),"QulityMetrix_"))
  df <- rbind(df,a)
}
head(df)
head(df$variant)
df$variant <- "INDEL"
df_indel <- df


df <- rbind(df_snp,df_indel)
head(df)
tail(df)


pdf(file = "~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/concordance.pdf", width = 4, height = 4)
df %>% #mutate(ID = str_replace(ID,"concordance_","")) %>%
  mutate(ID = ifelse(ID == "concordance_revio_nanopore","PacBio vs ONT",ifelse(ID == "concordance_wgs_revio","Illumina vs PacBio",'Illumina vs ONT'))) %>%
  mutate(ID = factor(ID,levels=c("Illumina vs PacBio","Illumina vs ONT","PacBio vs ONT"))) %>% 
  mutate(type = ifelse(type == "Accuracy",'ACC',ifelse(type == 'Precision','PPV','TPR'))) %>%
  mutate(variant = factor(variant, levels=c("SNP","INDEL"))) %>%
  ggplot(aes(x=type,y=Val,fill=type)) + 
  geom_boxplot() + 
  #theme_step1() + 
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
        ) + 
  facet_wrap(variant~ID)
dev.off()

df %>% #mutate(ID = str_replace(ID,"concordance_","")) %>%
  mutate(ID = ifelse(ID == "concordance_revio_nanopore","Revio vs PromethION",ifelse(ID == "concordance_wgs_revio","NovaSeq6000 vs Revio",'NovaSeq6000 vs PromethION'))) %>%
  mutate(ID = factor(ID,levels=c("NovaSeq6000 vs Revio","NovaSeq6000 vs PromethION","Revio vs PromethION"))) %>%
  mutate(variant = factor(variant, levels=c("SNP","INDEL"))) %>% #head()
  group_by(variant,ID,type) %>%  #head()
  summarise(mean = mean(Val)) %>% #head()
  pivot_wider(names_from = type, values_from = mean ) -> a
  


df %>% #mutate(ID = str_replace(ID,"concordance_","")) %>%
  mutate(ID = ifelse(ID == "concordance_revio_nanopore","PacBio vs ONT",ifelse(ID == "concordance_wgs_revio","Illumina vs PacBio",'Illumina vs ONT'))) %>%
  mutate(ID = factor(ID,levels=c("Illumina vs PacBio","Illumina vs ONT","PacBio vs ONT"))) %>% 
  mutate(type = ifelse(type == "Accuracy",'ACC',ifelse(type == 'Precision','PPV','TPR'))) %>%
  #mutate(variant = factor(variant, levels=c("SNP","INDEL"))) %>%
  filter(variant == "SNP") %>%
  ggplot(aes(x=type,y=Val,fill=type)) + 
  geom_boxplot() + 
#  labs(title="test") +
  ylim(c(0.980,1)) + 
  theme_step1() + 
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size=18)) +
  facet_wrap(~ID) -> p1

df %>% #mutate(ID = str_replace(ID,"concordance_","")) %>%
  mutate(ID = ifelse(ID == "concordance_revio_nanopore","PacBio vs ONT",ifelse(ID == "concordance_wgs_revio","Illumina vs PacBio",'Illumina vs ONT'))) %>%
  mutate(ID = factor(ID,levels=c("Illumina vs PacBio","Illumina vs ONT","PacBio vs ONT"))) %>% 
  mutate(type = ifelse(type == "Accuracy",'ACC',ifelse(type == 'Precision','PPV','TPR'))) %>%
  #mutate(variant = factor(variant, levels=c("SNP","INDEL"))) %>%
  filter(variant == "INDEL") %>%
  ggplot(aes(x=type,y=Val,fill=type)) + 
  geom_boxplot() + 
  ylim(c(0.980,1)) + 
  theme_step1() + 
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size=18)) + 
  facet_wrap(~ID) -> p2


pdf(file = "~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/concordance.pdf", width = 15, height = 3)
cowplot::plot_grid(p1,p2,labels = c("A","B"))
dev.off()
 