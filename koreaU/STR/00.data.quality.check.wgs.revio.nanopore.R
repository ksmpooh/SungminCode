## bam coverage
library(tidyverse)
setwd("~/Desktop/KCDC/pangenome/bam.stats/")
header = c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")
df_all <- NULL

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
df$batch <- "Nanopore(M)"
df_all <- rbind(df_all,df)

setwd("~/Desktop/KCDC/pangenome/bam.stats/2023.pro.KCHIP.Revio/bam.withunmapped/")

flist <- system("ls | grep coverage",intern = T)
df <- NULL
for (i in flist) {
  tmp <- read.table(paste0(i))
  #sample_id <- str_replace(unlist(strsplit(i, split = "\\.", perl=T))[4],"_mapped","")
  colnames(tmp)<-header
  #tmp$ID <- sample_id
  tmp$ID <- str_replace(str_replace(i,".pbmm2_hg38.bam.coverage",""),".merge","")
  df <- rbind(df,tmp)
}
df$batch <- "Revio(D)"
df_all <- rbind(df_all,df)
head(df_all)



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
df$batch <- "Revio(M)"
df_all <- rbind(df_all,df)
head(df_all)


setwd("~/Desktop/KCDC/pangenome/bam.stats/wgs/")
flist <- system("ls | grep coverage",intern = T)
flist
df <- NULL
for (i in flist) {
  tmp <- read.table(paste0(i))
  colnames(tmp)<-header
  tmp$ID <- str_replace(i,".bam.coverage","")
  df <- rbind(df,tmp)
}
df$batch <- "Illumina"
df_all <- rbind(df_all,df)
head(df_all)


head(df)

df_all %>% count(batch)


df_all %>% filter(rname != "chrM") %>% mutate(batch = ifelse(batch %in% c("Revio(D)","Revio(M)"),"Revio",batch)) %>% group_by(batch,ID) %>%
  summarise(numreads = sum(numreads),coverage=mean(coverage),meandepth=mean(meandepth),meanbaseq=mean(meanbaseq),meanmapq=mean(meanmapq)) -> df_3
head(df_3)


df_all %>% filter(rname != "chrM") %>% group_by(batch,ID) %>%
  summarise(numreads = sum(numreads),coverage=mean(coverage),meandepth=mean(meandepth),meanbaseq=mean(meanbaseq),meanmapq=mean(meanmapq)) -> df

df_all %>% filter(rname != "chrM") %>% group_by(batch) %>%
  summarise(numreads = sum(numreads),coverage=mean(coverage),meandepth=mean(meandepth),meanbaseq=mean(meanbaseq),meanmapq=mean(meanmapq)) -> df
df_all %>% filter(rname != "chrM") %>% group_by(ID,batch) %>% summarise(numreads = sum(numreads)) %>% group_by(batch) %>% summarise(numreads=mean(numreads))

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



df_all %>% filter(rname != "chrM") %>% #head()  
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


head(df_3)

df_all %>% filter(rname != "chrM") %>% mutate(batch = ifelse(batch %in% c("Revio(D)","Revio(M)"),"Revio",batch)) %>% group_by(batch) %>% 
  summarise(numreads = sum(numreads),coverage=mean(coverage),meandepth=mean(meandepth),meanbaseq=mean(meanbaseq),meanmapq=mean(meanmapq)) -> df
head(df)

df_3 %>% select(-numreads) %>%pivot_longer(3:6) %>% 
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

