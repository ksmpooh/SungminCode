
setwd("~/Desktop/KCDC/long_read/2022/short-read/01.bam.stats/")
library(tidyverse)
ls_dedup <- system("ls | grep dedup.bam.coverage", intern = T)
ls_nodup <- system("ls | grep sorted.bam.coverage", intern = T)


#setwd("~/Desktop/KCDC/HLA_seq/02.bam.stats/short-read/trimmed/")
#ls_trimmed <- system("ls ~/Desktop/KCDC/HLA_seq/02.bam.stats/short-read/trimmed/ | grep bam.coverage", intern = T)
#setwd("~/Desktop/KCDC/HLA_seq/02.bam.stats/short-read/trimmed.dedup/")
#ls_trimmed_dedup <- system("ls ~/Desktop/KCDC/HLA_seq/02.bam.stats/short-read/trimmed.dedup/ | grep bam.coverage", intern = T)

head(ls_nodup)

#sample_id <- unlist(strsplit(ls_dedup[1], split = "\\.", perl=T))[4]
#sample_id
#file_suffix <- ".subreads.gt-cl.ext.stats"
#file_suffix <- ".merge.stats"
###
###
#unlist(strsplit(i, split = "\\.", perl=T))[4]
df_dedup = data.frame()
header = c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")
for (i in ls_dedup) {
  tmp <- read.table(paste0(i))
  sample_id <- unlist(strsplit(i, split = "\\.", perl=T))[4]
  colnames(tmp)<-header
  tmp$ID <- sample_id
  df_dedup <- rbind(df_dedup,tmp)
}

df_nodup = data.frame()
header = c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")
for (i in ls_nodup) {
  tmp <- read.table(paste0(i))
  sample_id <- unlist(strsplit(i, split = "\\.", perl=T))[4]
  colnames(tmp)<-header
  tmp$ID <- sample_id
  df_nodup <- rbind(df_nodup,tmp)
}

head(df_dedup)
mean(df_nodup$meandepth)

df_dedup$type <- 'MarkDup'
df_nodup$type <- "BeforeRmDup"

df <- rbind(df_dedup,df_nodup)

head(df)
par(mfrow = c(2, 3))
dev.off()
for (i in 4:9) {
#  print(i)
  ggplot(df,aes(y=colnames(df)[i],x=type,fill=type)) + 
    geom_boxplot() +
    theme(legend.position = "none")
  
}


#rname	Reference name / chromosome
#startpos	Start position
#endpos	End position (or sequence length)
#numreads	Number reads aligned to the region (after filtering)
#covbases	Number of covered bases with depth >= 1
#coverage	Percentage of covered bases [0..100]
#meandepth	Mean depth of coverage
#meanbaseq	Mean baseQ in covered region
#meanmapq	Mean mapQ of selected reads

g1 <- ggplot(df,aes(y=numreads,x=type,fill=type)) + 
  geom_boxplot() +
  xlab("Number Reads aligned") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none")


g2<-ggplot(df,aes(y=meandepth,x=type,fill=type)) + 
  xlab("Mean depth of coverage") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none")+
  geom_boxplot()

g3<-ggplot(df,aes(y=meanmapq,x=type,fill=type)) + 
  xlab("Mean mapQ of selected reads") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none")+
  geom_boxplot()
g4 <- cowplot::get_legend(ggplot(df,aes(y=meandepth,x=type,fill=type)) + 
                            geom_boxplot())


legend <- get_legend(ggplot(df,aes(y=meandepth,x=type,fill=type)) + 
                       geom_boxplot() +  theme(legend.position = "top") + 
                       theme(legend.key = element_rect(fill = "white")))


grid.draw(g4)
head(df)
library(grid)
library(gridExtra)
grid.arrange(g1,g2, g3,g4,ncol = 4)
grid.arrange(arrangeGrob(g1,g2, g3, ncol=3),g4,ncol = 2)
grid.arrange(g1, g2, g3,legend, ncol=3, nrow = 2, 
             layout_matrix = rbind(c(1,2,3), c(4,4,4)),
             widths = c(2.7, 2.7, 2.7), heights = c(2, 0.3))


g1 <- ggplot(df,aes(y=numreads,x=type,fill=type)) + 
  geom_boxplot() +
  xlab("Number Reads aligned") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none")


g2<-ggplot(df,aes(y=meandepth,x=type,fill=type)) + 
  xlab("Mean depth of coverage") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none")+
  geom_boxplot()

g3<-ggplot(df,aes(y=covbases,x=type,fill=type)) + 
  xlab("Number of covered bases with depth") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none")+
  geom_boxplot()

g4<-ggplot(df,aes(y=coverage,x=type,fill=type)) + 
  xlab("Percentage of covered bases") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none")+
  geom_boxplot()

g5<-ggplot(df,aes(y=meanbaseq,x=type,fill=type)) + 
  xlab("Mean baseQ in covered region") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none")+
  geom_boxplot()

g6<-ggplot(df,aes(y=meanmapq,x=type,fill=type)) + 
  xlab("Mean mapQ of selected reads") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none")+
  geom_boxplot()


get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(ggplot(df,aes(y=meandepth,x=type,fill=type)) + 
                       geom_boxplot() +  theme(legend.position = "top") + 
                       theme(legend.title = element_blank()) + 
                       theme(legend.key = element_rect(fill = "white")))

blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()
grid.draw(g4)
head(df)
library(grid)
library(gridExtra)
library(ggplot2)


grid.arrange(g1, g2, g3,g4,g5,g6,legend, ncol=3, nrow = 3, 
             top = textGrob("Short-Read Bam statistics"),
             layout_matrix = rbind(c(1,2,3), c(4,5 ,6), c(7,7,7)),
             widths = c(2.7, 2.7,2.7), heights = c(2.5, 2.5,0.3))



#######
ggbarplot(gg_data, x = "count", y = "sample", fill = "len") +
  facet_grid(~type) + 
  theme(axis.text.x=element_text(size=7,angle=90), 
        axis.text.y=element_text(size=10), 
        panel.background=element_rect(fill="white",colour="Dark Blue", linetype="solid",size=1),
        panel.grid.major=element_line(color="grey"),
        axis.line.x=element_line(size=1,colour="Dark Blue",linetype="solid"),
        axis.line.y=element_line(size=1,colour="Dark Blue",linetype="solid"),
        axis.ticks=element_line(size=0.1,colour="white"),
        axis.title=element_text(size=14,colour="Black",face="bold"),
        title=element_text(size=14,colour="Black",face="bold"), 
        legend.position="top",strip.text.x=element_text(size=14,colour="Black",face="bold"),
        strip.text.y=element_text(size=14,colour="Black",face="bold"),
        axis.text=element_text(colour="Brown",face="bold",size=17),
        axis.text.x.bottom=element_text(colour="Brown",face="bold",size=12))


####3 trimmied

#setwd("~/Desktop/KCDC/HLA_seq/02.bam.stats/short-read/trimmed/")
ls_trimmed <- system("ls ~/Desktop/KCDC/HLA_seq/02.bam.stats/short-read/trimmed/ | grep bam.coverage", intern = T)
#setwd("~/Desktop/KCDC/HLA_seq/02.bam.stats/short-read/trimmed.dedup/")
ls_trimmed_dedup <- system("ls ~/Desktop/KCDC/HLA_seq/02.bam.stats/short-read/trimmed.dedup/ | grep bam.coverage", intern = T)

setwd("~/Desktop/KCDC/HLA_seq/02.bam.stats/short-read/trimmed/")
ls_trimmed <- system("ls ~/Desktop/KCDC/HLA_seq/02.bam.stats/short-read/trimmed/ | grep bam.coverage", intern = T)
df_trimmed = data.frame()
header = c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")
ls_trimmed
for (i in ls_trimmed) {
  tmp <- read.table(paste0(i))
  sample_id <- unlist(strsplit(i, split = "\\.", perl=T))[4]
  colnames(tmp)<-header
  tmp$ID <- sample_id
  df_trimmed <- rbind(df_trimmed,tmp)
}
tmp
head(df_trimmed)

setwd("~/Desktop/KCDC/HLA_seq/02.bam.stats/short-read/trimmed.dedup/")
ls_trimmed_dedup <- system("ls ~/Desktop/KCDC/HLA_seq/02.bam.stats/short-read/trimmed.dedup/ | grep bam.coverage", intern = T)
df_trimmed_dedup = data.frame()
header = c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")
for (i in ls_trimmed_dedup) {
  tmp <- read.table(paste0(i))
  sample_id <- unlist(strsplit(i, split = "\\.", perl=T))[4]
  colnames(tmp)<-header
  tmp$ID <- sample_id
  df_trimmed_dedup <- rbind(df_trimmed_dedup,tmp)
}

head(df_trimmed_dedup)



df_dedup$type <- 'RmDup'
df_nodup$type <- "Original"
df_trimmed$type <- "Trimmed"
df_trimmed_dedup$type <- "Trimmed_RmDup"
df <- rbind(df_dedup,df_nodup)
df <- rbind(df,df_trimmed)
df <- rbind(df,df_trimmed_dedup)

head(df)
par(mfrow = c(2, 3))
dev.off()
for (i in 4:9) {
  #  print(i)
  ggplot(df,aes(y=colnames(df)[i],x=type,fill=type)) + 
    geom_boxplot() +
    theme(legend.position = "none")
  
}
head(df)
df %>% group_by(type) %>%
  pivot_longer(colnames(df)[4:9],names_to = "category",values_to = "Value") %>% #head()
  ggplot(aes(y=Value,x=type,fill=type))+
  geom_boxplot() + 
  facet_wrap(~category)
#rname	Reference name / chromosome
#startpos	Start position
#endpos	End position (or sequence length)
#numreads	Number reads aligned to the region (after filtering)
#covbases	Number of covered bases with depth >= 1
#coverage	Percentage of covered bases [0..100]
#meandepth	Mean depth of coverage
#meanbaseq	Mean baseQ in covered region
#meanmapq	Mean mapQ of selected reads
library(grid)
library(gridExtra)


g1 <- ggplot(df,aes(y=numreads,x=type,fill=type)) + 
  geom_boxplot() +
  xlab("Number Reads aligned") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none")


g2<-ggplot(df,aes(y=meandepth,x=type,fill=type)) + 
  xlab("Mean depth of coverage") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none")+
  geom_boxplot()

g3<-ggplot(df,aes(y=covbases,x=type,fill=type)) + 
  xlab("Number of covered bases with depth") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none")+
  geom_boxplot()

g4<-ggplot(df,aes(y=coverage,x=type,fill=type)) + 
  xlab("Percentage of covered bases") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none")+
  geom_boxplot()

g5<-ggplot(df,aes(y=meanbaseq,x=type,fill=type)) + 
  xlab("Mean baseQ in covered region") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none")+
  geom_boxplot()

g6<-ggplot(df,aes(y=meanmapq,x=type,fill=type)) + 
  xlab("Mean mapQ of selected reads") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none")+
  geom_boxplot()


get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(ggplot(df,aes(y=meandepth,x=type,fill=type)) + 
                       geom_boxplot() +  theme(legend.position = "top") + 
                       theme(legend.title = element_blank()) + 
                       theme(legend.key = element_rect(fill = "white")))

blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

head(df)
library(grid)
library(gridExtra)
library(ggplot2)


grid.arrange(g1, g2, g3,g4,g5,g6,legend, ncol=3, nrow = 3, 
             top = textGrob("Short-Read Bam statistics"),
             layout_matrix = rbind(c(1,2,3), c(4,5 ,6), c(7,7,7)),
             widths = c(2.7, 2.7,2.7), heights = c(2.5, 2.5,0.3))


