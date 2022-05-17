
setwd("~/Desktop/KCDC/long_read/2022/short-read/01.bam.stats/")
library(tidyverse)
ls_dedup <- system("ls | grep dedup.bam.coverage", intern = T)
ls_nodup <- system("ls | grep sorted.bam.coverage", intern = T)
head(ls_nodup)

sample_id <- unlist(strsplit(ls_dedup[1], split = "\\.", perl=T))[4]
sample_id
#file_suffix <- ".subreads.gt-cl.ext.stats"
file_suffix <- ".merge.stats"
###
data1 <- read.table(paste0(sample_id, file_suffix), sep=" ", header=F)
data1$V3 <- abs(data1$V3)
data1 <- data1[order(data1$V1, data1$V5),]
data1$V7 <- sample_id
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
g1 <- ggplot(df,aes(y=numreads,x=type,fill=type)) + 
  geom_boxplot() +
  xlab("Number Reads aligned") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = "none")


g2<-ggplot(df,aes(y=meandepth,x=type,fill=type)) + 
  xlab("Mean depth") + 
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
grid.draw(g4)
head(df)
library(grid)
library(gridExtra)
grid.arrange(g1,g2, g3,g4,ncol = 4)
grid.arrange(arrangeGrob(g1,g2, g3, ncol=3),g4,ncol = 2)

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


