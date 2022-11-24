## long read bam box plot

#HLA.Longread.Seq.NIH19KT0745_mapped.Q20.bam.flagstats HLA.Longread.Seq.NIH19KT3809_mapped.Q20.bam.flagstats
#HLA.Longread.Seq.NIH19KT0745_mapped.Q20.bam.stats     HLA.Longread.Seq.NIH19KT3809_mapped.Q20.bam.stats
#HLA.Longread.Seq.NIH19KT0745_mapped.bam.coverage
setwd("~/Desktop/KCDC/long_read/2022/long-read/01.bam.stats/")
library(tidyverse)



ls_Q20 <- system("ls | grep Q20.bam.coverage", intern = T)
ls_noQ <- system("ls | grep mapped.bam.coverage", intern = T)
head(ls_Q20)

sample_id <- str_replace(unlist(strsplit(ls_Q20[1], split = "\\.", perl=T))[4],"_mapped","")
sample_id
#file_suffix <- ".subreads.gt-cl.ext.stats"
file_suffix <- ".merge.stats"
###
#unlist(strsplit(i, split = "\\.", perl=T))[4]
df_Q20 = data.frame()
header = c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")
for (i in ls_Q20) {
  tmp <- read.table(paste0(i))
  #sample_id <- str_replace(unlist(strsplit(i, split = "\\.", perl=T))[4],"_mapped","")
  colnames(tmp)<-header
  #tmp$ID <- sample_id
  tmp$ID <- str_replace(i,".coverage","")
  df_Q20 <- rbind(df_Q20,tmp)
}


head(df_Q20)

df_noQ = data.frame()
header = c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")
for (i in ls_noQ) {
  tmp <- read.table(paste0(i))
  #sample_id <- str_replace(unlist(strsplit(i, split = "\\.", perl=T))[4],"_mapped","")
  colnames(tmp)<-header
  #tmp$ID <- sample_id
  tmp$ID <- str_replace(i,".coverage","")
  df_noQ <- rbind(df_noQ,tmp)
}






df_Q20$type <- 'FilterQ20'
df_noQ$type <- "BeforeFilter"

df <- rbind(df_Q20,df_noQ)
head(df)

av_length <- read.table("../01.bam.average.length.txt")
av_length<-av_length %>% mutate(ID = str_replace(V1,".stats:","")) %>% #head()
  rename("mean_readLength" = V2) %>% select(ID,mean_readLength)

head(av_length)
head(df)
df <- df %>% inner_join(av_length)

ggplot(df,aes(x=mean_readLength,fill= type)) +
  geom_histogram( bins = 30) + 
  xlab("Mean Read Length") + 
  #theme(axis.title.y = element_blank()) +
  #theme(axis.text.x = element_blank()) + 
  theme(legend.position = "bottom")
  #geom_boxplot()



#header = c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")
ggplot(df,aes(y=numreads,x=type,fill=type)) + 
  geom_boxplot()
ggplot(df,aes(y=covbases,x=type,fill=type)) + 
  geom_boxplot()
ggplot(df,aes(y=coverage,x=type,fill=type)) + 
  geom_boxplot()
ggplot(df,aes(y=meandepth,x=type,fill=type)) + 
  geom_boxplot()
ggplot(df,aes(y=meanbaseq,x=type,fill=type)) + 
  geom_boxplot()
ggplot(df,aes(y=meanmapq,x=type,fill=type)) + 
  geom_boxplot()


ggplot(df,aes(y=colnames(df)[1],x=type,fill=type)) + 
  geom_boxplot()


ggplot(df,aes(y=colnames(df)[1],x=type,fill=type)) + 
  geom_boxplot()
colnames(df)

#[1] "rname"           "startpos"        "endpos"          "numreads"        "covbases"        "coverage"        "meandepth"      
#[8] "meanbaseq"       "meanmapq"        "ID"              "type"            "mean_readLength"


#rname	Reference name / chromosome
#startpos	Start position
#endpos	End position (or sequence length)
#numreads	Number reads aligned to the region (after filtering)
#covbases	Number of covered bases with depth >= 1
#coverage	Percentage of covered bases [0..100]
#meandepth	Mean depth of coverage
#meanbaseq	Mean baseQ in covered region
#meanmapq	Mean mapQ of selected reads


dev.off()
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
             top = textGrob("Long-Read Bam statistics"),
             layout_matrix = rbind(c(1,2,3), c(4,5 ,6), c(7,7,7)),
             widths = c(2.7, 2.7,2.7), heights = c(2.5, 2.5,0.3))





grid.arrange(g1, g2, g3,g4, ncol=3, nrow = 2, 
             layout_matrix = rbind(c(1,2,3), c(4,4,4)),
             widths = c(2.7, 2.7, 2.7), heights = c(2.5, 0.2))
layout_matrix

grid.arrange(g1,g2, g3,g4,ncol = 3,nrow=2,layout_matrix = rbind(c(1,3), c(4,4)))

layout_matrix = rbind(c(1,2), c(3,3))
layout_matrix
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


###### trimmed

setwd("~/Desktop/KCDC/long_read/2022/long-read/01.bam.stats/")
library(tidyverse)



ls_Q20 <- system("ls | grep Q20.bam.coverage", intern = T)
ls_noQ <- system("ls | grep mapped.bam.coverage", intern = T)
head(ls_Q20)

df_Q20 = data.frame()
header = c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")
for (i in ls_Q20) {
  tmp <- read.table(paste0(i))
  colnames(tmp)<-header
  tmp$ID <- str_replace(i,".coverage","")
  df_Q20 <- rbind(df_Q20,tmp)
}


head(df_Q20)

df_noQ = data.frame()
header = c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")
for (i in ls_noQ) {
  tmp <- read.table(paste0(i))
  colnames(tmp)<-header
  tmp$ID <- str_replace(i,".coverage","")
  df_noQ <- rbind(df_noQ,tmp)
}


setwd("~/Desktop/KCDC/HLA_seq/02.bam.stats/long-read/trimmed/")
ls_trimmed <- system("ls | grep bam.coverage", intern = T)
df = data.frame()
header = c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")
for (i in ls_trimmed) {
  tmp <- read.table(paste0(i))
  colnames(tmp)<-header
  tmp$ID <- str_replace(i,".coverage","")
  df <- rbind(df,tmp)
}
df_trimmed <- df

setwd("~/Desktop/KCDC/HLA_seq/02.bam.stats/long-read/trimmed.Q20/")
ls_trimmed_Q20 <- system("ls | grep bam.coverage", intern = T)
df = data.frame()
header = c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")
for (i in ls_trimmed_Q20) {
  tmp <- read.table(paste0(i))
  colnames(tmp)<-header
  tmp$ID <- str_replace(i,".coverage","")
  df <- rbind(df,tmp)
}
df_trimmed_Q20 <- df


df_Q20$type <- 'Q20'
df_noQ$type <- "Original"
df_trimmed$type <- 'Trimmed'
df_trimmed_Q20$type <- "Trimmed_Q20"



df <- rbind(df_noQ,df_Q20)
df <- rbind(df,df_trimmed)
df <- rbind(df,df_trimmed_Q20)
head(df)


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
             top = textGrob("Long-Read Bam statistics"),
             layout_matrix = rbind(c(1,2,3), c(4,5 ,6), c(7,7,7)),
             widths = c(2.7, 2.7,2.7), heights = c(2.5, 2.5,0.3))
head(df)
a<-df %>% 
  pivot_longer(cols = colnames(df)[4:9],values_to = "Value",names_to = "category") %>%
  group_by(type,category) %>%
  summarise(mean = mean(Value)) %>%
  pivot_wider(names_from = category,values_from = mean)
a  
df %>% group_by(type)
