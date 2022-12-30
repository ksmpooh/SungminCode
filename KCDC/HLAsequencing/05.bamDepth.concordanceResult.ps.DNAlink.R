#### concordance vs Depth
library(readxl)
library(tidyverse)
library(stringr)
library(grid)
library(gridExtra)
library(cowplot)
library(ggplot2)
library(ggbreak)
library(ggpmisc)
library(ggpubr)
library(patchwork) # To display 2 charts together
library(hrbrthemes)

setwd("~/Desktop/KCDC/HLA_seq/DNAlink/")
long <- read_excel("long_Exon.xlsx",skip = 1)
short <- read_excel("short_Exon.xlsx",skip = 1)
head(long)
head(short)
head(target)
target_region <- long %>% select(startpos,endpos) %>% unique()
dim(target_region)
table(long$startpos %in% short$startpos)
table(long$endpos %in% short$endpos)

#data <- read.table("/Volumes/DATA/HLA_seq/05.concordance/20221206_DV.VQSR.hard/concordance_result_3/QulityMetrix_ShortGATKhardfilter.LongDV.txt",header = T)
head(data)
df <- data %>% #inner_join(maf1) %>% 
  select(pos,Accuracy) #%>% 
#  pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val")# %>% head()
head(df)  
#"chr4:100,001-100,001",
df %>% mutate(hg19 = paste0("chr6:",pos,"-",pos)) %>% select(hg19)-> df1

#write.table(df1,"hg19.txt",row.names = F,quote = F,sep = "\t")
hg38 <- read.table("hglft_genome_267e2_11bd0.bed")
head(hg38)
dim(df)
dim(hg38)
head(target_region)

df1 <- df %>% cbind(hg38) %>% mutate(hg38 = str_split_fixed(V1,"-",2)[,2]) %>% select(hg38,Accuracy) %>% mutate(hg38 = as.double(hg38))
#df %>% 
df1 %>% filter(hg38 >= target_region$startpos && hg38 <= target_region$endpos)
head(df1)
#write.table(df1,"hg38_accuracy.txt",col.names =T,row.names = F, quote = F)
#write.table(target_region,"hg38_target.txt",col.names =T,row.names = F, quote = F)
str(df1)
str(long)
head(long)
df1[df1$hg38 >29942531 && df1$hg38 < 29942626,]

df%>% filter(hg38 >= target_region$startpos) %>% filter(hg38 <= target_region$endpos)
df1%>% filter(hg38 >= target_region$startpos) %>% filter(hg38 <= target_region$endpos)
head(df1)
head(long)


#ref <- read.table("hg38_target_accurac_regionSUM.txt",header = T)
ref <- read.table("hg38_target_accurac_regionMean.txt",header = T)
head(ref)

head(long)
short %>% select(startpos,endpos,meandepth,Gene,Exon) %>%
  group_by(Gene,startpos,endpos,Exon) %>%
  summarise(meandepth=mean(meandepth)) %>% mutate(type = "Short") %>% 
  inner_join(ref) -> short1

short %>% select(startpos,endpos,meandepth,Gene,Exon) %>%
  group_by(Gene,startpos,endpos,Exon) %>%
  summarise(meandepth=mean(meandepth)) %>% mutate(type = "Short") %>% #head()
  inner_join(ref) -> short1


long %>% select(startpos,endpos,meandepth,Gene,Exon) %>%
  group_by(Gene,startpos,endpos,Exon) %>%
  summarise(meandepth=mean(meandepth)) %>% mutate(type = "Long") %>%
  inner_join(ref) %>% 
  filter(Gene %in% c("HLA-B","HLA-C")) %>%
#  group_by(Gene) %>%
  ggplot(aes(x=meandepth,y=accuracy_sum,color= Gene)) +
  geom_line()


long %>% select(startpos,endpos,meandepth,Gene,Exon) %>%
  group_by(Gene,startpos,endpos,Exon) %>%
  summarise(meandepth=mean(meandepth)) %>% mutate(type = "Long") %>%
  inner_join(ref) -> long2


long %>% select(startpos,endpos,meandepth,Gene,Exon) %>%
  group_by(Gene,startpos,endpos,Exon) %>%
  summarise(meandepth=mean(meandepth)) %>% #head()
  mutate(type = "Long") %>%
  inner_join(ref) %>% rbind(short1) %>%
  filter(Gene %in% c("HLA-A","HLA-B","HLA-C","HLA-DRB1")) %>%
  ggplot(aes(x=meandepth,y=accuracy_sum,color=Gene)) +
  geom_point() + geom_line() + 
  facet_grid(~type)

head(long2)
long2 %>% 
  filter(Gene %in% c("HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1")) %>%
  ggplot(aes(x=meandepth,y=accuracy_sum,color=Gene)) +
  geom_point() + geom_line() + 
  geom_text(aes(x=meandepth, y=accuracy_sum, group=Exon, label=Exon), size = 3, hjust=0, vjust=-1) + 
  xlab(element_text("Mapping Depth (Mean)")) + ylab(element_text("Concordance Test Accuracy (Short vs Long)")) + 
  labs(title = "Long-read") + 
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none") ->a1

short1 %>% 
  filter(Gene %in% c("HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1")) %>%
  ggplot(aes(x=meandepth,y=accuracy_sum,color=Gene)) +
  geom_point() + geom_line() + 
  geom_text(aes(x=meandepth, y=accuracy_sum, group=Exon, label=Exon), size = 3, hjust=0, vjust=-1) + 
  xlab(element_text("Mapping Depth (Mean)")) + ylab(element_text("Concordance Test Accuracy (Short vs Long)")) + 
  labs(title = "Short-read") + 
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="right") ->a2

long2 %>% 
  filter(Gene %in% c("HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1")) %>%
  ggplot(aes(x=meandepth,y=accuracy_sum,color=Gene)) +
  geom_point() + geom_line() + 
  geom_text(aes(x=meandepth, y=accuracy_sum, group=Exon, label=Exon), size = 3, hjust=0, vjust=-1) + 
  xlab(element_text("Mapping Depth (Mean)")) + ylab(element_text("Concordance Test Accuracy (Short vs Long)")) + 
  labs(title = "Long-read") + 
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="bottom") ->a3
a3 <- get_legend(a3)

ggarrange(a1,a2,ncol = 2,nrow = 1,widths = c(5,6))
a1
a2
a3
grid.arrange(pA, pB, pC,pDRB,pDQA,pDQB,pDPA,pDPB,   widths = c(1,1,2),#labels = "AUTO",
             layout_matrix = rbind(c(1,1,1), c(2,2,2),c(3,3,3),c(4,4,4),c(5,5,6),c(7,8,8)))


grid.arrange(a1, a2, a3,  widths = c(5,5),
             layout_matrix = rbind(c(1,2), c(3,3)),)

head(ref)
long2 %>% 
  filter(Gene %in% c("HLA-B","HLA-C")) %>%
  ggplot(aes(x=meandepth,y=accuracy_mean,color=Gene)) +
  geom_point() + geom_line() + 
  geom_text(aes(x=meandepth, y=accuracy_mean, group=Exon, label=Exon), size = 5, hjust=1, vjust=1) + 
  xlab(element_text("Mapping Depth (Mean)")) + ylab(element_text("Concordance Test Accuracy (Short vs Long)")) + 
  labs(title = "Long-read") + 
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none") #->a1

short1 %>% 
  filter(Gene %in% c("HLA-B","HLA-C")) %>%
  ggplot(aes(x=meandepth,y=accuracy_mean,color=Gene)) +
  geom_point() + geom_line() + 
  geom_text(aes(x=meandepth, y=accuracy_mean, group=Exon, label=Exon), size = 5, hjust=1, vjust=1) + 
  xlab(element_text("Mapping Depth (Mean)")) + ylab(element_text("Concordance Test Accuracy (Short vs Long)")) + 
  labs(title = "Short-read") + 
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="right") ->a2

ggarrange(a1,a2,ncol = 2,nrow = 1,widths = c(5,6))


head(long2)
head(short1)
long2$type = "Long-read"
short1$type = "Short-read"
long2 %>% group_by(Gene) %>% summarise(n=sum(count))


long2 %>% filter(Gene=="HLA-B") %>%
  ggplot() + 
  geom_point(aes(x=meandepth,y=accuracy_mean,color=Gene)) + 
  geom_line(aes(x=meandepth,y=accuracy_mean,color=Gene)) + 
  geom_text(aes(x=meandepth, y=accuracy_mean, group=Exon, label=Exon), size = 5, hjust=1, vjust=1) +
  geom_point(aes(x=meandepth,y=count/20,color="# of variant")) + 
  scale_y_continuous(name = "Mean of concordance accuracy",sec.axis = sec_axis(~ .*20,name = '# of variant')) + 
  theme(legend.title = element_blank(),legend.position="bottom") + 
  scale_color_discrete(labels=c('# of variant', 'Accuracy'))

  
short1 %>% filter(Gene=="HLA-B") %>%
  ggplot() + 
  geom_point(aes(x=meandepth,y=accuracy_mean,color=Gene)) + 
  geom_line(aes(x=meandepth,y=accuracy_mean,color=Gene)) + 
  geom_text(aes(x=meandepth, y=accuracy_mean, group=Exon, label=Exon), size = 5, hjust=1, vjust=1) +
  geom_point(aes(x=meandepth,y=count/20,color="# of variant")) + 
  scale_y_continuous(name = "Mean of concordance accuracy",sec.axis = sec_axis(~ .*20,name = '# of variant')) + 
  theme(legend.title = element_blank(),legend.position="bottom") + 
  scale_color_discrete(labels=c('# of variant', 'Accuracy'))
  

  geom_bar(data = df1 %>% filter(!grepl("TS",type)),
           aes(x=factor(Method,levels=c("DV_unfiltered","DV_filtered","GATK_unfiltered","GATK_VSQR","GATK_Hardfiltering")),
               y=Val,fill = factor(type,levels=c("Multiallelic SNP site","Multiallelic site","INDEL","SNP"))),stat = "identity") + 
  geom_point(data = df1 %>% filter(grepl("TS",type)),aes(x=Method,y=Val*200000,color=type)) + 
  #scale_fill_manual(values=c("TS/TV"="#002955" ,"TS/TV(1st ALT)"="#074ca1")) +
  #geom_point(data = df1 %>% filter(grepl("TS",type)),aes(x=Method,y=Val*200000,color=c(rgb(0,0,0,0),rgb(0,0,1,0.5)))) + 
  scale_y_continuous(name = 'Count (Bar)',sec.axis = sec_axis(~ ./200000,name = 'TS/TV ratio (Point)'))+
  #scale_color_discrete_qualitative(palette = "Cold") + 
  #theme_ipsum() + 
  theme(legend.title = element_blank(),legend.position="bottom") +
  xlab(element_blank()) + ylab(element_text("# of variant"))+
  facet_grid(~Seqeuncing) #->a1


  
  

######################

setwd("~/Desktop/KCDC/HLA_seq/DNAlink/")
long <- read_excel("Long_Gene.xlsx",skip = 1)
short <- read_excel("Short_Gene.xlsx",skip = 1)

long <- read_excel("long_Exon.xlsx",skip = 1)
short <- read_excel("short_Exon.xlsx",skip = 1)
head(long)
head(short)


long %>% mutate(type = "Long") %>% #summarise(mean = mean(meandepth))
  filter(Gene %in% c("HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1")) %>% #head()
  #na.omit() %>%
  summarise(mean = mean(meandepth))

long %>% mutate(type = "Long") %>%
  filter(Gene %in% c("HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1")) %>%
  select(ID,Gene,meandepth,type) %>%
  ggplot(aes(x=Gene,y=meandepth,fill=Gene)) + 
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none") + 
  geom_boxplot() + 
  geom_hline(yintercept=62.4, linetype='dashed', color='black', size=1)


short %>% #mutate(type = "Long") %>% #summarise(mean = mean(meandepth))
  filter(Gene %in% c("HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1")) %>% #head()
  #na.omit() %>%
  summarise(mean = mean(meandepth))

short %>% mutate(type = "Short") %>%
  filter(Gene %in% c("HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1")) %>%
  select(ID,Gene,meandepth,type) %>%
  ggplot(aes(x=Gene,y=meandepth,fill=Gene)) + 
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none") + 
  geom_boxplot() + 
  geom_hline(yintercept=233, linetype='dashed', color='black', size=1)


long <- read_excel("long_Exon.xlsx",skip = 1)
short <- read_excel("short_Exon.xlsx",skip = 1)
head(long)
head(short)


short %>% filter(Gene %in% c("HLA-B")) %>% mutate(Exon = as.factor(Exon)) %>% summarise(n = mean(meandepth))
short %>% filter(Gene %in% c("HLA-B")) %>% mutate(Exon = as.factor(Exon)) %>%
  select(ID,meandepth,Exon) %>% #dim()
  group_by(Exon) %>% #head()
  mutate(meandepth = meandepth/226) %>% #head()
  ggplot(aes(x=Exon,y=meandepth,fill=Exon)) + 
  geom_boxplot() +
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none") +
  ylab(element_text("(Mean depth of an exon/mean Depth)*100"))





long %>% filter(Gene %in% c("HLA-B")) %>% mutate(Exon = as.factor(Exon)) %>% summarise(n = mean(meandepth))
long %>% filter(Gene %in% c("HLA-B")) %>% mutate(Exon = as.factor(Exon)) %>%
  select(ID,meandepth,Exon) %>% #dim()
  group_by(Exon) %>% #head()
  mutate(meandepth = meandepth/5.01) %>% #head()
  ggplot(aes(x=Exon,y=meandepth,fill=Exon)) + 
  geom_boxplot() +
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none") +
  ylab(element_text("(Mean depth of an exon/mean Depth)*100"))


short %>% filter(Gene %in% c("HLA-C")) %>% mutate(Exon = as.factor(Exon)) %>% summarise(n = mean(meandepth))
short %>% filter(Gene %in% c("HLA-C")) %>% mutate(Exon = as.factor(Exon)) %>%
  select(ID,meandepth,Exon) %>% #dim()
  group_by(Exon) %>% #head()
  mutate(meandepth = meandepth/129) %>% #head()
  ggplot(aes(x=Exon,y=meandepth,fill=Exon)) + 
  geom_boxplot() +
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none") +
  ylab(element_text("(Mean depth of an exon/mean Depth)*100"))


long %>% filter(Gene %in% c("HLA-C")) %>% mutate(Exon = as.factor(Exon)) %>% summarise(n = mean(meandepth))
long %>% filter(Gene %in% c("HLA-C")) %>% mutate(Exon = as.factor(Exon)) %>%
  select(ID,meandepth,Exon) %>% #dim()
  group_by(Exon) %>% #head()
  mutate(meandepth = meandepth/19) %>% #head()
  ggplot(aes(x=Exon,y=meandepth,fill=Exon)) + 
  geom_boxplot() +
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none") +
  ylab(element_text("(Mean depth of an exon/mean Depth)*100"))



long %>% mutate(type = "Long") %>% rbind(short %>% mutate(type = "Short")) -> df



df %>% filter(Gene %in% c("HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1")) %>%
  select(ID,Gene,meandepth,type) %>% #head()
  ggplot(aes(x=Gene,y=meandepth,fill=type)) + 
  geom_boxplot() +
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="bottom")

table(df$Exon)
df %>% filter(Gene %in% c("HLA-B")) %>% mutate(Exon = as.factor(Exon)) %>%
  select(ID,meandepth,type,Exon) %>%
  ggplot(aes(x=Exon,y=meandepth,fill=type)) + 
  geom_boxplot() +
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="bottom")

df %>% filter(type == "Long",Gene %in% c("HLA-B")) %>% mutate(Exon = as.factor(Exon)) %>%
  select(ID,meandepth,type,Exon) %>%
  ggplot(aes(x=Exon,y=meandepth,fill=type)) + 
  geom_boxplot() +
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none")


df %>% filter(Gene %in% c("HLA-C")) %>% mutate(Exon = as.factor(Exon)) %>%
  select(ID,meandepth,type,Exon) %>%
  ggplot(aes(x=Exon,y=meandepth,fill=type)) + 
  geom_boxplot() +
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="bottom")

df %>% filter(type == "Long",Gene %in% c("HLA-C")) %>% mutate(Exon = as.factor(Exon)) %>%
  select(ID,meandepth,type,Exon) %>%
  ggplot(aes(x=Exon,y=meandepth,fill=type)) + 
  geom_boxplot() +
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none")
