## HLA type haplotype frequency plot

library(tidyverse)
library(stringr)
library(alluvial)
library(ggalluvial)
library(RColorBrewer)


setwd("~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/")
df <- readxl::read_xlsx("HLA.typing.Final.result_529sample.xlsx")
dim(df)
ref <- read.table("../MakeReferencePanel/Final_520sample.index.txt",header = T)
head(ref)
dim(ref)
head(df)

data(majors)
head(majors)
majors %>% count(semester)
majors %>% count(curriculum)
majors %>% count(student)
'''
id	A	A	B	B	C	C	DPB1	DPB1	DQB1	DQB1	DRB1	DRB1
1	26:04	29:67	54:16	07:218	07:37	01:113	294:01	04:02:04	03:220	03:208	11:182	13:14:01
2	24:DXTW	02:570	14:DYBE	13:07	01:02:17	02:58	459:01	472:01	03:JHRJ	03:01:29	12:CVT	14:23:04
3	29:67	02:570	35:54	13:07	07:02:47	02:58	45:01	472:01	02:14	03:01:29	13:121	14:23:04
4	66:10	02:77	15:154	14:DYBE	12:130	16:02:12	526:01	20:01	05:101	06:158	16:30	16:30
5	29:36	03:217	35:54	15:01:30	14:02:09	05:52	09:EMP	88:01	03:82	06:158	14:23:04	11:95
6	03:93	02:454	35:54	15:316	12:13	08:67	288:01	504:01	06:01:04	03:01:26	04:52	08:18
7	30:13	24:02:55	39:27	18:19	12:30	01:113	504:01	205:01	03:01:17	03:122	14:23:04	13:116
8	30:13	03:93	18:12	18:50	01:11	01:11	26:01	456:01	05:52	03:04P	11:147	15:DUUW
9	68:94N	02:259	15:154	57:21	03:02:04	14:71	294:01	02:01:07	06:54N	03:220	13:116	14:73
'''
df %>% filter(ID %in% ref$KBAID) %>% #dim()
  select(-Sample) %>% pivot_longer(2:17) %>% mutate(value = paste0(str_split_fixed(value,":",4)[,1],":",str_split_fixed(value,":",4)[,2])) %>%
  mutate(value = str_split_fixed(value,"\\*",2)[,2]) %>% pivot_wider(names_from = "name",values_from = "value") -> df
colnames(df)
colnames(df) <- c("id","A","A","B","B","C","C","DQA1","DQA1","DQB1","DQB1","DPA1","DPA1","DPB1","DPB1","DRB1","DRB1")
head(df)  
df1_HLA <- df
######## ohter ref
df <-read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLAtype_check/other_panel/noman/Han.hg19.haplegendtovcf.modify.hlatype_fornomenclean_2field.chped")
head(df)

df %>% select(-V2,-V3,-V4,-V5,-V6) %>% pivot_longer(2:17) %>%
  mutate(value = paste0(str_split_fixed(value,":",4)[,1],":",str_split_fixed(value,":",4)[,2])) %>% #head()
  mutate(value = str_split_fixed(value,"\\*",2)[,2]) %>% pivot_wider(names_from = "name",values_from = "value") -> df
head(df)
#colnames(df) <- c("id","A","A","B","B","C","C","DQA1","DQA1","DQB1","DQB1","DPA1","DPA1","DPB1","DPB1","DRB1","DRB1")
colnames(df) <- c("id","A","A","B","B","C","C","DPA1","DPA1","DPB1","DPB1","DQA1","DQA1","DQB1","DQB1","DRB1","DRB1")
#write.table(df,"~/Tool/Hapl-o-Mat/examplePopulations/a/HanChinese_HLA.2field.HLAtype.forHapl-m-MAT.txt",col.names = T,row.names = F,quote = F,sep = "\t")


df <-read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLAtype_check/other_panel/noman/PanKor_merged.hg19.haplegendtovcf.modify.hlatype_fd.IDchange_fornomenclean_2field.chped")
head(df)

df %>% select(-V2,-V3,-V4,-V5,-V6) %>% pivot_longer(2:17) %>%
  mutate(value = paste0(str_split_fixed(value,":",4)[,1],":",str_split_fixed(value,":",4)[,2])) %>% #head()
  mutate(value = str_split_fixed(value,"\\*",2)[,2]) %>% pivot_wider(names_from = "name",values_from = "value") -> df
head(df)
#colnames(df) <- c("id","A","A","B","B","C","C","DQA1","DQA1","DQB1","DQB1","DPA1","DPA1","DPB1","DPB1","DRB1","DRB1")
colnames(df) <- c("id","A","A","B","B","C","C","DPA1","DPA1","DPB1","DPB1","DQA1","DQA1","DQB1","DQB1","DRB1","DRB1")
#write.table(df,"~/Tool/Hapl-o-Mat/examplePopulations/a/PanKor_HLA.2field.HLAtype.forHapl-m-MAT.txt",col.names = T,row.names = F,quote = F,sep = "\t")



###### 4 digit
#df <-read.table("~/Tool/Hapl-o-Mat/examplePopulations/a/run/genotypes.dat")
df <-read.table("~/Tool/Hapl-o-Mat/examplePopulations/a/run/hfs.dat")
df_eps <-read.table("~/Tool/Hapl-o-Mat/examplePopulations/a/run/epsilon.dat")
df_geno <-read.table("~/Tool/Hapl-o-Mat/examplePopulations/a/run/genotypes.dat")
head(df_eps)
head(df_geno)
table(df_geno$V3)
#df <-read.table("~/Tool/Hapl-o-Mat/examplePopulations/a/run/epsilon.dat")
head(df)

ref <- read.table("~/Desktop/KCDC/HLAimputation/HLA_type_frequency.txt",header = T)
#ref <- read.table("~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211")
head(ref)
ref %>% filter(digit == "4digit") %>% #mutate(check = ifelse(freq < 0.01,0,1)) %>% head()
  mutate(HLA_type = ifelse(freq < 0.01,"other",HLA_type)) -> ref1
###### 4 digit

ref %>% filter(digit == "4digit") %>% #mutate(check = ifelse(freq < 0.01,0,1)) %>% head()
  mutate(HLA_type = ifelse(freq < 0.01,"rare",HLA_type)) -> ref1



df %>% #head()
  mutate(count = ifelse(V2 > 0.005,"1","2")) %>% #head()#count(count) 
  separate( V1, c("A", "B", "C","DPA1","DPB1","DQA1","DQB1","DRB1"), sep = "~", remove = TRUE) %>%  #head()
  mutate(A = ifelse(A %in% ref1$HLA_type,A,'other')) %>% #count(A)
  mutate(B = ifelse(B %in% ref1$HLA_type,B,'other')) %>%
  mutate(C = ifelse(C %in% ref1$HLA_type,C,'other')) %>%
  mutate(DRB1 = ifelse(DRB1 %in% ref1$HLA_type,DRB1,'other')) %>%
  mutate(DPA1 = ifelse(DPA1 %in% ref1$HLA_type,DPA1,'other')) %>%
  mutate(DPB1 = ifelse(DPB1 %in% ref1$HLA_type,DPB1,'other')) %>%
  mutate(DQA1 = ifelse(DQA1 %in% ref1$HLA_type,DQA1,'other')) %>% 
  mutate(DQB1 = ifelse(DQB1 %in% ref1$HLA_type,DQB1,'other')) %>% #-> df1 #%>% #head()
  ggplot(aes(y=V2,
             axis1 = A,axis2=B,axis3 = C,axis4 = DRB1,axis5 = DQA1,axis6 = DQB1,axis7 = DPA1,axis8 = DPB1)) +
  #geom_alluvium(aes(fill=count),curve_type = "linear") +
  geom_alluvium(curve_type = "linear") + 
  scale_x_continuous(breaks = 1:8, labels = c("A", "B", "C","DRB1","DQA1","DQB1","DPA1","DPB1")) + 
  geom_stratum(width = 1/2) + 
  theme_minimal() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=3)
  #scale_fill_manual(values=c("#999999", "#56B4E9"))

c
head(df)
df %>% #head()
  #mutate(count = ifelse(V2 < 0.01,"0","1")) %>% #head()#count(count)
  mutate(count = rank(-V2)) %>% #head()#count(count) 
  mutate(count = ifelse(count < 6,count,"other")) %>%
  separate( V1, c("A", "B", "C","DPA1","DPB1","DQA1","DQB1","DRB1"), sep = "~", remove = TRUE) %>%  #head()
  mutate(A = ifelse(A %in% ref1$HLA_type,A,'other')) %>% #count(A)
  mutate(B = ifelse(B %in% ref1$HLA_type,B,'other')) %>%
  mutate(C = ifelse(C %in% ref1$HLA_type,C,'other')) %>%
  mutate(DRB1 = ifelse(DRB1 %in% ref1$HLA_type,DRB1,'other')) %>%
  mutate(DPA1 = ifelse(DPA1 %in% ref1$HLA_type,DPA1,'other')) %>%
  mutate(DPB1 = ifelse(DPB1 %in% ref1$HLA_type,DPB1,'other')) %>%
  mutate(DQA1 = ifelse(DQA1 %in% ref1$HLA_type,DQA1,'other')) %>% 
  mutate(DQB1 = ifelse(DQB1 %in% ref1$HLA_type,DQB1,'other')) %>% #-> df1 #%>% #head()
  mutate(subject = seq(1, n())) %>% #head()
  gather(key, value, -count , -subject, -V2) %>% #head()
  arrange(key, V2) %>% #head()
  #mutate(key = factor(key, levels = c("X1", "X2", "X3", "X4"))) %>% 
  ggplot(
  aes(x = factor(key,levels = c("A", "B", "C","DRB1","DQA1","DQB1","DPA1","DPB1")),
      y = V2,
      stratum = value, 
      alluvium = subject,
      label = value))+
  geom_flow(aes(fill = count),curve_type = "linear") +
  #geom_flow() +
  geom_stratum() +
  geom_text(stat = "stratum")+
  scale_fill_manual(values = c("red","green","blue","#d8b365","#5ab4ac","grey")) + 
  guides(fill=guide_legend(title="Top Rank of Haplotype",nrow = 1)) + 
  #scale_fill_manual(values=c("#BAB3B3EB", brewer.pal(5, "Set2"))) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.title=element_text(size=20),
        legend.text = element_text(size=15),
        axis.text.x = element_text(size=14),
        panel.background = element_blank())
 
  #theme_void()

df %>% #head()
  #mutate(count = ifelse(V2 < 0.01,"0","1")) %>% #head()#count(count)
  mutate(count = rank(-V2)) %>% #head()#count(count) 
  mutate(count = ifelse(count < 6,count,"Rare")) %>%
  separate( V1, c("A", "B", "C","DPA1","DPB1","DQA1","DQB1","DRB1"), sep = "~", remove = TRUE) %>%  #head()
  mutate(A = ifelse(A %in% ref1$HLA_type,A,'Rare')) %>% #count(A)
  mutate(B = ifelse(B %in% ref1$HLA_type,B,'Rare')) %>%
  mutate(C = ifelse(C %in% ref1$HLA_type,C,'Rare')) %>%
  mutate(DRB1 = ifelse(DRB1 %in% ref1$HLA_type,DRB1,'Rare')) %>%
  mutate(DPA1 = ifelse(DPA1 %in% ref1$HLA_type,DPA1,'Rare')) %>%
  mutate(DPB1 = ifelse(DPB1 %in% ref1$HLA_type,DPB1,'Rare')) %>%
  mutate(DQA1 = ifelse(DQA1 %in% ref1$HLA_type,DQA1,'Rare')) %>% 
  mutate(DQB1 = ifelse(DQB1 %in% ref1$HLA_type,DQB1,'Rare')) -> df1 #%>% #head()
  #mutate(subject = seq(1, n())) %>% #head()
  #gather(key, value, -count , -subject, -V2) %>% #head()
  #arrange(key, V2) -> df1
head(df1)
writexl::write_xlsx(df1,"~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/HLAhaploty_freqeuncy_2field_usingHapl-O-mat.xlsx")
library(ggalluvial)
head(vaccinations)
head(df1)
  #mutate(key = factor(key, levels = c("X1", "X2", "X3", "X4"))) %>% 
df1 %>% mutate(subject = seq(1, n())) %>% #head()
  arrange(A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1) %>%
ggplot(
       aes(axis1=A,
           axis2=B,
           axis3=C,
           axis4=DRB1,
           axis5=DPA1,
           axis6=DPB1,
           axis7=DQA1,
           axis8=DQB1,
           y=V2)) + 
  geom_alluvium(aes(fill = count)) +
  geom_stratum() + 
  geom_text(stat = "stratum",aes(label = after_stat(stratum)))+
  scale_fill_manual(values = c("red","green","blue","#d8b365","#5ab4ac","grey")) + 
  guides(fill=guide_legend(title="Top Rank of Haplotype",nrow = 1)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.title=element_text(size=20),
        legend.text = element_text(size=15),
        axis.text.x = element_text(size=14),
        panel.background = element_blank())


#https://r-charts.com/flow/ggalluvial/


df <- read.table("~/Desktop/")

head(mtcars)
head(df1)
df1 %>% #mutate(A=V2,B=V2,C=V2,DPA1=V2,DPB1=V2,DQA1=V2,DQB1=V2,DRB1=V2) %>% 
  select(-V2,-count) %>% as.matrix() -> ttt

cov.wt(ttt,df1$V2)

df1 %>% select(-count) %>% 
gather(key, value,-V2) -> df2
head(df1)
head(df2)

ggplot(df2,aes(x=key,y=key,fill=V2)) + 
  geom_tile()

df1 %>% select(-count) %>% head()
  gather(key, value,-V2) %>% #-> df2
  select(-value) %>%
  pivot_wider(names_from = key,values_from = V2)

df1 %>% select(-count) -> df2

#cov.wt(df2 %>% select(-V2) %>% as.matrix(),wt=V2)
tttt <- cov.wt(df2$key,wt=df$V2,cor=TRUE)

'''
head(df1)
rank(df1$V2)
df1 %>% mutate(rank = rank(-V2)) -> df1
alluvial(df1 %>% select(-count,-V2,-rank,A,B,C,DRB1,DPA1,DPB1,DQA1,DQB1),
         freq = df1$V2,border = NA,alpha = 0.5,
         cex=0.75,
         col = case_when(#df1$count == "1" ~ "blue",
                         df1$count == "2" ~ "grey",
                         df1$rank %in% seq(1:10) ~ 'red')
)
'''
### g

df <- read.table("~/Tool/Hapl-o-Mat/examplePopulations/a/run/g_hfs.dat")
ref <- read_xlsx("~/Desktop/KCDC/HLAimputation/99.forPaper/HLA.type.520samples.Ggroup.frequency.xlsx")
ref %>% mutate(HLA_type = ifelse(Frequency < 0.03,"other",value)) -> ref1
head(ref1)
head(df)
ref %>% mutate(check = str_split_fixed(value,":",3)[,3]) %>% group_by(gene) %>%
  count(gene,check)


display.brewer.all()
#c <- brewer.pal(5, "Set1")
c <- c("#BAB3B3EB","#FF3300","#FFFF00","#336600","#0033FF","#990099")

df %>% #group_by(V1) %>% #rank(df$V2)
  #mutate(Rank = dense_rank(desc(V2)))
  #head()
  mutate(count = ifelse(V2 < 0.01,"0","1")) %>% #head()#count(count) 
  separate( V1, c("A", "B", "C","DPA1","DPB1","DQA1","DQB1","DRB1"), sep = "~", remove = TRUE) %>%  #head()
  mutate(A = ifelse(A %in% ref1$HLA_type,A,'other')) %>% #count(A)
  mutate(B = ifelse(B %in% ref1$HLA_type,B,'other')) %>%
  mutate(C = ifelse(C %in% ref1$HLA_type,C,'other')) %>%
  mutate(DRB1 = ifelse(DRB1 %in% ref1$HLA_type,DRB1,'other')) %>%
  mutate(DPA1 = ifelse(DPA1 %in% ref1$HLA_type,DPA1,'other')) %>%
  mutate(DPB1 = ifelse(DPB1 %in% ref1$HLA_type,DPB1,'other')) %>%
  mutate(DQA1 = ifelse(DQA1 %in% ref1$HLA_type,DQA1,'other')) %>% 
  mutate(DQB1 = ifelse(DQB1 %in% ref1$HLA_type,DQB1,'other')) %>% #-> df1 #%>% #head()
  mutate(subject = seq(1, n())) %>% 
  #mutate(rank = ifelse(count == 1 & subject < 6,subject,0)) %>% #head()
  #mutate(freq = ifelse(V2 > 0.05,"0.05 ~ 1.00",ifelse(V2 <= 0.05 & V2>0.01,"0.01 ~ 0.05","0.00 ~ 0.01"))) %>%
  mutate(freq = ifelse(V2 > 0.05,"0.05 ~ 1.00",ifelse(V2 <= 0.05 & V2>0.01,"0.01 ~ 0.05","0.00 ~ 0.01"))) %>%
  gather(key, value, -freq, -count , -subject, -V2) %>% #head()
  mutate(value = ifelse(value == "other",value,str_split_fixed(value,"\\*",2)[,2])) %>% #head()
  arrange(key, V2) ->t


library(RColorBrewer)
#paletteName <- 'Set1' # 'Dark2'
brewer.pal(name=paletteName,n=3)
colorsPerCat <- brewer.pal(name=paletteName,n=3)
head(t)
t %>% arrange(-V2)

ref %>% filter(gene == "A") %>% arrange(-Frequency)
t %>%
  #mutate(key = factor(key, levels = c("X1", "X2", "X3", "X4"))) %>% 
  ggplot(
    aes(x = factor(key,levels = c("A", "B", "C","DRB1","DQA1","DQB1","DPA1","DPB1")),
        y = V2,
        stratum = value, 
        alluvium = subject,
        #fill = factor(rank),
        label = value)
        )+
  geom_flow(aes(fill = factor(freq)),
            curve_type = "linear",alpha = 1) +
  geom_stratum(width = 3/5) +
  geom_text(stat = "stratum",size=3)+
  scale_fill_manual(values=c("#BAB3B3EB","#E41A1C","#377EB8")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=15,face = "bold"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank())
#theme_void()



head(df1_HLA)
head(df1)


grep()
df1_HLA %>% select(grep("\\.1",colnames(df1_HLA))) %>% 
  gather(key,HLA1) %>% mutate(HLA1 = str_split_fixed(key,"str_s"))


head(df_geno)



#### freq changte

df <- readxl::read_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/HLAhaploty_freqeuncy_2field_usingHapl-O-mat.xlsx")
head(df)
ref <- read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLAtype_check/KMHC_HLAtype_frequency.txt",header = T)
head(ref)
out <- NULL
out <- df
for (i in c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1")) {
  ref1 <- ref %>% filter(Gene==i) %>% select(value,prop)
  colnames(ref1) <- c(i,paste0(i,"_freq"))
  out <- out %>% left_join(ref1)
}
head(out)
library(corrr)
out$prod <- out$A_freq * out$B_freq * out$C_freq * out$DPA1_freq * out$DPB1_freq * out$DQA1_freq * out$DQB1_freq * out$DRB1_freq
ggplot(out,aes(x=V2,y=prod)) +
  geom_point()

out %>% select(grep("_freq",colnames(out))) %>% as.matrix() -> a
#out %>% select(grep("_freq",colnames(out))) -> a
head(a)
a<- correlate(a)
#a<- cor(a)  
head(a)
a <- as.data.frame(a)
a[is.na(a)] <- 1
#a[is.na(a)] <- colnames(a)
#a[1,2] <- "A"
#a[1,2] <- "A"
#a[is.na(a),"A_freq"]


rownames(a) <- a$term
a %>% select(-term) %>% as.matrix() -> a
a

  
library(corrplot)
corrplot(a,
         method="color",
         type="upper",
         addCoef.col = "black",
         #order="hclust",
         col=brewer.pal(n=8, name="RdBu"))



##### haplotype frequency
kmhc <- readxl::read_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/HLAhaploty_freqeuncy_2field_usingHapl-O-mat.xlsx")
kmhc <- read.table("~/Tool/Hapl-o-Mat/examplePopulations/a/run/hfs.dat")
han <- read.table("~/Tool/Hapl-o-Mat/examplePopulations/a/run/han_hfs.dat")
pan <- read.table("~/Tool/Hapl-o-Mat/examplePopulations/a/run/PanKor_hfs.dat")
head(han)
head(kmhc)
head(pan)

#han$ref <- "Han"
#pan$ref <- "PanKor"
#kmhc$ref <- "KMHC"

colnames(han)[2] <- "Han"
colnames(pan)[2] <- "PanKor"
colnames(kmhc)[2] <- "KMHC"

#kmhc %>% rbind(han) %>% rbind(pan)
kmhc %>% left_join(han) %>% left_join(pan) %>% mutate(kmhc_rank = seq(1,nrow(kmhc))) %>%
  writexl::write_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/HLAhaplotype.freq.compare.withotherref.8gene.2field.xlsx")
  
kmhc %>% left_join(han) %>% mutate(kmhc_rank = seq(1,nrow(kmhc))) %>% #head()
  na.omit() %>%
  ggplot(aes(x=KMHC,y=Han,color=kmhc_rank)) + 
  geom_point()

