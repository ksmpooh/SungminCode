library(tidyverse)

setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_imgt3320/02.processing/")
setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_imgt3320/03.allele.matching/")
ref <- read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLA.4digit.520sample.nomenclean_withheader.chped",header = T)
df <-read.table("KBA.g1_5CV.SNP2HLAHLAimp_fd.cmp_Nomencleaner.txt",header = T)
#df <- read.table()
mt <- read.table("KBA.g1_5CV.SNP2HLAHLAimp_fd.cmp_Nomencleaner_missINFO.txt",header = T)
head(df)
head(ref)
head(mt)


ref %>% filter(!(FID %in% mt$ID)) %>% 
  pivot_longer(7:22) %>% select(value) %>% unique() -> c


head(df)
mt %>% mutate(check = ifelse(Imp %in% c$value,1,0)) %>% filter(check != 0) %>%  #head()
  pivot_wider(names_from = gene,values_from = Imp)



setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_imgt3320/03.allele.matching/")
b <- read.table("KBA.g1_5CV.SNP2HLAHLAimp_fd.cmp_Nomencleaner.txt",header = T)
b_mt <- read.table("KBA.g1_5CV.SNP2HLAHLAimp_fd.cmp_Nomencleaner_missINFO.txt",header = T)
a <- read.table("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/check/KBA.g1_5CV.SNP2HLAHLAimp_fd.cmp_Nomencleaner.txt",header = T)
a_mt <-read.table("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/check/KBA.g1_5CV.SNP2HLAHLAimp_fd.cmp_Nomencleaner_missINFO.txt",header = T)
head(a)
head(b)

setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_han/03.allele.matching/")
a <- read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/check/KBA.g1_5CV.SNP2HLAHLAimp_HanREF_fd.cmp_Nomencleaner.txt",header = T)
b <- read.table("KBA.g1_5CV.SNP2HLAHLAimp_HanREF_fd.cmp_Nomencleaner.txt",header = T)

#a <- read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/check/KBA.g1_5CV.SNP2HLAHLAimp_HanREF_fd.cmp_Nomencleaner.fdvstd.txt",header = T)
#b <- read.table("KBA.g1_5CV.SNP2HLAHLAimp_HanREF_fd.cmp_Nomencleaner.fdvstd.txt",header = T)
head(a)
head(b)

a %>% select(ID,empty_Sum,match_Sum,wrong_Sum) -> a
b %>% select(ID,empty_Sum,match_Sum,wrong_Sum) -> b

colnames(a) <- c("ID","empty_a","match_a","wrong_a")
colnames(b) <- c("ID","empty_b","match_b","wrong_b")

c <- merge(a,b)
head(c)
plot(c$empty_a,c$empty_b)
plot(c$match_a,c$match_b)
plot(c$wrong_a,c$wrong_b)

#c %>% filter(wrong_a == 2, wrong_b == 1)
a %>% filter(ID == "NIH19KT2254")
a_mt %>% filter(ID == "NIH19KT2254")
b_mt %>% filter(ID == "NIH19KT2254")




setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLAtype_check/other_panel/noman")

df <-read.table("Han.hg19.haplegendtovcf.modify.hlatype_fornomenclean_2field.chped")
head(df)
df %>% pivot_longer(7:22) %>% select(value) %>%
  unique()-> out

colnames(out) <- "Han"
dim(out)
#write.table(out,"../hlatype/Han.Nomen.2field",col.names = T,row.names = F,quote = F)


df <-read.table("PanKor_merged.hg19.haplegendtovcf.modify.hlatype_fd.IDchange_fornomenclean_2field.chped")
head(df)
df %>% pivot_longer(7:22) %>% select(value) %>%
  unique()-> out
dim(out)
colnames(out) <- "PanKor"
head(out)
#write.table(out,"../hlatype/PanKor.Nomen.2field",col.names = T,row.names = F,quote = F)

df <-read.table("HLAtyping.1000genomePhase3.rouph_forHLAPATAS_edit_fornomenclean_2field.chped")
head(df)
df %>% pivot_longer(c(7,8,9,10,11,12,19,20,21,22)) %>% select(value) %>%
  unique()-> out

colnames(out) <- "1KGP"
head(out)
#unique(out$`1KGP`)
dim(out)
#write.table(out,"../hlatype/1KGP.Nomen.2field",col.names = T,row.names = F,quote = F)


kmhc_type %>% unique() %>% dim()



### frequency check to make over 1%
ref <- read.table("/Users/ksmpooh/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/IMGT3320/HLA.typing.Final.result_520sample_IMGT3320convert_forReference_2field_withheader.chped",header = T)
head(ref)

## 1 field
ref %>% select(-IID,-FID,-pID,-mID,-SEX,-PHENO) %>% pivot_longer(1:16,names_to = 'Gene') %>% 
  mutate(Gene=str_split_fixed(Gene,"\\.",2)[,1]) %>% 
  mutate(Gene=str_split_fixed(Gene,"HLA_",2)[,2]) %>% 
  mutate(value=str_split_fixed(value,":",2)[,1]) %>% 
  filter(value !=0) %>% 
  count(Gene,value) %>% 
  group_by(Gene) %>% 
  mutate(prop = prop.table(n)) %>% #summarise(sum = sum(n))
  mutate(freq = ifelse(prop<0.01,"rare",ifelse(prop<0.05,"less common","common"))) -> hla_freq_1field
head(hla_freq_1field)
#writexl::write_xlsx(hla_freq_1field,"~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/KMHC.HLAtype.freq_1field.xlsx")
## 2 field
ref %>% select(-IID,-FID,-pID,-mID,-SEX,-PHENO) %>% pivot_longer(1:16,names_to = 'Gene') %>%
  mutate(Gene=str_split_fixed(Gene,"\\.",2)[,1]) %>% 
  mutate(Gene=str_split_fixed(Gene,"HLA_",2)[,2]) %>%
  filter(value !=0) %>%
  count(Gene,value) %>% 
  group_by(Gene) %>%
  mutate(prop = prop.table(n)) %>% 
  mutate(freq = ifelse(prop<0.01,"rare",ifelse(prop<0.05,"less common","common"))) -> hla_freq

table(hla_freq$freq)
head(hla_freq)


ggplot(hla_freq,aes(x=Gene,y=prop,color=freq)) + 
  geom_bar(stat='identity')

hla_freq %>% count(Gene,freq) %>% 
  mutate(freq = ifelse(freq == "rare","Rare: (0, 0.01)",ifelse(freq == "less common","Less common: [0.01, 0.05) ","Common: [0.05,1)"))) %>% #head()
  ggplot(aes(x=Gene,y=n,fill=freq)) +
  geom_bar(stat='identity', position = 'dodge') + 
  theme(legend.position = "bottom",
        legend.text=element_text(size=11),
        axis.title.x = element_text(size = 14,face = "bold"),
        axis.title.y = element_text(size = 14,face = "bold"))


head(hla_freq)

#write.table(hla_freq,"~/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLAtype_check/KMHC_HLAtype_frequency.txt",sep = "\t",quote = F,col.names = T,row.names = F)

ref <- read.table("/Users/ksmpooh/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/IMGT3320/HLA.typing.Final.result_520sample_IMGT3320convert_forReference_2field.chped",header = F)
h_ref <-read.table("/Users/ksmpooh/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/IMGT3320/HLA.typing.Final.result_520sample_IMGT3320convert_forReference_4field.chped",header = T)
head(ref)
head(h_ref)
head(hla_freq)

ref %>% mutate(across(c(V7:V22),~ifelse(. %in% hla_freq[hla_freq$freq!="rare",]$value,.,"0"))) -> out
head(out)
colnames(out) <- colnames(h_ref)
write.table(out,"/Users/ksmpooh/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/IMGT3320/HLA.typing.Final.result_520sample_IMGT3320convert_forReference_2field_notrare.chped",row.names = F,col.names = T,sep = "\t",quote = F)


ref %>% mutate(across(c(V7:V22),~ifelse(. %in% hla_freq[hla_freq$freq=="common",]$value,.,"0"))) -> out
colnames(out) <- colnames(h_ref)
head(out)
write.table(out,"/Users/ksmpooh/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/IMGT3320/HLA.typing.Final.result_520sample_IMGT3320convert_forReference_2field_common.chped",row.names = F,col.names = T,sep = "\t",quote = F)


ref %>% mutate(across(c(V7:V22),~ifelse(. %in% hla_freq[hla_freq$freq=="less common",]$value,.,"0"))) -> out

head(out)
colnames(out) <- colnames(h_ref)
write.table(out,"/Users/ksmpooh/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/IMGT3320/HLA.typing.Final.result_520sample_IMGT3320convert_forReference_2field_lesscommon.chped",row.names = F,col.names = T,sep = "\t",quote = F)

ref %>% mutate(across(c(V7:V22),~ifelse(. %in% hla_freq[hla_freq$freq=="rare",]$value,.,"0"))) -> out
colnames(out) <- colnames(h_ref)
head(out)
write.table(out,"/Users/ksmpooh/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/IMGT3320/HLA.typing.Final.result_520sample_IMGT3320convert_forReference_2field_rare.chped",row.names = F,col.names = T,sep = "\t",quote = F)


head(out)
head(h_ref)
head(ref)
head(hla_freq)

custom.col <- c("#FFDB6D", "#C4961A", "#F4EDCA", 
                "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")

df <- read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLAtype_check/KMHC_HLAtype_frequency.txt",header = T)
head(df)


df %>% mutate(freq = ifelse(n==1,"singleton",ifelse(n == 2,"doubleton",
                                                    ifelse(n > 2 & freq == "rare","doubleton ~ rare",
                                                           ifelse(freq == "less_common","less common",freq))))) %>% #$head
  group_by(Gene,freq) %>% count(Gene) %>% #writexl::write_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/KMHC.HLAtype.freq.count.xlsx")
  ggplot(aes(x=factor(freq,c('singleton','doubleton','doubleton ~ rare','less common','common')),y=n,fill=Gene))+
  geom_bar(stat="identity") + 
  scale_fill_manual(values = rev(custom.col)) + 
  ylab("Count") + 
  theme(axis.title.x = element_blank(),
        #axis.title.y = element_text(),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(size=10),
        axis.ticks.x = element_blank())
        #axis.ticks.y = element_blank()
        #panel.background = element_blank())


df %>% mutate(freq = ifelse(n==1,"singleton",ifelse(n == 2,"doubleton",
                                                      ifelse(n > 2 & freq == "rare","doubleton ~ rare",freq)))) %>% #$head
  group_by(Gene,freq) %>% count(Gene) %>% #head
  pivot_wider(names_from = Gene,values_from = n) %>%
  writexl::write_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/KMHC.HLAtype.freq.count.xlsx")
  


#######
ref <- readxl::read_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/KMHC.HLAtype.freq.count.xlsx")
ref <- read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLAtype_check/KMHC_HLAtype_frequency.txt",header = T)
#ref %>% writexl::write_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/KMHC.HLAtype.freq.xlsx")
head(ref)
ref %>% count(freq)
flist = grep(list.files("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_imgt3320/01.impresult/"),pattern = ".txt", invert=FALSE, value=TRUE)

df <- NULL
for (a in flist) {
  tmp <- read.table(paste0("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_imgt3320/01.impresult/",a))
  tmp$batch <- a
  df <- rbind(df,tmp)
}
head(df)
head(ref)
ref %>% select(Gene,value,freq) -> ref
head(ref)

df %>% mutate(batch = str_split_fixed(batch,"_",2)[,1]) %>% 
  mutate(value = str_split_fixed(V1,"_",2)[,2]) %>% 
  left_join(ref) %>% na.omit() %>% #dim()
  select(batch,value,Gene,freq,V4) %>% filter(V4 != 0) %>% #head()
  group_by(Gene,freq,value) %>%
  summarise(DS2 = mean(V4)) %>%
  ggplot(aes(x=freq,y=DS2,fill=Gene)) + 
  geom_boxplot() + 
  scale_fill_manual(values = rev(custom.col)) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15,face = "bold"),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.ticks.x = element_blank())
#axis.ticks.y = element_blank()
#panel.background = element_blank())


  
