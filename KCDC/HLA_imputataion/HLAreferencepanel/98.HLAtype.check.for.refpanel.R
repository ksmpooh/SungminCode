## HLA type compore 2023 03 20
## michigan vs Han vs Pan+korea vs KBA
library(tidyverse)
library(stringr)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(ggthemes)
library(reshape2)
setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLAtype_check/")

michigan <- read.table("michigan.hla.type.info")
michigan <- michigan[,1:2]
head(michigan)
michigan %>% mutate("Ref" ="Multi-ethnic","Gene" = str_split_fixed(michigan$V1,"\\*",2)[,1], "type" = str_split_fixed(michigan$V1,"_",2)[,2]) %>% #head()
  mutate(digit = ifelse(str_split_fixed(type,":",2)[,2] == "","td","fd")) %>% #head()
  select(-V1,-V2) -> michigan
#michigan$Ref <- "Multi-ethnic"
#michigan$Gene <- str_split_fixed(michigan$V1,"\\*",2)[,1]
#michigan$type <- str_split_fixed(michigan$V1,"_",2)[,2]
michigan %>% head()
'
michigan %>% filter(digit=="fd") %>% select(type) %>%
  rename("michigan" = type) %>% #head()
  write.table("other_panel/hlatype/michigan.Nomen.2field.txt",col.names = T,row.names = F,quote = F)
'

'
           Ref  Gene     type digit
1 Multi-ethnic HLA_A     A*01    td
2 Multi-ethnic HLA_A  A*01:01    fd
3 Multi-ethnic HLA_A  A*01:02    fd

'
table(michigan$Gene)
han <- read.table("Han.hg19.haplegendtovcf.modify.hlatype_fd.txt",header = T)
pan <- read.table("PanKor_merged.hg19.haplegendtovcf.modify.hlatype_fd.IDchange_fornomenclean.txt",header = F) %>%
  select(-V2,-V3,-V4,-V5,-V6)

#tmp <- readxl::read_xlsx("HLAtype.compare.IMGT3320_check.xlsx",sheet = 2)
#head(tmp)
#kmhc <- read.table("HLA.4digit.520sample.nomenclean.chped") %>%
kmhc <- read.table("~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/IMGT3320/HLA.typing.Final.result_520sample_IMGT3320convert_forReference_2field.chped") %>%
  select(-V2,-V3,-V4,-V5,-V6)

ref <- read.table('../Final_520sample.index.txt',header = T)
head(ref)
kmhc <- read.table("../HLA.type.result.8genes.merged.4digit_529sample_forMAKEreference.txt",header = T) %>%
  select(-IID,-pID,-mID,-SEX,-PHENO) %>% filter(FID %in% ref$KBAID)

kmhc <- read.table("~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/IMGT3320/HLA.typing.Final.result_520sample_IMGT3320convert_forReference_2field.chped") %>%
  select(-V2,-V3,-V4,-V5,-V6)
head(kmhc)

#kgp <- read.table("/Users/ksmpooh/Desktop/KCDC/HLAimputation/1000genome/HLAtyping.1000genomePhase3.SampleInfo.typeInfo.txt",header = T)
kgp <- read.table("/Users/ksmpooh/Desktop/KCDC/HLAimputation/1000genome/HLAtyping.1000genomePhase3.rouph_forHLAPATAS_edit_fornomenclean_2field.chped",header = F) %>%
  select(-V2,-V3,-V4,-V5,-V6)  

head(kgp)

head(kmhc)
head(han)
head(pan)
colnames(pan) <- colnames(han)
colnames(kmhc) <- colnames(han)
colnames(kgp) <- colnames(han)
#pivot_longer(cols = c(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DRB1,HLA_DQA1,HLA_DQB1,overall),names_to = 'Gene',values_to = 'Accuracy') %>%
han %>% pivot_longer(2:ncol(han),names_to = 'Gene',values_to = 'type') %>% #head()
  mutate(Gene = str_replace_all(Gene,"\\.1","")) %>% mutate(Gene = str_replace_all(Gene,"\\.2","")) %>% #count(Gene)#head()
  mutate(Ref = "Han Chinese") %>% select(-ID) %>% unique() %>% filter(type != "0") -> han_type

pan %>% pivot_longer(2:ncol(han),names_to = 'Gene',values_to = 'type') %>% #head()
  mutate(Gene = str_replace_all(Gene,"\\.1","")) %>% mutate(Gene = str_replace_all(Gene,"\\.2","")) %>% #count(Gene)#head()
  mutate(Ref = "Pan-Kor") %>% select(-ID) %>% unique() %>% filter(type != "0") -> pan_type

kmhc %>% pivot_longer(2:ncol(han),names_to = 'Gene',values_to = 'type') %>% #head()
  mutate(Gene = str_replace_all(Gene,"\\.1","")) %>% mutate(Gene = str_replace_all(Gene,"\\.2","")) %>% #count(Gene)#head()
  mutate(Ref = "KMHC") %>% select(-ID) %>% unique() %>% filter(type != "0") -> kmhc_type

kgp %>% #select(-Region,-Population,-Sample.ID) %>% #head()
  select(-ID,-HLA_DPA1.1,-HLA_DPA1.2,-HLA_DPB1.1,-HLA_DPB1.2,-HLA_DQA1.1,-HLA_DQA1.2) %>% #head()
  pivot_longer(1:10,names_to = 'Gene',values_to = 'type') %>% #head()
  mutate(Gene = str_replace_all(Gene,"\\.1","")) %>% mutate(Gene = str_replace_all(Gene,"\\.2","")) %>% #count(Gene)#head()
  mutate(Ref = "1KGP") %>% unique() %>% filter(type != "0") -> kgp_type


head(kgp_type)
head(han_type)
table(han_type$Gene)
michigan %>% filter(digit == 'fd') %>% select(-digit) %>%
  rbind(han_type) %>% rbind(pan_type) %>% rbind(kmhc_type) %>% mutate(Gene = str_replace_all(Gene,"HLA_","")) -> fd

michigan %>% filter(digit == 'fd') %>% select(-digit) %>%
  rbind(han_type) %>% rbind(pan_type) %>% rbind(kmhc_type) %>% rbind(kgp_type) %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) -> fd_withkgp


head(fd)
head(td)
td %>% count(Ref)
fd %>% mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique() -> td
fd_withkgp %>% mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique() -> td_withkgp
fd$digit <- "4 digit"
td$digit <- "2 digit"

fd$digit <- "Two-field"
td$digit <- "One-field"


fd_withkgp$digit <- "4 digit"
td_withkgp$digit <- "2 digit"

fd_withkgp$digit <- "Two-field"
td_withkgp$digit <- "One-field"


fd %>% rbind(td) %>%
  count(Gene,Ref,digit) -> test
  

#fd %>% rbind(td) %>% writexl::write_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/HLAtype_4panel.xlsx")

fd %>% rbind(td) %>%
  count(Gene,Ref,digit) %>% #head()
  ggplot(aes(x=Gene,y=n,fill = factor(Ref,levels = c("KMHC","Multi-ethnic","Han Chinese","Pan-Kor")))) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  xlab("HLA gene") + ylab("Number of HLA type") + 
  scale_fill_discrete(name="Ref.Panel") +
  theme(legend.title=element_text(size=12,face = "bold"),
        legend.text=element_text(size=10),
        axis.title.x = element_text(size = 13,face = "bold"),
        axis.title.y = element_text(size = 13,face = "bold"),
        legend.position = "right") +
  facet_grid(~digit) + 
  theme(strip.text.x = element_text(size = 15,face = "bold"))


fd_withkgp %>% rbind(td_withkgp) %>% 
  filter(Gene %in% c("A","B","C","DRB1","DQB1")) %>% #dim()
  count(Gene,Ref,digit) %>% #head()
  ggplot(aes(x=Gene,y=n,fill = factor(Ref,levels = c("KMHC","Multi-ethnic","Han Chinese","Pan-Kor","1KGP")))) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  xlab("HLA gene") + ylab("Number of HLA type") + 
  scale_fill_discrete(name="Ref.Panel") +
  theme(legend.title=element_text(size=12,face = "bold"),
        legend.text=element_text(size=10),
        axis.title.x = element_text(size = 13,face = "bold"),
        axis.title.y = element_text(size = 13,face = "bold"),
        legend.position = "bottom") +
  facet_grid(~digit) + 
  theme(strip.text.x = element_text(size = 13,face = "bold"))



michigan %>% filter(digit == 'fd') %>% select(-digit) %>%
  rbind(han_type) %>% rbind(pan_type) %>% rbind(kmhc_type) %>% mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% head()

VennDiagram::venn.diagram(x=list())


x <- list(KMHC = kmhc_type$type,`Multi-ethnic` = michigan[michigan$digit =='fd',]$type,`Han Chinese`  = han_type$type,`Pan-Kor` = pan_type$type)
x <- list(KMHC = td[td$Ref == "KMHC",]$type,`Multi-ethnic` = td[td$Ref == "Multi-ethnic",]$type,`Han Chinese`  = td[td$Ref == "Han Chinese",]$type,`Pan-Kor` = td[td$Ref == "Pan-Kor",]$type)
ggvenn::ggvenn(
  x,
  stroke_size = 0.5, set_name_size = 4
)
head(kmhc_type)
dim(kmhc_type)
kmhc_type %>% unique() %>% dim()
han_type %>% unique() %>% dim()
x <- list(KMHC = kmhc_type$type,`Multi-ethnic` = michigan[michigan$digit =='fd',]$type,`Han Chinese`  = han_type$type,`Pan-Kor` = pan_type$type,`1KGP` = kgp_type$type)
ggvenn::ggvenn(
  x,
  stroke_size = 0.5, set_name_size = 4
)

head(fd_withkgp)
table(fd_withkgp$Ref)
fd_withkgp <- fd_withkgp %>% filter(Gene %in% c("A","B","C","DRB1","DQB1"))
td_withkgp <- td_withkgp %>% filter(Gene %in% c("A","B","C","DRB1","DQB1"))
x <- list(KMHC = td_withkgp[td_withkgp$Ref=="KMHC",]$type,`Multi-ethnic` = td_withkgp[td_withkgp$Ref=="Multi-ethnic",]$type,`Han Chinese`  = td_withkgp[td_withkgp$Ref=="Han Chinese",]$type,`Pan-Kor` = td_withkgp[td_withkgp$Ref=="Pan-Kor",]$type,`1KGP` = td_withkgp[td_withkgp$Ref=="1KGP",]$type)
ggVennDiagram::ggVennDiagram(x) + 
  ggplot2::scale_fill_gradient(low="white",high = "blue")


x <- list(KMHC = fd_withkgp[fd_withkgp$Ref=="KMHC",]$type,`Multi-ethnic` = fd_withkgp[fd_withkgp$Ref=="Multi-ethnic",]$type,`Han Chinese`  = fd_withkgp[fd_withkgp$Ref=="Han Chinese",]$type,`Pan-Kor` = fd_withkgp[fd_withkgp$Ref=="Pan-Kor",]$type,`1KGP` = fd_withkgp[fd_withkgp$Ref=="1KGP",]$type)
ggVennDiagram::ggVennDiagram(x) + 
  ggplot2::scale_fill_gradient(low="white",high = "blue")

fd_withkgp %>% count(Ref,Gene)

head(kmhc_type)
head(michigan)
df_venn 
michigan %>% filter(digit=='fd') %>% select(-digit) %>%
  rbind(kmhc_type) %>% rbind(han_type) %>% rbind(pan_type) %>% #head()
  #filter(Gene == "HLA_DPA1") -> df_venn
  #filter(Gene == "HLA_DPB1") -> df_venn
  filter(Gene == "HLA_DQA1") -> df_venn
  #filter(Gene == "HLA_DQB1") -> df_venn
head(df_venn)
#table(df_venn
table(df_venn$Ref)
#x<- list(KMHC = df_venn[df_venn == 'KMHC',]$type,`Multi-ethnic` = df_venn[df_venn == 'Multi-ethnic',]$type,`Han Chinese`  = df_venn[df_venn == 'Han Chinese',]$type,`Pan-Kor` = df_venn[df_venn == 'Pan-Kor',]$type)
x<- list(KMHC = df_venn[df_venn$Ref == 'KMHC',]$type,`Pan-Kor` = df_venn[df_venn$Ref == 'Pan-Kor',]$type)
ggvenn::ggvenn(
  x,
  stroke_size = 0.5, set_name_size = 7,text_size = 6
)

ggvenn::ggvenn(
  x,
  stroke_size = 0.5, set_name_size = 7,text_size = 6,
)


x<- list(KMHC = df_venn[df_venn == 'KMHC',]$type,`Pan-Kor` = df_venn[df_venn == 'Pan-Kor',]$type)
ggvenn::ggvenn(
  x,
  stroke_size = 0.5, set_name_size = 4
)


  
# 2 digit

michigan %>% filter(digit == 'td') %>% select(-digit) %>% #head()
  rbind(han_type) %>% rbind(pan_type) %>% rbind(kmhc_type) %>% mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% 
  mutate(type = str_split_fixed(type,":",2)[,1]) %>% #head()
  unique() %>%
  count(Ref,Gene) %>% #head()
  group_by(Ref) %>%
  pivot_wider(names_from = Gene,values_from = n)

michigan %>% filter(digit == 'td') %>% select(-digit) %>% #head()
  rbind(han_type) %>% rbind(pan_type) %>% rbind(kmhc_type) %>% rbind(kgp_type) %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% 
  mutate(type = str_split_fixed(type,":",2)[,1]) %>% #head()
  unique() %>%
  count(Ref,Gene) %>% #head()
  group_by(Ref) %>%
  pivot_wider(names_from = Gene,values_from = n)

michigan %>% filter(digit != 'td') %>% select(-digit) %>% #head()
  rbind(han_type) %>% rbind(pan_type) %>% rbind(kmhc_type) %>% rbind(kgp_type) %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% 
  #mutate(type = str_split_fixed(type,":",2)[,1]) %>% #head()
  unique() %>%
  count(Ref,Gene) %>% #head()
  group_by(Ref) %>%
  pivot_wider(names_from = Gene,values_from = n)

head(kgp_type)
kgp_type %>% count(Gene)
kgp_type %>% dim()
kgp_type %>% unique %>% dim()

michigan %>% filter(digit == 'td') %>% select(-digit) %>% #head()
  rbind(han_type) %>% rbind(pan_type) %>% rbind(kmhc_type) %>% mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% 
  mutate(type = str_split_fixed(type,":",2)[,1]) %>% #head()
  unique() %>%
  count(Ref,Gene,type) %>% #head()
  select(-n,-Gene) -> a
a %>% na.omit() ->a
table(a$Ref)
head(a)

ggvenn::ggvenn(
  list(KMHC = a[a$Ref == "KMHC",]$type,
       `Multi-ethnic` = a[a$Ref == "Multi-ethnic",]$type,
      `Han Chinese` = a[a$Ref == "Han Chinese",]$type,
      `Pan-Kor`= a[a$Ref == "Pan-Kor",]$type),
  stroke_size = 0.5, set_name_size = 4
)

a %>% count(Ref)
head(a)
#####
han %>% pivot_longer(2:ncol(han),names_to = 'Gene',values_to = 'type') %>% #head()
  mutate(Gene = str_replace_all(Gene,"\\.1","")) %>% mutate(Gene = str_replace_all(Gene,"\\.2","")) %>% #count(Gene)#head()
  mutate(Ref = "Han Chinese") %>% select(-ID) ->han_freq
  
han_freq %>% na.omit() %>%
  mutate(Gene = str_split_fixed(Gene,"_",2)[,2]) %>%
  group_by(Ref,Gene) %>%
  count(type) %>% #head()
  mutate(Frequency = prop.table(n)) %>% 
  filter(type %in% c("C*08:41"," DRB1*14:141","DRB1*12:17","DQA1*05:06","DQA1*05:08",
                     "C*08:01","DRB1*14:03","DRB1*12:01","DQA1*05:03","DQA1*05:05"))

  

pan %>% pivot_longer(2:ncol(han),names_to = 'Gene',values_to = 'type') %>% #head()
  mutate(Gene = str_replace_all(Gene,"\\.1","")) %>% mutate(Gene = str_replace_all(Gene,"\\.2","")) %>% #count(Gene)#head()
  mutate(Ref = "Pan-Kor") %>% select(-ID) -> pan_freq

kmhc %>% pivot_longer(2:ncol(han),names_to = 'Gene',values_to = 'type') %>% #head()
  mutate(Gene = str_replace_all(Gene,"\\.1","")) %>% mutate(Gene = str_replace_all(Gene,"\\.2","")) %>% #count(Gene)#head()
  mutate(Ref = "KMHC") %>% select(-ID) ->kmhc_freq

head(kmhc_freq)

kmhc_freq %>% rbind(han_freq) %>% rbind(pan_freq) %>% #head()
  na.omit() %>%
  mutate(Gene = str_split_fixed(Gene,"_",2)[,2]) %>%
  group_by(Ref,Gene) %>%
  count(type) %>% #head()
  mutate(Frequency = prop.table(n)) %>% #->han_freq
  filter(Ref %in% c("Han Chinese","KMHC")) %>%
  select(-n) %>%
  pivot_wider(names_from = Ref,values_from = Frequency) %>% #head()
  na.omit() %>% #count(Gene) #%>% summarise(sum = sum(n)) #head()
  #corrr::correlate()
  ggplot(aes(x=KMHC,y=`Han Chinese`,color=Gene)) +
  geom_point() + 
  stat_smooth(method = 'lm', se=F, color='blue') +
  geom_abline(slope=1, intercept = 0,color = 'gray') + 
  geom_text(size = 8, x=0.3, y=0, label="Cor = 0.86",color ="black")+
  theme(legend.title=element_text(size=13,face = "bold"),
      #legend.position = "bottom",
      legend.text=element_text(size=11),
      axis.title.x = element_text(size = 14,face = "bold"),
      axis.title.y = element_text(size = 14,face = "bold"))
  

kmhc_freq %>% rbind(han_freq) %>% rbind(pan_freq) %>% #head()
  na.omit() %>%
  mutate(Gene = str_split_fixed(Gene,"_",2)[,2]) %>%
  group_by(Ref,Gene) %>%
  count(type) %>% #head()
  mutate(Frequency = prop.table(n)) %>% #->han_freq
  filter(Ref %in% c("Pan-Kor","KMHC")) %>%
  select(-n) %>%
  pivot_wider(names_from = Ref,values_from = Frequency) %>% #head()
  na.omit() %>% #count(Gene) #%>% summarise(sum = sum(n)) #head()
  #corrr::correlate()
  ggplot(aes(x=KMHC,y=`Pan-Kor`,color=Gene)) +
  geom_point() + 
  stat_smooth(method = 'lm', se=F, color='blue') +
  geom_abline(slope=1, intercept = 0,color = 'gray') + 
  geom_text(size = 8,x=0.3, y=0, label="Cor = 0.71",color ="black") + 
  theme(legend.title=element_text(size=13,face = "bold"),
        #legend.position = "bottom",
        legend.text=element_text(size=11),
        axis.title.x = element_text(size = 14,face = "bold"),
        axis.title.y = element_text(size = 14,face = "bold"))





  


### 5 vold type check

ref <- read.table('../Final_520sample.index.txt',header = T)
head(ref)
kmhc <- read.table("../HLA.type.result.8genes.merged.4digit_529sample_forMAKEreference.txt",header = T) %>%
  select(-IID,-pID,-mID,-SEX,-PHENO) %>% filter(FID %in% ref$KBAID)
head(kmhc)
colnames(kmhc) <- colnames(han)
ncol(kmhc)
ref %>% merge(kmhc,by.x = "KBAID",by.y="ID") %>% #head()
  pivot_longer(3:18,names_to = 'Gene',values_to = 'type') %>% #head()
  mutate(Gene = str_replace_all(Gene,"\\.1","")) %>% mutate(Gene = str_replace_all(Gene,"\\.2","")) %>% #count(Gene)#head()
  unique() -> a

table(a$group)
a[a$group == 1,]$type

x <- list(g1  = a[a$group == 1,]$type,
          g2  = a[a$group == 2,]$type,
          g3  = a[a$group == 3,]$type,
          g4  = a[a$group == 4,]$type,
          g5  = a[a$group == 5,]$type)
ggvenn::ggvenn(
  list(one_fold  = a[a$group == 1,]$type,
         two_fold  = a[a$group == 2,]$type,
         three_fold  = a[a$group == 3,]$type,
         four_fold  = a[a$group == 4,]$type,
         five_fold  = a[a$group == 5,]$type)
)

VennDiagram::venn.diagram(list(g1  = a[a$group == 1,]$type,
                               g2  = a[a$group == 2,]$type,
                               g3  = a[a$group == 3,]$type,
                  four_fold  = a[a$group == 4,]$type,
                  five_fold  = a[a$group == 5,]$type),
                  filename = "test.png"
             )

venn::venn(x,zcolor = "style")


###### fd td

head(fd)
head(td)
michigan %>% filter(digit == "fd") %>% select(-digit) -> michigan_type

fd %>% rbind(td) %>% select(-Gene) %>%
  pivot_wider(names_from = Ref,values_from = type) %>% #head()
  count()

head(td)
head(fd)
han_type
pan_type
michigan_type
kmhc_type
kgp_type

table(kmhc_type$type %in% michigan_type$type)
head(kmhc_type)
kmhc_type %>% mutate("count" = 1) %>%  #head()
  mutate(count = ifelse(type %in% michigan_type$type,count+1,count)) %>% #head()
  mutate(count = ifelse(type %in% han_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% pan_type$type,count+1,count)) %>% #head()
  count(Ref,count) -> kmhc_type_count

head(kmhc_type_count)

michigan_type %>% mutate("count" = 1) %>%
  mutate(count = ifelse(type %in% kmhc_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% han_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% pan_type$type,count+1,count)) %>% count(Ref,count) -> michigan_type_count

michigan_type_count

han_type %>% mutate("count" = 1) %>%
  mutate(count = ifelse(type %in% kmhc_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% michigan_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% pan_type$type,count+1,count)) %>% count(Ref,count) -> han_type_count

pan_type %>% mutate("count" = 1) %>%
  mutate(count = ifelse(type %in% kmhc_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% han_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% michigan_type$type,count+1,count)) %>% count(Ref,count) -> pan_type_count

kmhc_type_count %>% rbind(michigan_type_count) %>% rbind(han_type_count) %>% rbind(pan_type_count) %>% #head()
  mutate(count = as.factor(count),"digit" = "Two-field") -> fd_count

head(fd_count)
fd_count %>% count(Ref)
#kmhc_type2 <- kmhc_type
#michigan_type2 <- michigan_type
#han_type2 <- han_type
#pan_type2 <- pan_type

kmhc_type %>% mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique() -> kmhc_type2
michigan_type %>% mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique() -> michigan_type2
han_type %>% mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique() ->han_type2
pan_type %>% mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique() ->pan_type2

kmhc_type2 %>% #mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique() %>% head()
  mutate("count" = 1) %>% 
  mutate(count = ifelse(type %in% michigan_type2$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% han_type2$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% pan_type2$type,count+1,count)) %>% count(Ref,count) -> kmhc_type_count
#head(kmhc_type_count)

#michigan %>% filter(digit == "td") %>% select(-digit) %>%
#michigan_type %>% mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique() -> michigan_type
michigan_type2 %>%
  mutate("count" = 1) %>%
  mutate(count = ifelse(type %in% kmhc_type2$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% han_type2$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% pan_type2$type,count+1,count)) %>% count(Ref,count) -> michigan_type_count

#han_type %>% mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique() ->han_type
han_type2 %>% #mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique() %>%
  mutate("count" = 1) %>%
  mutate(count = ifelse(type %in% kmhc_type2$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% michigan_type2$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% pan_type2$type,count+1,count)) %>% count(Ref,count) -> han_type_count

#pan_type %>% mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique() ->pan_type
pan_type2 %>% #mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique()%>%
  mutate("count" = 1) %>%
  mutate(count = ifelse(type %in% kmhc_type2$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% han_type2$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% michigan_type2$type,count+1,count)) %>% count(Ref,count) -> pan_type_count

kmhc_type_count %>% rbind(michigan_type_count) %>% rbind(han_type_count) %>% rbind(pan_type_count) %>% #head()
  mutate(count = as.factor(count),"digit" = "One-field") -> td_count

td_count
fd_count

#head(han_type_count)
fd_count %>% rbind(td_count) %>% group_by(Ref) %>% #head()
  mutate(per=paste0(n,"(",round(n/sum(n)*100, 0), "%)")) %>%
  mutate(per=ifelse(count==1,per,NA)) %>%# head()
  ggplot(aes(x= factor(Ref,levels = c("KMHC","Multi-ethnic","Han Chinese","Pan-Kor")),y=n,fill=count)) +
  geom_bar(stat="identity") + 
  geom_text(size = 4,aes(label = per),hjust = 0.5, vjust = -1, position = "stack") +
  #update_geom_defaults("text", list(size = 12)) +
  scale_fill_discrete(name="# of overlapping HLA type") +
  ylim(c(0,450)) +
  ylab("# of HLA type") + 
  theme(legend.title=element_text(size=13,face = "bold"),
        legend.position = "bottom",
        legend.text=element_text(size=11),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14,face = "bold")) +
  facet_grid(~digit) + 
  theme(strip.text.x = element_text(size = 13,face = "bold"))



#### Michigan vs KMHC type

library("gridExtra") 

ggvenn::ggvenn(
  list(KMHC = fd[(fd$Ref == "KMHC") & (fd$Gene == "A"),]$type,
       `Multi-ethnic` = fd[(fd$Ref == "Multi-ethnic") & (fd$Gene == "A"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p1
p1
ggvenn::ggvenn(
  list(KMHC = fd[(fd$Ref == "KMHC") & (fd$Gene == "B"),]$type,
       `Multi-ethnic` = fd[(fd$Ref == "Multi-ethnic") & (fd$Gene == "B"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p2

ggvenn::ggvenn(
  list(KMHC = fd[(fd$Ref == "KMHC") & (fd$Gene == "C"),]$type,
       `Multi-ethnic` = fd[(fd$Ref == "Multi-ethnic") & (fd$Gene == "C"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p3

ggvenn::ggvenn(
  list(KMHC = fd[(fd$Ref == "KMHC") & (fd$Gene == "DPA1"),]$type,
       `Multi-ethnic` = fd[(fd$Ref == "Multi-ethnic") & (fd$Gene == "DPA1"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p4

ggvenn::ggvenn(
  list(KMHC = fd[(fd$Ref == "KMHC") & (fd$Gene == "DPB1"),]$type,
       `Multi-ethnic` = fd[(fd$Ref == "Multi-ethnic") & (fd$Gene == "DPB1"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p5

ggvenn::ggvenn(
  list(KMHC = fd[(fd$Ref == "KMHC") & (fd$Gene == "DQA1"),]$type,
       `Multi-ethnic` = fd[(fd$Ref == "Multi-ethnic") & (fd$Gene == "DQA1"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p6

ggvenn::ggvenn(
  list(KMHC = fd[(fd$Ref == "KMHC") & (fd$Gene == "DQB1"),]$type,
       `Multi-ethnic` = fd[(fd$Ref == "Multi-ethnic") & (fd$Gene == "DQB1"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p7

ggvenn::ggvenn(
  list(KMHC = fd[(fd$Ref == "KMHC") & (fd$Gene == "DRB1"),]$type,
       `Multi-ethnic` = fd[(fd$Ref == "Multi-ethnic") & (fd$Gene == "DRB1"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p8

library(cowplot)

plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,labels = c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1"),
          ncol = 3)


#### 1KGP vs KMHC type

library("gridExtra") 
head(fd_withkgp)
table(fd_withkgp$Ref)
ggvenn::ggvenn(
  list(KMHC = fd_withkgp[(fd_withkgp$Ref == "KMHC") & (fd_withkgp$Gene == "A"),]$type,
       `1KGP` = fd_withkgp[(fd_withkgp$Ref == "1KGP") & (fd_withkgp$Gene == "A"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p1
p1
ggvenn::ggvenn(
  list(KMHC = fd_withkgp[(fd_withkgp$Ref == "KMHC") & (fd_withkgp$Gene == "B"),]$type,
       `1KGP` = fd_withkgp[(fd_withkgp$Ref == "1KGP") & (fd_withkgp$Gene == "B"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p2
ggvenn::ggvenn(
  list(KMHC = fd_withkgp[(fd_withkgp$Ref == "KMHC") & (fd_withkgp$Gene == "C"),]$type,
       `1KGP` = fd_withkgp[(fd_withkgp$Ref == "1KGP") & (fd_withkgp$Gene == "C"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p3

ggvenn::ggvenn(
  list(KMHC = fd_withkgp[(fd_withkgp$Ref == "KMHC") & (fd_withkgp$Gene == "DRB1"),]$type,
       `1KGP` = fd_withkgp[(fd_withkgp$Ref == "1KGP") & (fd_withkgp$Gene == "DRB1"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p4

ggvenn::ggvenn(
  list(KMHC = fd_withkgp[(fd_withkgp$Ref == "KMHC") & (fd_withkgp$Gene == "DQB1"),]$type,
       `1KGP` = fd_withkgp[(fd_withkgp$Ref == "1KGP") & (fd_withkgp$Gene == "DQB1"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p5


library(cowplot)

plot_grid(p1,p2,p3,p4,p5,labels = c("A","B","C","DRB1","DQB1"),
          ncol = 2)

table(kgp_type$Gene)
p1
p2
p3
p4
p5
#### HLA type venn diagram by gene

head(fd)
table(fd$Ref)
ggvenn::ggvenn(
  list(KMHC = fd[(fd$Ref == "KMHC") & (fd$Gene == "A"),]$type,
       `Multi-ethnic` = fd[(fd$Ref == "Multi-ethnic") & (fd$Gene == "A"),]$type,
       `Han Chinese` = fd[(fd$Ref == "Han Chinese") & (fd$Gene == "A"),]$type,
       `Pan-Kor` = fd[(fd$Ref == "Pan-Kor") & (fd$Gene == "A"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p1

p1
ggvenn::ggvenn(
  list(KMHC = fd[(fd$Ref == "KMHC") & (fd$Gene == "B"),]$type,
       `Multi-ethnic` = fd[(fd$Ref == "Multi-ethnic") & (fd$Gene == "B"),]$type,
  `Han Chinese` = fd[(fd$Ref == "Han Chinese") & (fd$Gene == "B"),]$type,
  `Pan-Kor` = fd[(fd$Ref == "Pan-Kor") & (fd$Gene == "B"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p2

ggvenn::ggvenn(
  list(KMHC = fd[(fd$Ref == "KMHC") & (fd$Gene == "C"),]$type,
       `Multi-ethnic` = fd[(fd$Ref == "Multi-ethnic") & (fd$Gene == "C"),]$type,
  `Han Chinese` = fd[(fd$Ref == "Han Chinese") & (fd$Gene == "C"),]$type,
  `Pan-Kor` = fd[(fd$Ref == "Pan-Kor") & (fd$Gene == "C"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p3

ggvenn::ggvenn(
  list(KMHC = fd[(fd$Ref == "KMHC") & (fd$Gene == "DPA1"),]$type,
       `Multi-ethnic` = fd[(fd$Ref == "Multi-ethnic") & (fd$Gene == "DPA1"),]$type,
  `Han Chinese` = fd[(fd$Ref == "Han Chinese") & (fd$Gene == "DPA1"),]$type,
  `Pan-Kor` = fd[(fd$Ref == "Pan-Kor") & (fd$Gene == "DPA1"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p4

ggvenn::ggvenn(
  list(KMHC = fd[(fd$Ref == "KMHC") & (fd$Gene == "DPB1"),]$type,
       `Multi-ethnic` = fd[(fd$Ref == "Multi-ethnic") & (fd$Gene == "DPB1"),]$type,
  `Han Chinese` = fd[(fd$Ref == "Han Chinese") & (fd$Gene == "DPB1"),]$type,
  `Pan-Kor` = fd[(fd$Ref == "Pan-Kor") & (fd$Gene == "DPB1"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p5

ggvenn::ggvenn(
  list(KMHC = fd[(fd$Ref == "KMHC") & (fd$Gene == "DQA1"),]$type,
       `Multi-ethnic` = fd[(fd$Ref == "Multi-ethnic") & (fd$Gene == "DQA1"),]$type,
  `Han Chinese` = fd[(fd$Ref == "Han Chinese") & (fd$Gene == "DQA1"),]$type,
  `Pan-Kor` = fd[(fd$Ref == "Pan-Kor") & (fd$Gene == "DQA1"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p6

ggvenn::ggvenn(
  list(KMHC = fd[(fd$Ref == "KMHC") & (fd$Gene == "DQB1"),]$type,
       `Multi-ethnic` = fd[(fd$Ref == "Multi-ethnic") & (fd$Gene == "DQB1"),]$type,
       `Han Chinese` = fd[(fd$Ref == "Han Chinese") & (fd$Gene == "DQB1"),]$type,
       `Pan-Kor` = fd[(fd$Ref == "Pan-Kor") & (fd$Gene == "DQB1"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p7

ggvenn::ggvenn(
  list(KMHC = fd[(fd$Ref == "KMHC") & (fd$Gene == "DRB1"),]$type,
       `Multi-ethnic` = fd[(fd$Ref == "Multi-ethnic") & (fd$Gene == "DRB1"),]$type,
  `Han Chinese` = fd[(fd$Ref == "Han Chinese") & (fd$Gene == "DRB1"),]$type,
  `Pan-Kor` = fd[(fd$Ref == "Pan-Kor") & (fd$Gene == "DRB1"),]$type),
  text_size = 3,
  set_name_size = 4
) -> p8

library(cowplot)

plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,labels = c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1"),
          ncol = 2)

plot_grid(p1,p2,labels = c("A","B"),
          ncol = 2)

plot_grid(p3,p8,labels = c("C","DRB1"),
          ncol = 2)

plot_grid(p4,p5,labels = c("DPA1","DPB1"),
          ncol = 2)

plot_grid(p6,p7,labels = c("DQA1","DQB1"),
          ncol = 2)


### vs Park
#setwd("~/Desktop/KCDC/HLA_seq/")
park <- readxl::read_xlsx("~/Desktop/KCDC/HLA_seq/HLAtype.freq_park2016.xlsx",sheet = 3)
park %>% filter(HLAgene != "HLA-DQB1") %>%
  mutate(type = value,"Park" = Frequency/100) %>% #head()
  filter(value %in% kmhc_freq$type ) %>%
  mutate(HLAgene = gsub(x = HLAgene, pattern = "HLA-", replacement = "")) -> park

#kmhc_freq %>% na.omit() %>%
   #mutate('type' = value) -> kmhc_freq
head(kmhc_freq)
head(park)
park$type = park$value
ggvenn::ggvenn(
  list(KMHC = kmhc_freq$type,
       Park = park$type),
  text_size = 5,
  set_name_size = 4
)




head(park)
head(kmhc_freq)

kmhc_freq %>% group_by(Gene) %>% 
  na.omit() %>% #head()
  count(type) %>%
  mutate("KMHC" = prop.table(n),"HLAgene" = gsub(x = Gene, pattern = "HLA_", replacement = "")) %>% #head()
  filter(type %in% park$type) -> kmhc_freq1


head(park)
head(kmhc_freq1)
kmhc_freq1 %>% select(HLAgene,type,KMHC)
park %>% select(HLAgene,type,Park) %>% #head()
  inner_join(kmhc_freq1 %>% select(HLAgene,type,KMHC)) %>% #head()
  ggscatter(.,x='KMHC',y='Park',color='HLAgene',
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            xlab = "KMHC",
            ylab = "Park (2016)",
            #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
            cor.coeff.args = list(method = "pearson", label.sep = "\n")) + 
  theme(legend.position = "right")


head(kmhc_freq)
head(han_freq)
head(pan_freq)

kmhc_freq %>% group_by(Gene) %>% filter(type != 0) %>%
  na.omit() %>% #head()
  count(type) %>%
  mutate("KMHC" = prop.table(n),"HLAgene" = gsub(x = Gene, pattern = "HLA_", replacement = "")) %>% 
  select(HLAgene,type,KMHC)-> kmhc_freq1

han_freq %>% group_by(Gene) %>% filter(type != 0) %>%
  na.omit() %>% #head()
  count(type) %>%
  mutate("Han Chinese" = prop.table(n),"HLAgene" = gsub(x = Gene, pattern = "HLA_", replacement = "")) %>%
  select(HLAgene,type,`Han Chinese`)-> han_freq1

pan_freq %>% group_by(Gene) %>% filter(type != 0) %>%
  na.omit() %>% #head()
  count(type) %>%
  mutate("Pan-Kor" = prop.table(n),"HLAgene" = gsub(x = Gene, pattern = "HLA_", replacement = "")) %>% 
  select(HLAgene,type,`Pan-Kor`)-> pan_freq1


head(kmhc_freq1)
head(han_freq1)
head(pan_freq1)

kmhc_freq1 %>% group_by() %>% 
  mutate(common = ifelse(KMHC >= 0.05,"common",ifelse(KMHC < 0.01,"rare","less_common"))) %>%
  count(common)

head()
head(michigan)
head(kmhc_type)
head(kgp_type)
kmhc_type %>% rbind(han_type) %>% rbind(pan_type) %>% rbind(michigan %>% filter(digit == "fd") %>% select(-digit)) %>% 
  rbind(kgp_type) %>%
  group_by(Ref) %>% filter(type != "0") %>%
  count(Gene) %>% #head()
  pivot_wider(names_from = Gene,values_from = n) -> fd_type

head(fd_withkgp)
head(td_withkgp)

head(fd)
fd %>% count(Ref,Gene) %>%
  pivot_wider(names_from = Gene,values_from = n)

fd_type
head(han_type)
han_type %>% filter(type == "0")
kmhc_type %>% filter(type == "0")

kmhc_type %>% rbind(han_type) %>% rbind(pan_type) %>% rbind(michigan %>% filter(digit == "fd") %>% select(-digit)) %>% 
  rbind(kgp_type) %>%
  group_by(Ref) %>% filter(type != "0") %>%
  count(Gene) %>% #head()
  pivot_wider(names_from = Gene,values_from = n) -> fd_type

kmhc_type %>% rbind(han_type) %>% rbind(pan_type) %>% rbind(michigan %>% filter(digit == "fd") %>% select(-digit)) %>% 
  rbind(kgp_type) %>% #head()
  mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique() %>%
  group_by(Ref) %>% filter(type != "0") %>%
  count(Gene) %>% 
  pivot_wider(names_from = Gene,values_from = n) -> td_type
  
fd_type$digit <- "4"
td_type$digit <- "2"

fd_type %>% rbind(td_type) %>% select(Ref,HLA_A,HLA_B,HLA_C,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,HLA_DRB1,digit) %>%
  writexl::write_xlsx("HLAtype.compareRef.IMGT3320.xlsx")
  #=CONCATENATE(B7,"/",B2)

kmhc_type %>% group_by(Gene) %>%
  count(Gene) #%>% sum(n)

kmhc_freq1 %>% inner_join(han_freq1) %>%
  ggscatter(.,x='KMHC',y="Han Chinese",color='HLAgene',
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            xlab = "KMHC",
            ylab = "Han Chinese",
            xlim = c(0,0.5),
            ylim = c(0,0.5),
            #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
            cor.coeff.args = list(method = "pearson", label.sep = "\n")) + 
  theme(legend.position = "right")

kmhc_freq1 %>% inner_join(han_freq1) %>% #head()
  correlate()


#### venn bar
df <- readxl::read_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/HLAtype_4panel.xlsx")
head(df)
table(df$Ref)
df %>% filter(digit=="Two-field") %>% select(Ref,Gene,type) %>% filter(Ref == "Multi-ethnic") -> multi
df %>% filter(digit=="Two-field") %>% select(Ref,Gene,type) %>% filter(Ref == "KMHC") -> kmhc
df %>% filter(digit=="Two-field") %>% select(Ref,Gene,type) %>% filter(Ref == "Han Chinese") -> han
df %>% filter(digit=="Two-field") %>% select(Ref,Gene,type) %>% filter(Ref == "Pan-Kor") -> pan


  
head(kmhc) 
head(multi)
kmhc %>% mutate(venn=ifelse(type %in% multi$type,"0",Ref)) -> kmhc1
multi %>% filter(!(type %in% kmhc$type)) %>% mutate(venn=Ref) %>% rbind(kmhc1) %>%
  mutate(Ref = "KMHC vs Multi") -> multi1

kmhc %>% mutate(venn=ifelse(type %in% pan$type,"0",Ref)) -> kmhc1
pan %>% filter(!(type %in% kmhc$type)) %>% mutate(venn=Ref) %>% rbind(kmhc1) %>%
  mutate(Ref = "KMHC vs Pan") -> pan1

kmhc %>% mutate(venn=ifelse(type %in% han$type,"0",Ref)) -> kmhc1
han %>% filter(!(type %in% kmhc$type)) %>% mutate(venn=Ref) %>% rbind(kmhc1) %>%
  mutate(Ref = "KMHC vs Han") -> han1


kmhc %>% mutate(venn=ifelse(type %in% multi$type,'1','0')) -> kmhc1
multi %>% filter(!(type %in% kmhc$type)) %>% mutate(venn='2') %>% rbind(kmhc1) %>%
  mutate(Ref = "KMHC vs Multi") -> multi1

kmhc %>% mutate(venn=ifelse(type %in% pan$type,'1','0')) -> kmhc1
pan %>% filter(!(type %in% kmhc$type)) %>% mutate(venn='3') %>% rbind(kmhc1) %>%
  mutate(Ref = "KMHC vs Pan") -> pan1

kmhc %>% mutate(venn=ifelse(type %in% han$type,'1','0')) -> kmhc1
han %>% filter(!(type %in% kmhc$type)) %>% mutate(venn='4') %>% rbind(kmhc1) %>%
  mutate(Ref = "KMHC vs Han") -> han1



library(ggplot2)
#library(HH)

multi1 %>% rbind(pan1) %>% rbind(han1) %>%
  ggplot(aes(y=Ref,x=Ref,fill=venn)) +
  geom_bar(position="stack", stat="identity") +
  theme_fivethirtyeight() + 
  coord_flip()


multi1 %>% rbind(pan1) %>% rbind(han1) %>% 
HH::likert(Ref., positive.order = FALSE, main="Female Responses to know how students react to maths anxiety.", scales = list(y = list(cex = .65)))

head(df)

df %>% filter(digit == "Two-field") %>% dplyr::select(Ref,type) %>% #head()
  pivot_wider(names_from = Ref,values_from = type) -> df1

multi1 %>% rbind(pan1) %>% rbind(han1) %>% #head()
  count(Ref,venn) -> df1
head(df1)  
#writexl::write_xlsx(df1,"~/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLAtype_check/HLAtype_forvennbar.xlsx")

df1 <- readxl::read_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLAtype_check/HLAtype_forvennbar.xlsx",na = "NA",sheet = 2)

head(df1)
df1 %>% select(type,KMHC,`Multi-ethnic`,Intersection,`Han Chinese`,`Pan-Kor`) %>% #dplyr::select(-venn) %>% #head()
#  pivot_wider(names_from = Ref,values_from = n,values_fill = NA) %>%
  HH::likert( ordered=FALSE)


df1 %>% pivot_longer(cols = 2:6,values_to = "count") %>% #head()
  na.omit() %>% ##head()
  mutate(name = fct_relevel(name,"KMHC","Intersection","Multi-ethnic","Han Chinese","Pan-Kor"),name=fct_rev(name)) %>% #count(type)
  #mutate(name =)
  #mutate(name = fct_relevel(name,"KMHC","Inter","Mulit","Han","Pan")) %>% #count(type)

  ggplot(aes(x=factor(type,levels = c("KMHC vs Pan-Kor","KMHC vs Han Chinese","KMHC vs Multi-ethnic")),y=count,fill=name)) +
  #ggplot(aes(x=type,y=count,fill=name)) +
  #ggplot(aes(x=factor()type,y=count,fill=name)) +
  geom_col() + 
  geom_text(aes(label=count,fontface="bold"),position = position_stack(vjust = 0.5),color="white") + 
  coord_flip() + 
  #scale_fill_discrete(breaks = c("KMHC","Inter","Mulit","Han","Pan")) +
  #scale_x_discrete(breaks = c("KMHC vs Han","KMHC vs Pan","KMHC vs Multi")) + 
  #scale_fill_manual(values=c("#F8766D", "grey","#7CAE00","#00BFC4","#C77CFF"),guide=guide_legend(c("KMHC","Inter","Mulit","Han","Pan"))) +
  scale_fill_manual(values=c("#C77CFF","#00BFC4","#7CAE00","grey","#F8766D"),guide=guide_legend(c("KMHC","Intersection","Multi-ethnic","Han Chinese","Pan-Kor"))) +
  #scale_fill_manual(values=c("#C77CFF","#00BFC4","#7CAE00","grey","#F8766D"),guide=guide_legend(c("KMHC","Inter","Mulit","Han","Pan"))) +
  #scale_fill_manual(values=c("#F8766D", "grey","#7CAE00","#00BFC4","#C77CFF",guide=guide_legend(reverse = TRUE)))+
#  scale_fill_viridis_d() + 
  labs(x=NULL) + 
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "none") -> a

df1 %>% pivot_longer(cols = 2:6,values_to = "count") %>% #head()
  





  #scale_fill_manual(values=c("#C77CFF","#00BFC4","#7CAE00","grey","#F8766D"),guide=guide_legend(c("KMHC","Inter","Mulit","Han","Pan"))) +





'''
df %>% mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>%
  na.omit() %>% #head()
  group_by(HLAgene) %>% #head()
  mutate(Frequency = prop.table(n)) %>% #writexl::write_xlsx("HLAtype.freq.xlsx")
  mutate(type = "KMHC") %>% select(-n) %>% filter(HLAgene %in% ref2$HLAgene)-> a
'''
  

'''
c2_ngs %>% select(RecInfo,DonInfo,28:44) %>% mutate("type" = "NGS") %>% #head()
  rbind(c2_imp) %>%
  pivot_longer(3:19,names_to = "theme",values_to = "count") %>% #head()
  pivot_wider(names_from = type,values_from = count) %>% 
  filter(!theme %in% c("Total_eps","All_ABV_ClassII","allDRB","allDQB","allDQA","allDPB","allDPA")) %>% #count(theme)
  #mutate(theme = factor(theme,levels = c("Total_eps","All_ABV_ClassII","allDRB","allDQB","allDQA","allDPB","allDPA"))) %>%
  #ggplot(aes(x=NGS,y=HLAimp,color=theme)) + 
  ggscatter(.,x='NGS',y='HLAimp',color='theme',
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            facet.by = "theme",
            xlab = "NGS based HLA typing",
            ylab = "HLA imputation",nrow =2,scales='free',
            #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
            cor.coeff.args = list(method = "pearson", label.sep = "\n")) + 
  theme(legend.position = "none",
        strip.text.x = element_text(size = 12,face = "bold"))
'''