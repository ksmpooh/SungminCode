## HLA type compore 2023 03 20
## michigan vs Han vs Pan+korea vs KBA
library(tidyverse)
library(stringr)
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



'''
           Ref  Gene     type digit
1 Multi-ethnic HLA_A     A*01    td
2 Multi-ethnic HLA_A  A*01:01    fd
3 Multi-ethnic HLA_A  A*01:02    fd

'''
table(michigan$Gene)
han <- read.table("Han.hg19.haplegendtovcf.modify.hlatype_fd.txt",header = T)
pan <- read.table("PanKor_merged.hg19.haplegendtovcf.modify.hlatype_fd.IDchange_fornomenclean.txt",header = F) %>%
  select(-V2,-V3,-V4,-V5,-V6)

kmhc <- read.table("HLA.4digit.520sample.nomenclean.chped") %>%
  select(-V2,-V3,-V4,-V5,-V6)

ref <- read.table('../Final_520sample.index.txt',header = T)
head(ref)
kmhc <- read.table("../HLA.type.result.8genes.merged.4digit_529sample_forMAKEreference.txt",header = T) %>%
  select(-IID,-pID,-mID,-SEX,-PHENO) %>% filter(FID %in% ref$KBAID)
head(kmhc)
head(han)
head(pan)
colnames(pan) <- colnames(han)
colnames(kmhc) <- colnames(han)
#pivot_longer(cols = c(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DRB1,HLA_DQA1,HLA_DQB1,overall),names_to = 'Gene',values_to = 'Accuracy') %>%
han %>% pivot_longer(2:ncol(han),names_to = 'Gene',values_to = 'type') %>% #head()
  mutate(Gene = str_replace_all(Gene,"\\.1","")) %>% mutate(Gene = str_replace_all(Gene,"\\.2","")) %>% #count(Gene)#head()
  mutate(Ref = "Han Chinese") %>% select(-ID) %>% unique()-> han_type

pan %>% pivot_longer(2:ncol(han),names_to = 'Gene',values_to = 'type') %>% #head()
  mutate(Gene = str_replace_all(Gene,"\\.1","")) %>% mutate(Gene = str_replace_all(Gene,"\\.2","")) %>% #count(Gene)#head()
  mutate(Ref = "Pan-Kor") %>% select(-ID) %>% unique() -> pan_type

kmhc %>% pivot_longer(2:ncol(han),names_to = 'Gene',values_to = 'type') %>% #head()
  mutate(Gene = str_replace_all(Gene,"\\.1","")) %>% mutate(Gene = str_replace_all(Gene,"\\.2","")) %>% #count(Gene)#head()
  mutate(Ref = "KMHC") %>% select(-ID) %>% unique() -> kmhc_type


head(han_type)
table(han_type$Gene)
michigan %>% filter(digit == 'fd') %>% select(-digit) %>%
  rbind(han_type) %>% rbind(pan_type) %>% rbind(kmhc_type) %>% mutate(Gene = str_replace_all(Gene,"HLA_","")) -> fd

head(fd)
fd %>% mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique() -> td
fd$digit <- "4 digit"
td$digit <- "2 digit"

fd %>% rbind(td) %>%
  count(Gene,Ref,digit) %>% #head()
  ggplot(aes(x=Gene,y=n,fill = factor(Ref,levels = c("KMHC","Multi-ethnic","Han Chinese","Pan-Kor")))) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  xlab("HLA gene") + ylab("Number of HLA type") + 
  scale_fill_discrete(name="Ref.Panel") +
  theme(legend.title=element_text(size=12,face = "bold"),
        legend.text=element_text(size=10),
        axis.title.x = element_text(size = 13,face = "bold"),
        axis.title.y = element_text(size = 13,face = "bold")) +
  facet_grid(~digit) + 
  theme(strip.text.x = element_text(size = 20,face = "bold"))


michigan %>% filter(digit == 'fd') %>% select(-digit) %>%
  rbind(han_type) %>% rbind(pan_type) %>% rbind(kmhc_type) %>% mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% head()

VennDiagram::venn.diagram(x=list())

x <- list(KMHC = kmhc_type$type,`Multi-ethnic` = michigan[michigan$digit =='fd',]$type,`Han Chinese`  = han_type$type,`Pan-Kor` = pan_type$type)
ggvenn::ggvenn(
  x,
  stroke_size = 0.5, set_name_size = 4
)
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
han_type
pan_type
michigan_type
kmhc_type
table(kmhc_type$type %in% michigan_type$type)

kmhc_type %>% mutate("count" = 1) %>%
  mutate(count = ifelse(type %in% michigan_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% han_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% pan_type$type,count+1,count)) %>% count(Ref,count) -> kmhc_type_count
  
michigan_type %>% mutate("count" = 1) %>%
  mutate(count = ifelse(type %in% kmhc_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% han_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% pan_type$type,count+1,count)) %>% count(Ref,count) -> michigan_type_count

han_type %>% mutate("count" = 1) %>%
  mutate(count = ifelse(type %in% kmhc_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% michigan_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% pan_type$type,count+1,count)) %>% count(Ref,count) -> han_type_count

pan_type %>% mutate("count" = 1) %>%
  mutate(count = ifelse(type %in% kmhc_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% han_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% michigan_type$type,count+1,count)) %>% count(Ref,count) -> pan_type_count

kmhc_type_count %>% rbind(michigan_type_count) %>% rbind(han_type_count) %>% rbind(pan_type_count) %>% #head()
  mutate(count = as.factor(count),"digit" = "4 digit") -> fd_count
  

kmhc_type %>% mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique() -> kmhc_type
kmhc_type %>%  mutate("count" = 1) %>%
  mutate(count = ifelse(type %in% michigan_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% han_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% pan_type$type,count+1,count)) %>% count(Ref,count) -> kmhc_type_count

michigan_type %>% mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique() -> michigan_type
michigan_type %>%  mutate("count" = 1) %>%
  mutate(count = ifelse(type %in% kmhc_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% han_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% pan_type$type,count+1,count)) %>% count(Ref,count) -> michigan_type_count

han_type %>% mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique() ->han_type
han_type %>%  mutate("count" = 1) %>%
  mutate(count = ifelse(type %in% kmhc_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% michigan_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% pan_type$type,count+1,count)) %>% count(Ref,count) -> han_type_count

pan_type %>% mutate(type = str_split_fixed(type,":",2)[,1]) %>% unique() ->pan_type
pan_type %>%mutate("count" = 1) %>%
  mutate(count = ifelse(type %in% kmhc_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% han_type$type,count+1,count)) %>% 
  mutate(count = ifelse(type %in% michigan_type$type,count+1,count)) %>% count(Ref,count) -> pan_type_count

kmhc_type_count %>% rbind(michigan_type_count) %>% rbind(han_type_count) %>% rbind(pan_type_count) %>% #head()
  mutate(count = as.factor(count),"digit" = "2 digit") -> td_count

#head(han_type_count)
fd_count %>% rbind(td_count) %>% group_by(Ref) %>% #head()
  mutate(per=paste0(n,"(",round(n/sum(n)*100, 0), "%)")) %>%
  mutate(per=ifelse(count==1,per,NA)) %>%
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
