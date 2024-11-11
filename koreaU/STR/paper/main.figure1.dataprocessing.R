#### figure 작업
library(tidyverse)
## QC info
load(file="~/Desktop/KU/@research/STR/Rdata/shortread.eh.66sample.qc.info.RData")
load(file="~/Desktop/KU/@research/STR/Rdata/longread.tgrt.66sample.qc.info.RData")
trgt <- tgrt
tgrt <- NULL

complex_str <- read.table("~/Desktop/KU/@research/STR/figure/complex_STR.txt",header = T)
dup_str <- read.table("~/Desktop/KU/@research/STR/figure/rm_dup_str_list.txt",header = T)

head(complex_str)
head(dup_str)

complex_str %>% filter(STR_ID != ID) -> complex_onlyEH_rmlist
head(complex_onlyEH_rmlist)

head(eh)
eh %>% select(STR_ID,FILTER,ID) %>% filter(!(STR_ID %in% complex_onlyEH_rmlist$STR_ID)) %>%
  filter(!(STR_ID %in% dup_str$ID)) %>% count(ID,FILTER) -> eh_count

eh_count %>% group_by(ID) %>% #head()
  mutate(prop = prop.table(n)*100) %>% filter(FILTER == "PASS") %>%
  select(ID,prop) %>% mutate(type = "Short-read") %>% write.table("figure1.eh_pass_prop.txt",row.names = F,col.names = T,quote = F,sep = "\t")


head(trgt)
trgt %>% select(ID,TRID,AP) %>% filter(!(TRID %in% dup_str$ID)) %>%
  filter(!(AP %in% c(".","NA"))) %>% #select(-AM) %>% #dim()
  mutate(AP = as.numeric(AP)) %>% filter(AP >=0.8) %>% count(ID) -> trgt_count2

trgt %>% filter(MC != 0) %>% select(ID,TRID,AP) %>% filter(!(TRID %in% dup_str$ID)) %>%
  filter(!(AP %in% c(".","NA"))) %>% #select(-AM) %>% #dim()
  mutate(AP = as.numeric(AP)) -> trgt_count


head(trgt_count)
trgt_count  %>% group_by(ID) %>% summarise(AP_mean=mean(AP)) %>%
  mutate(type = "Long-read") %>% write.table("figure1.trgt_ap.mean_prop.txt",row.names = F,col.names = T,quote = F,sep = "\t")


trgt %>% select(TRID) %>% unique() -> trgt_type
dim(trgt_type) * 2
321298 * 2
trgt_count %>% filter(TRID %in% complex_str$ID)

trgt_count2 %>% mutate(prop=n/(321298*2)*100) %>%
  select(ID,prop) %>% mutate(type = "Long-read") %>% write.table("figure1.trgt_pass_prop.txt",row.names = F,col.names = T,quote = F,sep = "\t")


#trgt_pass_prop <- read.table("figure1.trgt_pass_prop.txt",header = T)
trgt_pass_prop <- read.table("figure1.trgt_ap.mean_prop.txt",header = T)
eh_pass_prop <- read.table("figure1.eh_pass_prop.txt",header = T)
head(trgt_pass_prop)
trgt_pass_prop %>% rbind(eh_pass_prop) %>% #head()
  ggplot(aes(x=type,y=prop,fill=type)) + 
  geom_violin() + 
#  theme_bw() + 
  theme_step1()



c("#F8766D", '#00BA38', '#619CFF')

eh_pass_prop %>% #head()
  ggplot(aes(x=type,y=prop,fill=type)) + 
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black") + 
  #  theme_bw() + 
  theme_step1() +
 # ylim(c(0.99,1)) + 
  
  labs(y="Proportion of PASS") + 
  scale_fill_manual(values = c("Short-read" = "#619CFF")) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) 


trgt_pass_prop %>% 
  ggplot(aes(y=AP_mean,x=type,fill="#F8766D")) + 
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black") + 
  theme_step1() + 
#  ylim(c(0.99,1)) + 
  labs(y="mean of AP") + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())


## venn
library(tidyverse)

setwd("~/Desktop/KU/@research/STR/figure/")
eh_ID <- read.table("~/Desktop/KU/@research/STR/figure/eh_pass_QC_STRlist_without_complexSTRotherID.txt",header = T)
trgt_ID <- read.table("~/Desktop/KU/@research/STR/figure/trgt_AP08_QC_STRlist.txt",header = T)
complex_str <- read.table("~/Desktop/KU/@research/STR/figure/complex_STR.txt",header = T)
dup_str <- read.table("~/Desktop/KU/@research/STR/figure/rm_dup_str_list.txt",header = T)

head(dup_str)
head(complex_str)
eh_ID %>% head()

complex_str %>% filter(STR_ID != ID) -> complex_onlyEH_rmlist

head(complex_onlyEH_rmlist)

ru <- read_table("~/Desktop/KU/@research/STR/eh.v5_w_gangstr.v13.polymorphic.JSONtoBED.bed",col_names = F)
ru %>% mutate(MOTIFS = str_split_fixed(X4,";",3)[,2]) %>% 
  mutate(MOTIFS = str_split_fixed(MOTIFS,"=",2)[,2]) %>% #head()
  mutate(ID = str_split_fixed(X4,";",3)[,1]) %>% 
  mutate(ID = str_split_fixed(ID,"=",2)[,2]) %>% #head()
  select(X1,X2,X3,MOTIFS,ID) %>% filter(!(ID %in%dup_str$ID)) -> ru


head(eh_ID)
head(ru)
dim(eh_ID) #314192
eh_ID %>% filter(STR_ID %in% ru$ID) %>% dim()  #314181
eh_ID %>% filter(!(STR_ID %in%dup_str$ID)) %>% dim() #314181
eh_ID %>% filter(!(STR_ID %in%dup_str$ID)) -> eh_ID

complex_str %>% select(ID) %>% unique() %>% rename(STR_ID = ID)
trgt_ID %>% dim() #306264
rbind(trgt_ID,complex_str %>% select(ID) %>% unique() %>% rename(STR_ID = ID)) %>% dim() #306270
rbind(trgt_ID,complex_str %>% select(ID) %>% unique() %>% rename(STR_ID = ID)) -> trgt_ID
trgt_ID %>% filter(!(STR_ID %in%dup_str$ID)) %>% dim() #306260
trgt_ID %>% filter(STR_ID %in%dup_str$ID) %>% dim() #10
trgt_ID %>% filter(!(STR_ID %in%dup_str$ID)) -> trgt_ID

head(trgt_ID)
table(eh_ID$STR_ID %in% trgt_ID$STR_ID)
#FALSE   TRUE 
# 11530 302662 
table(trgt_ID$STR_ID %in% eh_ID$STR_ID)
#FALSE   TRUE 
#3598 302662

trgt_ID %>% rbind(eh_ID) %>% unique() -> all_pass
table(all_pass$STR_ID %in% ru$ID)
#TRUE 
#317779 
table(ru$ID %in% all_pass$STR_ID)
#FALSE   TRUE 
#3519 317779 

library(VennDiagram) 
venn_data <- list(
  ru = ru$ID,
  trgt = trgt_ID$STR_ID,
  eh = eh_ID$STR_ID
)

venn.diagram(
  x = venn_data,
  category.names = c("DB", "Long-read", "Short-read"), # 각 집합에 이름을 붙입니다.
  filename = "test.png", # 파일 이름
  output = TRUE, # 결과를 화면에도 출력할지 결정
  imagetype = "png", # 이미지 파일 형식
  height = 480, # 이미지 높이
  width = 480, # 이미지 너비
  resolution = 300, # 이미지 해상도
  compression = "lzw" # 압축 방법
)

library(ggplot2)
library(ggvenn)
ggvenn(venn_data,c("DB", "Long-read", "Short-read"),
       auto_scale = T)

VennDiagram <- venn.diagram(x = venn_data, filename = NULL,
                            category.names = c("STR DB", "Long-read", "Short-read"),
                            col = c("#440154AA", "#21908DAA", "#FDE725AA"), # 원의 색상을 반투명하게 설정
                            fill = c("#440154AA", "#21908DAA", "#FDE725AA"), # 원 내부 색상도 반투명하게 설정
                            alpha = 0.3, # 원의 반투명도를 설정
                            cex = 1.9, # 카테고리 레이블 크기
                            cat.cex = 1.9, # 집합 이름의 글자 크기
                            cat.fontface = "bold", # 글자 폰트 스타일 설정
                            cat.fontfamily = "Arial")
cowplot::plot_grid(VennDiagram)


head(ru)
head(trgt_ID)
head(eh_ID)
#colnames(trgt_ID) <- "TRGT"
#colnames(eh_ID) <- "EH"
trgt_ID$TRGT <- trgt_ID$STR_ID
eh_ID$EH <- eh_ID$STR_ID


head(ru)
#ru %>% select(ID) %>% rename(STR_ID = ID) %>% left_join(trgt_ID) %>% left_join(eh_ID) %>%
  #rename(STR_DB = STR_ID) %>%
  #write.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_STR_ID_forVenn.txt",row.names = F,col.names = T,quote = F,sep = "\t")


####### STR unit check

library(tidyverse)

setwd("~/Desktop/KU/@research/STR/figure/")
eh_ID <- read.table("~/Desktop/KU/@research/STR/figure/eh_pass_QC_STRlist_without_complexSTRotherID.txt",header = T)
trgt_ID <- read.table("~/Desktop/KU/@research/STR/figure/trgt_AP08_QC_STRlist.txt",header = T)
complex_str <- read.table("~/Desktop/KU/@research/STR/figure/complex_STR.txt",header = T)
dup_str <- read.table("~/Desktop/KU/@research/STR/figure/rm_dup_str_list.txt",header = T)
head(dup_str)
ru <- read_table("~/Desktop/KU/@research/STR/eh.v5_w_gangstr.v13.polymorphic.JSONtoBED.bed",col_names = F)
head(ru)
ru %>% mutate(MOTIFS = str_split_fixed(X4,";",3)[,2]) %>% 
  mutate(MOTIFS = str_split_fixed(MOTIFS,"=",2)[,2]) %>% #head()
  mutate(ID = str_split_fixed(X4,";",3)[,1]) %>% 
  mutate(ID = str_split_fixed(ID,"=",2)[,2]) %>% #head()
  select(X1,X2,X3,MOTIFS,ID) -> ru

head(ru)
#head(common_str)
colnames(ru) <- c("chrom","start","end","MOTIFS","STR")
head(ru)
ru %>% mutate(RU.length=str_length(MOTIFS)) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) %>% write.table("STR_RU.length_GC.info.txt",col.names = T,row.names = F,quote = F,sep = "\t")
  
#ru %>% mutate(STR_region = paste0(chrom,"_",start,"_",end)) %>% filter(STR %in% common_str$TRID)-> ru#ref <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.pass.IDregion.txt")
#dim(ru)

venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_STR_ID_forVenn.txt",header = T)
head(ru);dim(ru)
#ref %>% filter
head(venn_rawdata)
dim(venn_rawdata)
venn_rawdata %>% filter(is.na(EH)) %>% filter(!is.na(TRGT))
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR

ru %>% filter(STR %in% common_STR$STR_DB) %>% select(MOTIFS) -> common_motif

common_motif %>% mutate(RU.length = ifelse(str_detect(MOTIFS,","),"CRU",str_length(MOTIFS))) %>%
  count(RU.length) -> length_common_motif

#write.table(length_common_motif,"figure1.length_common_motif.count.txt",col.names = T,row.names = F,quote = F,sep = "\t")
ru.scale = c(seq(2,14),"15+")
ru.scale
#table(length_common_motif$length)

length_common_motif %>% 
  mutate(RU.length = ifelse(RU.length %in% ru.scale,RU.length,"15+")) %>%
  ggplot(aes(x=factor(RU.length,levels=ru.scale),y=n)) + 
  geom_bar(stat = 'identity') + 
  theme_step1() + 
  scale_y_continuous(labels = scales::comma) + 
  labs(x="Length of Repeat Units",y="# of STRs") -> p1
p1

length_common_motif %>% #filter(is.na(length))
  mutate(RU.length = ifelse(RU.length %in% ru.scale,RU.length,"15+")) %>%
  ggplot(aes(x=factor(RU.length,levels=ru.scale),y=n)) + 
  geom_bar(stat = 'identity') + 
  theme_step1() + 
  scale_y_continuous(labels = scales::comma) + 
  xlim(c(seq(7,14),"15+")) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())-> p2
p2

combined_plot <- ggdraw() +
  draw_plot(p1, 0, 0, 1, 1) +           # p1을 원본 크기로 배치
  draw_plot(p2, 0.4, 0.4, 0.6, 0.6)     # p2를 오른쪽 위에 작게 배치

combined_plot

length_common_motif %>% mutate(pct = prop.table(n)) %>% #head()
  mutate(RU.length = factor(RU.length,levels=c(seq(2:30),"other"))) %>% arrange(RU.length) -> a


########## STR unit count

setwd("~/Desktop/KU/@research/STR/figure/")
complex_str <- read.table("~/Desktop/KU/@research/STR/figure/complex_STR.txt",header = T)
dup_str <- read.table("~/Desktop/KU/@research/STR/figure/rm_dup_str_list.txt",header = T)
head(dup_str)
ru <- read_table("~/Desktop/KU/@research/STR/eh.v5_w_gangstr.v13.polymorphic.JSONtoBED.bed",col_names = F)
head(ru)
ru %>% mutate(MOTIFS = str_split_fixed(X4,";",3)[,2]) %>% 
  mutate(MOTIFS = str_split_fixed(MOTIFS,"=",2)[,2]) %>% #head()
  mutate(ID = str_split_fixed(X4,";",3)[,1]) %>% 
  mutate(ID = str_split_fixed(ID,"=",2)[,2]) %>% #head()
  select(X1,X2,X3,MOTIFS,ID) -> ru

head(ru)
colnames(ru) <- c("chrom","start","end","MOTIFS","STR")

venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_STR_ID_forVenn.txt",header = T)
head(ru);dim(ru)
#ref %>% filter
head(venn_rawdata)
dim(venn_rawdata)
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR

ru %>% filter(STR %in% common_STR$STR_DB) %>% select(MOTIFS) -> common_motif

#common_motif %>% count(MOTIFS) %>% arrange(-n) %>% write.table("figure1.frequency_common_motif.txt",col.names = T,row.names = F,quote = F,sep = "\t")
common_motif %>% count(MOTIFS) %>% arrange(-n) %>% mutate(Rank = rank(-n)) %>% 
  mutate(new_n = ifelse(n %in% c(1:10),n,"11+")) %>% 
  count(new_n) -> common_motif_count_frequency

#write.table(common_motif_count_frequency,"figure1.frequency_common_motif.count.txt",col.names = T,row.names = F,quote = F,sep = "\t")

common_motif %>% count(MOTIFS) %>% arrange(-n) %>% mutate(Rank = rank(-n)) %>% 
  mutate(pct = prop.table(n)) %>%
  filter(Rank %in% c(1:20)) -> a

head(common_motif_count_frequency)
common_motif_count_frequency %>% mutate(pct = prop.table(n)) %>% #head()
  mutate(new_n = factor(new_n,levels=c(seq(1:10),"11+"))) %>% arrange(new_n) -> a

dim(venn_rawdata$STR_DB)
dim(venn_rawdata$TRGT)

venn_rawdata %>% select(TRGT) %>% na.omit() %>% dim() #306260
venn_rawdata %>% select(EH) %>% na.omit() %>% dim() #314181

  
  
common_motif %>% count(MOTIFS) %>% arrange(-n) %>% mutate(Rank = rank(-n)) %>% 
  mutate(new_n = ifelse(n %in% c(1:10),n,"11+")) %>% 
  count(new_n) %>% #head()
  ggplot(aes(x=factor(new_n,c("1","2","3","4","5","6","7","8","9","10","11+")),y=n)) + 
  geom_bar(stat='identity') + 
  labs(x="Frequency of Repeat Units",y="# of STRs") +
  theme_step1() + 
  scale_y_break(c(1000, 4800),ticklabels = c(1,2,3,4900,5000)) + 
  theme(legend.position = "none",
        axis.title.y = element_blank())
  
common_motif %>% count(MOTIFS) %>% arrange(-n) %>% mutate(Rank = rank(-n)) %>% 
  mutate(new_n = ifelse(n %in% c(1:10),n,"11+")) %>% 
  count(new_n) %>% #head()
  ggplot(aes(x=factor(new_n,c("1","2","3","4","5","6","7","8","9","10","11+")),y=n)) + 
  geom_bar(stat='identity') + 
  labs(x="Frequency of Repeat Units",y="# of STRs") +
  theme_step1() + 
  xlim(c("5","6","7","8","9","10")) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
  



p1

head(ru)
ru %>%
  select(STR,MOTIFS) %>% count(MOTIFS) %>% arrange(-n) %>% mutate(Rank = rank(-n)) -> str_RU_count

trans <- function(x){pmin(x,500) + 0.05*pmax(x-500,0)}
yticks <- c(0,100,200,400,600,800,4900,5000)

str_RU_count %>% mutate(new_n = ifelse(n %in% c(1:10),n,">10")) %>% 
  #group_by(RU) %>%
  count(new_n) %>% #head()
  ggplot(aes(x=factor(new_n,c("1","2","3","4","5","6","7","8","9","10",">10")),y=n,fill=new_n)) + 
  geom_bar(stat='identity') +
  scale_y_break(c(1000, 4700),ticklabels = c(1,2,3,4900,5000)) + 
  
  labs(x="# of RU types",y="# of STRs")  + 
  theme(#axis.title.x = element_blank(),
    legend.position = 'None',
    #axis.y.
    axis.text = element_text(size = 10)) -> p2
p2



