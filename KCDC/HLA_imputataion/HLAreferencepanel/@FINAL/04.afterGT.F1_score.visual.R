## F1 score
library(tidyverse)
library(cowplot)
setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/Final/F1_score/common")
ref_2 <- read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLAtype_check/KMHC_HLAtype_frequency.txt",header = T) %>% select(value,prop,freq)
head(ref_2)
colnames(ref_2)[1]<-c("HLAtype")
ref <-read.table("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLAtype_check/other_panel/hlatype/intersect.HLAtype.in4Ref.txt")
colnames(ref) <- "Gene"



flist = grep(list.files("./"),pattern = "after.GTmatchin", value=TRUE)
flist
df <- NULL
for (i in flist) {
  a <- read.table(i,header = T)
  if (grepl("michigan",i)) {
    a$Ref <- "Multi-ethnic"  
  } else if(grepl("PanKor",i)){
    a$Ref <- "Pan-Kor"  
  }else if(grepl("Han",i)){
    a$Ref <- "Han Chinese"  
  }else{
    a$Ref <- "KMHC"
  }
  
  a$filename <- i
  df <- rbind(df,a)
}
head(df)
df %>% mutate(CV = str_split_fixed(filename,"_",2)[,1]) %>% mutate(HLAtype = str_split_fixed(ID,"_",2)[,2],Gene = str_split_fixed(HLAtype,"\\*",2)[,1] ) %>% #head()
  mutate(F1_score = TP/(TP+(FP+FN)/2)) %>% select(Gene,HLAtype,Ref,CV,F1_score) %>% #head()
  filter(HLAtype %in% ref$Gene) %>% #inner_join(ref_2) %>%
  group_by(Ref,CV) %>% na.omit() %>%
  summarise(mean = mean(F1_score)) -> overall
head(overall)
overall$Gene <- "Overall"
colnames(overall)[1] <- "name"

df %>% mutate(CV = str_split_fixed(filename,"_",2)[,1]) %>% mutate(HLAtype = str_split_fixed(ID,"_",2)[,2],Gene = str_split_fixed(HLAtype,"\\*",2)[,1] ) %>% #head()
  mutate(F1_score = TP/(TP+(FP+FN)/2)) %>% select(Gene,HLAtype,Ref,CV,F1_score) %>% #head()
  filter(HLAtype %in% ref$Gene) %>% #inner_join(ref_2) %>%
  group_by(Gene,HLAtype,CV) %>%
  #pivot_wider(names_from = Ref,values_from = F1_score) %>% #na.omit() %>% nrow()
  pivot_wider(names_from = Ref,values_from = F1_score) %>% 
  pivot_longer(4:7) %>% #head()
  inner_join(ref_2) %>% #ungroup %>% count(name)
  filter(freq == "less_common") %>% na.omit() %>%
  group_by(name,Gene,CV) %>% 
  summarise(mean = mean(value)) %>% #head()
  rbind(overall) %>%
  mutate(name = factor(name, levels = c("KMHC","Multi-ethnic","Han Chinese", "Pan-Kor"), labels = c("KMHC", "Multi-ethnic","Han Chinese", "Pan-Kor"))) %>% 
  ggplot(aes(x=name,y=mean,fill=name)) +
  geom_boxplot() + 
  facet_grid(~Gene)

  

df %>% mutate(CV = str_split_fixed(filename,"_",2)[,1]) %>% mutate(HLAtype = str_split_fixed(ID,"_",2)[,2],Gene = str_split_fixed(HLAtype,"\\*",2)[,1] ) %>% #head()
  mutate(F1_score = TP/(TP+(FP+FN)/2)) %>% 
  select(Gene,HLAtype,Ref,CV,F1_score) %>% #head()
  filter(HLAtype %in% ref$Gene) %>% #inner_join(ref_2) %>%
  group_by(Gene,HLAtype,CV) %>% #head()
  pivot_wider(names_from = Ref,values_from = F1_score) %>% na.omit() %>% #head()
  pivot_longer(4:7) %>% 
  inner_join(ref_2) -> new_df
  
hla_gene <- new_df %>% ungroup %>% select(HLAtype,prop,freq) %>% unique()
hla_gene
setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/Final/F1_score/common/plot")

head(new_df)
new_df %>% filter(HLAtype == hla_gene[1,"HLAtype"]) %>% 
  mutate(name = factor(name, levels = c("KMHC","Multi-ethnic","Han Chinese", "Pan-Kor"), labels = c("KMHC", "Multi-ethnic","Han Chinese", "Pan-Kor"))) %>% 
  ggplot(aes(x=name,y=value,color=name)) +
  geom_point() +
  ggtitle(paste0(hla_gene[1,]$HLAtype,":",hla_gene[1,]$freq,",",hla_gene[1,]$prop)) + 
  theme(legend.position = "bottom")
  
hla_gene[1,"HLAtype"]
for (i in 1:nrow(hla_gene)) {
  new_df %>% filter(HLAtype == hla_gene[i,"HLAtype"]) %>% 
    ggplot(aes(x=name,y=value,color=name)) +
    geom_point() +
    ggtitle(paste0(hla_gene[i,]$HLAtype,":",hla_gene[i,]$freq,",",hla_gene[i,]$prop)) + 
    theme(legend.position = "bottom") -> a
  png(paste0(hla_gene[i,]$freq,":",hla_gene[i,]$HLAtype,".png"))
  plot(a)
  dev.off()
  
}

###
df %>% mutate(CV = str_split_fixed(filename,"_",2)[,1]) %>% mutate(HLAtype = str_split_fixed(ID,"_",2)[,2],Gene = str_split_fixed(HLAtype,"\\*",2)[,1] ) %>% #head()
  mutate(F1_score = TP/(TP+(FP+FN)/2)) %>% 
  select(Gene,HLAtype,Ref,CV,F1_score) %>% #head()
  filter(HLAtype %in% ref$Gene) %>% #inner_join(ref_2) %>%
  group_by(Gene,HLAtype,CV) %>% #head()
  mutate(Ref = paste0("F1_score_",Ref)) %>%
  pivot_wider(names_from = Ref,values_from = F1_score) %>% #na.omit() %>% #head()
  inner_join(ref_2) -> a
#pivot_longer(4:7) %>% 

df %>% mutate(CV = str_split_fixed(filename,"_",2)[,1]) %>% mutate(HLAtype = str_split_fixed(ID,"_",2)[,2],Gene = str_split_fixed(HLAtype,"\\*",2)[,1] ) %>% #head()
  mutate(F1_score = TP/(TP+(FP+FN)/2)) %>% 
  select(Gene,HLAtype,Ref,CV,TP) %>% #head()
  filter(HLAtype %in% ref$Gene) %>% #inner_join(ref_2) %>%
  group_by(Gene,HLAtype,CV) %>% #head()
  mutate(Ref = paste0("TP_",Ref)) %>%
  pivot_wider(names_from = Ref,values_from = TP) %>% #na.omit() %>% #head()
  inner_join(ref_2) -> b

head(a)
head(b)


a %>% merge(b) %>% #colnames()
  select(CV,Gene,HLAtype,prop,freq,TP_KMHC,F1_score_KMHC,'TP_Multi-ethinic','F1_score_Multi-ethinic','TP_Han Chinese','F1_score_Han Chinese',
         'TP_Pan-Kor','F1_score_Pan-Kor') %>% writexl::write_xlsx("../HLAtype.F1_score.TP.HLAimputation.4ref.compore.xlsx")



df %>% mutate(CV = str_split_fixed(filename,"_",2)[,1]) %>% mutate(HLAtype = str_split_fixed(ID,"_",2)[,2],Gene = str_split_fixed(HLAtype,"\\*",2)[,1] ) %>% #head()
  mutate(F1_score = TP/(TP+(FP+FN)/2)) %>% 
  select(Gene,HLAtype,Ref,CV,F1_score) %>% #head()
  filter(HLAtype %in% ref$Gene) %>% #inner_join(ref_2) %>%
  group_by(Gene,HLAtype,CV) -> a

df %>% mutate(CV = str_split_fixed(filename,"_",2)[,1]) %>% mutate(HLAtype = str_split_fixed(ID,"_",2)[,2],Gene = str_split_fixed(HLAtype,"\\*",2)[,1] ) %>% #head()
  mutate(F1_score = TP/(TP+(FP+FN)/2)) %>% 
  select(Gene,HLAtype,Ref,CV,TP) %>% #head()
  filter(HLAtype %in% ref$Gene) %>% #inner_join(ref_2) %>%
  group_by(Gene,HLAtype,CV) -> b
head(a)
head(b)

hla_gene <- new_df %>% ungroup %>% select(HLAtype,prop,freq) %>% unique()
hla_gene
setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/Final/F1_score/common/plot")
a %>% merge(b) %>%  left_join(hla_gene) -> c
head(c)
#c %>% left_join(hla_gene) %>% filter(is.na())
for (i in 1:nrow(hla_gene)) {
  c %>% filter(HLAtype == hla_gene[i,]$HLAtype) %>% na.omit() %>%
    ggplot(aes(x=Ref,y=F1_score,color=CV)) +
    geom_point() +
    ggtitle("F1_score") + 
    theme(legend.position = "bottom") -> p1
  c %>% filter(HLAtype == hla_gene[i,]$HLAtype) %>%  na.omit() %>%
    ggplot(aes(x=Ref,y=TP,color=CV)) +
    geom_point() +
    ggtitle("TP") +
    #ggtitle(paste0("CV\t",hla_gene[i,]$HLAtype,"\t",hla_gene[i,]$freq,",",hla_gene[i,]$prop)) + 
    theme(legend.position = "bottom") -> p2
  plot_row <- plot_grid(p1, p2)
  print(hla_gene[i,]$HLAtype)
  title <- ggdraw() + 
    draw_label(
      paste0(hla_gene[i,]$HLAtype,"\t",hla_gene[i,]$freq,",",hla_gene[i,]$prop),
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  plot_grid(
    title, plot_row,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  ) -> p

  png(paste0(hla_gene[i,]$freq,"_",round(hla_gene[i,]$prop,3),"_",hla_gene[i,]$HLAtype,".png"))
  plot(p)
  dev.off()
}

plot_row <- plot_grid(p1, p2)
hla_gene
# now add the title
c %>% filter(HLAtype == hla_gene[i,]$HLAtype) %>%  #head()
  ggplot(aes(x=Ref,y=F1_score,color=CV)) +
  geom_point() +
  #geom_boxplot()+
  ggtitle("F1_score") + 
  theme(legend.position = "bottom")
c %>% filter(HLAtype == hla_gene[i,"HLAtype"]) %>% 
  ggplot(aes(x=Ref,y=TP,color=CV)) +
  geom_point() +
  ggtitle("TP") +
  #ggtitle(paste0("CV\t",hla_gene[i,]$HLAtype,"\t",hla_gene[i,]$freq,",",hla_gene[i,]$prop)) + 
  theme(legend.position = "bottom") -> p2
plot_row <- plot_grid(p1, p2)
p1
