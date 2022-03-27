## ggplot usnig TCGA data
library(tidyverse)
library(cowplot)
setwd("~/Desktop/KU/2021_Fall/biovisual/ggplot_data/")

load('data.TCGA_LUAD_LUSC.gene_expression.Rdata')
# 환자 정보
colnames(d_luad) #adeno carcinoma
colnames(d_lusc) #scar # 담배 많이 핀사람 암종
# expression data
colnames(e_luad) 
colnames(e_lusc)

d_luad$age_at_diagnosis #일수로 되어 있음
hist(d_luad$age_at_diagnosis/365)

## lung cancer의 ad by gender?
## 성별에 따른 진단 나이

ggplot(d_luad,aes(x = gender, y = age_at_diagnosis/365)) +
  geom_boxplot()

ggplot(d_luad,aes(x = paper_Smoking.Status, y = age_at_diagnosis/365)) +
  geom_boxplot() + coord_flip()


ggplot(d_luad,aes(x = paper_Smoking.Status, y = age_at_diagnosis/365)) +
  geom_boxplot() + facet_wrap(~gender) +
#  scale_y_discrete(br"Current reformed smoker for > 15 years","Current reformed smoker for < or = 15 years") + 
  coord_flip()


table(d_luad$paper_Smoking.Status)


order_smoker = c('Current reformed smoker for > 15 years', 
                 'Lifelong Non-smoker', 
                 'Current reformed smoker for < or = 15 years',
                 'Current smoker')

d_luad %>% 
  filter(!is.na(paper_Smoking.Status) & paper_Smoking.Status %in% order_smoker) %>%
  mutate(paper_Smoking.Status = factor(paper_Smoking.Status, levels=order_smoker)) %>%
  ggplot(data = ., aes(paper_Smoking.Status, age_at_diagnosis/365)) + 
  geom_boxplot() + 
  facet_wrap(~gender) + # gender
  coord_flip() # y,x 변경


d_lusc %>%  select(paper_Smoking.Status) %>%
  table()

p1 <- d_luad %>% 
  filter(!is.na(paper_Smoking.Status) & paper_Smoking.Status %in% order_smoker) %>%
  mutate(paper_Smoking.Status = factor(paper_Smoking.Status, levels=order_smoker)) %>%
  ggplot(data = ., aes(paper_Smoking.Status, age_at_diagnosis/365)) + 
  geom_boxplot(aes(fill=paper_Smoking.Status)) + facet_wrap(~gender) + 
  geom_jitter() + 
  coord_flip() + 
  scale_fill_discrete(breaks = rev(order_smoker)) + 
  labs(x='', y='Age at Diagnosis (years)', 
       title='Smoking status of Lung cancer adenocarcinoma',
       fill='Smoking status')

p2 <- d_lusc %>% 
  filter(!is.na(paper_Smoking.Status) & paper_Smoking.Status %in% order_smoker) %>%
  mutate(paper_Smoking.Status = factor(paper_Smoking.Status, levels=order_smoker)) %>%
  ggplot(data = ., aes(paper_Smoking.Status, age_at_diagnosis/365)) + 
  geom_boxplot(aes(fill=paper_Smoking.Status)) + 
  facet_wrap(~gender) + # gender
  geom_jitter() + # 점이 추가 됨
  scale_fill_discrete(breaks = rev(order_smoker)) + # breaks는 숫서 
  coord_flip() + # y,x 변경 
  labs(x = ' ',y= 'Age at Diagnosis (years)',
       title = "Smoking history of Lung squamous cell carcinoma",
       fill = 'Smoking status')


p <- plot_grid(p1,p2,ncol = 1,labels = c("A","B"))

ggsave("plot.lung.cancer_smoking.pdf",p,height = 10, width = 16)


# violin -> density를 같이 표현이 가능함
d_lusc %>% 
  filter(!is.na(paper_Smoking.Status) & paper_Smoking.Status %in% order_smoker) %>%
  mutate(paper_Smoking.Status = factor(paper_Smoking.Status, levels=order_smoker)) %>%
  ggplot(data = ., aes(paper_Smoking.Status, age_at_diagnosis/365)) + 
  geom_violin(aes(fill=paper_Smoking.Status)) + 
  facet_wrap(~gender) + # gender
  geom_jitter() + # 점이 추가 됨
  scale_fill_discrete(breaks = rev(order_smoker)) + # breaks는 숫서 
  coord_flip() + # y,x 변경 
  labs(x = ' ',y= 'Age at Diagnosis (years)',
       title = "Smoking history of Lung squamous cell carcinoma",
       fill = 'Smoking status')


##

table(d_luad$ajcc_pathologic_stage)
d_luad$age_at_diagnosis
d_luad$years_smoked
## step1) Correlation between years of smoking and Age at diagnosis
## step2) By ajcc_pathologic_stage


d_luad %>% select(ajcc_pathologic_stage,age_at_diagnosis,years_smoked) %>%
  ggplot(data=.,aes(age_at_diagnosis/365,years_smoked)) +
  geom_point() + 
  stat_smooth(method = 'lm') +
  facet_wrap(~ajcc_pathologic_stage)
  
####gene expression

genes = c('ENSG00000008118', 'ENSG00000008128')
names(genes) = c('CAMK1G', 'CDK11A')
#'CAMK1G' 
#'CDK11A'

head(e_luad)

e_luad %>% as.data.frame() %>% 
  rownames_to_column('gene_id')

e_luad %>% as.data.frame() %>% 
  rownames_to_column('gene_id') %>%
  filter(gene_id %in% genes) %>% 
  gather(barcode, expr, -gene_id) ## convert to long format


dd = e_luad %>% as.data.frame() %>% rownames_to_column('gene_id') %>%
  filter(gene_id %in% genes) %>% gather(barcode, expr, -gene_id) ## -gene_id 제거를 위해 - 표시

dd1 = merge(d_luad, dd, by="barcode")

dd1 %>% ggplot(aes(gender,expr)) + geom_boxplot()


## row = geneID , col sample
dd1 = dd %>% pivot_wider(names_from = gene_id,values_from = expr)
head(dd1)
dd1 %>% ggplot(aes(log2(ENSG00000008118),log2(ENSG00000008128))) +
  geom_point() 


## correlation between two genes by ajcc_pathologic_stage

head(dd1)
head(d_luad)
d <- merge(as.data.frame(dd1),d_luad %>% select(barcode,ajcc_pathologic_stage),by='barcode')
d %>% ggplot(aes(log2(ENSG00000008118),log2(ENSG00000008128))) + 
  geom_point() + 
  facet_wrap(~ajcc_pathologic_stage) + 
  stat_smooth(method = 'lm') + 
  labs(x = 'CAMK1G', y = 'CDK11A',
       title = "Gene expression")
