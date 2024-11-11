## main.figure3 .data processing # patho # before chrX
library(tidyverse)


head(patho_complex)
patho_complex <- read_table("~/Desktop/KU/@research/STR/figure/complexSTR_prep.EH_TRGT_merged.txt")
eh_trgt_merge_simple_pass_intersect_forConcordance <- read_table("~/Desktop/KU/@research/STR/figure/figure2/eh_trgt_merge_simple_pass_intersect_forConcordance.txt")
head(eh_trgt_merge_simple_pass_intersect_forConcordance)
#eh_trgt_merge_simple_pass_intersect_forConcordance %>% filter(ID == 'NIH23J3904558')
eh_trgt_merge_simple_pass_intersect_forConcordance %>% select(ID:TRGT_STR2,RU:EH_STR2) %>% filter(str_detect(STR_ID,"chr")) %>%
  mutate(STR1 = TRGT_STR1-EH_STR1,STR2 = TRGT_STR2-EH_STR2) %>% 
  mutate(STR1 = ifelse(STR1 == 0,1,0),STR2 = ifelse(STR2 == 0,1,0)) %>% #head()
  mutate(concordance_point = STR1 + STR2) %>% group_by(STR_ID) %>% 
  summarise(concordance_point_sum = sum(concordance_point)) %>% 
  mutate(concordance_rate = concordance_point_sum/130) -> concordance_bySTR_simpleSTR_normalSTR
head(concordance_bySTR_simpleSTR_normalSTR)
dim(concordance_bySTR_simpleSTR_normalSTR)[1]
eh_trgt_merge_simple_pass_intersect_forConcordance %>% select(ID:TRGT_STR2,RU:EH_STR2) %>% filter(str_detect(STR_ID,"chr")) %>% 
  mutate(STR1 = TRGT_STR1-EH_STR1,STR2 = TRGT_STR2-EH_STR2) %>% 
  mutate(STR1 = ifelse(STR1 == 0,1,0),STR2 = ifelse(STR2 == 0,1,0)) %>% #head()
  mutate(concordance_point = (STR1 + STR2)) %>% group_by(ID) %>% #head()
  summarise(concordance_rate = sum(concordance_point)/(dim(concordance_bySTR_simpleSTR_normalSTR)[1]*2)) -> concordance_byID_simpleSTR_normalSTR
  #mutate(concordance_rate = concordance_point_sum/(311542*2)) -> concordance_byID 
#311542*2
head(concordance_byID_simpleSTR_normalSTR)
head(concordance_bySTR_simpleSTR_normalSTR)
#eh_trgt_merge_simple_pass_intersect_forConcordance %>% filter(!str_detect(STR_ID,"chr")) -> eh_trgt_merge_simple_pass_intersect_forConcordance_onlysimplepath
#eh_trgt_merge_simple_pass_intersect_forConcordance_onlysimplepath %>% write.table("~/Desktop/KU/@research/STR/figure/simplepathogenic_STR_prep.EH_TRGT_merged.txt",col.names = T,row.names = F,quote = F,sep = "\t")

patho_simple <- read_table("~/Desktop/KU/@research/STR/figure/simplepathogenic_STR_prep.EH_TRGT_merged.txt")
patho_complex$STR_type <- "complex"
patho_simple$STR_type <- "simple"

head(patho_complex)
head(patho_simple)

patho_simple %>% mutate(new_ID = STR_ID) %>% select(-RU) %>% rbind(patho_complex) %>% #head()
  select(ID:TRGT_STR2,EH_STR1:new_ID) %>%
  mutate(STR1 = TRGT_STR1-EH_STR1,STR2 = TRGT_STR2-EH_STR2) %>% 
  mutate(STR1 = ifelse(STR1 == 0,1,0),STR2 = ifelse(STR2 == 0,1,0)) %>% #head()
  mutate(concordance_point = STR1 + STR2) -> patho_merge_concordance_point 
  
patho_merge_concordance_point %>% #filter(STR_type == "complex") %>%
  group_by(new_ID) %>% summarise(concordance_point_sum = sum(concordance_point)) %>% 
  mutate(concordance_rate = concordance_point_sum/130)


patho_merge_concordance_point %>% filter(STR_type == "complex") %>% count(STR_ID) %>% mutate(n=n*2) %>% rename(allele_count = n)-> complex_allele_count
patho_merge_concordance_point %>% filter(STR_type == "complex") %>% select(STR_ID,new_ID) %>% unique() %>% count(STR_ID) %>% rename(str_number = n) -> complex_strcount_bySTR_ID
head(complex_strcount_bySTR_ID)
head(complex_allele_count)

patho_merge_concordance_point %>% filter(STR_type == "complex") %>% count(STR_ID,concordance_point) %>% #head()
  left_join(complex_allele_count) %>% #head()
  group_by(STR_ID) %>% #count(new_ID)
  reframe(concordance_rate = sum(concordance_point*n)/allele_count) %>% unique() -> complexSTR_concordance_bySTR_ID
head(complexSTR_concordance_bySTR_ID)

patho_merge_concordance_point %>% filter(STR_type == "complex") %>% group_by(ID,STR_ID) %>%
  summarise(STR1_sum = sum(STR1),STR2_sum = sum(STR2)) %>% left_join(complex_strcount_bySTR_ID) %>%
  mutate(concordance_point = ifelse(STR1_sum == str_number & STR2_sum == str_number,2,ifelse(STR1_sum != str_number & STR2_sum != str_number,0,1))) %>%
  group_by(STR_ID) %>% summarise(concordance_point_sum = sum(concordance_point)) %>%
  mutate(concordance_rate = concordance_point_sum/130)
unique(complexSTR_concordance_bySTR_ID)
patho_merge_concordance_point %>% filter(STR_type == "complex") %>% count(STR_ID,concordance_point) %>% head()
  
#concordance_byID_simpleSTR_normalSTR
patho_merge_concordance_point %>% filter(STR_type == "complex") %>% group_by(ID,STR_ID) %>% head()
patho_merge_concordance_point %>% filter(STR_type == "complex") %>% group_by(ID,STR_ID) %>% summarise(mean_STR1 = mean(STR1),mean_STR2 = mean(STR2)) %>% #head()
  mutate(mean_STR1 = ifelse(mean_STR1 == 1,mean_STR1,0),mean_STR2 = ifelse(mean_STR2 == 1,mean_STR2,0)) %>% 
  mutate(concordance_point = mean_STR1 + mean_STR2) -> concordance_point_complexSTR
head(concordance_point_complexSTR)
concordance_point_complexSTR$STR_ID %>% unique() %>% length() -> number_of_complex_STR
concordance_point_complexSTR %>% group_by(ID) %>% summarise(concordance_rate = sum(concordance_point)/(number_of_complex_STR*2)) -> concordance_byID_complexSTR_pathoSTR



patho_merge_concordance_point %>% filter(STR_type != "complex") %>% select(STR_ID) %>% unique() %>% dim() -> number_of_simple_STR_patho
head(number_of_simple_STR_patho)
patho_merge_concordance_point %>% filter(STR_type != "complex") %>% group_by(ID) %>% 
  summarise(concordance_rate = sum(concordance_point)/(number_of_simple_STR_patho[1]*2)) -> concordance_byID_simpleSTR_pathoSTR


patho_merge_concordance_point %>% filter(STR_type != "complex") %>% select(ID,STR_ID,concordance_point)

# path : simple + complex by ID
patho_merge_concordance_point %>% filter(STR_type == "complex") %>% group_by(ID,STR_ID) %>%
  summarise(STR1_sum = sum(STR1),STR2_sum = sum(STR2)) %>% left_join(complex_strcount_bySTR_ID) %>%
  mutate(concordance_point = ifelse(STR1_sum == str_number & STR2_sum == str_number,2,ifelse(STR1_sum != str_number & STR2_sum != str_number,0,1))) %>% 
  select(ID,STR_ID,concordance_point) %>%
  rbind(patho_merge_concordance_point %>% filter(STR_type != "complex") %>% select(ID,STR_ID,concordance_point)) -> patho_merge_concordance_point_prep

head(patho_merge_concordance_point_prep)
patho_merge_concordance_point_prep %>% ungroup()%>%select(STR_ID) %>% unique() %>% dim() -> number_of_patho_STR

patho_merge_concordance_point_prep  %>% group_by(ID) %>%
  summarise(concordance_rate = sum(concordance_point)/(number_of_patho_STR[1]*2)) -> concordance_byID_all_pathoSTR



head(concordance_byID_complexSTR_pathoSTR)
head(concordance_byID_simpleSTR_pathoSTR)
head(concordance_byID_all_pathoSTR)
head(concordance_byID_simpleSTR_normalSTR)
  
concordance_byID_complexSTR_pathoSTR$g = "Pathogenic"
concordance_byID_simpleSTR_pathoSTR$g = "Pathogenic"
concordance_byID_all_pathoSTR$g = "Pathogenic"
concordance_byID_simpleSTR_normalSTR$g = "non-Pathogenic"

concordance_byID_complexSTR_pathoSTR$str_type = "(complexSTR)"
concordance_byID_simpleSTR_pathoSTR$str_type = "(simpleSTR)"
concordance_byID_all_pathoSTR$str_type = "(Overall)"
concordance_byID_simpleSTR_normalSTR$str_type = "(simpleSTR)"


concordance_byID_complexSTR_pathoSTR %>% 
  rbind(concordance_byID_simpleSTR_pathoSTR) %>%
  rbind(concordance_byID_all_pathoSTR) %>%
  rbind(concordance_byID_simpleSTR_normalSTR) %>% #write.table("~/Desktop/KU/@research/STR/figure/figure3/concordance.byID.nonpatho.patho.simpleSTR.complexSTR.txt",col.names = T,row.names = F,quote = F,sep = "\t")
  mutate(new_g = paste0(g,"\n",str_type)) %>%
  ggplot(aes(x=factor(new_g,levels=c("non-Pathogenic\n(simpleSTR)","Pathogenic\n(Overall)","Pathogenic\n(simpleSTR)","Pathogenic\n(complexSTR)")),y=concordance_rate,fill=factor(new_g,levels=c("non-Pathogenic\n(simpleSTR)","Pathogenic\n(Overall)","Pathogenic\n(simpleSTR)","Pathogenic\n(complexSTR)")))) + 
  geom_boxplot() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black",alpha=0.8) + 
  theme_bw() +
  theme_step1() + 
  theme(legend.position = 'none',
        axis.title.x=element_blank())


#####





  