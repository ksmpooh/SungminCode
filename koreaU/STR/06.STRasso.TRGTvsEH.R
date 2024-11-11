library(tidyverse)
library(data.table)
library(ggpubr)
library(ggbreak)

#### STR asso
################ ASSO data pro

lc <- read_table("~/Desktop/KU/@research/STR/02.compare/EH.common.pass.QC.info_LC_type.process.txt")
str_match_score_freqeuncy <- read_table("~/Desktop/KU/@research/STR/02.compare/STR_type/str_match_score_frequency.txt")
df <- read_table("~/Desktop/KU/@research/STR/02.compare/STR.TRGT_0.8upper.EH_pass.common.merge.onlyReaptnumber.with_dfiff.txt")

head(str_match_score_freqeuncy)
head(df)
head(lc)

df %>% left_join(lc) -> df
df %>% mutate(TRGT_Allele12_repeatcount_diff = abs(TRGT_STR1 - TRGT_STR2)) -> df

df %>% select(ID,STR,MOTIFS,EH_STR1,TRGT_STR1,STR1,EH_type1,TRGT_Allele12_repeatcount_diff) -> df1
df %>% select(ID,STR,MOTIFS,EH_STR2,TRGT_STR2,STR2,EH_type2,TRGT_Allele12_repeatcount_diff) -> df2
head(df1)
head(df2)
df1$allele <- "A1"
df2$allele <- "A2"

head(df1)
head(df2)
colnames(df1) <- c("ID","STR","MOTIFS","EH","TRGT","STR_repeat_diff","EH_type","TRGT_Allele12_repeatcount_diff","Allele")
colnames(df2) <- c("ID","STR","MOTIFS","EH","TRGT","STR_repeat_diff","EH_type","TRGT_Allele12_repeatcount_diff","Allele")

df1 %>% rbind(df2) %>% mutate(RU.length = str_length(MOTIFS)) %>% #head()
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> df


#write.table(df,"~/Desktop/KU/@research/STR/02.compare/STR_forASSO.txt",col.names = T,row.names = F,quote = F,sep = "\t")

########

df <- read_table("~/Desktop/KU/@research/STR/02.compare/STR_forASSO.txt")
str_match_score_freqeuncy <- read_table("~/Desktop/KU/@research/STR/02.compare/STR_type/str_match_score_frequency.txt")
head(df)
head(str_match_score_freqeuncy)

ref_db <- readxl::read_xlsx("~/Desktop/KCDC/paper/STR/2023_guo_sup_2.xlsx",sheet = 1)
head(ref_db)
ref_db %>% select(chr,`STR start position (GRCh38)`,`STR end position (GRCh38)`,motif,`Reference tract length`,`Genomic annotation`,`Distance to nearest TSS`) -> ref_db
colnames(ref_db) <- c("chrom","start","end","MOTIFS","Referencetractlength","anotation","DistancetoTSS")
ref_db %>% mutate(STR_region = paste0(chrom,"_",start,"_",end)) -> ref_db

head(str_match_score_freqeuncy)
head(ref_db)
ru <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.pass.IDregion.txt")
ru %>% left_join(ref_db) %>% mutate(anotation = ifelse(is.na(anotation),"chrX",anotation)) %>% unique()-> ru

head(ru)

head(df)

ru %>% left_join(str_match_score_freqeuncy) %>% 
  mutate(RU.length=str_length(MOTIFS)) %>% 
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) %>%  
  rename(annotation = anotation)-> df
head(df)  

model <- lm(concordance~anotation,data=df)
summary(model)

model <- lm(concordance~RU.length,data=df)
summary(model)

model <- lm(concordance~GC,data=df)
summary(model)

model <- lm(concordance~annotation+RU.length+GC,data=df)
summary(model)

head(lc)

df <- read_table("~/Desktop/KU/@research/STR/02.compare/STR_forASSO.txt")

head(df)
head(ru)
df %>% mutate(new_ID = paste0(ID,"_",Allele)) %>% #head()
  select(new_ID,STR,MOTIFS,EH,TRGT,RU.length,GC) -> df

head(df)
df %>% mutate(match = ifelse(EH == TRGT,1,0)) -> df
toy <- df %>% filter(STR == 'AFF2')
head(toy)
head(str_match_score_freqeuncy$)
result <- NULL
for(i in str_match_score_freqeuncy ){
  if (colnames(c1_c2_eplet_forAsso)[i] != "theme") {
    c1_c2_eplet_forAsso[,i] <- scale(as.numeric(c1_c2_eplet_forAsso[,i]))
    temp <- glm(paste("rej_tot ~ ",colnames(c1_c2_eplet_forAsso)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group",sep=""), data=c1_c2_eplet_forAsso %>% filter(theme == "EMS"), family="binomial")
    result <- rbind(result, c("EMS",colnames(c1_c2_eplet_forAsso)[i],as.vector(summary(temp)$coefficients[2,])))
    temp <- glm(paste("rej_tot ~ ",colnames(c1_c2_eplet_forAsso)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group",sep=""), data=c1_c2_eplet_forAsso %>% filter(theme != "EMS"), family="binomial")
    result <- rbind(result, c("SMMS",colnames(c1_c2_eplet_forAsso)[i],as.vector(summary(temp)$coefficients[2,])))
  }  
}

####
result <- NULL
for(i in 3:27){
  if (colnames(c1_c2_eplet_forAsso)[i] != "theme") {
    c1_c2_eplet_forAsso[,i] <- scale(as.numeric(c1_c2_eplet_forAsso[,i]))
    temp <- glm(paste("rej_tot ~ ",colnames(c1_c2_eplet_forAsso)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group",sep=""), data=c1_c2_eplet_forAsso %>% filter(theme == "EMS"), family="binomial")
    result <- rbind(result, c("EMS",colnames(c1_c2_eplet_forAsso)[i],as.vector(summary(temp)$coefficients[2,])))
    temp <- glm(paste("rej_tot ~ ",colnames(c1_c2_eplet_forAsso)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group",sep=""), data=c1_c2_eplet_forAsso %>% filter(theme != "EMS"), family="binomial")
    result <- rbind(result, c("SMMS",colnames(c1_c2_eplet_forAsso)[i],as.vector(summary(temp)$coefficients[2,])))
  }  
}
head(result)
as.vector(summary(temp)$coefficients[2,])
summary(temp)

colnames(result) <- c("theme","ID","Estimate","SE","Z","P")

result <- as.data.frame(result)
for (i in 3:6) {
  result[,i] <- as.numeric(result[,i])
  
}


###