# aggR
library(tidyverse)

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/kmhc/v1/")

flist = grep(list.files("./"),pattern = "aggR", value=TRUE)
flist
df <- NULL
#c("Bin.Aggregated.by.MAF","Average.MAF","Variants","Imputation.R2","Gold.MAF","Imputed.MAF")

head(df)
for (i in flist) {
  tmp <- read.table(i)
  tmp$fold <- i 
  df <- rbind(df,tmp)
}
head(df)
colnames(df) <- c("Bin.Aggregated.by.MAF","Average.MAF","Variants","Imputation.R2","Gold.MAF","Imputed.MAF","fold")
tmp$V1
colnames(df)
c("(0.002000,0.005000]","(0.005000,0.010000]","(0.010000,0.015000]","(0.015000,0.020000]","(0.020000,0.035000]","(0.035000,0.050000]","(0.050000,0.100000]",
  "(0.100000,0.200000]","(0.200000,0.300000]","(0.300000,0.400000]","(0.400000,0.500000]")

head(df)



df %>% count(Bin.Aggregated.by.MAF)
df %>% select(Bin.Aggregated.by.MAF,fold,Imputation.R2) %>% #head()
  ggplot(aes(x=Bin.Aggregated.by.MAF,y=Imputation.R2)) +
  geom_boxplot()



df %>% #select(Bin.Aggregated.by.MAF,fold,Imputation.R2) %>% #head()
  ggplot(aes(x=Bin.Aggregated.by.MAF,y=Variants)) +
  geom_boxplot()


##### compare
id_ref <-read.table("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/HLAID_HLAtype_posinfo.txt")
colnames(id_ref) <- c("SNP.ID","HLAtype")
id_ref %>% mutate(HLAtype = str_split_fixed(HLAtype,"_",2)[,2]) -> id_ref
freq_ref <- readxl::read_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/KMHC.HLAtype.freq.xlsx") %>% select(-n,-prop)
head(freq_ref)
head(id_ref)
colnames(freq_ref)[2] <- "HLAtype"

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/kmhc/v1/")
#setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/kmhc/v2/")
flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
kmhc <- NULL
for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  kmhc <- rbind(kmhc,df)
}

head(kmhc)

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/han/v1/")
#setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/han/v1/")
flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
han <- NULL
for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  han <- rbind(han,df)
}

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/pan/v1/")
#setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/pan/v1/")
flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
pan <- NULL
for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  pan <- rbind(pan,df)
}

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/multi/v1/")
#setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/multi/v1/")
flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
multi <- NULL
#flist
####
head(id_ref)
multi_ref <- read.table("check/check.g1.vcf.gz")
head(multi_ref)
multi_ref$SNP.ID <- paste0("6:",multi_ref$V2,":A:T")
multi_ref %>% mutate(V1 = str_split_fixed(V1,"_",2)[,2]) %>% #head()
  inner_join(id_ref) #%>% #dim()
####
head(df)
for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  multi <- rbind(multi,df)
}

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/1kgp/v1")
#setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/1kgp/v1/")
flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
kgp <- NULL
for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  kgp <- rbind(kgp,df)
}

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/")


head(kmhc)
head(multi)
head(pan)
head(han)
head(kgp)

kmhc %>% rbind(multi) %>% rbind(pan) %>% rbind(han) %>% rbind(kgp) %>% 
  mutate(Ref = str_split_fixed(fold,"\\.",4)[,1]) %>%
  mutate(fold = str_split_fixed(fold,"\\.",4)[,2]) %>% #head()
  filter(Imputation.R2 != 0,Ref != "1KGP") %>% #head()
  mutate(Ref = factor(Ref,c("KMHC","Multi","Han","PanKor"),labels = c("KMHC","Multi-ehthnic","Han Chinese","Pan-Kor"))) -> df

head(df)
table(df$Ref)
df %>%
  group_by(Ref,fold,Gene,freq) %>%
  summarise(Aggr.R2 = mean(Imputation.R2)) -> df1

df %>% group_by(Ref,fold,freq) %>% #head()
  summarise(Aggr.R2 = mean(Imputation.R2)) %>% #dim()
  mutate(Gene = "Overall") -> df2
head(df1)
head(df2)

df1 %>% ungroup() %>% count(Ref)

df1 %>% rbind(df2) %>%
  ggplot(aes(x=freq,y=Aggr.R2,fill=Ref)) +
  geom_boxplot() +
  facet_wrap(~Gene,nrow = 3)
  #scale_fill_manual(values=c("KMHC"="#F8766D","#C77CFF","#00BFC4","#7CAE00"),guide=guide_legend(c("KMHC","Multi-ethnic","Han Chinese","Pan-Kor")))
  #scale_fill_manual(values=c("#C77CFF","#00BFC4","#7CAE00","#F8766D"),guide=guide_legend(c("KMHC","Multi-ethnic","Han Chinese","Pan-Kor")))
df1 %>% rbind(df2) -> v1

df %>% #head()
  group_by(Ref,Gene,HLAtype,freq) %>%
  summarise(Aggr.R2 = mean(Imputation.R2)) -> df1

df %>% 
  group_by(Ref,HLAtype,freq) %>%
  summarise(Aggr.R2 = mean(Imputation.R2)) %>% #dim()
  mutate(Gene = "Overall") -> df2

df1 %>% filter(Gene == "A" & Ref == "Han Chinese")
df1 %>% filter(Gene == "A" & Ref == "Multi-ehthnic") %>% ungroup() %>% count(freq)


df1 %>% rbind(df2) %>%
  ggplot(aes(x=freq,y=Aggr.R2,fill=Ref)) +
  geom_boxplot() +
  facet_wrap(~Gene,nrow = 3)


freq_ref <- readxl::read_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/KMHC.HLAtype.freq.xlsx")
colnames(freq_ref)[2] <- "HLAtype"
freq_ref 
head(id_ref)
freq_ref %>% left_join(id_ref) %>% select(HLAtype,SNP.ID,prop) %>% write.table("/Volumes/DATA/HLAreferencePanel/aggR/HLAtype.af.txt",col.names = T,row.names = F,quote = F,sep = "\t")
##### all varaint
id_ref <-read.table("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/HLAID_HLAtype_posinfo.txt")
colnames(id_ref) <- c("SNP.ID","HLAtype")
id_ref %>% mutate(HLAtype = str_split_fixed(HLAtype,"_",2)[,2]) -> id_ref
freq_ref <- readxl::read_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/KMHC.HLAtype.freq.xlsx") %>% select(-n,-prop)
head(freq_ref)
head(id_ref)
colnames(freq_ref)[2] <- "HLAtype"

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/kmhc/v2/")
flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
kmhc <- NULL
for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  kmhc <- rbind(kmhc,df)
}

head(kmhc)

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/han/v2")

flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
han <- NULL
for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  han <- rbind(han,df)
}

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/pan/v2")
#setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/pan/v1/")
flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
pan <- NULL
for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  pan <- rbind(pan,df)
}

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/multi/v2")
#setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/multi/v1/")
flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
multi <- NULL

####

for (i in flist) {
    df <-read.table(i)
    colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
    df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
    multi <- rbind(multi,df)
  }

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/1kgp/v2")
#setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/1kgp/v1/")
flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
kgp <- NULL
for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  kgp <- rbind(kgp,df)
}

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/")


head(kmhc)
head(multi)
head(pan)
head(han)
head(kgp)

kmhc %>% rbind(multi) %>% rbind(pan) %>% rbind(han) %>% rbind(kgp) %>% 
  mutate(Ref = str_split_fixed(fold,"\\.",4)[,1]) %>%
  mutate(fold = str_split_fixed(fold,"\\.",4)[,2]) %>% #head()
  filter(Imputation.R2 != 0,Ref != "1KGP") %>% #head()
  mutate(Ref = factor(Ref,c("KMHC","Multi","Han","PanKor"),labels = c("KMHC","Multi-ehthnic","Han Chinese","Pan-Kor"))) -> df

head(df)
table(df$Ref)
df %>% na.omit() %>% filter(Ref == "Multi-ehthnic") -> multi_check
df %>% na.omit() %>%
  group_by(Ref,fold,Gene,freq) %>% #head()
  summarise(Aggr.R2 = mean(Imputation.R2)) -> df1


df %>% na.omit() %>%
  group_by(Ref,fold,freq) %>%
  summarise(Aggr.R2 = mean(Imputation.R2)) %>%
  mutate(Gene = "Overall") -> df2


df1 %>% rbind(df2) %>% 
  ggplot(aes(x=freq,y=Aggr.R2,fill=Ref)) +
  geom_boxplot() +
  facet_wrap(~Gene,nrow = 3)


df %>% na.omit() %>% 
  group_by(Ref,Gene,HLAtype,freq) %>% #dim()
  summarise(Aggr.R2 = mean(Imputation.R2)) -> df1


df %>% na.omit() %>% #head()
  group_by(Ref,HLAtype,freq) %>% #dim()
  summarise(Aggr.R2 = mean(Imputation.R2)) %>% #dim()
  mutate(Gene = "Overall") -> df2

df1 %>% ungroup() %>% count(Ref)

df1 %>% rbind(df2) %>% 
  ggplot(aes(x=freq,y=Aggr.R2,fill=Ref)) +
  geom_boxplot() +
  facet_wrap(~Gene,nrow = 3)


df %>% #na.omit() %>%
  group_by(Ref,fold,Gene,freq) %>% #head()
  summarise(Aggr.R2 = mean(Imputation.R2)) -> df1


df %>% #na.omit() %>%
  group_by(Ref,fold,freq) %>%
  summarise(Aggr.R2 = mean(Imputation.R2)) %>%
  mutate(Gene = "Overall") -> df2


df1 %>% rbind(df2) %>% 
  ggplot(aes(x=freq,y=Aggr.R2,fill=Ref)) +
  geom_boxplot()# +
  #facet_wrap(~Gene,nrow = 3)





df1 %>% rbind(df2) %>%  #head()
  filter(Gene=='DPA1')

df %>% na.omit() %>% #head()
  group_by(Ref,Gene,HLAtype,freq) %>% #head()
  summarise(Aggr.R2 = mean(Imputation.R2)) %>% #head()
  ggplot(aes(x=freq,y=Aggr.R2,fill=Ref)) +
  geom_boxplot() +
  facet_wrap(~Gene,nrow = 3)


v1 %>%
  ggplot(aes(x=freq,y=Aggr.R2,fill=Ref)) +
  geom_boxplot() +
  facet_wrap(~Gene,nrow = 3)


df1 %>% rbind(df2) %>% 
  ggplot(aes(x=freq,y=Aggr.R2,fill=Ref)) +
  geom_boxplot() +
  facet_wrap(~Gene,nrow = 3)

head(han)
han %>% na.omit()


##### 1 field 
id_ref <-read.table("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/HLAtype_pos.txt")
colnames(id_ref) <- c("SNP.ID","HLAtype")
id_ref %>% mutate(HLAtype = str_split_fixed(HLAtype,"_",2)[,2]) %>% 
  filter(!grepl(":",HLAtype))-> id_ref

#freq_ref <- readxl::read_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/KMHC.HLAtype.freq.xlsx") %>% select(-n,-prop)
freq_ref <- readxl::read_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/KMHC.HLAtype.freq_1field.xlsx") %>% select(-n,-prop)
head(freq_ref)
head(id_ref)
colnames(freq_ref)[2] <- "HLAtype"

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/kmhc/v2/")
flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
kmhc <- NULL
for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  kmhc <- rbind(kmhc,df)
}

head(kmhc)

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/han/v2")

flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
han <- NULL
for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  han <- rbind(han,df)
}

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/pan/v2")
#setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/pan/v1/")
flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
pan <- NULL
for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  pan <- rbind(pan,df)
}

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/multi/v2")
#setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/multi/v1/")
flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
multi <- NULL

####

for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  multi <- rbind(multi,df)
}

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/1kgp/v2")
#setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/1kgp/v1/")
flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
kgp <- NULL
for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  kgp <- rbind(kgp,df)
}

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/")


head(kmhc)
head(multi)
head(pan)
head(han)
head(kgp)

kmhc %>% rbind(multi) %>% rbind(pan) %>% rbind(han) %>% rbind(kgp) %>% 
  mutate(Ref = str_split_fixed(fold,"\\.",4)[,1]) %>%
  mutate(fold = str_split_fixed(fold,"\\.",4)[,2]) %>% #head()
  filter(Imputation.R2 != 0,Ref != "1KGP") %>% #head()
  mutate(Ref = factor(Ref,c("KMHC","Multi","Han","PanKor"),labels = c("KMHC","Multi-ehthnic","Han Chinese","Pan-Kor"))) -> df

df %>% na.omit() %>%
  group_by(Ref,fold,Gene,freq) %>% #head()
  summarise(Aggr.R2 = mean(Imputation.R2)) -> df1
head(df1)

df %>% na.omit() %>%
  group_by(Ref,fold,freq) %>%
  summarise(Aggr.R2 = mean(Imputation.R2)) %>%
  mutate(Gene = "Overall") -> df2
head(df2)

df1 %>% rbind(df2) %>% 
  ggplot(aes(x=freq,y=Aggr.R2,fill=Ref)) +
  geom_boxplot() +
  facet_wrap(~Gene,nrow = 3)


df %>% na.omit() %>% 
  group_by(Ref,Gene,HLAtype,freq) %>% #dim()
  summarise(Aggr.R2 = mean(Imputation.R2)) -> df1


df %>% na.omit() %>% #head()
  group_by(Ref,HLAtype,freq) %>% #dim()
  summarise(Aggr.R2 = mean(Imputation.R2)) %>% #dim()
  mutate(Gene = "Overall") -> df2

df1 %>% ungroup() %>% count(Ref)

df1 %>% rbind(df2) %>% 
  ggplot(aes(x=freq,y=Aggr.R2,fill=Ref)) +
  geom_boxplot() +
  facet_wrap(~Gene,nrow = 3)

####### v4
id_ref <-read.table("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/HLAtype_pos.txt")
colnames(id_ref) <- c("SNP.ID","HLAtype")
id_ref %>% mutate(HLAtype = str_split_fixed(HLAtype,"_",2)[,2]) %>% 
  filter(grepl(":",HLAtype))-> id_ref

freq_ref <- readxl::read_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/KMHC.HLAtype.freq.xlsx") %>% select(-n,-prop)
#freq_ref <- readxl::read_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/KMHC.HLAtype.freq_1field.xlsx") %>% select(-n,-prop)
head(freq_ref)
head(id_ref)
colnames(freq_ref)[2] <- "HLAtype"

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/kmhc/v4/")
flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
kmhc <- NULL
for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  kmhc <- rbind(kmhc,df)
}

head(kmhc)

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/han/v4")

flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
han <- NULL
for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  han <- rbind(han,df)
}

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/pan/v4")
#setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/pan/v1/")
flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
pan <- NULL
for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  pan <- rbind(pan,df)
}

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/multi/v4")
#setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/multi/v1/")
flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
multi <- NULL

####

for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  multi <- rbind(multi,df)
}

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/1kgp/v4")
#setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/1kgp/v1/")
flist <- grep(list.files("./"),pattern = "aggR.RSquare", value=TRUE)
kgp <- NULL
for (i in flist) {
  df <-read.table(i)
  colnames(df) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")
  df %>% left_join(id_ref) %>% left_join(freq_ref) %>% select(Gene,HLAtype,Imputation.R2,freq) %>% mutate(fold = i) -> df
  kgp <- rbind(kgp,df)
}

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/")


head(kmhc)
head(multi)
head(pan)
head(han)
head(kgp)

kmhc %>% rbind(multi) %>% rbind(pan) %>% rbind(han) %>% rbind(kgp) %>% 
  mutate(Ref = str_split_fixed(fold,"\\.",4)[,1]) %>%
  mutate(fold = str_split_fixed(fold,"\\.",4)[,2]) %>% #head()
  filter(Imputation.R2 != 0,Ref != "1KGP") %>% #head()
  mutate(Ref = factor(Ref,c("KMHC","Multi","Han","PanKor"),labels = c("KMHC","Multi-ehthnic","Han Chinese","Pan-Kor"))) -> df

head(kmhc)
kmhc %>% na.omit()

df %>% na.omit()
df %>% na.omit() %>%
  group_by(Ref,fold,Gene,freq) %>% #head()
  summarise(Aggr.R2 = mean(Imputation.R2)) -> df1
head(df1)

df %>% na.omit() %>%
  group_by(Ref,fold,freq) %>%
  summarise(Aggr.R2 = mean(Imputation.R2)) %>%
  mutate(Gene = "Overall") -> df2
head(df2)

df1 %>% rbind(df2) %>% 
  ggplot(aes(x=freq,y=Aggr.R2,fill=Ref)) +
  geom_boxplot() +
  facet_wrap(~Gene,nrow = 3)


df %>% na.omit() %>% 
  group_by(Ref,Gene,HLAtype,freq) %>% #dim()
  summarise(Aggr.R2 = mean(Imputation.R2)) -> df1


df %>% na.omit() %>% #head()
  group_by(Ref,HLAtype,freq) %>% #dim()
  summarise(Aggr.R2 = mean(Imputation.R2)) %>% #dim()
  mutate(Gene = "Overall") -> df2

df1 %>% ungroup() %>% count(Ref)

df1 %>% rbind(df2) %>% 
  ggplot(aes(x=freq,y=Aggr.R2,fill=Ref)) +
  geom_boxplot() +
  facet_wrap(~Gene,nrow = 3)



multi %>% na.omit() %>% filter(Imputation.R2 != 0) %>% head()
head(df)
table(df$Ref)
df %>% na.omit() %>% filter(Imputation.R2 != 0,Ref != "1KGP") %>% filter(Ref == "Multi-ehthnic") %>% filter(Gene == "A")

#### multi check

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/multi/v1/")
flist
flist = grep(list.files("./"),pattern = "R2", value=TRUE)
df <- NULL
for (i in flist) {
  tmp <- read.table(i)
  tmp$g <- i
  df <- rbind(df,tmp)
}
head(tmp)
head(df)
ref <- readxl::read_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/KMHC.HLAtype.freq.xlsx")
head(ref)
ref %>% mutate(HLAtype=value) %>% select(Gene,HLAtype,freq) -> ref

df %>% mutate(HLAtype = str_split_fixed(V1,"HLA_",2)[,2]) %>%
  mutate(Gene = str_split_fixed(HLAtype,"\\*",2)[,1]) %>% 
  mutate(fold = str_split_fixed(g,"\\.",2)[,1],R2=V2) %>% #head()
  select(fold,Gene,HLAtype,R2) %>% #head()
  #filter(HLAtype %in% ref$value) %>% 
  inner_join(ref) %>% #head()
  group_by(Gene,HLAtype,freq) %>%
  summarise(meanR2 = mean(R2)) %>% #head() 
  ggplot(aes(x=freq,y=meanR2,fill=Gene)) +
  geom_boxplot()
  
