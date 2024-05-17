library(writexl)
library(tidyverse)
library(data.table)
library(ggpubr)
library(ggthemes)
library(reshape2)
library(ggplot2)

nano <- readxl::read_xlsx("~/Desktop/KCDC/pangenome/00.datacheck/2023_pro_long.xlsx",sheet = 1)
revio48 <- readxl::read_xlsx("~/Desktop/KCDC/pangenome/00.datacheck/2023_pro_long.xlsx",sheet = 3)

head(nano)
head(revio48)
nano %>% select(1,2,3,4,5) -> nano
revio48 %>% select(1,2,3,7,4) -> revio48
colnames(revio48) <- colnames(nano)

nano$platform <- "Nanopore"
revio48$platform <- "Revio"

nano %>% rbind(revio48) -> df
#df$`Read N50 (bp)` <- as.numeric(df$`Read N50 (bp)`)

#df$`Average Read Length (bp)` <- as.numeric(df$`Average Read Length (bp)`)
head(df)
ggplot(df,aes(y=`Total Bases (bp)`,x=platform,fill=platform)) +
  geom_boxplot() + 
  theme(legend.position = 'none') -> a
a
ggplot(df,aes(y=`Total Reads`,x=platform,fill=platform)) +
  geom_boxplot() + 
  theme(legend.position = 'none') -> b
b
ggplot(df,aes(y=`Read N50 (bp)`,x=platform,fill=platform)) +
  geom_boxplot() + 
  theme(legend.position = 'none') -> c
c
ggplot(df,aes(y=`Average Read Length (bp)`,x=platform,fill=platform)) +
  geom_boxplot() + 
  theme(legend.position = 'none')-> d
d
ggarrange(b,d,c,a, nrow = 1)
c
head(df) 

df %>% pivot_longer(2:5) %>% group_by(platform,name) %>%
  summarise(mean=mean(value), sd=sd(value)) -> a

###
ref <- read.table("/Users/ksmpooh/Desktop/KCDC/pangenome/00.datacheck/KBA.Long_Revio_Nanopore_short.IDmatchinagtable.txt",header = T) #%>% na.omit()
head(ref)
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/00.datacheck/01.vcf.stat")
df <-read.table("filter.info.txt",header = T)
wgs <- read.table("wgs/auto/filter.info.txt",header = T)
head(df)
head(wgs)

wgs %>% mutate(ID = str_remove(file,".autosomal.vcf.gz")) %>% mutate(platform = "WGS") %>% left_join(ref %>% mutate(ID = Illumina)) %>% na.omit() %>%
  mutate(ID = KBAv1) %>% select(-colnames(ref)) -> wgs
  
df %>% mutate(ID = str_split_fixed(file,"\\.",5)[,2]) %>% filter(grepl("autosomal",file)) %>% filter(!grepl("merge",file)) %>% #dim()
  mutate(platform = str_split_fixed(file,"\\.",5)[,1]) %>%
  mutate(platform = str_split_fixed(platform,"_",2)[,1])-> df
head(df)
colnames(df)

df %>% rbind(wgs) -> df
#df %>% left_join(wgs) -> df

df %>% filter(sample == 1) %>% select(ID,platform,SNPs,indels) %>% 
  pivot_longer(3:4) %>% 
  ggplot(aes(x=ID,y=value,fill=name)) +
  geom_bar(stat='identity') +
  #xlab(element_blank())
  theme(axis.text.x = element_blank()) +
  ylab("Count") + xlab("Sample") + 
  facet_grid(~platform) 

df %>% filter(sample == 1) %>% select(ID,platform,records,SNPs,indels) %>%
  #pivot_longer(3:5) %>% #head()
  ggplot(aes(x=value,fill=platform))+
  geom_histogram() +
  facet_grid(~name)
  
df %>% filter(sample == 1) %>% #select(ID,platform,records,SNPs,indels) %>%
  ggplot(aes(x=records,fill=platform))+
  geom_histogram() -> a

df %>% filter(sample == 1) %>% #select(ID,platform,records,SNPs,indels) %>%
  ggplot(aes(x=SNPs,fill=platform))+
  geom_histogram() -> b


df %>% filter(sample == 1) %>% #select(ID,platform,records,SNPs,indels) %>%
  ggplot(aes(x=indels,fill=platform))+
  geom_histogram() -> c



head(df)
df %>% 
  select(platform,ID,records,SNPs,indels,others,multiallelic_sites,multiallelic_SNP_sites,ts,tv,ts_tv,ts_tv_1alt) %>% #head()
  pivot_longer(3:12) %>%
  group_by(platform,name) %>% #head()
  summarise(average = mean(value)) %>% 
  pivot_wider(names_from = name,values_from = average) -> a

a

#### AF ºñ±³
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/00.datacheck/01.vcf.stat")
nano <- fread("Nanopore_repanel.merge.rmINFO_autosomal.norm_setID.onlySNP.vcf.gz.AF_stat.txt.gz")
revio <- fread("Revio_kchip.merge.deep_PASS_autosomal.norm_setID.onlySNP.vcf.gz.AF_stat.txt.gz")
head(nano)
nano <- nano[,c(1:6)]
revio <- revio[,c(1:6)]

colnames(nano) <- c("ID","F_MISSING","NS","AN","AF","AC")
colnames(revio) <- c("ID","F_MISSING","NS","AN","AF","AC")


head(nano)
head(revio)
nano$platform <- "Nanopore"
revio$platform <- "Revio"
dim(nano)
dim(revio)
nano %>% filter(F_MISSING <= 0.1) %>% dim()
nano %>% filter(F_MISSING <= 0.2) %>% dim()
nano %>% filter(F_MISSING <= 0.3) %>% dim()


revio %>% filter(F_MISSING <= 0.1) %>% dim()
revio %>% filter(F_MISSING <= 0.2) %>% dim()
revio %>% filter(F_MISSING <= 0.3) %>% dim()



nano <- nano %>% filter(ID %in% revio$ID)
revio <- revio %>% filter(ID %in% nano$ID)

#nano <- nano %>% select(ID,F_MISSING,AC,platform)
#revio <- revio %>% select(ID,F_MISSING,AC,platform)
nano <- nano %>% select(ID,F_MISSING,AC)
revio <- revio %>% select(ID,F_MISSING,AC)

colnames(nano) <- c("ID","Nanopore_Missing","Nanopore_AC")
colnames(revio) <- c("ID","Revio_Missing","Revio_AC")

head(revio)
revio %>% left_join(nano) -> m

head(m)
dim(m)
m %>% filter(Revio_Missing <= 0.5,Nanopore_Missing <=0.5) -> m50
m %>% filter(Revio_Missing <= 0.9,Nanopore_Missing <=0.9) -> m90
m %>% filter(Revio_Missing <= 0.7,Nanopore_Missing <=0.7) -> m30

dim(m50)
m %>% filter(Revio_Missing <= 0.3,Nanopore_Missing <=0.3) -> m30
m %>% filter(Revio_Missing <= 0.2,Nanopore_Missing <=0.2) -> m20
m %>% filter(Revio_Missing <= 0.1,Nanopore_Missing <=0.1) -> m10


head(m10)
m10 %>% 
  ggscatter(.,x='Revio_AC',y="Nanopore_AC",
          add = "reg.line",
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          xlab = "Revio_AC",
          ylab = "Nanopore_AC",
          title = "Call rate > 0.9",
          #xlim = c(0,40),cor.coef.size = 10,
          #ylim = c(0,40),
          #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
          cor.coeff.args = list(method = "pearson", label.sep = "\t")) + 
  theme(legend.position = "right") -> p10
m20 %>% 
  ggscatter(.,x='Revio_AC',y="Nanopore_AC",
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            xlab = "Revio_AC",
            ylab = "Nanopore_AC",
            title = "Call rate > 0.8",
            #xlim = c(0,40),cor.coef.size = 10,
            #ylim = c(0,40),
            #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
            cor.coeff.args = list(method = "pearson", label.sep = "\t")) + 
  theme(legend.position = "right") -> p20

m30 %>% 
  ggscatter(.,x='Revio_AC',y="Nanopore_AC",
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            xlab = "Revio_AC",
            ylab = "Nanopore_AC",
            title = "Call rate > 0.7",
            #xlim = c(0,40),cor.coef.size = 10,
            #ylim = c(0,40),
            #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
            cor.coeff.args = list(method = "pearson", label.sep = "\t")) + 
  theme(legend.position = "right") -> p30

p30


p20
m90 %>% 
  ggscatter(.,x='Revio_AC',y="Nanopore_AC",
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            xlab = "Revio_AC",
            ylab = "Nanopore_AC",
            title = "Call rate > 0.1",
            #xlim = c(0,40),cor.coef.size = 10,
            #ylim = c(0,40),
            #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
            cor.coeff.args = list(method = "pearson", label.sep = "\t")) + 
  theme(legend.position = "right") -> p90
p10
### concordance

setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/00.datacheck/02.concordance/03")
#setwd("~/")

file_list <- list.files(pattern = "by_sample")
file_list

df <- NULL
for (i in file_list) {
  data <- read.table(i,header = T)
  data <- cbind(data, TP=(data$REF.REF + data$ALT_1.ALT_1 + data$ALT_2.ALT_2))
  data <- cbind(data, FN=(data$ALT_1.REF + data$ALT_2.REF + data$ALT_2.ALT_1))
  data <- cbind(data, FP=(data$REF.ALT_1 + data$REF.ALT_2 + data$ALT_1.ALT_2))
  
  data$Sensitivity <- data$TP/(data$TP+data$FN)
  data$Precision <- data$TP/(data$TP+data$FP)
  data$Accuracy <- data$TP/(data$TP+data$FP+data$FN)
  
  a <- data %>% #inner_join(maf1) %>% 
    select(sample,Sensitivity,Precision,Accuracy) %>% 
    pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% 
    mutate(title = str_remove(str_remove(i,".by_sample.txt"),"concordance_"))
  df <- rbind(df,a)
}
#head(data)
head(df)
table(df$title)

df %>% mutate(title = recode(title,"0000_0001"="KBA vs Revio","0000_0002" ="KBA vs Nanopore",
                             "0001_0002" ="Revio vs Nanopore",
                             "0002_0003" ="Revio vs Nanopore (ALL)"
                             )) %>%
  ggplot(aes(x=type,y=Val,color=type)) + 
  geom_boxplot(outlier.size = 0.1) + 
  facet_wrap(~title,ncol=4) + 
  labs(x= NULL,y=NULL,colour=NULL,title="Concordance Test Result by sample") +
  theme(plot.title = element_text(hjust = 0.5))

df %>% mutate(title = recode(title,"0000_0001"="KBA vs Revio","0000_0002" ="KBA vs Nanopore",
                             "0001_0002" ="Revio vs Nanopore",
                             "0002_0003" ="Revio vs Nanopore (ALL)")) %>% #head()
  filter(type == "Accuracy") %>% group_by(title) %>%
  summarise(mean = mean(Val))
  

#### wgs vs revio vs nanopore  
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/00.datacheck/02.concordance/newcall_withWGS/inter3all/")
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/00.datacheck/02.concordance/with_WGS/filter2_inter3all/")
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/00.datacheck/02.concordance/with_WGS/filter3_inter3all/")

file_list <- list.files(pattern = "by_sample")
file_list

df <- NULL
for (i in file_list) {
  data <- read.table(i,header = T)
  data <- cbind(data, TP=(data$REF.REF + data$ALT_1.ALT_1 + data$ALT_2.ALT_2))
  data <- cbind(data, FN=(data$ALT_1.REF + data$ALT_2.REF + data$ALT_2.ALT_1))
  data <- cbind(data, FP=(data$REF.ALT_1 + data$REF.ALT_2 + data$ALT_1.ALT_2))
  
  data$Sensitivity <- data$TP/(data$TP+data$FN)
  data$Precision <- data$TP/(data$TP+data$FP)
  data$Accuracy <- data$TP/(data$TP+data$FP+data$FN)
  
  a <- data %>% #inner_join(maf1) %>% 
    select(sample,Sensitivity,Precision,Accuracy) %>% 
    pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% 
    mutate(title = str_remove(str_remove(i,".by_sample.txt"),"concordance_"))
  df <- rbind(df,a)
}

setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/00.datacheck/02.concordance/with_WGS/inter_wgs_nanopore/")
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/00.datacheck/02.concordance/newcall_withWGS/inter_wgs_nanopore/")
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/00.datacheck/02.concordance/with_WGS/filter2_inter_wgs_nanopore/")
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/00.datacheck/02.concordance/with_WGS/filter3_inter_wgs_nanopore/")


file_list <- list.files(pattern = "by_sample")
file_list
for (i in file_list) {
  data <- read.table(i,header = T)
  data <- cbind(data, TP=(data$REF.REF + data$ALT_1.ALT_1 + data$ALT_2.ALT_2))
  data <- cbind(data, FN=(data$ALT_1.REF + data$ALT_2.REF + data$ALT_2.ALT_1))
  data <- cbind(data, FP=(data$REF.ALT_1 + data$REF.ALT_2 + data$ALT_1.ALT_2))
  
  data$Sensitivity <- data$TP/(data$TP+data$FN)
  data$Precision <- data$TP/(data$TP+data$FP)
  data$Accuracy <- data$TP/(data$TP+data$FP+data$FN)
  
  a <- data %>% #inner_join(maf1) %>% 
    select(sample,Sensitivity,Precision,Accuracy) %>% 
    pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% 
    #mutate(title = str_remove(str_remove(i,".by_sample.txt"),"concordance_"))
    mutate(title = "WGS_Nanopore")
  df <- rbind(df,a)
}
head(df)

setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/00.datacheck/02.concordance/newcall_withWGS/inter_wgs_revio/")
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/00.datacheck/02.concordance/with_WGS/filter2_inter_wgs_revio/")
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/00.datacheck/02.concordance/with_WGS/filter3_inter_wgs_revio/")
file_list <- list.files(pattern = "by_sample")
file_list
for (i in file_list) {
  data <- read.table(i,header = T)
  data <- cbind(data, TP=(data$REF.REF + data$ALT_1.ALT_1 + data$ALT_2.ALT_2))
  data <- cbind(data, FN=(data$ALT_1.REF + data$ALT_2.REF + data$ALT_2.ALT_1))
  data <- cbind(data, FP=(data$REF.ALT_1 + data$REF.ALT_2 + data$ALT_1.ALT_2))
  
  data$Sensitivity <- data$TP/(data$TP+data$FN)
  data$Precision <- data$TP/(data$TP+data$FP)
  data$Accuracy <- data$TP/(data$TP+data$FP+data$FN)
  
  a <- data %>% #inner_join(maf1) %>% 
    select(sample,Sensitivity,Precision,Accuracy) %>% 
    pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% 
    #mutate(title = str_remove(str_remove(i,".by_sample.txt"),"concordance_"))
    mutate(title = "WGS_Revio")
  df <- rbind(df,a)
}
head(df)


table(df$title)

df %>% mutate(title = recode(title,"0000_0001"="#WGS vs Revio","0000_0002" ="#WGS vs Nanopore",
                             "0001_0002" ="#Revio vs Nanopore",
                             "WGS_Revio" ="WGS vs Revio (2C)",
                             "WGS_Nanopore" ="WGS vs Nanopore (2C)")) %>%
  ggplot(aes(x=type,y=Val,color=type)) + 
  geom_boxplot(outlier.size = 0.1) + 
  facet_wrap(~title,ncol=5) + 
  labs(x= NULL,y=NULL,colour=NULL,title="Concordance Test Result by sample") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank())

df %>% mutate(title = recode(title,"0000_0001"="#WGS vs Revio","0000_0002" ="#WGS vs Nanopore",
                             "0001_0002" ="#Revio vs Nanopore",
                             "WGS_Revio" ="WGS vs Revio (2C)",
                             "WGS_Nanopore" ="WGS vs Nanopore (2C)" ))%>%
  filter(type == "Accuracy") %>% group_by(title) %>%
  summarise(mean = mean(Val)) -> a

#### ID venn
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/00.datacheck/02.concordance/with_WGS/")

nano <- read_table("Nanopore_repanel.merge.rmINFO_autosomal.norm_setID.onlySNP.03snpqc.vcf.gz.ID",col_names = F)
revio <- read_table("Revio_kchip.merge.deep_PASS_autosomal.norm_setID.onlySNP.03snpqc.vcf.gz.ID",col_names = F)
wgs <- read_table("WGS_8062_joint_normMulti_recal_rmVQSR_LA_annotated_FIL_final_rmLQ_Phased_rmCHR_rmAC_17sample_reheaderid.rechr.hg19to38.setID.filter.vcf.gz.ID",col_names = F) %>% unique()
wgs <- read_table("WGS_8062_joint_normMulti_recal_rmVQSR_LA_annotated_FIL_final_rmLQ_Phased_rmCHR_rmAC_17sample_reheaderid.rechr.hg19to38.setID.filter2.vcf.gz.ID",col_names = F) %>% unique()
wgs <- read_table("WGS_8062_joint_normMulti_recal_rmVQSR_LA_annotated_FIL_final_rmLQ_Phased_rmCHR_rmAC_17sample_reheaderid.rechr.hg19to38.setID.filter3.vcf.gz.ID",col_names = F) %>% unique()

setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/00.datacheck/02.concordance/newcall_withWGS/")
wgs <- read_table("WGS.merge_17samples_reheader.autosomal.norm_setID.onlySNP.03snpqc.vcf.gz.ID",col_names = F) %>% unique()
head(nano)
head(revio)
head(wgs)
wgs %>% unique() %>%
table(is.na(wgs$X1))
wgs %>% mutate(chr = str_split_fixed(X1,"_",3)[,1]) %>% count(chr) -> a
sum(a$n)
library(ggVennDiagram)

x <- list(Nanopore=nano$X1,Revio=revio$X1,WGS=wgs$X1)
x
library(ggvenn)
ggvenn(x)
