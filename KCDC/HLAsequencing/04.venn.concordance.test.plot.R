### short vs han 이거 꼭 저장해야함!!!1 KOGO 용

##### HLA seqeuncing
library(ggplot2)
library(ggbreak)
library(tidyverse) 

setwd("~/Desktop/KCDC/long_read/2022/VCF/onlySNP/DV/")
setwd("~/")

long <-read.table("longread.SNP.txt")
short <-read.table("shortread.SNP.txt")

head(long)

ref <- read.table("~/Desktop/KCDC/HLAimputation/IMPUTE4/Han.ref/ref.allele.with.ALT.NoHLA.txt",header = T)
head(ref)

gnomad <- read.table("~/Desktop/KCDC/long_read/analysis/03.check/gnomad/gnomad.genome_exome_HLA_AF.allele_set.txt")
gnomad_eas <- read.table("~/Desktop/KCDC/long_read/analysis/03.check/gnomad/gnomad.genome_exome_HLA_EAS.allele_set.txt")

colnames(gnomad) <- c("CHROM","POS","REF","ALT")
colnames(gnomad_eas) <- c("CHROM","POS","REF","ALT")

gnomad$ID <- paste0(gnomad$CHR,":",gnomad$POS,"_",gnomad$REF,"/",gnomad$ALT)
gnomad_eas$ID <- paste0(gnomad_eas$CHR,":",gnomad_eas$POS,"_",gnomad_eas$REF,"/",gnomad_eas$ALT)

long$ID <- paste0("6:",long$V2,"_",long$V3,"/",long$V4)
short$ID <- paste0("6:",short$V2,"_",short$V3,"/",short$V4)

ref$ID <- paste0(ref$chr,":",ref$hg19,"_",ref$ref,"/",ref$alt)


table(long$ID %in% ref$ID)

dim(ref)

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
library(VennDiagram)

head(long)
head(ref)
dim(ref)
head(gnomad)
head(gnomad_eas)
dim(gnomad)
dim(gnomad_eas)


table(long$ID %in% ref$ID)
#FALSE  TRUE 
#52070 21098 
nrow(long)/nrow(short)
nrow(long) - nrow(ref)
(nrow(long)/nrow(ref))
(nrow(short)/nrow(ref))

table(long$ID %in% gnomad$ID)
#FALSE  TRUE 
#20413 52755 

table(short$ID %in% ref$ID) 
#FALSE  TRUE 
#42192 20889 
20889/nrow(ref)
table(short$ID %in% gnomad_eas$ID)
#FALSE  TRUE 
#17379 45702 
table(short$ID %in% gnomad$ID)
#FALSE  TRUE 
#15798 47283
52755/47283

short1 <- short[!(short$ID %in% ref$ID),]
table(short1$ID %in% gnomad_eas$ID)
#FALSE  TRUE 
#17224 24968 
short1$gnomad_eas_check <- 'No'
head(short1)
short1[(short1$ID %in% gnomad_eas$ID),]$gnomad_eas_check <- 'Yes'
dim(short1)
head(short1)
colnames(short1)<- c("contig","pos","ref","alt","ID","gnomad_eas_check")
#write.table(short1,"~/Desktop/KCDC/long_read/2022/forKOGO2022/shortredad_genoad_eas.check.txt",col.names = T,row.names = F,quote = F,sep = "\t")
head(gnomad_eas)

setwd("~/Desktop/KCDC/long_read/2022/figure/")
head(short)
head(long)

venn.diagram(
  x=list(long$ID,short$ID,ref$ID),
  #category.names = c("1.HLAsequencing","2.HanREF","3.gnomAD_EAS"),
  category.names = c("Longread","Shortread","HanREF"),
  filename = "long.short.HanREF.DV.png",
  output = TRUE,
  
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 400,
  compression = "lzw",
  
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  cex = .4,
  #fontface = "bold",
  fontfamily = "sans",
  
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 18, 140),
  cat.dist = c(0.055, 0.055, 0.05),
  cat.fontfamily = "sans",
  rotation = 1
)

venn.diagram(
  x=list(long$ID,short$ID,gnomad$ID),
  #category.names = c("1.HLAsequencing","2.HanREF","3.gnomAD_EAS"),
  category.names = c("Longread","Shortread","gnomAD"),
  filename = "long.short.gnomAD.DV.png",
  output = TRUE,
  
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 400,
  compression = "lzw",
  
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  cex = .4,
  #fontface = "bold",
  fontfamily = "sans",
  
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 18, 140),
  cat.dist = c(0.055, 0.055, 0.05),
  cat.fontfamily = "sans",
  rotation = 1
)

venn.diagram(
  x=list(long$ID,short$ID,gnomad_eas$ID),
  #category.names = c("1.HLAsequencing","2.HanREF","3.gnomAD_EAS"),
  category.names = c("Longread","Shortread","gnomAD_EAS"),
  filename = "long.short.gnomAD_EAS.DV.png",
  output = TRUE,
  
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 400,
  compression = "lzw",
  
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  cex = .4,
  #fontface = "bold",
  fontfamily = "sans",
  
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 18, 140),
  cat.dist = c(0.055, 0.055, 0.05),
  cat.fontfamily = "sans",
  rotation = 1
)





venn.diagram(
  x=list(short$ID,ref$ID,gnomad_eas$ID),
  #category.names = c("1.HLAsequencing","2.HanREF","3.gnomAD_EAS"),
  category.names = c("KMHC","HanREF","gnomAD_EAS"),
  filename = "gnomAD-EAS_types_calling.png",
  output = TRUE,
  
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 400,
  compression = "lzw",
  
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  cex = .4,
  #fontface = "bold",
  fontfamily = "sans",
  
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 18, 140),
  cat.dist = c(0.055, 0.055, 0.05),
  cat.fontfamily = "sans",
  rotation = 1
)

venn.diagram(
  x=list(short$ID,ref$ID,gnomad$ID),
  category.names = c("KMHC","HanREF","gnomAD"),
  filename = "gnomAD_types_calling.png",
  output = TRUE,
  
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 400,
  compression = "lzw",
  
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  cex = .4,
  #fontface = "bold",
  fontfamily = "sans",
  
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 18, 140),
  cat.dist = c(0.055, 0.055, 0.05),
  cat.fontfamily = "sans",
  rotation = 1
)


19778 + 123085
19832 + 123192

136222 - 19778
136222 - 19787

exome <- read.table("~/Desktop/KCDC/long_read/analysis/03.check/gnomad/gnomad.exomes.HLA.AF.withEAS_withoutAFzero_withoutpoint_withoutEASzero_withoutpoint.txt")
genome <- read.table("~/Desktop/KCDC/long_read/analysis/03.check/gnomad/gnomad.genomes.HLA.AF.withEAS_withoutAFzero_withoutpoint_withoutEASzero_withoutpoint.txt")

#exome <- read.table("~/Desktop/KCDC/long_read/analysis/03.check/gnomad/gnomad.exomes.HLA.AF.withEAS.txt")
#genome <- read.table("~/Desktop/KCDC/long_read/analysis/03.check/gnomad/gnomad.genomes.HLA.AF.withEAS.txt")
head(exome)
head(genome)
table(genome$V2 %in% exome$V2)

colnames(exome) <- c("CHROM","POS","REF","ALT","AF","EAS")
colnames(genome) <- c("CHROM","POS","REF","ALT","AF","EAS")

exome$ID <- paste0(exome$CHR,":",exome$POS,"_",exome$REF,"/",exome$ALT)
genome$ID <- paste0(genome$CHR,":",genome$POS,"_",genome$REF,"/",genome$ALT)
table(genome$ID %in% exome$ID)


df <- rbind(genome,exome[!(exome$ID %in% genome$ID),])


head(df)
head(short1)
table(short1$ID %in% df$ID)
table(short1$gnomad_eas_check)

out <- merge(short1,df[,c("ID","EAS")])
head(out)
short_af <- read.table("~/Desktop/KCDC/long_read/2022/forKOGO2022/short.DV.AF.txt")
head(short_af)
short_af$ID <- paste0("6:",short_af$V2,"_",short_af$V3,"/",short_af$V4)

af_check <- merge(out[,c("ID","EAS","gnomad_eas_check")],short_af[,c("ID","V5")])
head(af_check)

library(tidyverse)
library(ggpmisc)
library(ggpubr)


af_check <- af_check %>%  filter(gnomad_eas_check == "Yes")
cor(af_check$EAS,af_check$V5) # [1] 0.8531091
head(af_check)
lm(af_check$V5~af_check$EAS)

ggplot(af_check, aes(y=EAS,x=V5)) +
  geom_point(size=0.3,shape=23,colour=rgb(0,0,1,0.6)) +
  geom_abline(slope=1, intercept=0,colour='red',size=2) + 
  labs(title = "Scatter plot of allele frequency",
       y= "gnomAD EAS",x="KMHC") + 
  theme(plot.title = element_text(hjust = 0.5,size = 20),panel.background = element_blank())
#  stat_poly_line()

library(ggpmisc)

dev.off()

## after annotation

setwd("/Users/ksmpooh/Desktop/KCDC/long_read/2022/forKOGO2022/anno")

a <- read.csv("HLA.shortread.DV_VEP.vcf_summary_dataprocessing.txt")
b <- read.csv("HLA.shortread.DV.noHanREF_VEP.vcf_summary_dataprocessing.txt")


a <- a %>% filter(grepl("most",index)) %>% mutate(type = str_split_fixed(index,'-',2)[,2]) 
b<- b %>% filter(grepl("most",index)) %>% mutate(type = str_split_fixed(index,'-',2)[,2])
head(a)
head(b)
anno_out <- merge(a[,c(3,2)],b[,c(3,2)])
head(anno_out)
str(anno_out)
anno_out$noInHanREF_count <- as.numeric(anno_out$noInHanREF_count)
anno_out$ALL_count <- as.numeric(anno_out$ALL_count)
anno_out$ALL_per <- anno_out$ALL_count/sum(anno_out$ALL_count) * 100
anno_out$noInHanREF_per <- anno_out$noInHanREF_count/sum(anno_out$noInHanREF_count) * 100
anno_out <- anno_out[,c(1,2,4,3,5)]
#write.csv(anno_out,"anno_Result.csv")
anno_out %>% pivot_longer(ALL_count,)

ggplot(anno_out,aes(x=type,color=noInHanREF_count)) +
  geom_histogram()



#### concordance test

#setwd("/Users/ksmpooh/Desktop/KCDC/long_read/2022/VCF/Final.VCF/onlySNP/DV/concordance")
#setwd("/Volumes/DATA/HLA_seq_2021/forKOGO/chip_short_concordance/inter/")
setwd("~/Desktop/KCDC/long_read/2022/VCF/Final.VCF/concordance/concordance_result/")
#setwd("~/")

file_list <- list.files(pattern = ".txt")
file_list
data <- read.table(file_list[1],header = T)
head(data)
#data <- read.table("QulityMetricx_forKOGO.txt", header = T)


data <- cbind(data, TP=(data$REF.REF + data$ALT_1.ALT_1 + data$ALT_2.ALT_2))
data <- cbind(data, FN=(data$ALT_1.REF + data$ALT_2.REF + data$ALT_2.ALT_1))
data <- cbind(data, FP=(data$REF.ALT_1 + data$REF.ALT_2 + data$ALT_1.ALT_2))

library(tidyverse)
data$Sensitivity <- data$TP/(data$TP+data$FN)
data$Precision <- data$TP/(data$TP+data$FP)
data$Accuracy <- data$TP/(data$TP+data$FP+data$FN)


data %>% #inner_join(maf1) %>% 
  select(Sensitivity,Precision,Accuracy) %>% 
  pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>%
  ggplot(aes(x=type,y=Val,fill=type)) +
  geom_boxplot()

data %>% #inner_join(maf1) %>% 
  select(Sensitivity,Precision,Accuracy) %>% 
  pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% 
  mutate(title = "hi")
  
str_remove(str_remove(file_list[1],".txt"),"QulityMetricx_")
 

data %>% #inner_join(maf1) %>%  
  select(Sensitivity,Precision,Accuracy) %>% #head()
  summarise(Sensitivity=mean(Sensitivity))

  #pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% head()
  group_by(type) %>%
  summarise(mean=mean(Val))


df <- data %>% #inner_join(maf1) %>% 
  select(pos,Sensitivity,Precision,Accuracy) %>% 
  pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% 
  mutate(title = str_remove(str_remove(file_list[1],".txt"),"QulityMetricx_"))
head(df)
file_list[-1]

for (i in file_list[-1]) {
  data <- read.table(i,header = T)
  data <- cbind(data, TP=(data$REF.REF + data$ALT_1.ALT_1 + data$ALT_2.ALT_2))
  data <- cbind(data, FN=(data$ALT_1.REF + data$ALT_2.REF + data$ALT_2.ALT_1))
  data <- cbind(data, FP=(data$REF.ALT_1 + data$REF.ALT_2 + data$ALT_1.ALT_2))

  data$Sensitivity <- data$TP/(data$TP+data$FN)
  data$Precision <- data$TP/(data$TP+data$FP)
  data$Accuracy <- data$TP/(data$TP+data$FP+data$FN)
  
  a <- data %>% #inner_join(maf1) %>% 
    select(pos,Sensitivity,Precision,Accuracy) %>% 
    pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% 
    mutate(title = str_remove(str_remove(i,".txt"),"QulityMetricx_"))
  df <- rbind(df,a)
}
df %>% mutate(title = recode(title,"DV_kbalong"="DV : KBA vs Long","DV_kbashort" ="DV : KBA vs Short",
                        "DV_longshort" = "DV : Long vs Short",
                        "GATK_kbalong" = "GATK : KBA vs Long",
                        "GATK_kbashort" = "GATK : KBA vs Short",
                        "GATK_longshort" = "GATK : Long vs Short",
                        "longread" = "Longread : DV vs GATK",
                        "shortread" = "Shortread : DV vs GATK")) %>%
  ggplot(aes(x=type,y=Val,color=type)) + 
  geom_boxplot(outlier.size = 0.1) + 
  facet_wrap(~title,ncol=4) + 
  labs(x= NULL,y=NULL,colour=NULL,title="Concordance Test Result by SNP") +
  theme(plot.title = element_text(hjust = 0.5))

head(df)
table(df$title)/3


setwd("~/Desktop/KCDC/long_read/2022/VCF/Final.VCF/concordance/")
file_list <- list.files(pattern = "by.sample.txt")
data <- read.table(file_list[1],header = T)

data <- cbind(data, TP=(data$REF.REF + data$ALT_1.ALT_1 + data$ALT_2.ALT_2))
data <- cbind(data, FN=(data$ALT_1.REF + data$ALT_2.REF + data$ALT_2.ALT_1))
data <- cbind(data, FP=(data$REF.ALT_1 + data$REF.ALT_2 + data$ALT_1.ALT_2))


data$Sensitivity <- data$TP/(data$TP+data$FN)
data$Precision <- data$TP/(data$TP+data$FP)
data$Accuracy <- data$TP/(data$TP+data$FP+data$FN)

head(data)


data %>% #inner_join(maf1) %>% 
  select(sample,Sensitivity,Precision,Accuracy) %>% 
  pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% 
  mutate(title = str_remove(str_remove(file_list[1],".by_sample.txt"),"QulityMetricx_")) %>%
  ggplot(aes(x=type,y=Val,fil=title)) +
  geom_boxplot()
  
df <- data %>% #inner_join(maf1) %>% 
  select(sample,Sensitivity,Precision,Accuracy) %>% 
  pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% 
  mutate(title = str_remove(str_remove(file_list[1],".by_sample.txt"),"QulityMetricx_"))

file_list[-1]
for (i in file_list[-1]) {
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
    mutate(title = str_remove(str_remove(i,".by_sample.txt"),"QulityMetricx_"))
  df <- rbind(df,a)
}
head(df)
table(df$title)


means <- aggregate(Val ~  type, df, mean)

df %>% mutate(title = recode(title,"DV_kbalong"="DV : KBA vs Long","DV_kbashort" ="DV : KBA vs Short",
                             "DV_longshort" = "DV : Long vs Short",
                             "GATK_kbalong" = "GATK : KBA vs Long",
                             "GATK_kbashort" = "GATK : KBA vs Short",
                             "GATK_longshort" = "GATK : Long vs Short",
                             "longread" = "Longread : DV vs GATK",
                             "shortread" = "Shortread : DV vs GATK")) %>%
  ggplot(aes(x=type,y=Val,color=type)) + 
  geom_boxplot() + 
  facet_wrap(~title,ncol=4) + 
  labs(x= NULL,y=NULL,colour=NULL,title="Concordance Test Result by Sample") +
  theme(plot.title = element_text(hjust = 0.5))

  
a<- df %>% mutate(title = recode(title,"DV_kbalong"="DV : KBA vs Long","DV_kbashort" ="DV : KBA vs Short",
                             "DV_longshort" = "DV : Long vs Short",
                             "GATK_kbalong" = "GATK : KBA vs Long",
                             "GATK_kbashort" = "GATK : KBA vs Short",
                             "GATK_longshort" = "GATK : Long vs Short",
                             "longread" = "Longread : DV vs GATK",
                             "shortread" = "Shortread : DV vs GATK")) %>%
  group_by(title,type) %>%
  summarise(mean=mean(Val)) %>%
  pivot_wider(names_from = type,values_from = mean)
head(df)
head(a)

