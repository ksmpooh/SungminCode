### short vs han 이거 꼭 저장해야함!!!1 KOGO 용

##### HLA seqeuncing
library(ggplot2)
library(ggbreak) 

setwd("~/Desktop/KCDC/long_read/2022/VCF/onlySNP/DV/")

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



library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
library(VennDiagram)

head(long)
head(ref)
head(gnomad)
head(gnomad_eas)

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


venn.diagram(
  x=list(df2$ID,ref$ID,gnomad_eas$ID),
  category.names = c("1.HLAsequencing","2.HanREF","3.gnomAD_EAS"),
  filename = "gnomAD-EAS_1types_calling.png",
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

af_check <- af_check %>%  filter(gnomad_eas_check == "Yes")
cor(af_check$EAS,af_check$V5) # [1] 0.8531091


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

#write.csv(anno_out,"anno_Result.csv")
anno_out %>% pivot_longer(ALL_count,)

ggplot(anno_out,aes(x=type,color=noInHanREF_count)) +
  geom_histogram()
