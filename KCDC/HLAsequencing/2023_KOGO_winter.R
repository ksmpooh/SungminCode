################ 2023 KOGO winter Poster 

setwd("/Volumes/DATA/HLA_seq/winterKOGO_2023/")
setwd("~/")
library(ggplot2)
library(ggbreak)
library(tidyverse) 
library(ggpmisc)
library(ggpubr)
library(VennDiagram)
ggvenn::geom_venn()
long <- read.table("HLA.Longread.Seq.Q20.Deepvariant_Variantcalling.GLnexus_Jointcalling.sampleIDcheck_normQC.updateID.setID.variINFO.txt")
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

ref$ID <- paste0(ref$chr,":",ref$hg19,"_",ref$ref,"/",ref$alt)

head(ref)
head(long)
head(gnomad_eas)
head(gnomad)

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
library(VennDiagram)


venn.diagram(
  x=list(long$ID,gnomad$ID,gnomad_eas$ID),
  #category.names = c("1.HLAsequencing","2.HanREF","3.gnomAD_EAS"),
  category.names = c("KMHC","gnomAD","gnomAD_EAS"),
  filename = "long.gnomAD.DV.png",
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
head(long)
head(ref)
myCol <- brewer.pal(3, "Pastel2")
myCol[1:2]

venn.diagram(
  x=list(long$ID,ref$ID),
  #category.names = c("1.HLAsequencing","2.HanREF","3.gnomAD_EAS"),
  category.names = c("KMHC","HanREF"),
  filename = "long.Hanref.png",
  output = TRUE,
  
  imagetype="png",
  height = 600, 
  width = 600, 
  resolution = 400,
  compression = "lzw",
  
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],
  
  cex = .4,
  #fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 18),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
  #rotation = 1
)




exome <- read.table("~/Desktop/KCDC/long_read/analysis/03.check/gnomad/gnomad.exomes.HLA.AF.withEAS_withoutAFzero_withoutpoint_withoutEASzero_withoutpoint.txt")
genome <- read.table("~/Desktop/KCDC/long_read/analysis/03.check/gnomad/gnomad.genomes.HLA.AF.withEAS_withoutAFzero_withoutpoint_withoutEASzero_withoutpoint.txt")

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
head(long)
table(long$ID %in% df$ID)
table(long$gnomad_eas_check)
long %>% select(ID,V5) -> long_af
  
out <- merge(long_af,df[,c("ID","EAS")])
head(out)

head(out)
nrow(out)
lm(out$V5~out$EAS)
cor(out$V5,out$EAS)
table(out$V5,out$EAS)

ggplot(out, aes(y=EAS,x=V5)) +
  geom_point(size=0.3,shape=23,colour=rgb(0,0,1,0.6)) +
  geom_abline(slope=1, intercept=0,colour='red',size=2) + 
  labs(title = "Scatter plot of allele frequency",
       y= "gnomAD EAS",x="KMHC") + 
  theme(plot.title = element_text(hjust = 0.5,size = 20),panel.background = element_blank())
#  stat_poly_line()


#####concordance Test
# long vs KBA


#setwd("~/Desktop/KCDC/long_read/2022/VCF/Final.VCF/concordance/concordance_result/")
#setwd("~/")

#data <- read.table("/Volumes/DATA/HLA_seq/winterKOGO_2023/concordance_result/QulityMetrix_LongDV.KBA.txt", header = T)
data <- read.table("~/Desktop/KCDC/HLA_seq/winterKOGO_2023/QulityMetrix_LongDVnorm.KBA.txt", header = T)

dim(data)
data <- cbind(data, TP=(data$REF.REF + data$ALT_1.ALT_1 + data$ALT_2.ALT_2))
data <- cbind(data, FN=(data$ALT_1.REF + data$ALT_2.REF + data$ALT_2.ALT_1))
data <- cbind(data, FP=(data$REF.ALT_1 + data$REF.ALT_2 + data$ALT_1.ALT_2))

library(tidyverse)
data$Sensitivity <- data$TP/(data$TP+data$FN)
data$Precision <- data$TP/(data$TP+data$FP)
data$Accuracy <- data$TP/(data$TP+data$FP+data$FN)


data <- data %>% mutate(ID = paste0(chr,":",pos,"_",ref,"/",alt)) 
head(data)
head(long)

fun_mean <- function(x){
  return(data.frame(y=mean(x),label=mean(x,na.rm=T)))}

head(long)
data %>% inner_join(long) %>% 
  select(Sensitivity,Precision,Accuracy,V5) %>%
  arrange(V5) %>%
  mutate(MAF = ifelse(V5>=0.5,"0.5 <= MAF",ifelse(V5 < 0.5 & V5>=0.1,"0.1 <= MAF < 0.5","MAF < 0.1"))) %>%
  pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>%
  ggplot(aes(x=MAF,y=Val,fill=type)) +
  #geom_boxplot()
  geom_violin()
#s
  

data %>% inner_join(long) %>%
  select(Accuracy,V5) %>%
  mutate(MAF = ifelse(V5>=0.5,"0.5 <= MAF",ifelse(V5 < 0.5 & V5>=0.1,"0.1 <= MAF < 0.5","MAF < 0.1"))) %>%
  pivot_longer(cols = Accuracy,names_to = 'type',values_to = "Val") %>%
  ggplot(aes(x=factor(MAF,level =c("MAF < 0.1","0.1 <= MAF < 0.5","0.5 <= MAF")) ,y=Val,fill=MAF)) +
  geom_boxplot() + 
  geom_point(aes(col=MAF)) + 
  labs(x="MAF",y="",title="Accuracy") + 
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none",axis.text.x=element_blank()) -> p1


data %>% inner_join(long) %>% 
  select(Sensitivity,V5) %>%
  mutate(MAF = ifelse(V5>=0.5,"0.5 <= MAF",ifelse(V5 < 0.5 & V5>=0.1,"0.1 <= MAF < 0.5","MAF < 0.1"))) %>%
  pivot_longer(cols = Sensitivity,names_to = 'type',values_to = "Val") %>%
  ggplot(aes(x=factor(MAF,level =c("MAF < 0.1","0.1 <= MAF < 0.5","0.5 <= MAF")) ,y=Val,fill=MAF)) +
  geom_boxplot() + 
  geom_point(aes(col=MAF)) + 
  labs(x="MAF",y="",title="Sensitivity") + 
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none",axis.text.x=element_blank()) -> p2

data %>% inner_join(long) %>% na.omit() %>%
  select(Precision,V5) %>%
  mutate(MAF = ifelse(V5>=0.5,"0.5 <= MAF",ifelse(V5 < 0.5 & V5>=0.1,"0.1 <= MAF < 0.5","MAF < 0.1"))) %>%
  pivot_longer(cols = Precision,names_to = 'type',values_to = "Val") %>%
  ggplot(aes(x=factor(MAF,level =c("MAF < 0.1","0.1 <= MAF < 0.5","0.5 <= MAF")) ,y=Val,fill=MAF)) +
  geom_boxplot() + 
  geom_point(aes(col=MAF)) + 
  labs(x="MAF",y="",title="Precision") + 
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none",axis.text.x=element_blank()) -> p3




data %>% inner_join(long) %>%
  #data %>% inner_join(maf1) %>% 
  select(Precision,V5) %>%
  arrange(V5) %>%
  mutate(MAF = ifelse(V5>=0.5,"0.5 <= MAF",ifelse(V5 < 0.5 & V5>=0.1,"0.1 <= MAF < 0.5","MAF < 0.1"))) %>%
  pivot_longer(cols = Precision,names_to = 'type',values_to = "Val") %>%
  ggplot(aes(x=factor(MAF,level =c("MAF < 0.1","0.1 <= MAF < 0.5","0.5 <= MAF")) ,y=Val,fill=MAF)) +
  geom_boxplot() + 
  scale_fill_discrete(breaks=c("MAF < 0.1","0.1 <= MAF < 0.5","0.5 <= MAF")) ->  p4




p4 <- get_legend(p4)
p4
cowplot::plot_grid(p1,p2,p3,p4,ncol=4,labels = c("A","B","C"),rel_widths = c(1,1,1,0.5))



data %>% inner_join(long) %>% na.omit() %>%
  select(Sensitivity,Precision,Accuracy,V5) %>%
  mutate(MAF = ifelse(V5>=0.5,"0.5 <= MAF",ifelse(V5 < 0.5 & V5>=0.1,"0.1 <= MAF < 0.5","MAF < 0.1"))) %>%
  group_by(MAF) %>%
  summarise(Sensitivity = mean(Sensitivity),Precision = mean(Precision),Accuracy = mean(Accuracy)) #-> a

data %>% inner_join(long) %>% na.omit() %>%
  select(Sensitivity,Precision,Accuracy,V5) %>%
  mutate(MAF = ifelse(V5>=0.5,"0.5 <= MAF",ifelse(V5 < 0.5 & V5>=0.1,"0.1 <= MAF < 0.5","MAF < 0.1"))) %>%
  #group_by(MAF) %>%
  summarise(Sensitivity = mean(Sensitivity),Precision = mean(Precision),Accuracy = mean(Accuracy)) -> b
b$MAF = 'ALL'
head(b)
a %>%  rbind(b) -> b

head(ref)
nrow(long)/nrow(ref)

### VCF 
'
      1 ##fileformat=VCFv4.1
      2 #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
      3 6       28478882        6:28478882_T/G  T       G       .       .       .
      4 6       28478884        6:28478884_G/T  G       T       .       .       .
      5 6       28479657        6:28479657_G/C  G       C       .       .       .
      6 6       28479730        6:28479730_C/A  C       A       .       .       .
      7 6       28479796        6:28479796_A/T  A       T       .       .       .
'
long %>% select(V1,V2,ID,V3,V4) %>% 
  rename("CHROM" = V1,"POS" = V2,"REF" = V3,"ALT"=V4) %>%
  mutate(QUAL=".",FILTER=".",INFO=".") %>% #dim()
  write.table("~/Desktop/KCDC/HLA_seq/winterKOGO_2023/longDV.vari.INFO.forAnno.vcf",quote=F,row.names=F,col.names = T,sep = "\t")


long %>% select(V1,V2,ID,V3,V4) %>% 
  rename("CHROM" = V1,"POS" = V2,"REF" = V3,"ALT"=V4) %>%
  mutate(QUAL=".",FILTER=".",INFO=".") %>% 
  filter(!(ID %in% ref$ID)) %>% #dim()
  write.table("~/Desktop/KCDC/HLA_seq/winterKOGO_2023/longDV_noinHanREF.vari.INFO.forAnno.vcf",quote=F,row.names=F,col.names = T,sep = "\t")


## annotation result


library(tidyverse)
library(writexl)
setwd("~/Desktop/KCDC/HLA_seq/winterKOGO_2023/")

a <- read.csv("longDV.vari.INFO.forAnno_VEP.vcf_summary_dataprocessing.txt")
b <- read.csv("longDV_noinHanREF.vari.INFO.forAnno_VEP.vcf_summary_dataprocessing.txt")

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
writexl::write_xlsx(anno_out,"anno_Result.xlsx")
anno_out %>% pivot_longer(ALL_count,)

ggplot(anno_out,aes(x=type,color=noInHanREF_count)) +
  geom_histogram()
