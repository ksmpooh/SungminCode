#/BDATA/smkim/imputation.tool.check/OUTPUTs/info_maf

#IMPUTE4.chr1.ID_info_maf.buffer.txt     

#impute5.chr1.ID_info_maf.eagle_vcftovcf_pasing.usingKRG1KGP.txt
#impute5.chr1.ID_info_maf.eagle_plinktohap_phasing.usingKRG1KGP.txt
#impute5.chr1.ID_info_maf.shapeit4_pasing.usingKRG1KGP.txt

#minimac4.eagle_phasing.default_window.KRG1KGP.txt
#minimac4.eagle_phasing.use_window.KRG1KGP.txt

#minimac4.shapeit4_phasing.Defualt_window.KRG1KGP.txt
#minimac4.shapeit4_phasing.use_window.KRG1KGP.txt

library(ggplot2)
setwd("~/Desktop/KCDC/imputation/info_maf/")

head(impute5_1)
impute5_1<- read.table("impute5.chr1.ID_info_maf.eagle_vcftovcf_pasing.usingKRG1KGP.txt",header = T)
impute5_2<- read.table("impute5.chr1.ID_info_maf.eagle_plinktohap_phasing.usingKRG1KGP.txt",header = T)
impute5_3<- read.table("impute5.chr1.ID_info_maf.shapeit4_pasing.usingKRG1KGP.txt",header = T)
summary(impute5_1)
plot(impute5_1)

ggplot(data = impute5_1,axis(x='MAF',y='INFO',side = 3))+
  geom_line()
plot.new()
