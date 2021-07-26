library(RColorBrewer)
library(VennDiagram)

####header 
setwd("/Users/ksmpooh/Desktop/KCDC/long_read/analysis/03.check/01.multi/ID")
setwd("/Users/ksmpooh/Desktop/KCDC/long_read/analysis/03.check/01.multi/compare_3tpyes/ID/")
front <- 28477797
head(df1)
head(df2)
head(df3)
head(gnomad)
df1 <- read.table("01.2020HLAseq.pass_vcfmerge_devidedMultiAllelic_ID.txt")
df2 <- read.table("02.2020HLAseq_bcftoolsMerge_PASS_devidedMultiAllelic_ID.txt")
df3 <- read.table("03.merge_usingGLnexus_devidedMultiAllelic_ID.txt")
gnomad <- read.table("~/Desktop/KCDC/long_read/analysis/03.check/gnomad/gnomad.ATL.allele.frequency.txt",header = T)

colnames(df1) <- c("CHROM","POS","REF","ALT")
colnames(df2) <- c("CHROM","POS","REF","ALT")
colnames(df3) <- c("CHROM","POS","REF","ALT")

df1$Real.POS <- front + df1$POS
df2$Real.POS <- front + df2$POS
df3$Real.POS <- front + df3$POS

df1$ID <- paste0(df1$CHR,":",df1$Real.POS,"_",df1$REF,"/",df1$ALT)
df2$ID <- paste0(df2$CHR,":",df2$Real.POS,"_",df2$REF,"/",df2$ALT)
df3$ID <- paste0(df3$CHR,":",df3$Real.POS,"_",df3$REF,"/",df3$ALT)
gnomad$ID <- paste0(gnomad$CHR,":",gnomad$POS,"_",gnomad$REF,"/",gnomad$ALT)
###




myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x=list(df1$ID,df2$ID,df3$ID),
  category.names = c("Variant.Calling_v1","Variant.Calling_v2","Multi-sample.Calling"),
  filename = "3types_compare_allmarker.png",
  output = TRUE,
  
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw",
  
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  cex = .4,
  fontface = "bold",
  fontfamily = "sans",
  
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 20, 120),
  cat.dist = c(0.055, 0.05, 0.05),
  cat.fontfamily = "sans",
  rotation = 1
)




venn.diagram(
  x=list(df1$ID,df3$ID,gnomad$ID),
  category.names = c("Variant.Calling","Multi-sample.Calling","gnomAD"),
  filename = "2types_compare_allmarker_withGnomad.png",
  output = TRUE,
  
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  cex = .4,
  fontface = "bold",
  fontfamily = "sans",
  
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 20, 120),
  cat.dist = c(0.055, 0.05, 0.05),
  cat.fontfamily = "sans",
  rotation = 1
)



myCol <- brewer.pal(4, "Pastel2")
venn.diagram(
  x=list(df1$ID,df2$ID,df3$ID,gnomad$ID),
  category.names = c("Variant.Calling_v1","Variant.Calling_v2","Multi-sample.Calling","gnomAD"),
  filename = "3types_compare_allmarker_withGnomad.png",
  output = TRUE,
  
  imagetype="png" ,
  height = 600 , 
  width = 600, 
  resolution = 300,
  compression = "lzw",
  
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  cex = .4,
  fontface = "bold",
  fontfamily = "sans",
  
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0, 0,0),
  cat.dist = c(-0.3, -0.3, 0.1,0.1),
  cat.fontfamily = "sans"
  #rotation = 1
)




##### gnomAD with Han
setwd("/Users/ksmpooh/Desktop/KCDC/long_read/analysis/03.check/01.multi/compare_3tpyes/ID/")
front <- 28477797
head(df1)
head(df2)
head(gnomad)
df1 <- read.table("01.2020HLAseq.pass_vcfmerge_devidedMultiAllelic_ID.txt")
df2 <- read.table("03.merge_usingGLnexus_devidedMultiAllelic_ID.txt")
#gnomad <- read.table("~/Desktop/KCDC/long_read/analysis/03.check/gnomad/gnomad.ATL.allele.frequency.txt",header = T)
gnomad <- read.table("~/Desktop/KCDC/long_read/analysis/03.check/gnomad/gnomad.genome_exome_HLA_AF.allele_set.txt")
gnomad_eas <- read.table("~/Desktop/KCDC/long_read/analysis/03.check/gnomad/gnomad.genome_exome_HLA_EAS.allele_set.txt")

colnames(df1) <- c("CHROM","POS","REF","ALT")
colnames(df2) <- c("CHROM","POS","REF","ALT")
colnames(gnomad) <- c("CHROM","POS","REF","ALT")
colnames(gnomad_eas) <- c("CHROM","POS","REF","ALT")

df1$Real.POS <- front + df1$POS
df2$Real.POS <- front + df2$POS

df1$ID <- paste0(df1$CHR,":",df1$Real.POS,"_",df1$REF,"/",df1$ALT)
df2$ID <- paste0(df2$CHR,":",df2$Real.POS,"_",df2$REF,"/",df2$ALT)
gnomad$ID <- paste0(gnomad$CHR,":",gnomad$POS,"_",gnomad$REF,"/",gnomad$ALT)
gnomad_eas$ID <- paste0(gnomad_eas$CHR,":",gnomad_eas$POS,"_",gnomad_eas$REF,"/",gnomad_eas$ALT)

setwd("/Users/ksmpooh/Desktop/KCDC/long_read/analysis/03.check/gnomad/")
#gnomad <- read.table("~/Desktop/KCDC/long_read/analysis/03.check/gnomad/gnomad.ATL.allele.frequency.txt",header = T)
ref <- read.table("~/Desktop/KCDC/HLAimputation/IMPUTE4/Han.ref/ref.allele.with.ALT.NoHLA.txt",header = T)
ref$ID <- paste0(ref$chr,":",ref$hg19,"_",ref$ref,"/",ref$alt)

head(gnomad)
head(gnomad_eas)
head(ref)


myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x=list(ref$ID,gnomad$ID,gnomad_eas$ID),
  category.names = c("Han ref.","gnomAD","gnomAD_EAS"),
  filename = "Han_gnomAD_gnomAD-EAS.png",
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
  cat.pos = c(-27, 20, 120),
  cat.dist = c(0.055, 0.05, 0.05),
  cat.fontfamily = "sans",
  rotation = 1
)


venn.diagram(
  x=list(df1$ID,df2$ID,gnomad$ID),
  category.names = c("Variant Calling","Multi-sample Calling","gnomAD"),
  filename = "gnomAD-All_2types_calling.png",
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
  cat.pos = c(-27, 18, 140),
  cat.dist = c(0.055, 0.055, 0.05),
  cat.fontfamily = "sans",
  rotation = 1
)

venn.diagram(
  x=list(df1$ID,df2$ID,gnomad_eas$ID),
  category.names = c("Variant Calling","Multi-sample Calling","gnomAD_EAS"),
  filename = "gnomAD-EAS_2types_calling.png",
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
  cat.pos = c(-27, 18, 140),
  cat.dist = c(0.055, 0.055, 0.05),
  cat.fontfamily = "sans",
  rotation = 1
)


## with han 
myCol <- brewer.pal(4, "Pastel2")
venn.diagram(
  x=list(df1$ID,df2$ID,gnomad$ID,gnomad_eas$ID),
  category.names = c("Variant.Calling","Multi-sample.Calling","gnomAD","gnomAD_EAS"),
  filename = "2types.calling_2type.gnomAD.png",
  output = TRUE,
  
  imagetype="png" ,
  height = 600 , 
  width = 600, 
  resolution = 300,
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
  cat.pos = c(0, 3, 0,0),
  cat.dist = c(-0.3, -0.3, 0.1,0.1),
  cat.fontfamily = "sans"
  #rotation = 1
)

venn.diagram(
  x=list(df1$ID,df2$ID,gnomad$ID,ref$ID),
  category.names = c("Variant.Calling","Multi-sample.Calling","gnomAD","Han ref"),
  filename = "2types.calling_gnomAD-ALL_Han.png",
  output = TRUE,
  
  imagetype="png" ,
  height = 600 , 
  width = 600, 
  resolution = 300,
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
  cat.pos = c(0, 3, 0,0),
  cat.dist = c(-0.3, -0.3, 0.1,0.1),
  cat.fontfamily = "sans"
  #rotation = 1
)

venn.diagram(
  x=list(df1$ID,df2$ID,gnomad_eas$ID,ref$ID),
  category.names = c("Variant.Calling","Multi-sample.Calling","gnomAD_EAS","Han ref"),
  filename = "types.calling_gnomAD-EAS_Han.png",
  output = TRUE,
  
  imagetype="png" ,
  height = 600 , 
  width = 600, 
  resolution = 300,
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
  cat.pos = c(0, 3, 0,0),
  cat.dist = c(-0.3, -0.3, 0.1,0.1),
  cat.fontfamily = "sans"
  #rotation = 1
)
