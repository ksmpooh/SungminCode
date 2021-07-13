library(RColorBrewer)
library(VennDiagram)

####header 
setwd("/Users/ksmpooh/Desktop/KCDC/long_read/analysis/03.check/01.multi/ID")
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
