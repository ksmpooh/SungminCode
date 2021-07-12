##### gnoad 비교, #### freq

##### gnoad 비교
setwd("~/Desktop/KCDC/long_read/analysis/01.variant.calling/")
front <- 28477797
df <- read.table("each/uniq.list.txt",header = T)
deepvariant <- read.table("HLA.merge_usingGLnexus_marker.QUAL.txt",header = T)

df$Real.POS <- front + df$POS
deepvariant$Real.POS <- front + deepvariant$pos
df$ID <- paste0(df$CHR,":",df$Real.POS,"_",df$REF,"/",df$ALT)
deepvariant$ID <- paste0(deepvariant$chr,":",deepvariant$Real.POS,"_",deepvariant$ref,"/",deepvariant$alt)

ref <- read.table("~/Desktop/KCDC/HLAimputation/IMPUTE4/Han.ref/ref.allele.with.ALT.NoHLA.txt",header = T)
ref$ID <- paste0(ref$chr,":",ref$hg19,"_",ref$ref,"/",ref$alt)

gnomad <- read.table("gnomad/gnomad.ATL.allele.frequency.txt",header = T)
gnomad$ID <- paste0(gnomad$CHROM,":",gnomad$POS,"_",gnomad$REF,"/",gnomad$ALT)


head(ref)
head(df)
head(deepvariant)
head(gnomad)


library(RColorBrewer)
myCol <- brewer.pal(4, "Pastel2")
library(VennDiagram)
venn.diagram(
  x=list(df$ID,deepvariant$ID,ref$ID,gnomad$ID),
  category.names = c("Variant calling","Multi-sample calling","Han ref","gnomAD"),
  filename = "Marker_compare_withGnomad.png",
  output = TRUE,
  #output = FALSE,
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
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
  cat.pos = c(0,5,0,0),
  cat.dist = c(-0.3, -0.3, 0.1,0.1),
  cat.fontfamily = "sans"
  #rotation = 1
)

myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x=list(df$ID,deepvariant$ID,gnomad$ID),
  category.names = c("Variant calling","Multi-sample calling","gnomAD"),
  filename = "Marker_compare_with_noHAN.Gnomad.png",
  output = TRUE,
  #output = FALSE,
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
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
  cat.pos = c(-10,10,0),
  cat.dist = c(0, 0, -0.4),
  cat.fontfamily = "sans",
  rotation = 1
)


#### freq
library(stringr)
setwd("~/Desktop/KCDC/long_read/analysis/01.variant.calling/vcf/freq/")
variant.call <- read.table("../../each/uniq.list.txt",header = T)
multi.call <- read.table("HLA.merge_usingGLnexus_frq.frq",header = T)

head(multi.call)
head(variant.call)

head(str_split_fixed(multi.call$SNP,"_",4)[,2])
multi.call$POS <- str_split_fixed(multi.call$SNP,"_",4)[,2]
front <- 28477797
#df$Real.POS <- front + df$POS
variant.call$Real.POS <-  front + variant.call$POS
multi.call$Real.POS <-  front + as.integer(multi.call$POS)

hist(multi.call$NCHROBS)
