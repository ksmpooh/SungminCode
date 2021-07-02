##### HLA seqeuncing
library(ggplot2)
setwd("~/Desktop/KCDC/long_read/analysis/01.variant.calling/")
front <- 28477797
df <- read.table("each/uniq.list.txt",header = T)
head(df)
df$Real.POS <- front + df$POS

hist(df$N,breaks = 11,col = rgb(0,0,1,0.3),xlab = "# of duplicates SNP",
     main = "Deepvariant Variant calling marker count (90310)")




ref <- read.table("~/Desktop/KCDC/HLAimputation/IMPUTE4/Han.ref/ref.allele.txt",header = T)
head(ref)


deepvariant <- read.table("HLA.merge_usingGLnexus_marker.QUAL.txt",header = T)
gatk <- read.table("HLA.GATK.GenotypeGVCFs.12sample.Using.GATK.marker.info.txt",header = T)

head(deepvariant)
tail(deepvariant)
deepvariant$Real.POS <- front + deepvariant$pos
head(gatk)
gatk$Real.POS <- front + gatk$pos
plot(deepvariant$QUAL)

par(mfrow=c(2,1))
dev.off()
deepvariant[deepvariant$Real.POS %in% df$Real.POS,]
plot(x = deepvariant$Real.POS,y = deepvariant$QUAL,col=rgb(0,0,0,0.3), cex=0.5, pch=20)
plot(x = gatk$Real.POS,y=gatk$QUAL,col=rgb(0,0,0,0.3), cex=0.5, pch=20)
points(gatk[gatk$Real.POS %in% df$Real.POS,]$Real.POS,
       gatk[gatk$Real.POS %in% df$Real.POS,]$QUAL,
       col = rgb(1,0,0,0.1), cex = 0.5 , pch = 20)

#DP :Approximate read depth; some reads may have been filtered
#MQ : RMS Mapping Quality
#QD : Variant Confidence/Quality by Depth"

plot(gatk$DP,gatk$MQ)

plot(gatk[,5:ncol(gatk)-1])
