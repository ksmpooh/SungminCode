##### HLA seqeuncing
library(ggplot2)
library(ggbreak) 

setwd("~/Desktop/KCDC/long_read/analysis/01.variant.calling/")
front <- 28477797
df <- read.table("each/uniq.list.txt",header = T)
deepvariant <- read.table("HLA.merge_usingGLnexus_marker.QUAL.txt",header = T)

head(df)
head(deepvariant)
df$Real.POS <- front + df$POS
deepvariant$Real.POS <- front + deepvariant$pos
df$ID <- paste0(df$CHR,":",df$Real.POS,"_",df$REF,"/",df$ALT)
deepvariant$ID <- paste0(deepvariant$chr,":",deepvariant$Real.POS,"_",deepvariant$ref,"/",deepvariant$alt)

ref <- read.table("~/Desktop/KCDC/HLAimputation/IMPUTE4/Han.ref/ref.allele.with.ALT.NoHLA.txt",header = T)
head(ref)
##############
hist(df$N,breaks = 12,col = rgb(0,0,1,0.3),xlab = "# of duplicates SNP",
     main = "Deepvariant Variant calling marker count (90310)")

gatk <- read.table("HLA.GATK.GenotypeGVCFs.12sample.Using.GATK.marker.info.txt",header = T)

deepvariant$Real.POS <- front + deepvariant$pos
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



########marker check

head(deepvariant)
head(df)
df$new_N <- df$N
df[df$new_N == 1,]$new_N <- 0.9
hist(df$new_N,col = rgb(0,0,1,0.3),xlab = "# of duplicates SNP",
#     xlim=c(0,12), 
     breaks = 12,
     main = "Deepvariant Variant calling marker count (PASS : 90,310)")
abline(h=10000, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=5000, col=rgb(0,0,1,0.5), lty=3, lwd=2)

head(deepvariant)
head(df)
table(df$N)
hist(df$Real.POS)
hist(deepvariant$Real.POS)

table(deepvariant$ID %in% df$ID)
table(df$ID %in% deepvariant$ID)
df2 <- df[df$ID %in% deepvariant$ID,]
df2$new_N <- df2$N
df2[df2$new_N == 1,]$new_N <- 0.9
table(df2$new_N)
hist(df2$new_N,col = rgb(0,0,1,0.3),xlab = "# of duplicates SNP",
     #     xlim=c(0,12), 
     breaks = 12,
     main = "Multi-sample Calling marker in variant calling marker (65,908)")
abline(h=10000, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=5000, col=rgb(0,0,1,0.5), lty=3, lwd=2)

#####
head(df)
head(df2)
df$type <- "Variant.call"
df2$type <- "Multi.call"
out <-rbind(df,df2)
hist(out$N)
str(out)
ggplot(out,aes(x=N,color=type)) +
        geom_bar(fill = 'white')
        #geom_histogram(fill = "white",binwidth = 1,main = "# of duplicated Marker",xlim=c(0,12))
barplot(table(out$type,out$N), main = "Multi.call vs Varinat.call in Deepvariant",
        xlab = "# of duplicated Sample",
        ylab = "# of Marker",
#        breaks = 12,
        beside = T,col=c(rgb(0,1,0,1),rgb(0,0,1,1)))
legend(19,20000,legend = c("Multi calling","Variant calling"),
       col = c(rgb(0,1,0,1),rgb(0,0,1,1)),lty = 1,cex = 1.2,
       box.lty = 0)
######



### han ref 와 비교

head(ref)

ref$ID <- paste0(ref$chr,":",ref$hg19,"_",ref$ref,"/",ref$alt)


table(df$ID %in% ref$ID)
table(deepvariant$ID %in% ref$ID)


library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
library(VennDiagram)
venn.diagram(
        x=list(df$ID,deepvariant$ID,ref$ID),
        category.names = c("Variant calling","Multi-sample calling","Han ref"),
        filename = "Marker_compare.png",
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
