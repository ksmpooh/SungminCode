#### variant type compare

setwd("~/Desktop/KCDC/long_read/analysis/03.check/04.withoutMulti.allele/inter/anno/")

type1 <- read.table("0000.variant.type.txt",header = F)
type2 <- read.table("0001.variant.type.txt",header = F)

inter <- read.table("0002.variant.type.txt",header = F)
head(inter)
colnames(type1) <-c("CHROM","POS","REF","ALT","Type")
colnames(type2) <-c("CHROM","POS","REF","ALT","Type")
colnames(inter) <-c("CHROM","POS","REF","ALT","Type")


table(type1$Type)
table(type2$Type)

table(inter$Type)

#out <- data.frame()
a <- as.data.frame(t(data.frame(table(inter$Type))))
colnames(a) <- a["Var1",]
a$type <- "Intersect"

out <- a["Freq",]

a <- as.data.frame(t(data.frame(table(type1$Type))))
colnames(a) <- a["Var1",]
a$type <- "Variant.call"

out <- merge(out,a["Freq",],all = T)

a <- as.data.frame(t(data.frame(table(type2$Type))))
colnames(a) <- a["Var1",]
a$type <- "Multi.call"

out <- merge(out,a["Freq",],all = T)

head(a)
head(out)
colnames(out)

out <- out[,c(15,1:14,16:ncol(out))]
out <- out[c(2,3,1),]


write.table(out,"Variant.type.freq.txt",col.names = T,row.names = F,quote = F,sep = "\t")
