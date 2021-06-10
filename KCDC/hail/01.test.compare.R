setwd("~/Desktop/KCDC/hail/")

d1 <- read.table("gaws.output.tsv",header = T)
head(d1)
d1 <- d1[,c(1,3,4,5)]
summary(d1)

d2 <- read.table("JG.KR.ESRD.association.sub_Total.b.firth.chr12.132000001_133841815.epacts",header = T)
head(d2)
d2$locus <- paste0(d2$CHROM,":",d2$BEGIN)
d2 <- d2[,c("locus","BETA","CHISQ","PVALUE")]
colnames(d1) <- colnames(d2)
head(d1)
head(d2)

library(stringr)

d1$POS <- str_split_fixed(d1$locus,":",2)[,2]
head(d1)
d2$POS <- str_split_fixed(d2$locus,":",2)[,2]
head(d2)


plot(main = "Beta scatter plot",d1$BETA,d2$BETA,xlab = "Hail",ylab = "epacts")
cor(d1$BETA,d2$BETA)
cor.test(d1$BETA,d2$BETA)

par(mfrow=c(1,2))
hist(main = "Beta histogram - hail",d1$BETA,ylim = c(0,2500),col = adjustcolor("Red",alpha = 0.3))
hist(main = "Beta histogram - epacts",d2$BETA,ylim = c(0,2500),col = adjustcolor("Blue",alpha = 0.3))

dev.off()
summary(d1$BETA)
sd(d1$BETA)

summary(d2$BETA)
sd(d2$BETA)

plot(log10(d1$PVALUE),log10(d2$PVALUE),col = adjustcolor("blue",alpha = 0.5),main = "log10(p-value)",xlab= "hail p-value",ylab = "epacts p-value")
