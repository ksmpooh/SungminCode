library(stringr)
setwd("c:/Users/user/Desktop/KCDC/hail/test/")
d1 <- read.table("hail/gaws.output.T2D.tsv",header = T)
head(d1)
d1 <- d1[,c(1,3,4,5)]
d2 <- read.table("epacts/test.chr22.merge.list.txt",header = T)
head(d2)
d2 <- d2[,c(9,10,12)]
head(d1)
colnames(d1) <- c("locus","hail.BETA","hail.CHISQ","hail.PVALUE")
head(d2)
colnames(d2) <- c("epacts.PVALUE","epacts.BETA","epacts.CHISQ")


d <- cbind(d1,d2)
head(d)
d$p.check = 0
d[d$epacts.PVALUE <0.05,]$p.check = 1
#d <- unique.data.frame(d)
library(ggplot2)
ggplot(main = "association Result : hail.beta vs epacts.beta",d,aes(x=hail.BETA,y=epacts.BETA,color = ifelse(p.check==1,"red","darkgray"))) +
  geom_point() +
  scale_color_manual(guide=FALSE, values=c("red","darkgray"))


ggplot(d,aes(x=-log10(hail.PVALUE),y=-log10(epacts.PVALUE),color = ifelse(p.check==1,"red","darkgray"))) +
  ggtitle( "Association Result Compare : hail.p-value vs epacts.p-value") +
  geom_point() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(guide=FALSE, values=c("red","darkgray"))
