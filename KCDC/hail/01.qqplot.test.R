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
d$p.check <- as.factor(d$p.check)
#d <- unique.data.frame(d)
library(ggplot2)
ggplot(d,aes(x=hail.BETA,y=epacts.BETA,fill = p.check,color = ifelse(p.check==1,"red","darkgray"))) +
#ggplot(d,aes(x=hail.BETA,y=epacs.BETA,fill = p.check)) +
  ggtitle("association Result : hail.beta vs epacts.beta") +
  geom_point() + 
  theme(plot.title = element_text(hjust = 0.5),legend.title = element_text()) +
  #guides(fill = guide_legend(title = "epacts P.value")) +
  scale_fill_discrete(name = "epacts P.value",label = c("<0.05",">=0.05")) + 
#  scale_fill_manual(values = c("#FF0000","#a9a9a9")) 
#  scale_fill_gradientn(colours = c("red","darkgrey")) +
  scale_color_manual(guide=FALSE, values=c("red","darkgray"))

cor.test(d$hail.BETA,d$epacts.BETA)

ggplot(d,aes(x=-log10(hail.PVALUE),y=-log10(epacts.PVALUE),color = ifelse(p.check==1,"red","darkgray"))) +
  ggtitle( "Association Result Compare : hail.p-value vs epacts.p-value") +
  geom_point() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(guide=FALSE, values=c("red","darkgray"))




ggplot(d,aes(x=hail.BETA,y=epacts.BETA,fill = p.check,color = ifelse(p.check=="1","red","darkgray"))) +
#  scale_color_manual(guide=FALSE, values=c("darkgray","red")) + 
  scale_color_manual(guide=FALSE, values=c("darkgray","red")) + 
  #ggplot(d,aes(x=hail.BETA,y=epacs.BETA,fill = p.check)) +
  ggtitle("association Result : hail.beta vs epacts.beta") +
  geom_point()
  theme(plot.title = element_text(hjust = 0.5),legend.title = element_text()) +
  labs(fill = "epacts.P.value") + 
  scale_fill_discrete(labels = c(">=0.05","<0.05")) 
#  scale_colour_manual(values = c("#FF0000","#a9a9a9")) 
#  scale_colour_manual(values = c("darkgrey","red")) 

# guides(fill = guide_legend(title = "epacts P.value"))
  scale_fill_discrete(name = "epacts P.value",label = c("<0.05",">=0.05"))
  scale_color_manual(values = c("#FF0000","#a9a9a9")) 
  #  scale_fill_gradientn(colours = c("red","darkgrey")) +
 scale_color_manual(guide=FALSE, values=c("red","darkgray"))

 
 
 
 
summary(d$p.check)
plot(d$hail.BETA,d$epacts.BETA,main = "association Result : hail.beta vs epacts.beta"
     ,xlab = "hail.Beta",ylab = "epact.Beta"
     ,col = ifelse(d$p.check==0,rgb(0.5,0.5,0.5,0.2),rgb(1,0,0,0.5))
     #,col = rgb(0.6,0.6,0.6,0.3)
     ,pch =16)
legend("bottomright",title = "epacts p.value",legend = c("<0.05",">=0.05"),col = c("red","darkgray"),pch = 16,box.lty = 0)


hist(abs(d$epacts.BETA - d$hail.BETA))

plot(-log10(d$hail.PVALUE),-log10(d$epacts.PVALUE),main = "association Result : hail.p value vs epacts.p value"
     ,xlab = "-log10(hail p-value)",ylab = "-log10(epacts p-value)"
     ,col = ifelse(d$p.check==0,rgb(0.5,0.5,0.5,0.5),rgb(1,0,0,0.5)))
legend("bottomright",title = "epacts p.value",legend = c("<0.05",">=0.05"),col = c("red","darkgray"),pch = 16,box.lty = 0)
