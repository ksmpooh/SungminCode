setwd("c:/Users/user/Desktop/KCDC/HLAimputation/")


df <- read.csv("IMPUTE4/Han.ref/Result/compare.IMPvsNGS.all.gene.4digit.csv",header = T)
df <- read.csv("20200731/Han/compare.IMPvsNGS.all.gene.4digit.csv",header = T)
df <- read.csv("20200731/Pan/compare.IMPvsNGS.all.gene.4digit.csv",header = T)

df <- read.csv("IMPUTE4/Han.ref/Result/compare.IMPvsNGS.all.gene.2digit.csv",header = T)
df <- read.csv("20200731/Han/compare.IMPvsNGS.all.gene.2digit.csv",header = T)
df <- read.csv("20200731/Pan/compare.IMPvsNGS.all.gene.2digit.csv",header = T)


colnames(df)
df <-df[df$YSample != 'CDC015',]

df$match <- df$A.match + df$B.match + df$DRB1.match
df$wrong <- df$A.wrong + df$B.wrong + df$DRB1.wrong
df$accuracy <- df$match/(df$match+df$wrong)

boxplot(df$accuracy)

han.impute4.4digit <- df
han.cookHLA.4digit <- df
pan.cookHLA.4digit <- df

out1 <- han.cookHLA.4digit[,c("IID","accuracy")]
out2 <- pan.cookHLA.4digit[,c("IID","accuracy")]
out3 <- han.impute4.4digit[,c("IID","accuracy")]

out <- merge(out1,out2,by="IID")
head(out)
out <- merge(out,out3,by="IID")
head(out)
colnames(out)[2:4]<-c("han.cookHLA.4digit","pan.cookHLA.4digit","han.impute4.4digit")

#colnames(out)[2:4]<-c("han.cookHLA.2digit","pan.cookHLA.2digit","han.impute4.2digit")

#write.csv(out,"HLA.accuracy.A.B.DRB1.sampleAccuracy.csv",col.names = T,row.names = F,quote = F)
write.csv(out,"HLA.accuracy.A.B.DRB1.sampleAccuracy.2digit.csv",col.names = T,row.names = F,quote = F)

par(mfrow=c(1,3))
boxplot(out$han.cookHLA.4digit,main = "han.cookHLA.4digit",ylim = c(0.3,1))
boxplot(out$pan.cookHLA.4digit,main = "pan.cookHLA.4digit",ylim = c(0.3,1))
boxplot(out$han.impute4.4digit,main = "han.impute4.4digit",ylim = c(0.3,1))
dev.off()


par(mfrow=c(1,3))
boxplot(out$han.cookHLA.2digit,main = "han.cookHLA.2digit")
boxplot(out$pan.cookHLA.2digit,main = "pan.cookHLA.2digit")
boxplot(out$han.impute4.2digit,main = "han.impute4.2digit")
dev.off()


####0.6 이하 및 HLA 용역 ID 추가

out <- read.csv("HLA.accuracy.A.B.DRB1.sampleAccuracy.csv")
ref <- read.csv("NGS/HLA_NGS_typing_255samples_results_202002.csv")
rownames(ref) <-ref$KID
ref <- ref[,c(1,2)]
head(ref)
head(out)
a <- out[out$han.cookHLA.4digit <= 0.6 | out$pan.cookHLA.4digit <= 0.6 | out$han.impute4.4digit <= 0.6,]
a

a <- merge(a,ref,by.x = "IID",by.y = "KID",all.x = T)
a[,c(5,1,2,3,4)]
write.csv(a[,c(5,1,2,3,4)],"HLA.accuracy.A.B.DRB1.sampleAccuracy.0.6.csv",row.names = F,quote = F)

