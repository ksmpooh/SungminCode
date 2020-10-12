setwd("c:/Users/user/Desktop/KCDC/HLAimputation/")


df <- read.csv("IMPUTE4/Han.ref/Result/compare.IMPvsNGS.all.gene.4digit.csv",header = T)
#df <- read.csv("IMPUTE4/Han.ref/Result/compare.IMPvsNGS.all.gene.2digit.csv",header = T)

df <-df[df$YSample != 'CDC015',]

df$match <- df$A.match + df$B.match + df$DRB1.match
df$wrong <- df$A.wrong + df$B.wrong + df$DRB1.wrong
df$accuracy <- df$match/(df$match+df$wrong)

han.impute4.4digit <- df
out3 <- han.impute4.4digit[,c("IID","accuracy")]


df <- read.csv("IMPUTE4/Pan.ref/Result/compare.IMPvsNGS.all.gene.4digit.csv")
#df <- read.csv("IMPUTE4/Pan.ref/Result/compare.IMPvsNGS.all.gene.2digit.csv",header = T)

df <-df[df$YSample != 'CDC015',]
df$match <- df$A.match + df$B.match + df$DRB1.match
df$wrong <- df$A.wrong + df$B.wrong + df$DRB1.wrong
df$accuracy <- df$match/(df$match+df$wrong)
pan.impute4.4digit <- df
out4 <- pan.impute4.4digit[,c("IID","accuracy")]

df <- read.csv("20200731/Han/compare.IMPvsNGS.all.gene.4digit.csv",header = T)
#df <- read.csv("20200731/Han/compare.IMPvsNGS.all.gene.2digit.csv",header = T)

df <-df[df$YSample != 'CDC015',]

df$match <- df$A.match + df$B.match + df$DRB1.match
df$wrong <- df$A.wrong + df$B.wrong + df$DRB1.wrong
df$accuracy <- df$match/(df$match+df$wrong)
han.cookHLA.4digit <- df
out1 <- han.cookHLA.4digit[,c("IID","accuracy")]

df <- read.csv("20200731/Pan/compare.IMPvsNGS.all.gene.4digit.csv",header = T)
#df <- read.csv("20200731/Pan/compare.IMPvsNGS.all.gene.2digit.csv",header = T)
df <-df[df$YSample != 'CDC015',]

df$match <- df$A.match + df$B.match + df$DRB1.match
df$wrong <- df$A.wrong + df$B.wrong + df$DRB1.wrong
df$accuracy <- df$match/(df$match+df$wrong)

pan.cookHLA.4digit <- df
out2 <- pan.cookHLA.4digit[,c("IID","accuracy")]



boxplot(df$accuracy)

sum(table(out$Pan.impute4.2digit))

out <- merge(out1,out2,by="IID")
#colnames(out)[2:3]<-c()
head(out)
out <- merge(out,out3,by="IID")
out <- merge(out,out4,by="IID")
head(out)
colnames(out)[2:5]<-c("Han.cookHLA.4digit","Pan.cookHLA.4digit","Han.impute4.4digit","Pan.impute4.4digit")
#colnames(out)[2:5]<-c("Han.cookHLA.2digit","Pan.cookHLA.2digit","Han.impute4.2digit","Pan.impute4.2digit")
#colnames(out)[2:5]<-c("Han.cookHLA","Pan.cookHLA","Han.impute4","Pan.impute4")

str(out)
#write.csv(out,"HLA.accuracy.A.B.DRB1.sampleAccuracy.4digit.without.cdc015(254sample).csv",col.names = T,row.names = F,quote = F)
#write.csv(out,"HLA.accuracy.A.B.DRB1.sampleAccuracy.2digit.csv",col.names = T,row.names = F,quote = F)


########

out <- read.csv("HLA.accuracy.A.B.DRB1.sampleAccuracy.4digit.csv",header = T)
#out <- read.csv("HLA.accuracy.A.B.DRB1.sampleAccuracy.2digit.csv",header = T)
summary(out)
table(out$Han.impute4)
colnames(out)[2:5]<-c("Han.cookHLA","Pan.cookHLA","Han.impute4","Pan.impute4")
mean(out$Han.impute4)

boxplot(out[2:5],main = "HLA imputation Sample Accuracy : 4digit"
        ,xlab = "Method",ylab = "Accuracy",col = c("gold","darkgreen")
        )
dev.off()



par(mfrow=c(2,2))
hist(out$Han.cookHLA,main ="Han.cookHLA",xlab = "Accuracy")
hist(out$Pan.cookHLA,main ="Pan.cookHLA",xlab = "Accuracy")
hist(out$Han.impute4,main ="Han.impute4",xlab = "Accuracy")
hist(out$Pan.impute4,main ="Pan.impute4",xlab = "Accuracy")
####0.6 이하 및 HLA 용역 ID 추가

#out <- read.csv("HLA.accuracy.A.B.DRB1.sampleAccuracy.csv")
out <- read.csv("HLA.accuracy.A.B.DRB1.sampleAccuracy.4digit.csv")
ref <- read.csv("NGS/HLA_NGS_typing_255samples_results_202002.csv")
rownames(ref) <-ref$KID
ref <- ref[,c(1,2)]
head(ref)
head(out)
a <- out[out$Han.cookHLA.4digit <= 0.6 | out$Pan.cookHLA.4digit <= 0.6 | out$Han.impute4.4digit <= 0.6 |out$Pan.impute4.4digit <= 0.6,]
#a <- out[out$han.cookHLA <= 0.6 | out$pan.cookHLA <= 0.6 | out$han.impute4 <= 0.6 | out$pan.impute4 <= 0.6,]
a

a <- merge(a,ref,by.x = "IID",by.y = "KID",all.x = T)
a[,c(5,1,2,3,4)]
write.csv(a[,c(5,1,2,3,4)],"HLA.accuracy.A.B.DRB1.sampleAccuracy.0.6.csv",row.names = F,quote = F)



### A B DRB1 without CDC015

setwd("c:/Users/user/Desktop/KCDC/HLAimputation/")


#df <- read.csv("IMPUTE4/Han.ref/Result/compare.IMPvsNGS.all.gene.4digit.csv",header = T)
df <- read.csv("IMPUTE4/Han.ref/Result/compare.IMPvsNGS.all.gene.2digit.csv",header = T)
#df <-df[df$YSample != 'CDC015',]

df$match <- df$A.match + df$B.match + df$DRB1.match
df$wrong <- df$A.wrong + df$B.wrong + df$DRB1.wrong
df$accuracy <- df$match/(df$match+df$wrong)

han.impute4.4digit <- df
out3 <- han.impute4.4digit[,c("IID","accuracy")]


#df <- read.csv("IMPUTE4/Pan.ref/Result/compare.IMPvsNGS.all.gene.4digit.csv")
df <- read.csv("IMPUTE4/Pan.ref/Result/compare.IMPvsNGS.all.gene.2digit.csv",header = T)


#df <-df[df$YSample != 'CDC015',]
df$match <- df$A.match + df$B.match + df$DRB1.match
df$wrong <- df$A.wrong + df$B.wrong + df$DRB1.wrong
df$accuracy <- df$match/(df$match+df$wrong)
pan.impute4.4digit <- df
out4 <- pan.impute4.4digit[,c("IID","accuracy")]

#df <- read.csv("20200731/Han/compare.IMPvsNGS.all.gene.4digit.csv",header = T)
df <- read.csv("20200731/Han/compare.IMPvsNGS.all.gene.2digit.csv",header = T)

#df <-df[df$YSample != 'CDC015',]

df$match <- df$A.match + df$B.match + df$DRB1.match
df$wrong <- df$A.wrong + df$B.wrong + df$DRB1.wrong
df$accuracy <- df$match/(df$match+df$wrong)
han.cookHLA.4digit <- df
out1 <- han.cookHLA.4digit[,c("IID","accuracy")]

#df <- read.csv("20200731/Pan/compare.IMPvsNGS.all.gene.4digit.csv",header = T)
df <- read.csv("20200731/Pan/compare.IMPvsNGS.all.gene.2digit.csv",header = T)
#df <-df[df$YSample != 'CDC015',]

df$match <- df$A.match + df$B.match + df$DRB1.match
df$wrong <- df$A.wrong + df$B.wrong + df$DRB1.wrong
df$accuracy <- df$match/(df$match+df$wrong)

pan.cookHLA.4digit <- df
out2 <- pan.cookHLA.4digit[,c("IID","accuracy")]

out <- merge(out1,out2,by="IID")
#colnames(out)[2:3]<-c()
head(out)
out <- merge(out,out3,by="IID")
out <- merge(out,out4,by="IID")
head(out)
#colnames(out)[2:5]<-c("Han.cookHLA.4digit","Pan.cookHLA.4digit","Han.impute4.4digit","Pan.impute4.4digit")
colnames(out)[2:5]<-c("Han.cookHLA.2digit","Pan.cookHLA.2digit","Han.impute4.2digit","Pan.impute4.2digit")

write.csv(out,"HLA.accuracy.A.B.DRB1.sampleAccuracy.4digit.with.CDC015(255samples).csv",col.names = T,row.names = F,quote = F)
write.csv(out,"HLA.accuracy.A.B.DRB1.sampleAccuracy.2digit.with.CDC015(255samples).csv",col.names = T,row.names = F,quote = F)


fd <- out
td <- out
summary(fd)
summary(td)


############################################################
setwd("c:/Users/user/Desktop/KCDC/HLAimputation/")
df <- read.csv("cookHLAvsIMPUTE4.compare.Result.csv",header = T)
df <- df[c(1:4,11:14),]
head(df)
df$type
df$type2<-c("Pan.impute4","Han.impute4","Pan.cookHLA","Han.cookHLA","Pan.impute4","Han.impute4","Pan.cookHLA","Han.cookHLA")
library(stringr)

