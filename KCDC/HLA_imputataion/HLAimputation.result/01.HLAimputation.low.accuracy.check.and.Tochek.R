###20201021 
## HLA A-B-DRB1 결과낮은 sample 결과 확인

setwd("c:/Users/user/Desktop/KCDC/HLAimputation/")

df <- read.csv("HLA.accuracy.A.B.DRB1.sampleAccuracy.4digit.with.CDC015(255samples).csv",header = T)
head(df)
df <- read.csv("HLA.accuracy.A.B.DRB1.sampleAccuracy.2digit.with.CDC015(255samples).csv",header = T)



ref <- read.csv("../transplantation/HLAtyping/20200828/HLAtyping.alle.gene.4digit.csv",header = T)
ref1 <- read.csv("../transplantation/HLAtyping/HLAtyping.allGene.result.csv",header = T)
head(ref)
head(ref1)
a <- df[df$Han.cookHLA.4digit <= 0.6 | df$Pan.cookHLA.4digit <= 0.6 | df$Han.impute4.4digit <= 0.6 |df$Pan.impute4.4digit <= 0.6,]
a <- df[df$Han.cookHLA.4digit <= 0.7 | df$Pan.cookHLA.4digit <= 0.7 | df$Han.impute4.4digit <= 0.7 |df$Pan.impute4.4digit <= 0.7,]
a




b <- ref[ref$KID %in% a$IID ,]
b

c <- merge(a,b,by.x = "IID",by.y = "KID")
c


a <- ref1[ref1$X %in% b$YSample,]
a

v <- read.csv("NGS/HLA_NGS_typing_255samples_results_202002_modify_alltype_processing.4digit.csv")





##############20201026
setwd("c:/Users/user/Desktop/KCDC/")
df <- read.csv("HLAimputation/NGS/HLA_NGS_typing_255samples_results_202002_modify_alltype_processing.4digit.csv")
#df <- read.csv("HLAimputation/NGS/HLA_NGS_typing_255samples_results_202002_modify_alltype.csv")
head(df)

allgene <- read.csv("transplantation/HLAtyping/20200828/HLAtyping.alle.gene.4digit.csv")
head(allgene)
allgene <- allgene[,c(2,3,4,5,6,9,10)]
#allgene <- allgene[,c(1,3,4,5,6,9,10)]
colnames(df)<-colnames(allgene)


head(df)
head(allgene)

a <- merge(df[,c("KID","NGS_A.1","NGS_A.2")],allgene[,c("KID","NGS_A.1","NGS_A.2")],by = c("KID","NGS_A.1","NGS_A.2"))
b <- merge(df[,c("KID","NGS_B.1","NGS_B.2")],allgene[,c("KID","NGS_B.1","NGS_B.2")],by = c("KID","NGS_B.1","NGS_B.2"))
drb <- merge(df[,c("KID","NGS_DRB1.1","NGS_DRB1.2")],allgene[,c("KID","NGS_DRB1.1","NGS_DRB1.2")],by = c("KID","NGS_DRB1.1","NGS_DRB1.2"))

c <- !(df$KID %in% a$KID)





allgene <- read.csv("transplantation/HLAtyping/20200828/HLAtyping.alle.gene.4digit.csv")
allgene <- allgene[,c(1,2,3,4,5,6,9,10)]

head(df)
head(allgene)
colnames(df)[2:ncol(df)]<-c("before_A.1","before_A.2","before_B.1","before_B.2","before_DRB1.1","before_DRB1.2")
colnames(allgene)[3:ncol(allgene)]<-c("allgene_A.1","allgene_A.2","allgene_B.1","allgene_B.2","allgene_DRB1.1","allgene_DRB1.2") 

#allgene2 <- allgene
df2 <- df

colnames(df2)[2:ncol(df2)]<-c("before_A.2","before_A.1","before_B.2","before_B.1","before_DRB1.2","before_DRB1.1")
#colnames(allgene2)[3:ncol(allgene)]<-c("allgene_A.1","allgene_A.2","allgene_B.1","allgene_B.2","allgene_DRB1.1","allgene_DRB1.2") 


out_a <- merge(df[!(df$KID %in% a$KID),1:3],allgene[!(allgene$KID %in% a$KID),1:4],by = "KID")
out_a1 <- merge(df2[!(df2$KID %in% a$KID),1:3],allgene[!(allgene$KID %in% a$KID),1:4],by = "KID")

out_b <- merge(df[!(df$KID %in% b$KID),c(1,4,5)],allgene[!(allgene$KID %in% b$KID),c(1,2,5,6)],by = "KID")
out_drb1 <- merge(df[!(df$KID %in% drb$KID),c(1,6,7)],allgene[!(allgene$KID %in% drb$KID),c(1,2,7,8)],by = "KID")



out_a
out_a1

#out <- unique(out_a[duplicated(out_a,c(2,3)),c(3,2)])
out_a <- out_a[c(1,2,4,5,6,7,8,9),]

out_b
out_b <- out_b[c(2,5,6,7,10),]
 

out_drb1
out_drb1 <- out_drb1[c(1,3,4,5,6,9,10,12,13,14,15,16,17,18,19,20,25,26,27,28),]


head(ori_df)
ori_df <- read.csv("HLAimputation/NGS/HLA_NGS_typing_255samples_results_202002_modify_alltype.csv")
ori_allgene <- read.csv("transplantation/HLAtyping/HLAtyping.allGene.result.csv")
ori_allgene <- ori_allgene[,c("X","A_1","A_2","B_1","B_2","DRB1_1","DRB1_2")]

head(ori_df)
head(ori_allgene)

colnames(ori_df)[3:ncol(ori_df)] <-c("modify.A_1","modify.A_2","modify.B_1","modify.B_2","modify.DRB1_1","modify.DRB1_2")
colnames(ori_allgene)<- c("Sample","allgene.A_1","allgene.A_2","allgene.B_1","allgene.B_2","allgene.DRB1_1","allgene.DRB1_2")



ori_a <-merge(ori_df[,c(1,2,3,4)],ori_allgene[,c(1,2,3)],by="Sample")
ori_a <- ori_a[ori_a$KID %in% out_a$KID,]

ori_b <-merge(ori_df[,c(1,2,5,6)],ori_allgene[,c(1,4,5)],by="Sample")
ori_b <- ori_b[ori_b$KID %in% out_b$KID,]

ori_drb1 <-merge(ori_df[,c(1,2,7,8)],ori_allgene[,c(1,6,7)],by="Sample")
ori_drb1 <- ori_drb1[ori_drb1$KID %in% out_drb1$KID,]


#ori_ori <- read.csv("transplantation/HLAtyping/HLA_NGS_typing_255samples_results_202002.csv")
#colnames(ori_ori)<- c("Sample","KID","original.A_1","original.A_2","original.B_1","original.B_2","allgene.DRB1_1","allgene.DRB1_2")
#head(ori_ori)

write.csv(ori_a,"transplantation/HLAtyping/20201026_checklist/HLAtyping.A.gene.tocheck.csv",row.names = F)
write.csv(ori_b,"transplantation/HLAtyping/20201026_checklist/HLAtyping.B.gene.tocheck.csv",row.names = F)
write.csv(ori_drb1,"transplantation/HLAtyping/20201026_checklist/HLAtyping.DRB1.gene.tocheck.csv",row.names = F)


han <- read.csv("HLAimputation/IMPUTE4/Han.ref/Result/HLAimputation.all.gene.4digit.Result.csv")
pan <- read.csv("HLAimputation/IMPUTE4/Pan.ref/Result/HLAimputation.all.gene.4digit.Result.csv")

han <- han[,1:3]
pan <- pan[,1:3]
colnames(han)<-c("KID","Han_A.1","Han_A.2")
colnames(pan)<-c("KID","Pan_A.1","Pan_A.2")

out <- ori_a
out <- merge(out,han,by="KID",all.x = T)
out <- merge(out,pan,by="KID",all.x = T)


out <- out[,c(4,1,2,3,5,6,7,8,9,10)]

write.csv(out,"transplantation/HLAtyping/20201026_checklist/HLAtyping.A.gene.tocheck_with.IMPUTE.Result.csv",row.names = F)


head(han)
han <- read.csv("HLAimputation/IMPUTE4/Han.ref/Result/HLAimputation.all.gene.4digit.Result.csv")
pan <- read.csv("HLAimputation/IMPUTE4/Pan.ref/Result/HLAimputation.all.gene.4digit.Result.csv")

han <- han[,c(1,4,5)]
pan <- pan[,c(1,4,5)]
colnames(han)<-c("KID","Han_B.1","Han_B.2")
colnames(pan)<-c("KID","Pan_B.1","Pan_B.2")

out <- ori_b
out <- merge(out,han,by="KID",all.x = T)
out <- merge(out,pan,by="KID",all.x = T)


out <- out[,c(4,1,2,3,5,6,7,8,9,10)]
write.csv(out,"transplantation/HLAtyping/20201026_checklist/HLAtyping.B.gene.tocheck_with.IMPUTE.Result.csv",row.names = F)




han <- read.csv("HLAimputation/IMPUTE4/Han.ref/Result/HLAimputation.all.gene.4digit.Result.csv")
pan <- read.csv("HLAimputation/IMPUTE4/Pan.ref/Result/HLAimputation.all.gene.4digit.Result.csv")
head(han)
han <- han[,c(1,8,9)]
pan <- pan[,c(1,8,9)]
colnames(han)<-c("KID","Han_DRB1.1","Han_DRB1.2")
colnames(pan)<-c("KID","Pan_DRB1.1","Pan_DRB1.2")

out <- ori_drb1
out <- merge(out,han,by="KID",all.x = T)
out <- merge(out,pan,by="KID",all.x = T)

head(out)
out <- out[,c(4,1,2,3,5,6,7,8,9,10)]

write.csv(out,"transplantation/HLAtyping/20201026_checklist/HLAtyping.DRB1.gene.tocheck_with.IMPUTE.Result.csv",row.names = F)


##### check list 20201027
setwd("c:/Users/user/Desktop/KCDC/HLAimputation/")

impute4.han <- read.csv("20201026/impute4/han/compare.IMPvsNGS.A.B.DRB1.4digit.csv")
impute4.pan <- read.csv("20201026/impute4/pan/compare.IMPvsNGS.A.B.DRB1.4digit.csv")
cookHLA.han <- read.csv("20201026/cookHLA/han/compare.IMPvsNGS.A.B.DRB1.4digit.csv")
cookHLA.pan <- read.csv("20201026/cookHLA/pan/compare.IMPvsNGS.A.B.DRB1.4digit.csv")

a <- read.csv("../transplantation/HLAtyping/20201026_checklist/HLAtyping.A.gene.tocheck.csv")
b <- read.csv("../transplantation/HLAtyping/20201026_checklist/HLAtyping.B.gene.tocheck.csv")
drb <- read.csv("../transplantation/HLAtyping/20201026_checklist/HLAtyping.DRB1.gene.tocheck.csv")


df <- merge(impute4.han[,c(1,3,4)],cookHLA.han[,c(1,3,4)],by = "IID")
df <- merge(df,impute4.pan[,c(1,3,4)],by = "IID")
df <- merge(df,cookHLA.pan[,c(1,3,4)],by = "IID")
head(df)
colnames(df) <- c("KID","impute4.han.A_1","impute4.han.A_2","cookHLA.han.A_1","cookHLA.han.A_2","impute4.pan.A_1","impute4.pan.A_2","cookHLA.pan.A_1","cookHLA.pan.A_2")

a <- merge(a,df,by = "KID",all.x = T)


df <- merge(impute4.han[,c(1,10,11)],cookHLA.han[,c(1,10,11)],by = "IID")
df <- merge(df,impute4.pan[,c(1,10,11)],by = "IID")
df <- merge(df,cookHLA.pan[,c(1,10,11)],by = "IID")
head(df)
colnames(df) <- c("KID","impute4.han.A_1","impute4.han.A_2","cookHLA.han.A_1","cookHLA.han.A_2","impute4.pan.A_1","impute4.pan.A_2","cookHLA.pan.A_1","cookHLA.pan.A_2")

b <- merge(b,df,by = "KID",all.x = T)



df <- merge(impute4.han[,c(1,17,18)],cookHLA.han[,c(1,17,18)],by = "IID")
df <- merge(df,impute4.pan[,c(1,17,18)],by = "IID")
df <- merge(df,cookHLA.pan[,c(1,17,18)],by = "IID")
head(df)
colnames(df) <- c("KID","impute4.han.A_1","impute4.han.A_2","cookHLA.han.A_1","cookHLA.han.A_2","impute4.pan.A_1","impute4.pan.A_2","cookHLA.pan.A_1","cookHLA.pan.A_2")

drb <- merge(drb,df,by = "KID",all.x = T)


write.csv(a,"../transplantation/HLAtyping/20201026_checklist/HLAtyping.A.gene.tocheck_with.IMPUTE.Result.csv",row.names = F,quote = F)
write.csv(b,"../transplantation/HLAtyping/20201026_checklist/HLAtyping.B.gene.tocheck_with.IMPUTE.Result.csv",row.names = F,quote = F)
write.csv(drb,"../transplantation/HLAtyping/20201026_checklist/HLAtyping.DRB1.gene.tocheck_with.IMPUTE.Result.csv",row.names = F,quote = F)




a
b
drb


df <- read.csv("HLAimputation/20201026/HLA.accuracy.A.B.DRB1.sampleAccuracy.4digit.csv")
df1 <- read.csv("HLAimputation/HLA.accuracy.A.B.DRB1.sampleAccuracy.4digit.csv")
df <- df[df$IID %in% a$KID,]
df1 <- df1[df1$IID %in% a$KID,]

df <- df[df$IID %in% b$KID,]
df1 <- df1[df1$IID %in% b$KID,]

head(df1)
head(df)
out <- merge(df,df1,by = "IID")
out
