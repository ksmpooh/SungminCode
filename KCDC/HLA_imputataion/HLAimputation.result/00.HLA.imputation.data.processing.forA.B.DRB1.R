##### 20201026
##### HLA A B DRB1 using alltype processing NGS result
setwd("c:/Users/user/Desktop/KCDC/")

### 2digit
ref <- read.csv("HLAimputation/NGS/HLA_NGS_typing_255samples_results_202002_modify_alltype_processing.2digit.csv")
ref2 <- read.csv("transplantation/HLAtyping/20200828/HLAtyping.alle.gene.4digit.csv")
head(ref)
head(ref2)
colnames(ref) <- c('KID','NGS_A.1','NGS_A.2','NGS_B.1','NGS_B.2','NGS_DRB1.1','NGS_DRB1.2')
ref <- merge(ref2[,c(1:2)],ref,by='KID')
##cookHLA
cookHLA.han <- read.csv("HLAimputation/20200731/Han/HLAimputation.all.gene.2digit.Result.csv")
cookHLA.pan <- read.csv("HLAimputation/20200731/Pan/HLAimputation.all.gene.2digit.Result.csv")
impute4.han <- read.csv("HLAimputation/IMPUTE4/Han.ref/Result/HLAimputation.all.gene.2digit.Result.csv")
impute4.pan <- read.csv("HLAimputation/IMPUTE4/Pan.ref/Result/HLAimputation.all.gene.2digit.Result.csv")

head(cookHLA.han)

cookHLA.han<-merge(cookHLA.han[,c(1:5,8,9)],ref,by.x = "IID",by.y = "KID")[,c(1,8,2:7,9:14)]
cookHLA.pan<-merge(cookHLA.pan[,c(1:5,8,9)],ref,by.x = "IID",by.y = "KID")[,c(1,8,2:7,9:14)]
impute4.han<-merge(impute4.han[,c(1:5,8,9)],ref,by.x = "IID",by.y = "KID")[,c(1,8,2:7,9:14)]
impute4.pan<-merge(impute4.pan[,c(1:5,8,9)],ref,by.x = "IID",by.y = "KID")[,c(1,8,2:7,9:14)]

write.csv(cookHLA.han,"HLAimputation/20201026/cookHLA/han/MERGE.impResult.hlatyping.A_B_DRB1.gene.2digit.csv",row.names = F,quote = F)
write.csv(cookHLA.pan,"HLAimputation/20201026/cookHLA/pan/MERGE.impResult.hlatyping.A_B_DRB1.gene.2digit.csv",row.names = F,quote = F)
write.csv(impute4.han,"HLAimputation/20201026/impute4/han/MERGE.impResult.hlatyping.A_B_DRB1.gene.2digit.csv",row.names = F,quote = F)
write.csv(impute4.pan,"HLAimputation/20201026/impute4/pan/MERGE.impResult.hlatyping.A_B_DRB1.gene.2digit.csv",row.names = F,quote = F)


### 4digit

ref <- read.csv("HLAimputation/NGS/HLA_NGS_typing_255samples_results_202002_modify_alltype_processing.4digit.csv")
ref2 <- read.csv("transplantation/HLAtyping/20200828/HLAtyping.alle.gene.4digit.csv")
head(ref)
head(ref2)
colnames(ref) <- c('KID','NGS_A.1','NGS_A.2','NGS_B.1','NGS_B.2','NGS_DRB1.1','NGS_DRB1.2')
ref <- merge(ref2[,c(1:2)],ref,by='KID')
##cookHLA
cookHLA.han <- read.csv("HLAimputation/20200731/Han/HLAimputation.all.gene.4digit.Result.csv")
cookHLA.pan <- read.csv("HLAimputation/20200731/Pan/HLAimputation.all.gene.4digit.Result.csv")
impute4.han <- read.csv("HLAimputation/IMPUTE4/Han.ref/Result/HLAimputation.all.gene.4digit.Result.csv")
impute4.pan <- read.csv("HLAimputation/IMPUTE4/Pan.ref/Result/HLAimputation.all.gene.4digit.Result.csv")

cookHLA.han<-merge(cookHLA.han[,c(1:5,8,9)],ref,by.x = "IID",by.y = "KID")[,c(1,8,2:7,9:14)]
cookHLA.pan<-merge(cookHLA.pan[,c(1:5,8,9)],ref,by.x = "IID",by.y = "KID")[,c(1,8,2:7,9:14)]
impute4.han<-merge(impute4.han[,c(1:5,8,9)],ref,by.x = "IID",by.y = "KID")[,c(1,8,2:7,9:14)]
impute4.pan<-merge(impute4.pan[,c(1:5,8,9)],ref,by.x = "IID",by.y = "KID")[,c(1,8,2:7,9:14)]

write.csv(cookHLA.han,"HLAimputation/20201026/cookHLA/han/MERGE.impResult.hlatyping.A_B_DRB1.gene.4digit.csv",row.names = F,quote = F)
write.csv(cookHLA.pan,"HLAimputation/20201026/cookHLA/pan/MERGE.impResult.hlatyping.A_B_DRB1.gene.4digit.csv",row.names = F,quote = F)
write.csv(impute4.han,"HLAimputation/20201026/impute4/han/MERGE.impResult.hlatyping.A_B_DRB1.gene.4digit.csv",row.names = F,quote = F)
write.csv(impute4.pan,"HLAimputation/20201026/impute4/pan/MERGE.impResult.hlatyping.A_B_DRB1.gene.4digit.csv",row.names = F,quote = F)




### merge NGS which is modified result using alltype NGS
