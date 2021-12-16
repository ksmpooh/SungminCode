#setwd("~/Desktop/KCDC/HLAimputation/IMPUTE4/gen.calling.test/RESULTs/test/")
# recipient donor
#df <- read.table("~/Desktop/KCDC/HLAimputation/00.KOTRY_HLAimp/JG.ESRD.KR_2019/KR/JG.ESRD.KR.2019.HLAimputation_2digit_raw.raw",header = T)
#df <- read.table("~/Desktop/KCDC/HLAimputation/00.KOTRY_HLAimp/JG.ESRD.KR_2019/KR/JG.ESRD.KR.2019.HLAimputation_4digit_raw.raw",header = T)

df <- read.table("~/Desktop/KCDC/HLAimputation/00.KOTRY_HLAimp/KD_2019_2020/plink/JG.KD_with_rep.2digit_afterHardcall_raw.raw",header = T)
#df <- read.table("~/Desktop/KCDC/HLAimputation/00.KOTRY_HLAimp/KD_2019_2020/plink/JG.KD_with_rep.4digit_afterHardcall_raw.raw",header = T)

head(df)
colnames(df)
ncol(df)
n = nchar(colnames(df)[ncol(df)])
check <- substr(colnames(df)[ncol(df)],n-5,n)
if (check == ".A.P_P") {
  ori_col <- colnames(df) %>% str_replace_all("X6.","") %>% str_replace_all(".A.P","")
}else{
  ori_col <- colnames(df) %>% str_replace_all("X6.","") %>% str_replace_all(".P.A","")
  
}
ori_col


#'X6.HLA_DPB1_10.A.P_P'
#'X6.HLA_DQB1_06.P.A_P'


colnames(df) <- ori_col
df[is.na(df)]<- -1
colnames(df)
colnames(header)



library(stringr)
head(A)
A <- df[,c(1,2,grep("HLA_A_*",colnames(df)))]
B <- df[,c(1,2,grep("HLA_B_*",colnames(df)))]
C <- df[,c(1,2,grep("HLA_C_*",colnames(df)))]
DRB1 <- df[,c(1,2,grep("*DRB1*",colnames(df)))]
#DRB3 <- df[,c(1,2,grep("*DRB3*",colnames(df)))]

DPB1 <- df[,c(1,2,grep("*DPB1*",colnames(df)))]
DPA1 <- df[,c(1,2,grep("*DPA1*",colnames(df)))]

DQB1 <- df[,c(1,2,grep("*DQB1*",colnames(df)))]
DQA1 <- df[,c(1,2,grep("*DQA1*",colnames(df)))]




#DRB <- hla.subset(DRB,3)
hla.find<-function(df,concept){
  colcount <- ncol(df)
  print(colcount)
  rowcount <- nrow(df)
  print(rowcount)
  for(i in (1:rowcount)){
    n = 0
    for(j in (3:colcount)){
      if(df[i,j] == 2){
        if(n == 0){
          df[i,paste("IMP_",concept,'.1',sep = "")] <- as.integer(str_split_fixed(colnames(df)[j],"_",4)[3])
          df[i,paste("IMP_",concept,'.2',sep = "")] <- as.integer(str_split_fixed(colnames(df)[j],"_",4)[3])
          n = n + 2  
        }else if(n == 1){
          df[i,paste("IMP_",concept,'.2',sep = "")] <- as.integer(str_split_fixed(colnames(df)[j],"_",4)[3])
          df[i,paste("IMP_",concept,'.3',sep = "")] <- as.integer(str_split_fixed(colnames(df)[j],"_",4)[3])
          n = n + 2  
        }
      }else if(df[i,j] == 1){
        if(n == 0){
          df[i,paste("IMP_",concept,'.1',sep = "")] <- as.integer(str_split_fixed(colnames(df)[j],"_",4)[3])
          n = n + 1
        }else if(n == 1){
          df[i,paste("IMP_",concept,'.2',sep = "")] <- as.integer(str_split_fixed(colnames(df)[j],"_",4)[3])
          n = n + 1
        }else{
          df[i,paste("IMP_",concept,'.3',sep = "")]  <- as.integer(str_split_fixed(colnames(df)[j],"_",4)[3])
        }
      }
    }
  }
  print(head(df))
  return(df)
}


hla.find<-function(df,concept){
  colcount <- ncol(df)
  print(colcount)
  rowcount <- nrow(df)
  print(rowcount)
  for(i in (1:rowcount)){
    n = 0
    for(j in (3:colcount)){
      if(df[i,j] == 2){
        if(n == 0){
          df[i,paste("IMP_",concept,'.1',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          df[i,paste("IMP_",concept,'.2',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          n = n + 2  
        }else if(n == 1){
          df[i,paste("IMP_",concept,'.2',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          df[i,paste("IMP_",concept,'.3',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          n = n + 2  
        }
      }else if(df[i,j] == 1){
        if(n == 0){
          df[i,paste("IMP_",concept,'.1',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          n = n + 1
        }else if(n == 1){
          df[i,paste("IMP_",concept,'.2',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          n = n + 1
        }else{
          df[i,paste("IMP_",concept,'.3',sep = "")]  <- str_split_fixed(colnames(df)[j],"_",4)[3]
        }
      }
    }
  }
  print(head(df))
  return(df)
}



DRB1 <- hla.find(DRB1,"DRB1")
table(DRB1$DRB3)


DQA1 <- hla.find(DQA1,"DQA1")
DQB1 <- hla.find(DQB1,"DQB1")


DPA1 <- hla.find(DPA1,"DPA1")
DPB1 <- hla.find(DPB1,"DPB1")


#grep("*N_P",colnames(A))
#colnames(A)[6]<-"HLA_A_0122_P"
#colnames(A)[19]<-"HLA_A_0253_P"

A <- hla.find(A,"A")

B <- hla.find(B,"B")

C <- hla.find(C,"C")

head(DQA1)
table(DQA1$IMP_DQA1.3)
table(DQB1$IMP_DQB1.3)
table(DPA1$IMP_DPA1.3)
table(DPB1$IMP_DPB1.3)
table(DRB1$IMP_DRB1.3)
table(C$IMP_C.3)
table(B$IMP_B.3)
table(A$IMP_A.3)
head(A)
head(B)


out <- merge(A[,c('IID','IMP_A.1','IMP_A.2')],B[,c('IID','IMP_B.1','IMP_B.2')],by = 'IID')
out <- merge(out,C[,c('IID','IMP_C.1','IMP_C.2')],by = 'IID')
out <- merge(out,DRB1[,c('IID','IMP_DRB1.1','IMP_DRB1.2')],by = 'IID')
out <- merge(out,DPA1[,c('IID','IMP_DPA1.1','IMP_DPA1.2')],by = 'IID')
out <- merge(out,DPB1[,c('IID','IMP_DPB1.1','IMP_DPB1.2')],by = 'IID')
out <- merge(out,DQA1[,c('IID','IMP_DQA1.1','IMP_DQA1.2')],by = 'IID')
out <- merge(out,DQB1[,c('IID','IMP_DQB1.1','IMP_DQB1.2')],by = 'IID')


concept = "KR"
concept = "KD"

new_col <- c("KBA_ID")
for (i in c("A","B","C","DRB1","DPA1","DPB1","DQA1","DQB1")) {
  new_col <- c(new_col,paste0(concept,"_",i,".1"),paste0(concept,"_",i,"_B.1"))
}
new_col
colnames(out) <-new_col


KR <- out
KD <- out
#####merge NGS
Ref <- read_excel("~/Desktop/KCDC/HLAimputation/00.KOTRY_HLAimp/ALL(2019and2020).sampleID.xlsx") %>% filter(Prod == 2019) %>% 
  select(KBA_ID,OriID,ref)
head(Ref)


out <- merge(Ref%>%merge(KR,by="KBA_ID") %>% rename(KR_KBA_ID = 'KBA_ID',KR_OriID = 'OriID'),
             Ref %>% merge(KD,by="KBA_ID") %>% rename(KD_KBA_ID = 'KBA_ID',KD_OriID = 'OriID'),by='ref') 
colnames(out)
organ = "K"
paste0(organ,"R_A.1")
colnames(out)[4:19] <- c(paste0(organ,"R_A.1"),paste0(organ,"R_A.2"),paste0(organ,"R_B.1"),paste0(organ,"R_B.2"),paste0(organ,"R_C.1"),paste0(organ,"R_C.2"),paste0(organ,"R_DRB1.1"),paste0(organ,"R_DRB1.2"),paste0(organ,"R_DPA1.1"),paste0(organ,"R_DPA1.2"),paste0(organ,"R_DPB1.1"),paste0(organ,"R_DPB1.2"),paste0(organ,"R_DQA1.1"),paste0(organ,"R_DQA1.2"),paste0(organ,"R_DQB1.1"),paste0(organ,"R_DQB1.2"))
colnames(out)[22:37] <- c(paste0(organ,"D_A.1"),paste0(organ,"D_A.2"),paste0(organ,"D_B.1"),paste0(organ,"D_B.2"),paste0(organ,"D_C.1"),paste0(organ,"D_C.2"),paste0(organ,"D_DRB1.1"),paste0(organ,"D_DRB1.2"),paste0(organ,"D_DPA1.1"),paste0(organ,"D_DPA1.2"),paste0(organ,"D_DPB1.1"),paste0(organ,"D_DPB1.2"),paste0(organ,"D_DQA1.1"),paste0(organ,"D_DQA1.2"),paste0(organ,"D_DQB1.1"),paste0(organ,"D_DQB1.2"))
out <- out[,c(1,2,3,20,21,4:19,22:37)]
head(out)
colnames(out)
#write.csv(out,"~/Desktop/KCDC/HLAimputation/03.HLA_matching/2019_Kidney_HLAmathing_table_2digit.csv",row.names = F,quote = F)
#write.csv(out,"~/Desktop/KCDC/HLAimputation/03.HLA_matching/2019_Kidney_HLAmathing_table_4digit.csv",row.names = F,quote = F)


write.table(out,"~/Desktop/KCDC/HLAimputation/03.HLA_matching/2019_Kidney_HLAmathing_table_2digit.txt",row.names = F,quote = F,sep = "\t")
write.table(out,"~/Desktop/KCDC/HLAimputation/03.HLA_matching/2019_Kidney_HLAmathing_table_4digit.txt",row.names = F,quote = F,sep = "\t")


str(out)
write_csv(out,"~/Desktop/KCDC/HLAimputation/03.HLA_matching/2019_Kidney_HLAmathing_table_4digit.csv",row.names = F,quote = F)
write_csv(out,"~/Desktop/KCDC/HLAimputation/03.HLA_matching/2019_Kidney_HLAmathing_table_4digit.csv",row.names = F,quote = F)





