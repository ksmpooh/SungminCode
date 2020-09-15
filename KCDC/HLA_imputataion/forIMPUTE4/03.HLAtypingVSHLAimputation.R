
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/IMPUTE4/Han.ref/Result/")
setwd("c:/Users/user/Desktop/KCDC/HLAimputation/IMPUTE4/Pan.ref/Result/")
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/20200731/Pan/")
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/255sample/01.pan/")
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/255sample/02.han/")
#df <- read.table("JG.HLA.imputation_RAW.raw",header = T)
df <- read.table("test.HLA_raw.raw",header = T)
#df <- read.table("JG.HLA.imputation_RAW.raw",header = T)
head(df)
colnames(df)
ncol(df)
#grep("*DRB1*",colnames(df))
write.csv(t(colnames(df)),"header.txt",row.names = F)
#write.csv((colnames(df)),"header2.txt",row.names = F)
#t(colnames(df))
header <- read.csv("last.header.txt")
head(header)
ncol(header)
colnames(df) <- colnames(header)
df[is.na(df)]<- -1
#colnames(df) <- c("FID","IID","PAT","MAT","SEX","PHENOTYPE","HLA_A_01_P","HLA_A_0101_P","HLA_A_0103_P","HLA_A_0122N_P","HLA_A_02_P","HLA_A_0201_P","HLA_A_0202_P","HLA_A_0203_P","HLA_A_0205_P","HLA_A_0206_P","HLA_A_0207_P","HLA_A_0209_P","HLA_A_0210_P","HLA_A_0211_P","HLA_A_0217_P","HLA_A_0248_P","HLA_A_0253N_P","HLA_A_03_P","HLA_A_0301_P","HLA_A_0302_P","HLA_A_11_P","HLA_A_1101_P","HLA_A_1102_P","HLA_A_1103_P","HLA_A_11110_P","HLA_A_1112_P","HLA_A_1177_P","HLA_A_23_P","HLA_A_2301_P","HLA_A_24_P","HLA_A_2402_P","HLA_A_2403_P","HLA_A_2404_P","HLA_A_2407_P","HLA_A_2408_P","HLA_A_2410_P","HLA_A_2420_P","HLA_A_25_P","HLA_A_2501_P","HLA_A_26_P","HLA_A_2601_P","HLA_A_2602_P","HLA_A_2603_P","HLA_A_29_P","HLA_A_2901_P","HLA_A_2902_P","HLA_A_30_P","HLA_A_3001_P","HLA_A_3002_P","HLA_A_3004_P","HLA_A_3011_P","HLA_A_31_P","HLA_A_3101_P","HLA_A_3102_P","HLA_A_32_P","HLA_A_3201_P","HLA_A_33_P","HLA_A_3301_P","HLA_A_3303_P","HLA_A_34_P","HLA_A_3401_P","HLA_A_66_P","HLA_A_6601_P","HLA_A_68_P","HLA_A_6801_P","HLA_A_6802_P","HLA_A_6871_P","HLA_A_69_P","HLA_A_6901_P","HLA_A_74_P","HLA_A_7401_P","HLA_A_7402_P","HLA_C_01_P","HLA_C_0102_P","HLA_C_0103_P","HLA_C_0106_P","HLA_C_0108_P","HLA_C_0114_P","HLA_C_0130_P","HLA_C_0173_P","HLA_C_02_P","HLA_C_0202_P","HLA_C_03_P","HLA_C_0302_P","HLA_C_0303_P","HLA_C_0304_P","HLA_C_03100_P","HLA_C_0340_P","HLA_C_04_P","HLA_C_0401_P","HLA_C_0403_P","HLA_C_0406_P","HLA_C_0482_P","HLA_C_05_P","HLA_C_0501_P","HLA_C_0525_P","HLA_C_06_P","HLA_C_0602_P","HLA_C_0606_P","HLA_C_07_P","HLA_C_0701_P","HLA_C_0702_P","HLA_C_0704_P","HLA_C_0706_P","HLA_C_0718_P","HLA_C_0766_P","HLA_C_0767_P","HLA_C_08_P","HLA_C_0801_P","HLA_C_0802_P","HLA_C_0803_P","HLA_C_0822_P","HLA_C_0841_P","HLA_C_12_P","HLA_C_1202_P","HLA_C_1203_P","HLA_C_14_P","HLA_C_1402_P","HLA_C_1403_P","HLA_C_15_P","HLA_C_1502_P","HLA_C_1504_P","HLA_C_1505_P","HLA_C_1513_P","HLA_C_1543_P","HLA_C_16_P","HLA_C_1602_P","HLA_C_1604_P","HLA_C_17_P","HLA_C_1701_P","HLA_B_07_P","HLA_B_0702_P","HLA_B_0705_P","HLA_B_0706_P","HLA_B_08_P","HLA_B_0801_P","HLA_B_13_P","HLA_B_1301_P","HLA_B_1302_P","HLA_B_14_P","HLA_B_1401_P","HLA_B_1402_P","HLA_B_15_P","HLA_B_1501_P","HLA_B_1502_P","HLA_B_1503_P","HLA_B_1505_P","HLA_B_1507_P","HLA_B_1508_P","HLA_B_1511_P","HLA_B_1512_P","HLA_B_1513_P","HLA_B_1515_P","HLA_B_1517_P","HLA_B_1518_P","HLA_B_1519_P","HLA_B_1520_P","HLA_B_15220_P","HLA_B_1525_P","HLA_B_1527_P","HLA_B_1532_P","HLA_B_1546_P","HLA_B_1550_P","HLA_B_1558_P","HLA_B_18_P","HLA_B_1801_P","HLA_B_1802_P","HLA_B_27_P","HLA_B_2702_P","HLA_B_2704_P","HLA_B_2705_P","HLA_B_2706_P","HLA_B_2707_P","HLA_B_2711_P","HLA_B_2724_P","HLA_B_35_P","HLA_B_3501_P","HLA_B_3502_P","HLA_B_3503_P","HLA_B_3505_P","HLA_B_3508_P","HLA_B_3511_P","HLA_B_3514_P","HLA_B_3515_P","HLA_B_3520_P","HLA_B_3531_P","HLA_B_3543_P","HLA_B_37_P","HLA_B_3701_P","HLA_B_38_P","HLA_B_3801_P","HLA_B_3802_P","HLA_B_39_P","HLA_B_3901_P","HLA_B_3903_P","HLA_B_3905_P","HLA_B_3909_P","HLA_B_40_P","HLA_B_4001_P","HLA_B_4002_P","HLA_B_4003_P","HLA_B_4005_P","HLA_B_4006_P","HLA_B_4040_P","HLA_B_41_P","HLA_B_4101_P","HLA_B_4102_P","HLA_B_42_P","HLA_B_44_P","HLA_B_4402_P","HLA_B_4403_P","HLA_B_4409_P","HLA_B_4446_P","HLA_B_45_P","HLA_B_4501_P","HLA_B_46_P","HLA_B_4601_P","HLA_B_4603_P","HLA_B_47_P","HLA_B_4701_P","HLA_B_48_P","HLA_B_4801_P","HLA_B_4802_P","HLA_B_4803_P","HLA_B_49_P","HLA_B_4901_P","HLA_B_50_P","HLA_B_5001_P","HLA_B_51_P","HLA_B_5101_P","HLA_B_5102_P","HLA_B_5104_P","HLA_B_5107_P","HLA_B_5108_P","HLA_B_5142_P","HLA_B_52_P","HLA_B_5201_P","HLA_B_53_P","HLA_B_5301_P","HLA_B_54_P","HLA_B_5401_P","HLA_B_55_P","HLA_B_5501_P","HLA_B_5502_P","HLA_B_5504_P","HLA_B_5512_P","HLA_B_5524_P","HLA_B_56_P","HLA_B_5601_P","HLA_B_5603_P","HLA_B_5604_P","HLA_B_57_P","HLA_B_5701_P","HLA_B_58_P","HLA_B_5801_P","HLA_B_59_P","HLA_B_5901_P","HLA_B_67_P","HLA_B_6701_P","HLA_B_78_P","HLA_B_7802_P","HLA_B_81_P","HLA_B_8101_P","HLA_DRB1_01_P","HLA_DRB1_0101_P","HLA_DRB1_0102_P","HLA_DRB1_03_P","HLA_DRB1_0301_P","HLA_DRB1_0307_P","HLA_DRB1_0317_P","HLA_DRB1_0324_P","HLA_DRB1_04_P","HLA_DRB1_0401_P","HLA_DRB1_0402_P","HLA_DRB1_0403_P","HLA_DRB1_0404_P","HLA_DRB1_0405_P","HLA_DRB1_0406_P","HLA_DRB1_0407_P","HLA_DRB1_0408_P","HLA_DRB1_0410_P","HLA_DRB1_07_P","HLA_DRB1_0701_P","HLA_DRB1_08_P","HLA_DRB1_0801_P","HLA_DRB1_0802_P","HLA_DRB1_0803_P","HLA_DRB1_0809_P","HLA_DRB1_09_P","HLA_DRB1_0901_P","HLA_DRB1_10_P","HLA_DRB1_1001_P","HLA_DRB1_11_P","HLA_DRB1_1101_P","HLA_DRB1_1103_P","HLA_DRB1_1104_P","HLA_DRB1_1106_P","HLA_DRB1_1111_P","HLA_DRB1_12_P","HLA_DRB1_1201_P","HLA_DRB1_1202_P","HLA_DRB1_1210_P","HLA_DRB1_1217_P","HLA_DRB1_13_P","HLA_DRB1_1301_P","HLA_DRB1_1302_P","HLA_DRB1_1312_P","HLA_DRB1_1327_P","HLA_DRB1_1367_P","HLA_DRB1_14_P","HLA_DRB1_1402_P","HLA_DRB1_1403_P","HLA_DRB1_1404_P","HLA_DRB1_1405_P","HLA_DRB1_1407_P","HLA_DRB1_14141_P","HLA_DRB1_1454_P","HLA_DRB1_15_P","HLA_DRB1_1501_P","HLA_DRB1_1502_P","HLA_DRB1_1504_P","HLA_DRB1_16_P","HLA_DRB1_1601_P","HLA_DRB1_1602_P","HLA_DQA1_01_P","HLA_DQA1_0101_P","HLA_DQA1_0102_P","HLA_DQA1_0103_P","HLA_DQA1_0104_P","HLA_DQA1_0105_P","HLA_DQA1_02_P","HLA_DQA1_0201_P","HLA_DQA1_03_P","HLA_DQA1_0301_P","HLA_DQA1_0302_P","HLA_DQA1_0303_P","HLA_DQA1_04_P","HLA_DQA1_0401_P","HLA_DQA1_05_P","HLA_DQA1_0501_P","HLA_DQA1_0503_P","HLA_DQA1_0505_P","HLA_DQA1_0506_P","HLA_DQA1_0508_P","HLA_DQA1_0509_P","HLA_DQA1_06_P","HLA_DQA1_0601_P","HLA_DQB1_02_P","HLA_DQB1_0201_P","HLA_DQB1_0202_P","HLA_DQB1_03_P","HLA_DQB1_0301_P","HLA_DQB1_0302_P","HLA_DQB1_0303_P","HLA_DQB1_0305_P","HLA_DQB1_03100_P","HLA_DQB1_0319_P","HLA_DQB1_0329_P","HLA_DQB1_04_P","HLA_DQB1_0401_P","HLA_DQB1_0402_P","HLA_DQB1_05_P","HLA_DQB1_0501_P","HLA_DQB1_0502_P","HLA_DQB1_0503_P","HLA_DQB1_06_P","HLA_DQB1_0601_P","HLA_DQB1_0602_P","HLA_DQB1_0603_P","HLA_DQB1_0604_P","HLA_DQB1_0609_P","HLA_DQB1_0610_P","HLA_DPA1_01_P","HLA_DPA1_0103_P","HLA_DPA1_02_A","HLA_DPA1_0201_P","HLA_DPA1_0202_P","HLA_DPA1_03_P","HLA_DPA1_0301_P","HLA_DPA1_04_P","HLA_DPA1_0401_P","HLA_DPB1_01_P","HLA_DPB1_0101_P","HLA_DPB1_02_P","HLA_DPB1_0201_P","HLA_DPB1_0202_P","HLA_DPB1_03_P","HLA_DPB1_0301_P","HLA_DPB1_04_P","HLA_DPB1_0401_P","HLA_DPB1_0402_P","HLA_DPB1_05_P","HLA_DPB1_0501_P","HLA_DPB1_06_P","HLA_DPB1_0601_P","HLA_DPB1_09_P","HLA_DPB1_0901_P","HLA_DPB1_10_P","HLA_DPB1_1001_P","HLA_DPB1_10401_P","HLA_DPB1_10501_P","HLA_DPB1_10601_P","HLA_DPB1_10701_P","HLA_DPB1_11_P","HLA_DPB1_1101_P","HLA_DPB1_13_P","HLA_DPB1_1301_P","HLA_DPB1_13501_P","HLA_DPB1_13801_P","HLA_DPB1_14_P","HLA_DPB1_1401_P","HLA_DPB1_15_P","HLA_DPB1_1501_P","HLA_DPB1_16_P","HLA_DPB1_1601_P","HLA_DPB1_17_P","HLA_DPB1_1701_P","HLA_DPB1_19_P","HLA_DPB1_1901_P","HLA_DPB1_21_P","HLA_DPB1_2101_P","HLA_DPB1_23_P","HLA_DPB1_2301_P","HLA_DPB1_26_P","HLA_DPB1_2601_P","HLA_DPB1_27_P","HLA_DPB1_2701_P","HLA_DPB1_28_P","HLA_DPB1_2801_P","HLA_DPB1_45_P","HLA_DPB1_4501_P")

colnames(df)



library(stringr)

A <- df[,c(1,2,grep("HLA_A_*",colnames(df)))]
B <- df[,c(1,2,grep("HLA_B_*",colnames(df)))]
C <- df[,c(1,2,grep("HLA_C_*",colnames(df)))]
DRB1 <- df[,c(1,2,grep("*DRB1*",colnames(df)))]
#DRB3 <- df[,c(1,2,grep("*DRB3*",colnames(df)))]

DPB1 <- df[,c(1,2,grep("*DPB1*",colnames(df)))]
DPA1 <- df[,c(1,2,grep("*DPA1*",colnames(df)))]

DQB1 <- df[,c(1,2,grep("*DQB1*",colnames(df)))]
DQA1 <- df[,c(1,2,grep("*DQA1*",colnames(df)))]


hla.subset <- function(df,n){
  if( (n ==2) || (n == 4)){
    i = paste(n,"(nvalue) is OK")
    print(i)
    temp <- df[,1:2]
    #    print(temp)
    if(n == 2){
      for (i in 3:ncol(df)){
        #print(as.integer(str_split_fixed(colnames(df)[i],"_",4)[3]))
        if(as.integer(str_split_fixed(colnames(df)[i],"_",4)[3]) < 100){
          b <- df[,c(1,2,i)]
          temp <- merge(temp,b)
        }
      }
    }else{
      for (i in 3:ncol(df)){
        if(as.integer(str_split_fixed(colnames(df)[i],"_",4)[3]) >= 100){
          b <- df[,c(1,2,i)]
          temp <- merge(temp,b)
        }
      }
    }
    return(temp)
  }else{
    a = paste(n,"(nvalue) is wrong, only for 2,4")
    print(a)
    return(0)
  }
}
#DRB_td <- hla.subset(DRB,3)
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


#hla.subset(df,digit)
#hla.find(df_digit,"concept")



###################### 2digit
###################### 2digit

DRB1_td <- hla.subset(DRB1,2)
head(DRB1_td)
str(DRB1_td)
DRB1_td <- hla.find(DRB1_td,"DRB1")
table(DRB1_td$DRB3)


DQA1_td <- hla.subset(DQA1,2)
DQA1_td <- hla.find(DQA1_td,"DQA1")
DQB1_td <- hla.subset(DQB1,2)
DQB1_td <- hla.find(DQB1_td,"DQB1")


DPA1_td <- hla.subset(DPA1,2)
DPA1_td <- hla.find(DPA1_td,"DPA1")
DPB1_td <- hla.subset(DPB1,2)
DPB1_td <- hla.find(DPB1_td,"DPB1")


grep("*N_P",colnames(A))
colnames(A)[6]<-"HLA_A_0122_P"
colnames(A)[19]<-"HLA_A_0253_P"
A_td <- hla.subset(A,2)
A_td <- hla.find(A_td,"A")

B_td <- hla.subset(B,2)
B_td <- hla.find(B_td,"B")

C_td <- hla.subset(C,2)
C_td <- hla.find(C_td,"C")


table(DQA1_td$IMP_DQA1.3)
table(DQB1_td$IMP_DQB1.3)
table(DPA1_td$IMP_DPA1.3)
table(DPB1_td$IMP_DPB1.3)
table(DRB1_td$IMP_DRB1.3)
table(C_td$IMP_C.3)
table(B_td$IMP_B.3)
table(A_td$IMP_A.3)
head(A_td)

##HAN 2digit
DQA1_td[!is.na(DQA1_td$IMP_DQA1.3),]$IID
##PAN 2digit
DPA1_td[!is.na(DPA1_td$IMP_DPA1.3),]$IID

out <- merge(A_td[,c('IID','IMP_A.1','IMP_A.2')],B_td[,c('IID','IMP_B.1','IMP_B.2')],by = 'IID')
out <- merge(out,C_td[,c('IID','IMP_C.1','IMP_C.2')],by = 'IID')
out <- merge(out,DRB1_td[,c('IID','IMP_DRB1.1','IMP_DRB1.2')],by = 'IID')
out <- merge(out,DPA1_td[,c('IID','IMP_DPA1.1','IMP_DPA1.2')],by = 'IID')
out <- merge(out,DPB1_td[,c('IID','IMP_DPB1.1','IMP_DPB1.2')],by = 'IID')
out <- merge(out,DQA1_td[,c('IID','IMP_DQA1.1','IMP_DQA1.2')],by = 'IID')
out <- merge(out,DQB1_td[,c('IID','IMP_DQB1.1','IMP_DQB1.2')],by = 'IID')


write.csv(out,"HLAimputation.all.gene.2digit.Result.csv",row.names = F,quote = F)

head(out)
#####merge NGS

out <- read.csv("HLAimputation.all.gene.2digit.Result.csv")
ngs <- read.csv("../../../../transplantation/HLAtyping/20200828/HLAtyping.alle.gene.2digit.csv")
head(ngs)
ncol(ngs)
out <- merge(out,ngs,by.x = "IID",by.y = "KID")
ncol(out)
head(out)
out <-out[,c(1,18,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)]

row.names(out) <- out$IID
#out[DQA1_td[!is.na(DQA1_td$IMP_DQA1.3),]$IID,]
#out[DPA1_td[!is.na(DPA1_td$IMP_DPA1.3),]$IID,]
head(out)
nrow(out)
write.csv(out,"MERGE.impResult.hlatyping.all.gene.2digit.csv",row.names = F,quote = F)




cmp.result.2digit <- read.csv("compare.IMPvsNGS.all.gene.2digit.csv")

rownames(cmp.result.2digit) <- cmp.result.2digit$IID
nrow(cmp.result.2digit[DQA1_td[!is.na(DQA1_td$IMP_DQA1.3),]$IID,c("DQA1.match","DQA1.wrong","DQA1.empty")])






####################4 digit

DRB1_fd <- hla.subset(DRB1,4)
head(DRB1_fd)
DRB1_fd <- hla.find(DRB1_fd,"DRB1")
table(DRB1_fd$DRB3)


DQA1_fd <- hla.subset(DQA1,4)
DQA1_fd <- hla.find(DQA1_fd,"DQA1")
DQB1_fd <- hla.subset(DQB1,4)
DQB1_fd <- hla.find(DQB1_fd,"DQB1")


DPA1_fd <- hla.subset(DPA1,4)
DPA1_fd <- hla.find(DPA1_fd,"DPA1")

DPB1_fd <- hla.subset(DPB1,4)
DPB1_fd <- hla.find(DPB1_fd,"DPB1")


grep("*N_P",colnames(A))
colnames(A)[6]<-"HLA_A_0122_P"
colnames(A)[19]<-"HLA_A_0253_P"

A_fd <- hla.subset(A,4)
A_fd <- hla.find(A_fd,"A")

B_fd <- hla.subset(B,4)
B_fd <- hla.find(B_fd,"B")

C_fd <- hla.subset(C,4)
C_fd <- hla.find(C_fd,"C")


table(DQA1_fd$IMP_DQA1.3)
table(DQB1_fd$IMP_DQB1.3)
table(DPA1_fd$IMP_DPA1.3)
table(DPB1_fd$IMP_DPB1.3)
table(DRB1_fd$IMP_DRB1.3)
table(C_fd$IMP_C.3)
table(B_fd$IMP_B.3)
table(A_fd$IMP_A.3)
head(A_fd)


out <- merge(A_fd[,c('IID','IMP_A.1','IMP_A.2')],B_fd[,c('IID','IMP_B.1','IMP_B.2')],by = 'IID')
out <- merge(out,C_fd[,c('IID','IMP_C.1','IMP_C.2')],by = 'IID')
out <- merge(out,DRB1_fd[,c('IID','IMP_DRB1.1','IMP_DRB1.2')],by = 'IID')
out <- merge(out,DPA1_fd[,c('IID','IMP_DPA1.1','IMP_DPA1.2')],by = 'IID')
out <- merge(out,DPB1_fd[,c('IID','IMP_DPB1.1','IMP_DPB1.2')],by = 'IID')
out <- merge(out,DQA1_fd[,c('IID','IMP_DQA1.1','IMP_DQA1.2')],by = 'IID')
out <- merge(out,DQB1_fd[,c('IID','IMP_DQB1.1','IMP_DQB1.2')],by = 'IID')


write.csv(out,"HLAimputation.all.gene.4digit.Result.csv",row.names = F,quote = F)

head(out)
#####merge NGS
out <- read.csv("HLAimputation.all.gene.4digit.Result.csv")
ngs <- read.csv("../../../../transplantation/HLAtyping/20200828/HLAtyping.alle.gene.4digit.csv")
head(ngs)
ncol(ngs)
out <- merge(out,ngs,by.x = "IID",by.y = "KID")
ncol(out)
head(out)
out <-out[,c(1,18,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)]
write.csv(out,"MERGE.impResult.hlatyping.all.gene.4digit.csv",row.names = F,quote = F)



