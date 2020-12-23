### merge 2019 2020
# HLA 8 gene
# (A, B, C, DRB1, DQA1, DQB2, DPA1, DPB2)
# 2019 : 254
# 2020 : 265
# all : 519


setwd("c:/Users/user/Desktop/KCDC/HLAimputation/HLAtyping/")

ref <- read.csv("2019/HLA_NGS_typing_255samples_results_202002_modify_alltype_A.B.DRB1.csv")
hla19 <- read.csv("2019/HLAtyping.allGene.result.csv")
hla20 <- read.csv("2020/2020_HLAtyping_all.csv")
head(hla19)
head(hla20)
head(ref)

ref <- ref[,c(1,2)]
colnames(ref) <- c("YSample","KID")
hla19 <- hla19[hla19$X %in% ref$Sample,]
colnames(hla19) <- c("YSample","NGS_A.1","NGS_A.2","NGS_B.1","NGS_B.2","NGS_C.1","NGS_C.2"
                     ,"NGS_DRB1.1","NGS_DRB1.2","NGS_DRB3.1","NGS_DRB3.2" 
                     ,"NGS_DQA1.1","NGS_DQA1.2","NGS_DQB1.1","NGS_DQB1.2","NGS_DPA1.1","NGS_DPA1.2","NGS_DPB1.1","NGS_DPB1.2")
hla19 <- hla19[hla19$YSample != "CDC015",]

hla19 <- merge(ref,hla19,by = "YSample")



colnames(hla20)[1:2]<-c("KID","YSample")
hla20 <- hla20[,c("KID","YSample","NGS_A.1","NGS_A.2","NGS_B.1","NGS_B.2","NGS_C.1","NGS_C.2"
                ,"NGS_DRB1.1","NGS_DRB1.2","NGS_DRB3.1","NGS_DRB3.2" 
                ,"NGS_DQA1.1","NGS_DQA1.2","NGS_DQB1.1","NGS_DQB1.2","NGS_DPA1.1","NGS_DPA1.2","NGS_DPB1.1","NGS_DPB1.2")]

head(hla20)
tail(hla20)
hla.all <- rbind(hla19,hla20)

hla.all <- hla.all[,c("KID","YSample","NGS_A.1","NGS_A.2","NGS_B.1","NGS_B.2","NGS_C.1","NGS_C.2"
                      ,"NGS_DRB1.1","NGS_DRB1.2","NGS_DQA1.1","NGS_DQA1.2","NGS_DQB1.1","NGS_DQB1.2"
                      ,"NGS_DPA1.1","NGS_DPA1.2","NGS_DPB1.1","NGS_DPB1.2")]

write.csv(hla.all,"all/HLA.type.result.8genes.merged(2019.2020).csv",row.names = F,quote = F)
