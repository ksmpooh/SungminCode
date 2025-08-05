library(tidyverse)
setwd("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/KBA_QC/2nd/")

df <- read.table("PCA.txt",header = T)
df <- read.table("extra/PCA.txt",header = T)
ref <- readxl::read_xlsx("~/Desktop/KCDC/transplantation/00.sampleInfo/KOTRY.KBAv1.1????????????????????????????????????_20250311_v2.xlsx",sheet = 2)
ref2 <- readxl::read_xlsx("~/Desktop/KCDC/transplantation/00.sampleInfo/ALL(2019and2020).sampleID.withbCODE.xlsx")
yslist <- read.table("ys.list.txt")
head(yslist)
head(df)
head(ref2)
head(ref)
dim(ref)
dim(ref2)
head(ref2)
ref2 %>% filter(KBA_ID %in% c("NIH19KT63742","NIH19KT65602","NIH19KT65722"))
ref2[(ref2$bCODE %in% ref$bCODE),]$bCODE -> qcin_list
head(ref2)
ref2 %>% mutate(QC = ifelse(bCODE %in% qcin_list,1,0)) %>% filter(type %in% c("KR","KD")) %>%
  filter(QC == 1) -> qcin

ref2 %>% mutate(QC = ifelse(bCODE %in% qcin_list,1,0)) %>% filter(type %in% c("KR","KD")) %>%
  filter(QC == 0) -> qcout

head(qcout)
head(qcin)
dim(qcin)
dim(ref)


head(df)
df %>% mutate(ID = str_remove_all(str_split_fixed(FID,"_",6)[,6],".CEL")) %>% filter(str_detect(ID,"NIH")) -> df_nih
head(df_nih)
head(qcin)
df_nih %>% filter(PC1 <= -0.1) %>% count(ID %in% qcin$KBA_ID)
df_nih %>% filter(PC1 <= -0.1) %>% count(ID %in% qcout$KBA_ID)

#df_nih %>% 

head(ref2)
head(df)
df %>% mutate(ID = str_remove_all(str_split_fixed(FID,"_",6)[,6],".CEL")) %>% mutate(group = ifelse(str_detect(ID,"NIH"),"NIH","YS")) %>% #-> check
  left_join(ref2 %>% select(KBA_ID,type) %>% rename(ID = KBA_ID)) %>% #head()
   mutate(qc = ifelse(ID %in% qcin$KBA_ID,"QC in","QC out")) %>% #count(group,qc)
  #mutate(qc = ifelse(ID %in% qcin$KBA_ID,"QC in","QC out"))
  #mutate(qc = ifelse(ID %in% qcout$KBA_ID,"0","1")) %>% count(group,qc) 
  mutate(qc = ifelse(group == "YS","YS",qc)) %>% #count(group,qc)
  ggplot(aes(x=PC1,y=PC2,color=qc)) +
    #geom_point(aes(shape = as.factor(type))) +   # type을 factor로 변환하여 shape으로 사용
  geom_point(size = 1,alpha=0.8) + 
  geom_hline(yintercept = 0.1, linetype = "dashed",color = "purple") +  # y = 0.1 수평선
  geom_hline(yintercept = -0.1, linetype = "dashed",color = "purple") +  # y = -0.1 수평선
  geom_vline(xintercept = 0.1, linetype = "dashed",color = "purple") +  # x = 0.1 수직선
  geom_vline(xintercept = -0.05, linetype = "dashed",color = "purple")  # x = -0.1 수직선
  #scale_x_continuous(limits = c(-1, 1),breaks = seq(-1, 1, by = 0.1)) + 
  #scale_y_continuous(limits = c(-1, 1),breaks = seq(-1, 1, by = 0.1)) +  
  theme_bw()


df %>% mutate(ID = str_remove_all(str_split_fixed(FID,"_",6)[,6],".CEL")) %>% mutate(group = ifelse(str_detect(ID,"NIH"),"NIH","YS")) %>% #-> check
  mutate(ID = ifelse(str_detect(ID,"_2"),str_remove_all(ID,"_2"),ID)) %>%
    left_join(ref2 %>% select(KBA_ID,type) %>% rename(ID = KBA_ID)) %>% #head()
    mutate(qc = ifelse(ID %in% qcin$KBA_ID,"QC in","QC out")) %>% filter(group == "NIH") -> df_qc

head(df)
df %>% filter(FID %in% yslist$V1)
head(yslist)
yslist %>% filter(!(V1 %in% df$FID))
yslist %>% filter(str_detect(V1,"NIH")) -> yslist_NIH
dim(yslist_NIH)
head(df_qc)
getwd()
head(qcin)
df %>% mutate(ID = str_remove_all(str_split_fixed(FID,"_",6)[,6],".CEL")) %>% #mutate(group = ifelse(str_detect(ID,"NIH"),"NIH","YS")) %>% #-> check
  mutate(ID = ifelse(str_detect(ID,"_2"),str_remove_all(ID,"_2"),ID)) %>%
  filter(PC1 < 0.1 & PC1 > -0.1 & PC2 < 0.1 & PC2 > -0.1) %>% 
  filter(ID %in% qcin$KBA_ID) %>% 
  count(FID %in% yslist_NIH$V1)




df_qc %>% filter(PC1 < 0.1 & PC1 > -0.1 & PC2 < 0.1 & PC2 > -0.1) %>% 
  filter(ID %in% qcin$KBA_ID) %>% #head()
  select(FID,ID) -> df_qcin_NIH
qcout
head(df_qcin_NIH)
head(yslist)
df %>% mutate(ID = str_remove_all(str_split_fixed(FID,"_",6)[,6],".CEL")) %>% #mutate(group = ifelse(str_detect(ID,"NIH"),"NIH","YS")) %>% #-> check
  mutate(ID = ifelse(str_detect(ID,"_2"),str_remove_all(ID,"_2"),ID)) %>%
  filter(FID %in% yslist$V1) %>% filter(!(FID %in% df_qcin_NIH$FID)) %>% #dim()
  filter(!(ID %in% qcout$KBA_ID)) %>%
  select(FID,ID)-> extra_ys_list

head(extra_ys_list)



df_qcin_NIH %>% rbind(extra_ys_list) %>% mutate(ys_ID = ifelse(FID %in% yslist$V1,yslist$V2,ID)) %>% 
    write.table("change.ID.finalQCin.list.txt",col.names = F,quote = F,sep = "\t",row.names = F)

df_qc %>% filter(PC1 < 0.1 & PC1 > -0.1 & PC2 < 0.1 & PC2 > -0.1) %>% 
  select(ID,type,qc) %>% rename(KBA_ID = ID)-> df_qc1

df_qc1 %>% filter(type == "KR") %>% left_join(ref2 %>% select(KBA_ID,ref)) %>% #head()
  merge(df_qc1 %>% filter(type == "KD") %>% left_join(ref2 %>% select(KBA_ID,ref)),all = T,by='ref') %>% #head()
  count(type.x,type.y)
  count(type.x,type.y,qc.x,qc.y)

df_qc %>% filter(PC1 < 0.1 & PC1 > -0.05 & PC2 < 0.1 & PC2 > -0.1) %>% 
  select(ID,type,qc) %>% rename(KBA_ID = ID)-> df_qc2
df_qc2 %>% count(type)
#df_qc2 %>% filter(type == "KR") %>% left_join(ref2 %>% select(KBA_ID,ref)) %>% #head()
  merge(df_qc2 %>% filter(type == "KD") %>% left_join(ref2 %>% select(KBA_ID,ref)),all = T,by='ref') %>% #head()
  count(type.x,type.y)
  count(type.x,type.y,qc.x,qc.y)



############################## tranplantation AR QC
setwd("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/KBA_QC/1st/")
miss <-read.table("MISS.imiss",header = T)
het <- read.table("HET.het", header = T)


miss <- cbind(miss, CR=((1 - miss$F_MISS)*100))
het <- cbind(het, HET=((het$N.NM. - het$O.HOM.)/het$N.NM.)*100)

lowSample <- merge(miss, het, by="FID")

#pdf("2nd.miss_het.pdf",height=7,width=10)
plot(lowSample$HET, lowSample$F_MISS, 
     xlim=c(13,22), ylim=c(0,0.1), 
     xlab="heterozygosity rate",ylab="missing rate", 
     main="Missing vs heterozygosity", col=rgb(0,0,1,0.3), pch=16)
abline(v=15, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(v=18, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(h=0.03, col=rgb(1,0,0,1), lty=3, lwd=2)
points(lowSample[lowSample$HET < 15 | 18 < lowSample$HET | 0.03 < lowSample$F_MISS,]$HET,
       lowSample[lowSample$HET < 15 | 18 < lowSample$HET | 0.03 < lowSample$F_MISS,]$F_MISS,
       col=rgb(1,0,0,0.3), pch=16)
#dev.off()

rmList <- lowSample[0.03 < lowSample$F_MISS | lowSample$HET < 15 | 18 < lowSample$HET,]
dim(rmList)
rmList %>% mutate(ID = str_remove_all(str_split_fixed(FID,"_",6)[,6],".CEL")) %>% filter(ID %in% qcin$KBA_ID)


head(pca)
pca <- read.table("PCA.txt", header=T) 
pca %>% mutate(ID = str_remove_all(str_split_fixed(FID,"_",6)[,6],".CEL")) %>% mutate(group = ifelse(str_detect(ID,"NIH"),"NIH","YS")) -> pca
#pca <- read.table("PCA2.txt", header=T)

#pdf("./PREG.2nd.QC_PCA2.pdf", height = 10, width = 10)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.5),#shape=pca$gruop,
#     xlim=c(-0.5, 0.5), ylim=c(-0.5,0.5),
     xlab="PC1", ylab="PC2", main="PCA", cex=0.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)



points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.7), cex=0.5, pch=16)

dev.off()
pca %>% filter()


#### 2nd
setwd("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/KBA_QC/2nd/")
miss <-read.table("MISS.imiss",header = T)
het <- read.table("HET.het", header = T)


miss <- cbind(miss, CR=((1 - miss$F_MISS)*100))
het <- cbind(het, HET=((het$N.NM. - het$O.HOM.)/het$N.NM.)*100)

lowSample <- merge(miss, het, by="FID")

#pdf("2nd.miss_het.pdf",height=7,width=10)
plot(lowSample$HET, lowSample$F_MISS, 
     xlim=c(13,22), ylim=c(0,0.1), 
     xlab="heterozygosity rate",ylab="missing rate", 
     main="Missing vs heterozygosity", col=rgb(0,0,1,0.3), pch=16)
abline(v=15, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(v=18, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(h=0.03, col=rgb(1,0,0,1), lty=3, lwd=2)
points(lowSample[lowSample$HET < 15 | 18 < lowSample$HET | 0.03 < lowSample$F_MISS,]$HET,
       lowSample[lowSample$HET < 15 | 18 < lowSample$HET | 0.03 < lowSample$F_MISS,]$F_MISS,
       col=rgb(1,0,0,0.3), pch=16)
#dev.off()

rmList <- lowSample[0.03 < lowSample$F_MISS | lowSample$HET < 15 | 18 < lowSample$HET,]
dim(rmList)
rmList %>% mutate(ID = str_remove_all(str_split_fixed(FID,"_",6)[,6],".CEL")) %>% filter(ID %in% qcin$KBA_ID)


head(pca)
pca <- read.table("PCA.txt", header=T)
pca <- read.table("~/Desktop/KCDC/transplantation/QCKD_2019/01.1stQC/PCA.txt", header=T)
pca %>% mutate(ID = str_remove_all(str_split_fixed(FID,"_",6)[,6],".CEL")) %>% mutate(group = ifelse(str_detect(ID,"NIH"),"NIH","YS")) -> pca
#pca <- read.table("PCA2.txt", header=T)

#pdf("./PREG.2nd.QC_PCA2.pdf", height = 10, width = 10)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.5),#shape=pca$gruop,
          xlim=c(-0.5, 0.5), ylim=c(-0.5,0.5),
     xlab="PC1", ylab="PC2", main="PCA", cex=0.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)



points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.7), cex=0.5, pch=16)

dev.off()
pca %>% filter()






