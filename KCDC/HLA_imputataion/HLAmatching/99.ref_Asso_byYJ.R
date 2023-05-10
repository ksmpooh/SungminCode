### HLA eplet association

Covariates: age, sex, dm, cvd, retransplant, cmv_igg_reci, hcv_ab_reci, hbsag_reci, hla_ms_ab, hla_ms_dr, ind_atg, d_age, d_sex, dgf

### Make phenotype
"/Users/youngjinkim/Desktop/JOB2022/논문진행 중/TRANS/Immune_signal"
ref <- read.table("../rand_KID.txt", header=T)
pheno <- read.csv("Rejection_phenotype_coded_20230323.csv", header=T)
ref2 <- ref[!is.na(ref[,1]),]
pheno2 <- pheno[!is.na(pheno$random_id),]
out <- merge(ref2, pheno2, by.x="Random_ID", by.y="random_id")


### Without NA (in all vars)
> summary(as.factor(out2$rej_tot))
0    1
1167  473


### Association!
data <- read.table("KR.KD.immune.cell.co-signal_targetGene_RLpair.alleleMatching01.Score_Sum.txt", header=T)
pheno <- read.table("Rejection_phenotype_coded_KID_20230323.txt", header=T)
data2 <- merge(pheno, data, by.x="KCHIP_ID", by.y="KBA_ID.KR")

KCHIP_ID rej_tot AGE SEX dm cvd cmv_igg_reci hbsag_reci hcv_ab_reci
1 NIH19KT0023       1  61   1  0   0            1          0           0
2 NIH19KT0024       0  40   0  0   0            0          0           0
3 NIH19KT0025       0  47   0  0   0            1          0           0
4 NIH19KT0026       0  24   1  0   0            1          0           0
5 NIH19KT0029       0  52   0  0   0            1          0           0
6 NIH19KT0030       0  51   0  0   0            1          0           0
ind_atg D_AGE D_SEX dgf hla_ms_dr hla_ms_ab


#### Adjust #1
result <- NULL
for(i in 17:224){
  temp <- glm(paste("rej_tot ~ ",colnames(data2)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+hla_ms_dr+hla_ms_ab",sep=""), data=data2, family="binomial")
  result <- rbind(result, c(colnames(data2)[i],as.vector(summary(temp)$coefficients[2,])))
}
colnames(result) <- c("ID","Estimate","SE","Z","P")
write.table(result, "Pair_asso_adj1.txt", col.names=T, row.names=F, sep="\t", quote=F)

#### Adjust #2
result <- NULL
for(i in 17:224){
  temp <- glm(paste("rej_tot ~ ",colnames(data2)[i], "+AGE+SEX+dm+cvd+ind_atg+D_AGE+D_SEX+dgf+hla_ms_dr+hla_ms_ab",sep=""), data=data2, family="binomial")
  result <- rbind(result, c(colnames(data2)[i],as.vector(summary(temp)$coefficients[2,])))
}
colnames(result) <- c("ID","Estimate","SE","Z","P")
write.table(result, "Pair_asso_adj2.txt", col.names=T, row.names=F, sep="\t", quote=F)




### Main signal
data <- read.table("02.KR.KD.immune.cell.Main_signal_targetGene.alleleMatching01.Score_Sum.txt", header=T)
pheno <- read.table("Rejection_phenotype_coded_KID_20230323.txt", header=T)
data2 <- merge(pheno, data, by.x="KCHIP_ID", by.y="KBA_ID.KR")

#### Adjust #1
result <- NULL
for(i in 17:120){
  temp <- glm(paste("rej_tot ~ ",colnames(data2)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+hla_ms_dr+hla_ms_ab",sep=""), data=data2, family="binomial")
  result <- rbind(result, c(colnames(data2)[i],as.vector(summary(temp)$coefficients[2,])))
}
colnames(result) <- c("ID","Estimate","SE","Z","P")
write.table(result, "MainSig_asso_adj1.txt", col.names=T, row.names=F, sep="\t", quote=F)

#### Adjust #2
result <- NULL
for(i in 17:120){
  temp <- glm(paste("rej_tot ~ ",colnames(data2)[i], "+AGE+SEX+dm+cvd+ind_atg+D_AGE+D_SEX+dgf+hla_ms_dr+hla_ms_ab",sep=""), data=data2, family="binomial")
  result <- rbind(result, c(colnames(data2)[i],as.vector(summary(temp)$coefficients[2,])))
}
colnames(result) <- c("ID","Estimate","SE","Z","P")
write.table(result, "MainSig_asso_adj2.txt", col.names=T, row.names=F, sep="\t", quote=F)



### Co-signal
data <- read.table("02.KR.KD.immune.cell.co-signal_targetGene.alleleMatching01.Score_Sum.txt", header=T)
pheno <- read.table("Rejection_phenotype_coded_KID_20230323.txt", header=T)
data2 <- merge(pheno, data, by.x="KCHIP_ID", by.y="KBA_ID.KR")


#### Adjust #1
result <- NULL
for(i in 17:155){
  temp <- glm(paste("rej_tot ~ ",colnames(data2)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+hla_ms_dr+hla_ms_ab",sep=""), data=data2, family="binomial")
  result <- rbind(result, c(colnames(data2)[i],as.vector(summary(temp)$coefficients[2,])))
}
colnames(result) <- c("ID","Estimate","SE","Z","P")
write.table(result, "CoSig_asso_adj1.txt", col.names=T, row.names=F, sep="\t", quote=F)

#### Adjust #2
result <- NULL
for(i in 17:155){
  temp <- glm(paste("rej_tot ~ ",colnames(data2)[i], "+AGE+SEX+dm+cvd+ind_atg+D_AGE+D_SEX+dgf+hla_ms_dr+hla_ms_ab",sep=""), data=data2, family="binomial")
  result <- rbind(result, c(colnames(data2)[i],as.vector(summary(temp)$coefficients[2,])))
}
colnames(result) <- c("ID","Estimate","SE","Z","P")
write.table(result, "CoSig_asso_adj2.txt", col.names=T, row.names=F, sep="\t", quote=F)
