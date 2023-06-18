##### make genetic map for 1000 genome progect

library(readr)
library(tidyverse)
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/")
setwd("~/Desktop/KCDC/HLAimputation/")
#df <- readr::read_table2("1000genome/20181129_HLA_types_full_1000_Genomes_Project_panel.txt")
df <- read.table("1000genome/20181129_HLA_types_full_1000_Genomes_Project_panel.txt",header = T,sep = "\t",na.strings = "")
#df <- read.csv("1000genome/20181129_HLA_types_full_1000_Genomes_Project_panel.csv",header = T,na.strings = "none")
head(df)
sapply(colnames(df), function(x) grep("None", df[,x]))
sapply(colnames(df), function(x) grep("*", df[,x]))


#NA
df[df$Sample.ID == "HG01051",]

# *
df[df$Sample.ID == "NA19038",]

table(df$Sample.ID %in% fam$V1)
head(df$Sample.ID[99:150])


fam <- read.table("1000genome/1kgp.phase3.chr6.MHC.fam")
head(fam)

df <- df[df$Sample.ID %in% fam$V1,]
df[df == "None"] <- 0
df[is.na(df)] <- 0



df$FID <- df$Sample.ID
df$pID <- 0
df$mID <- 0
df$SEX <- 0
df$PHENO <- 0
df$HLA.DPA1.1 <- 0
df$HLA.DPA1.2 <- 0
df$HLA.DPB1.1 <- 0
df$HLA.DPB1.2 <- 0
df$HLA.DQA1.1 <- 0
df$HLA.DQA1.2 <- 0




colnames(df)

df <- df[,c("FID","FID","pID","mID","SEX","PHENO","HLA.A.1","HLA.A.2",
            "HLA.B.1","HLA.B.2","HLA.C.1","HLA.C.2","HLA.DPA1.1","HLA.DPA1.2","HLA.DPB1.1","HLA.DPB1.2",
            "HLA.DQA1.1","HLA.DQA1.2","HLA.DQB1.1","HLA.DQB1.2","HLA.DRB1.1","HLA.DRB1.2")]
head(df)


write.table(df,"1000genome/HLAtyping.1000genomePhase3.rouph.txt",sep = "\t",col.names = T,row.names = F,quote = F)


##### 20230327 for HLA-TAPAS

##### make genetic map for 1000 genome progect

library(readr)
library(tidyverse)
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/")
setwd("~/Desktop/KCDC/HLAimputation/")
#df <- readr::read_table2("1000genome/20181129_HLA_types_full_1000_Genomes_Project_panel.txt")
df <- read.table("1000genome/20181129_HLA_types_full_1000_Genomes_Project_panel.txt",header = T,sep = "\t",na.strings = "")
#df <- read.csv("1000genome/20181129_HLA_types_full_1000_Genomes_Project_panel.csv",header = T,na.strings = "none")
head(df)
sapply(colnames(df), function(x) grep("None", df[,x]))
sapply(colnames(df), function(x) grep("*", df[,x]))


#NA
df[df$Sample.ID == "HG01051",]

# *
df[df$Sample.ID == "NA19038",]

table(df$Sample.ID %in% fam$V1)
head(df$Sample.ID[99:150])


fam <- read.table("1000genome/1kgp.phase3.chr6.MHC.fam")
head(fam)

df <- df[df$Sample.ID %in% fam$V1,]
df[df == "None"] <- 0
df[is.na(df)] <- 0
head(df)

df[,4]
gene = str_split_fixed(str_replace(colnames(df)[4],"HLA.",""),"\\.",2)[,1]
gene
for (i in 4:ncol(df)) {
  gene = str_split_fixed(str_replace(colnames(df)[i],"HLA.",""),"\\.",2)[,1]
  #df %>% mutate(colnames(df)[i] = paste0(gene,"*",colnames(df)[i]))
  df[,i] <- paste0(gene,"*",df[,i])
}
head(df)
grepl("HLA",colnames(df))
colnames(df)[grepl("HLA",colnames(df))]
colnames(df)[grepl("HLA",colnames(df))] <-  str_replace_all(colnames(df)[grepl("HLA",colnames(df))],"HLA.","HLA_")
head(df)
write.table(df,"1000genome/HLAtyping.1000genomePhase3.SampleInfo.typeInfo.txt",col.names = T,row.names = F,quote = F,sep = "\t")


df$FID <- df$Sample.ID
df$IID <- df$Sample.ID
df$pID <- 0
df$mID <- 0
df$SEX <- 0
df$PHENO <- 0
df$HLA_DPA1.1 <- 0
df$HLA_DPA1.2 <- 0
df$HLA_DPB1.1 <- 0
df$HLA_DPB1.2 <- 0
df$HLA_DQA1.1 <- 0
df$HLA_DQA1.2 <- 0




colnames(df)

df <- df[,c("IID","FID","pID","mID","SEX","PHENO","HLA_A.1","HLA_A.2",
            "HLA_B.1","HLA_B.2","HLA_C.1","HLA_C.2","HLA_DPA1.1","HLA_DPA1.2","HLA_DPB1.1","HLA_DPB1.2",
            "HLA_DQA1.1","HLA_DQA1.2","HLA_DQB1.1","HLA_DQB1.2","HLA_DRB1.1","HLA_DRB1.2")]
head(df)


write.table(df,"1000genome/HLAtyping.1000genomePhase3.rouph.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(df,"1000genome/HLAtyping.1000genomePhase3.rouph_forHLAPATAS.txt",sep = "\t",col.names = F,row.names = F,quote = F)


####
##### make genetic map for 1000 genome progect

library(readr)
library(tidyverse)
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/")
setwd("~/Desktop/KCDC/HLAimputation/")
#df <- readr::read_table2("1000genome/20181129_HLA_types_full_1000_Genomes_Project_panel.txt")
df <- read.table("1000genome/20181129_HLA_types_full_1000_Genomes_Project_panel.txt",header = T,sep = "\t",na.strings = "")
#df <- read.csv("1000genome/20181129_HLA_types_full_1000_Genomes_Project_panel.csv",header = T,na.strings = "none")
head(df)
sapply(colnames(df), function(x) grep("None", df[,x]))
sapply(colnames(df), function(x) grep("*", df[,x]))


#NA
df[df$Sample.ID == "HG01051",]

# *
df[df$Sample.ID == "NA19038",]

table(df$Sample.ID %in% fam$V1)
head(df$Sample.ID[99:150])


fam <- read.table("1000genome/1kgp.phase3.chr6.MHC.fam")
head(fam)

df <- df[df$Sample.ID %in% fam$V1,]
df[df == "None"] <- 0
df[is.na(df)] <- 0



df$FID <- df$Sample.ID
df$pID <- 0
df$mID <- 0
df$SEX <- 0
df$PHENO <- 0
df$HLA.DPA1.1 <- 0
df$HLA.DPA1.2 <- 0
df$HLA.DPB1.1 <- 0
df$HLA.DPB1.2 <- 0
df$HLA.DQA1.1 <- 0
df$HLA.DQA1.2 <- 0




colnames(df)

df <- df[,c("FID","FID","pID","mID","SEX","PHENO","HLA.A.1","HLA.A.2",
            "HLA.B.1","HLA.B.2","HLA.C.1","HLA.C.2","HLA.DPA1.1","HLA.DPA1.2","HLA.DPB1.1","HLA.DPB1.2",
            "HLA.DQA1.1","HLA.DQA1.2","HLA.DQB1.1","HLA.DQB1.2","HLA.DRB1.1","HLA.DRB1.2")]
head(df)


write.table(df,"1000genome/HLAtyping.1000genomePhase3.rouph.txt",sep = "\t",col.names = T,row.names = F,quote = F)


##### 20230327 for HLA-TAPAS

##### make genetic map for 1000 genome progect

library(readr)
library(tidyverse)

df <- read.table("20181129_HLA_types_full_1000_Genomes_Project_panel.txt",header = T,sep = "\t",na.strings = "")

head(df)
sapply(colnames(df), function(x) grep("None", df[,x]))
sapply(colnames(df), function(x) grep("*", df[,x]))

fam <- read.table("1kgp.phase3.chr6.MHC.fam")
df <- df[df$Sample.ID %in% fam$V1,]
df[df == "None"] <- 0

df[is.na(df)] <- 0
gene = str_split_fixed(str_replace(colnames(df)[4],"HLA.",""),"\\.",2)[,1]
gene
for (i in 4:ncol(df)) {
  gene = str_split_fixed(str_replace(colnames(df)[i],"HLA.",""),"\\.",2)[,1]
  df[,i] <- paste0(gene,"*",df[,i])
}
colnames(df)[grepl("HLA",colnames(df))] <-  str_replace_all(colnames(df)[grepl("HLA",colnames(df))],"HLA.","HLA_")
write.table(df,"1000genome/HLAtyping.1000genomePhase3.SampleInfo.typeInfo.txt",col.names = T,row.names = F,quote = F,sep = "\t")


df$FID <- df$Sample.ID
df$IID <- df$Sample.ID
df$pID <- 0
df$mID <- 0
df$SEX <- 0
df$PHENO <- 0
df$HLA_DPA1.1 <- 0
df$HLA_DPA1.2 <- 0
df$HLA_DPB1.1 <- 0
df$HLA_DPB1.2 <- 0
df$HLA_DQA1.1 <- 0
df$HLA_DQA1.2 <- 0

df <- df[,c("IID","FID","pID","mID","SEX","PHENO","HLA_A.1","HLA_A.2",
            "HLA_B.1","HLA_B.2","HLA_C.1","HLA_C.2","HLA_DPA1.1","HLA_DPA1.2","HLA_DPB1.1","HLA_DPB1.2",
            "HLA_DQA1.1","HLA_DQA1.2","HLA_DQB1.1","HLA_DQB1.2","HLA_DRB1.1","HLA_DRB1.2")]

write.table(df,"1000genome/HLAtyping.1000genomePhase3.rouph_forHLAPATAS.txt",sep = "\t",col.names = F,row.names = F,quote = F)

