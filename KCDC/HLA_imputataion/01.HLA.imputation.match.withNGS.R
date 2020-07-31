setwd("c:/Users/user/Desktop/KCDC/transplantation/HLAtyping/")
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/NGS/")
td <- read.table("HLA_JG_2DGT_imputed.txt",header = T)
fd <- read.table("HLA_JG_4DGT_imputed.txt",header = T)
head(td)

df <- read.csv("HLA_NGS_typing_255samples_results_202002_modify_alltype",header = T)
head(df)               
#write.table(df,"HLA_NGS_sep.tab.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#df <- read.table("HLA_NGS_sep.tab.txt",header = T,na.strings = "0")
colnames(df)
head(df)
library(stringr)
#ngs1<-read.table("NGS.A.gene.txt",header = T)
ngs <-df
head(ngs)
#ngs1$A1 <- str_split_fixed()
head(str_split_fixed(ngs1$A.Allele1,":",2)[,1])
ngs$A1 <- str_split_fixed(str_split_fixed(ngs$A.Allele1,":",2)[,1],"\\*",2)[,2]
ngs$A2 <- str_split_fixed(str_split_fixed(ngs$A.Allele2,":",2)[,1],"\\*",2)[,2]
ngs$B1 <- str_split_fixed(str_split_fixed(ngs$B.Allele1,":",2)[,1],"\\*",2)[,2]
ngs$B2 <- str_split_fixed(str_split_fixed(ngs$B.Allele2,":",2)[,1],"\\*",2)[,2]
ngs$DRB1 <- str_split_fixed(str_split_fixed(ngs$BRB1.Allele1,":",2)[,1],"\\*",2)[,2]
ngs$DRB2 <- str_split_fixed(str_split_fixed(ngs$DRB1.Allele2,":",2)[,1],"\\*",2)[,2]

head(ngs)

ngs.2d <- subset(ngs,select =c("KID","A1","A2","B1","B2","DRB1","DRB2"))

features <- colnames(ngs.2d)[2:ncol(ngs.2d)]
features

for (i in features){
  ngs.2d[,i] <- strtoi(ngs.2d[,i])
}

head(ngs.2d)

td.ngs <- merge(td,ngs.2d,by.x = "ID",by.y = "KID")
head(td.ngs)        

for (i in features){
  for (j in 1:nrow(td.ngs)){
    if (td.ngs[j,paste0(i,".x")] == td.ngs[j,paste0(i,".y")]){
      td.ngs[j,paste(i,".match")] <- 1
    } 
  }
}


i = "A1"
j = 3
td.ngs[j,paste0(i,".x")]
paste0(i,".x")
td.ngs[j,paste0(i,".x")] == td.ngs[j,paste0(i,".y")]
if ((td.ngs$A1.x == td.ngs$A1.y) &&(td.ngs$A2.x == td.ngs$A2.y) ) {
  td.ngs$A1.match <- 2
}
head(td.ngs)
str(td.ngs)
