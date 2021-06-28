### HLa typing °á°ú ºóµµ

setwd("c:/Users/user/Desktop/KCDC/HLAimputation/HLAtyping/all/")

df <- read.csv("HLAtyping.alle.gene.4digit_2019.with.2020.csv")
ref <- read.csv("../../all/impute4.Han/compare.IMPvsNGS.all.gene.2digit.csv")
head(ref)
head(df)
df <- df[df$KID %in% ref$IID,]
head(df)

head(df)
table(df$NGS_A.1) + table(df$NGS_A.2)


a = df[,c('KID','NGS_A.1')]
colnames(a) <- c("ID","A")
b = df[,c('KID','NGS_A.2')]
colnames(b) <- c("ID","A")

A <- rbind(a,b)
A <- as.data.frame(table(A$A))
A$pec <- A$Freq /1022 * 100
A <- A[order(A$Freq,decreasing = TRUE),]
head(A)
hist(A$A)

for (gene in c("A","B","C","DRB1","DPA1","DPB1","DQA1","DQB1")) {
  print(gene)
  a <- df[,c('KID',paste0("NGS_",gene,".1"))]
  colnames(a) <- c("ID",'A')
  b <- df[,c('KID',paste0("NGS_",gene,".2"))]
  colnames(b) <- c("ID",'A')
  A <- rbind(a,b)
  A <- as.data.frame(table(A$A))
  A$pec <- A$Freq /sum(A$Freq) * 100
  A <- A[order(A$Freq,decreasing = TRUE),]
  
  write.csv(A,paste0("freq/HLA.",gene,"_4digit_freq.csv"),row.names = F,quote = F)
  
}

df = data.frame()
for (gene in c("A","B","C","DRB1","DPA1","DPB1","DQA1","DQB1")) {
  a = read.csv(paste0("freq/HLA.",gene,"_4digit_freq.csv"))
  a$HLA_gene <- gene
  df = rbind(df,a)
}


head(df)
write.csv(df,'freq/HLA.all.freq.csv',row.names = F)

ggplot(data=df,aes(y=pec,x=HLA_gene,fill=Var1)) +
  geom_bar("identity") +
  theme_minimal()


ggplot(df,aes(y=pec,x=HLA_gene)) +
  geom_bar("identity") +
  theme_minimal()


write.csv(df,"freq/HLA.all.freq.csv")
sum(df[df$HLA_gene == 'DPA1',]$pec)
