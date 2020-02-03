setwd("C:/Users/user/Desktop/KCDC/transplantation/")

type <- read.csv("sample_info/KchipJG_pairTable_20190908.csv",header = T)
t2 <- na.omit(type)
head(t2)
dim(t2)

t3 <- subset(type,type$Match == "NO")
dim(t3)

liver <- t2[grep("^LR",t2$OriID),]
kidney <- t2[grep("^KR",t2$OriID),]
dim(liver)
dim(kidney)

dim(liver)[1] + dim(kidney)[1]

head(type)
dim(type)
df <- read.table("sample_info/sample.info.with.type.txt",header = T)
head(df)
dim(df)

dim(df)[1] - dim(liver)[1]*2 - dim(kidney)[1] * 2

table(df$type)


ori <- read.table("summary.info.txt")
dim(ori)
head(ori)
ori <- ori[,c("V1","V5")]
colnames(ori) <- c("NewID","sex")

a <- merge(ori,df,by = "NewID")
head(a)
dim(a)
table(a$type)
#########################2Â÷ QC
second <- read.table("2nd/JG.2nd.QC_snpolisher_miss-het.het",header = T)
second <- subset(second,select = c("FID"))
df <- read.table("sample_info/sample.info.with.type.txt",header = T)
dim(df)
head(df)
head(second)
df <- merge(df,second,by.x = 'NewID',by.y = 'FID')
head(df)
dim(df)
table(df$type)
pca <- read.table("ethical/case_ID.txt")
head(pca)
dim(pca)
pca <- merge(pca,df,by.x = 'V1',by.y = 'NewID')
table(pca$type)
