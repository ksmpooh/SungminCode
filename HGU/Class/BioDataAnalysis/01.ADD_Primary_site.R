###1. 파일 읽어 오기
#저장할 변수 이름 <-read.csv("파일경로",sep='구분자',header = colname 유무)
TCGA<-read.csv("C:/Users/sungmin/Downloads/gene_expression_RNAseq_RSEM_fpkm#1.csv",sep = ',', header = T)
#TCGA<-read.csv("/home/tjahn/BioDataAnalysis/BioMarker/phenotype/gene_expression_RNAseq_RSEM_fpkm#1.csv",sep = ',', header = T)

CESC<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/CESC phenotype.csv",sep = ',', header = T)
COAD<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/COAD phenotype.csv",sep = ',', header = T)
PAAD<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/PAAD phenotype.csv",sep = ',', header = T)
STAD<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/STAD phenotype.csv",sep = ',', header = T)
UCEC<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/UCEC phenotype.csv",sep = ',', header = T)
UCS<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/UCS phenotype.csv",sep = ',', header = T)

#CESC<-read.csv("/home/tjahn/BioDataAnalysis/BioMarker/CESC phenotype.csv",sep = ',', header = T)
#COAD<-read.csv("/home/tjahn/BioDataAnalysis/BioMarker/COAD phenotype.csv",sep = ',', header = T)
#PAAD<-read.csv("/home/tjahn/BioDataAnalysis/BioMarker/PAAD phenotype.csv",sep = ',', header = T)
#STAD<-read.csv("/home/tjahn/BioDataAnalysis/BioMarker/STAD phenotype.csv",sep = ',', header = T)
#UCEC<-read.csv("/home/tjahn/BioDataAnalysis/BioMarker/UCEC phenotype.csv",sep = ',', header = T)
#UCS<-read.csv("/home/tjahn/BioDataAnalysis/BioMarker/UCS phenotype.csv",sep = ',', header = T)

#2. TCAG파일(gene 다 합친거, 정리 작업)
gene_name <- TCGA$sample #sample이란 column에 gene이름이 저장되어 있는 것을 따로 gene_name변수를 만들어 저장
tTCGA<-t(TCGA)   #original file이 row에는 gene, col에는 사람으로 되어 있어서 그것을 switch하는 작업, 
                 #row에 주로 object 즉 사람과 같은 객체를 저장, col에는 변수, 즉 특징같은 것을 저장 작업함
tTCGA<-tTCGA[2:nrow(tTCGA),] #1번 row는 sample이므로 그것을 빼고 저장
tTCGA<-as.data.frame(tTCGA)  #data.frame형식으로 변환
colnames(tTCGA)<-gene_name   #변수(colname)을 이전에 저장한 gene_name으로 변환

#CESC<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/CESC phenotype.csv",sep = ',', header = T)
#COAD<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/COAD phenotype.csv",sep = ',', header = T)
#PAAD<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/PAAD phenotype.csv",sep = ',', header = T)
#STAD<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/STAD phenotype.csv",sep = ',', header = T)
#UCEC<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/UCEC phenotype.csv",sep = ',', header = T)
#UCS<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/UCS phenotype.csv",sep = ',', header = T)
#filelist<-c(CESC,COAD,PAAD,STAD,UCEC,UCS)


###3.파일별 TCGA이름 저장 및 rowname을 숫자(1,2,3,4,...)에서 TCGA(TCGA_a10 /....)이름으로 변환
rownames(CESC)<-CESC$X_INTEGRATION
rownames(COAD)<-COAD$X_INTEGRATION
rownames(PAAD)<-PAAD$X_INTEGRATION
rownames(STAD)<-STAD$X_INTEGRATION
rownames(UCEC)<-UCEC$X_INTEGRATION
rownames(UCS)<-UCS$X_INTEGRATION

CESC_names<-CESC$X_INTEGRATION
COAD_names<-COAD$X_INTEGRATION
PAAD_names<-PAAD$X_INTEGRATION
STAD_names<-STAD$X_INTEGRATION
UCEC_names<-UCEC$X_INTEGRATION
UCS_names<-UCS$X_INTEGRATION


####4. tTCGA file에서 각 파일별 name을 row에서 불어와서 그것과 매칭되는 object 값을 따로 저장
#tTCGA[CESC_names,"Primary_site"] <- CESC$X_primary_site
a<-tTCGA[CESC_names,]
b<-tTCGA[COAD_names,]
c<-tTCGA[PAAD_names,]
d<-tTCGA[STAD_names,]
e<-tTCGA[UCEC_names,]
f<-tTCGA[UCS_names,]

####5. 따로 저장한 것에 "Primary_site"라는 variable를 만들어 col에 저장하고, 각 파일에서 얻는 site 값을 새로 만든 col에 저장 
a[,"Primary_site"]<-CESC$X_primary_site
b[,"Primary_site"]<-COAD$X_primary_site
c[,"Primary_site"]<-PAAD$X_primary_site
d[,"Primary_site"]<-STAD$X_primary_site
e[,"Primary_site"]<-UCEC$X_primary_site
f[,"Primary_site"]<-UCS$X_primary_site

#### 6. 나눈 파일들을 다시 하나로 만들기
#list <- c("b","c","d","e","f")
a<-rbind(a,b)
a<-rbind(a,c)
a<-rbind(a,d)
a<-rbind(a,e)
a<-rbind(a,f)


#### 7. 결과 저장 write.csv(저장된 변수이름,"저장할 파일 경로 및 이름")
write.csv(a,"/home/tjahn/BioDataAnalysis/BioMarker/result.csv",row.names = T)

CESC_names
COAD_names
PAAD_names
STAD_names
UCEC_names
UCS_names

a <- as.factor(CESC_names)
b <- as.factor(COAD_names)
c <- as.factor(PAAD_names)
d <- STAD_names
e<- UCEC_names
f<-UCS_names
name<-colnames(TCGA)
name<-gsub("\\.","\\-",name[2:1535])
intersect(name,CESC_names)
