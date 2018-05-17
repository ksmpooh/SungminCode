#TCGA<-read.csv("C:/Users/sungmin/Downloads/gene_expression_RNAseq_RSEM_fpkm#1.csv",sep = ',', header = T)
TCGA<-read.csv("/home/tjahn/BioDataAnalysis/BioMarker/phenotype/gene_expression_RNAseq_RSEM_fpkm#1.csv",sep = ',', header = T)

gene_name <- TCGA$sample
tTCGA<-t(TCGA)
tTCGA<-tTCGA[2:nrow(tTCGA),]
tTCGA<-as.data.frame(tTCGA)
colnames(tTCGA)<-gene_name

#CESC<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/CESC phenotype.csv",sep = ',', header = T)
#COAD<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/COAD phenotype.csv",sep = ',', header = T)
#PAAD<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/PAAD phenotype.csv",sep = ',', header = T)
#STAD<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/STAD phenotype.csv",sep = ',', header = T)
#UCEC<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/UCEC phenotype.csv",sep = ',', header = T)
#UCS<-read.csv("C:/Users/sungmin/Downloads/drive-download-20180517T041248Z-001/UCS phenotype.csv",sep = ',', header = T)
#filelist<-c(CESC,COAD,PAAD,STAD,UCEC,UCS)

CESC<-read.csv("/home/tjahn/BioDataAnalysis/BioMarker/CESC phenotype.csv",sep = ',', header = T)
COAD<-read.csv("/home/tjahn/BioDataAnalysis/BioMarker/COAD phenotype.csv",sep = ',', header = T)
PAAD<-read.csv("/home/tjahn/BioDataAnalysis/BioMarker/PAAD phenotype.csv",sep = ',', header = T)
STAD<-read.csv("/home/tjahn/BioDataAnalysis/BioMarker/STAD phenotype.csv",sep = ',', header = T)
UCEC<-read.csv("/home/tjahn/BioDataAnalysis/BioMarker/UCEC phenotype.csv",sep = ',', header = T)
UCS<-read.csv("/home/tjahn/BioDataAnalysis/BioMarker/UCS phenotype.csv",sep = ',', header = T)


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

tTCGA[CESC_names,"Primary_site"] <- CESC$X_primary_site
a<-tTCGA[CESC_names,]
b<-tTCGA[COAD_names,]
c<-tTCGA[PAAD_names,]
d<-tTCGA[STAD_names,]
e<-tTCGA[UCEC_names,]
f<-tTCGA[UCS_names,]
a[,"Primary_site"]<-CESC$X_primary_site
b[,"Primary_site"]<-COAD$X_primary_site
c[,"Primary_site"]<-PAAD$X_primary_site
d[,"Primary_site"]<-STAD$X_primary_site
e[,"Primary_site"]<-UCEC$X_primary_site
f[,"Primary_site"]<-UCS$X_primary_site

#list <- c("b","c","d","e","f")
a<-rbind(a,b)
a<-rbind(a,c)
a<-rbind(a,d)
a<-rbind(a,e)
a<-rbind(a,f)

write.csv(a,"/home/tjahn/BioDataAnalysis/BioMarker/result.csv")