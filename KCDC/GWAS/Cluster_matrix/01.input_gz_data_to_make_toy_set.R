setwd("C:/Users/user/Desktop/cluster_matrix/")
getwd()

file_list <- read.table("list_window.txt")
file_list <-unlist(file_list)

##############read file


for (f in file_list){
  toy<-data.frame()
  toy <- read.table(gzfile(paste0(f,"/Ps.performance.uniquevar.txt.gz"),"rt"),header = T )
  #assign(paste0("model_",f),subset(paste0("model_",f),select = c("probeset_id","ConversionType","H.W.p.Value","MinorAlleleFrequency")))
  #assign(paste0("model_",f),subset(paste0("model_",f),select = c("probeset_id","ConversionType")))
  toy <- subset(toy,select = c(probeset_id,ConversionType))
  write.table(toy,paste0("new_",f,".txt"),quote = F,row.names = F,col.names = T)
  on.exit(close(toy))
}

for (f in file_list[1]){
  assign(paste0(f),read.table(paste0("new_set/new_",f,".txt"),header = T))
}

for (f in file_list[1]){
  assign(paste0(f),subset(paste0(f),(ConversionType == "Other" | ConversionType == "OTV" | ConversionType =="callRateBelowThreshold")))
}
a<- subset(AS_V1_B1_2nd,(ConversionType == "Other" | ConversionType == "OTV" | ConversionType =="callRateBelowThreshold"))

c <- merge(AS_V1_B1_2nd,AS_V1_B2_2nd,by = "probeset_id")
length(intersect(AS_V1_B1_2nd$probeset_id,AS_V1_B2_2nd$probeset_id))

colnames(model_AS_V1_B1_2nd)

##############read file
for (f in file_list[1:1]){
  toy <- data.frame()
  #assign(paste0("model_",f), read.table(gzfile(paste0(f,"/Ps.performance.uniquevar.txt.gz"),"rt"),header = T ))
  
  
  toy <- read.table(gzfile(paste0(f,"/Ps.performance.uniquevar.txt.gz"),"rt"),header = T )
  #assign(toy,toy[sample(nrow(toy),10000),])
  #toy <- toy[sample(nrow(toy),10000),]
  #toy <- toy[1:10000,]
  #write.table(toy,paste0("Toyset/Toy_",f,".txt"),col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
  on.exit(close(toy))
  #on.exit(close(paste0("model_",f)))
}


for (f in file_list[1:4]){
  assign(paste0(f),read.table(paste0("Toyset/Toy_",f,".txt"),header = T,sep = "\t"))
}

levels(model_DS_V1_B2_2nd$ConversionType)
 
#rownames(model_DS_V1_B2_2nd)
Cluster_matrix <- data.frame(row.names = rownames(model_DS_V1_B2_2nd))
#Cluster_matrix <-cbind(Cluster_matrix,model_DS_V1_B2_2nd$ConversionType)

#50% >
length(which(model_DS_V1_B2_2nd$ConversionType =='Other'))
length(which(model_DS_V1_B2_2nd$ConversionType =='OTV'))
length(which(model_DS_V1_B2_2nd$ConversionType =='CallRateBelowThreshold'))
length(which(model_DS_V1_B2_2nd$ConversionType =='Other' || 'CallRateBelowThreshold' || 'OTV'))

a <- table(model_DS_V1_B2_2nd$ConversionType)
sum(a[c('Other' , 'OTV' , 'CallRateBelowThreshold')])
#######################
name1<-c("kim","lee","ha")
v <- c("1","2","3")
a <- data.frame(v,name1)
rownames(a)<-name1
v1 <- c("1","3","4")
name<-c("kim","ha","jui")
b<-data.frame(name,v1)
rownames(b)<-name
c <- merge(a,b,by = 0,all=T)



a<- subset(AS_V1_B1_2nd,(AS_V1_B1_2nd$ConversionType == "Other" |AS_V1_B1_2nd$ConversionType == "OTV" |AS_V1_B1_2nd$ConversionType =="callRateBelowThreshold"))
b<=