setwd("C:/Users/user/Desktop/cluster_matrix/")
getwd()

file_list <- read.table("list_window.txt")
file_list <-unlist(file_list)

##############read file
for (f in file_list[1:2]){
  toy <- data.frame()
  #assign(paste0("model_",f), read.table(gzfile(paste0(f,"/Ps.performance.uniquevar.txt.gz"),"rt"),header = T ))
  #assign(f, read.table(gzfile(paste0(f,"/Ps.performance.uniquevar.txt.gz"),"rt"),header = T ))
  #assign(toy, read.table(gzfile(paste0(f,"/Ps.performance.uniquevar.txt.gz"),"rt"),header = T ))
  toy <- read.table(gzfile(paste0(f,"/Ps.performance.uniquevar.txt.gz"),"rt"),header = T )
  #assign(toy,toy[sample(nrow(toy),10000),])
  toy <- toy[sample(nrow(toy),10000),]
  write.table(toy,paste0("Toyset/Toy_",f,".txt"),col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
  on.exit(close(toy))
}


for (f in file_list[1:2]){
  assign(paste0(f),read.table(paste0("Toyset/Toy_",f,".txt"),header = T,sep = "\t"))
}

length(intersect(AS_V1_B1_2nd$probeset_id,AS_V1_B2_2nd$probeset_id))
a <- subset(AS_V1_B1_2nd,select = c(probeset_id,ConversionType))
colnames(a)
