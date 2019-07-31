setwd("C:/Users/user/Desktop/cluster_matrix/")
getwd()

file_list <- read.table("list_window.txt")
file_list <-unlist(file_list)
  

for (f in file_list[1:2]){
  #assign(paste0("model_",f), read.table(gzfile(paste0(f,"/Ps.performance.uniquevar.txt.gz"),"rt"),header = T ))
  #assign(f, read.table(gzfile(paste0(f,"/Ps.performance.uniquevar.txt.gz"),"rt"),header = T ))
  assign(toy, read.table(gzfile(paste0(f,"/Ps.performance.uniquevar.txt.gz"),"rt"),header = T ))
  assign(toy,toy[sample(nrow(toy),10000),])
  write.csv(toy,paste0("Toyset/Toy_",f,".csv"),col.names = TRUE,row.names = FALSE,quote = FALSE)
}

Toy_AS_V1_B1_2nd <- AS_V1_B1_2nd[sample(nrow(AS_V1_B1_2nd),10000),]
Toy_AS_V1_B2_2nd <- AS_V1_B2_2nd[sample(nrow(AS_V1_B2_2nd),10000),]

#########################
zz = gzfile('AS_V1_B1_2nd/Ps.performance.uniquevar.txt.gz',"rt")
a <- read.table(zz,header = T)
close()
head(a)
??close()

