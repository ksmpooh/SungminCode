library(dplyr)
library(stringr)
setwd("C:/Users/user/Desktop/cluster_matrix/")
getwd()

file_list <- read.table("list_window.txt")
file_list <-unlist(file_list)

##############read file
# 한번에 주석 처리할때는 ctrl + shift + c
# 
# for (f in file_list){
#   toy<-data.frame()
#   toy <- read.table(gzfile(paste0(f,"/Ps.performance.uniquevar.txt.gz"),"rt"),header = T )
#   #assign(paste0("model_",f),subset(paste0("model_",f),select = c("probeset_id","ConversionType","H.W.p.Value","MinorAlleleFrequency")))
#   #assign(paste0("model_",f),subset(paste0("model_",f),select = c("probeset_id","ConversionType")))
#   toy <- subset(toy,select = c(probeset_id,ConversionType))
#   write.table(toy,paste0("new_",f,".txt"),quote = F,row.names = F,col.names = T)
#   on.exit(close(toy))
# }
toy <- data.frame()
# for (f in file_list[1]){
#   assign(paste0(f),read.table(paste0("new_set/new_",f,".txt"),header = T))
# }
for (f in file_list){
  toy <- read.table(paste0("new_set/v1tov2/new_",f,".txt"),header = T)
  toy <- subset(toy,(ConversionType == "Other" | ConversionType == "OTV" | ConversionType =="callRateBelowThreshold"))
  write.table(toy,paste0("new_set/only_error_set/only_error_",f,".txt"),quote = F, row.names = F, col.names = T)
}


#####################
venn <-data.frame()
for (f in file_list){
  assign(paste0(f),read.table(paste0("new_set/only_error_set/only_error_",f,".txt"),header = T))
}



