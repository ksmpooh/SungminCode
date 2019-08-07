setwd("C:/Users/user/Desktop/KCDC/Gastric/")
cell_list <-read.table("cel_file_list.txt",header = T)
rm_list <- read.table("merge_rmSample_1stQC.txt")

cell_list <- data.frame(do.call("rbind",strsplit(as.character(cell_list$cel_files),'_',fixed = T)))

c2 <-data.frame(cell_list$X6)
c3 <-data.frame(do.call("rbind",strsplit(as.character(c2$cell_list.X6),'.',fixed = T)))

rownames(cell_list)<-as.factor(c3$X1)
rownames(rm_list)<-rm_list$V1

#library(dplyr)
rownames(rm_list)
rownames(cell_list)
rmID <- row.names(rm_list)

rows_to_remove <- which(row.names(cell_list) %in% rmID)

Final_list <- cell_list_frame[-rows_to_remove,]

######################3
cell_list <-read.table("cel_file_list.txt",header = T)
rm_list <- read.table("merge_rmSample_1stQC.txt")

cell_list_frame <- data.frame(do.call("rbind",strsplit(as.character(cell_list$cel_files),'_',fixed = T)))

c2 <-data.frame(cell_list_frame$X6)
c3 <-data.frame(do.call("rbind",strsplit(as.character(c2$cell_list_frame.X6),'.',fixed = T)))

rownames(cell_list)<-as.factor(c3$X1)
rownames(rm_list)<-rm_list$V1

#library(dplyr)
rownames(rm_list)
rownames(cell_list)
rmID <- row.names(rm_list)

rows_to_remove <- which(row.names(cell_list) %in% rmID)

Final_list <- as.data.frame(cell_list[-rows_to_remove,])
colnames(Final_list) <-"cel_files"
#Final_list <-gsub("\"\"","",Final_list$cel_files)

write.table(Final_list,"2nd_QC_cel_file_list.txt",col.names = T,row.names = F,quote = F)
