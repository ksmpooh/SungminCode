

library(dplyr)

cluster_matrix <- read.table("C:/Users/user/Desktop/cluster_matrix/df.txt",header = T,sep = "\t",na.strings = c("","NA"))
cluster_matrix$na_count = apply(cluster_matrix, 1, function(x) sum(is.na(x)))
numberOfbatch <- ncol(cluster_matrix)-2
cluster_matrix$Batch_failure <- (numberOfbatch - cluster_matrix$na_count)/numberOfbatch *100


cluster_matrix$Other <- apply(cluster_matrix,1,function(x) length(which(x == "Other")))
cluster_matrix$OTV <- apply(cluster_matrix,1,function(x) length(which(x == "OTV")))
cluster_matrix$CallRateBelowThreshold <- apply(cluster_matrix,1,function(x) length(which(x == "CallRateBelowThreshold")))




#hi
a <- 1
a
