

library(dplyr)
setwd("C:/Users/user/Desktop/cluster_matrix/")
cluster_matrix <- read.table("df.txt",header = T,sep = "\t",na.strings = c("","NA"))
cluster_matrix$na_count = apply(cluster_matrix, 1, function(x) sum(is.na(x)))
numberOfbatch <- ncol(cluster_matrix)-2
cluster_matrix$Batch_failure <- (numberOfbatch - cluster_matrix$na_count)/numberOfbatch *100


cluster_matrix$Other <- apply(cluster_matrix,1,function(x) length(which(x == "Other")))
cluster_matrix$OTV <- apply(cluster_matrix,1,function(x) length(which(x == "OTV")))
cluster_matrix$CallRateBelowThreshold <- apply(cluster_matrix,1,function(x) length(which(x == "CallRateBelowThreshold")))



##############################probe count
marker_list <- cluster_matrix$probeset_id

v1_probeset <- read.table("Axiom_KORV1.0.na34.annot.probe_count_change_v1.1ID.txt",header = T)

length(intersect(v1_probeset$Probe_Set_ID,marker_list))
#v1_probeset <- read.table("Axiom_KORV1.0.na34.annot.probe_count.txt",header= T)
rownames(v1_probeset) <- v1_probeset$Probe_Set_ID
v1_probeset<-v1_probeset[marker_list,]


v2_probeset <-read.table("Axiom_KORV1.1.na35.annot.probe_count.txt",header = T)
length(intersect(v2_probeset$Probe_Set_ID,marker_list))
rownames(v2_probeset) <- v2_probeset$Probe_Set_ID
v2_probeset<-v2_probeset[marker_list,]


