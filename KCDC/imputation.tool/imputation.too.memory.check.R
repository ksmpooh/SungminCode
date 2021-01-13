setwd("c:/Users/user/Desktop/KCDC/imputation.tool/")


#################TEST
loglist <- read.table("memory.check/impute4/log.list.txt",header = T)
head(loglist)

df <- read.table("memory.check/impute4/100000001_105000000.log",header = T)
head(df)
table(df$COMMAND)
process_list <- data.frame(table(df$COMMAND))$Var1

out <- matrix(nrow = 1,ncol = 5)
out <- as.data.frame(out)
colnames(out) <-c("COMMAND","X.CPU","X.MEM","VSZ","RSS")


for (i in process_list) {
  ref = df[df$COMMAND == i,]    
  ref = ref[ref$RSS == max(ref$RSS),c("COMMAND","X.CPU","X.MEM","VSZ","RSS")]
  ref = ref[1,]
  #print(ref)
#  print(max(df[df$COMMAND == i,]$RSS))
  out = rbind(out,ref)
}
out
ref

##########g20201230 all
setwd("c:/Users/user/Desktop/KCDC/imputation.tool/")
loglist <- read.table("memory.check/impute4/log.list.txt",header = T)

out <- matrix(nrow = 1,ncol = 6)
out <- as.data.frame(out)
colnames(out) <-c("chunk","COMMAND","X.CPU","X.MEM","VSZ","RSS")


for (i in loglist$log) {
  df <- read.table(paste0("memory.check/impute4/",i,".log"),header = T)
  process_list <- data.frame(table(df$COMMAND))$Var1
  for (j in process_list) {
    ref = df[df$COMMAND == j,]
    ref = ref[ref$RSS == max(ref$RSS),c("COMMAND","X.CPU","X.MEM","VSZ","RSS")]
    ref = ref[1,]
    ref$chunk = i
    out = rbind(out,ref)
  }
}
head(out)
row.names(out) <- NULL
out <- out[2:nrow(out),]









