setwd("~/Desktop/KCDC/HLAimputation/")

ref <- read.table("NGS/255sampleID.txt")
head(ref)


df <- read.csv("patent/split_result/1M/compare/compare.IMPvsNGS.splitimp.A.B.DRB1.2digit.csv")
df <- read.csv("patent/split_result/1M/compare/Pruning.compare.IMPvsNGS.splitimp.A.B.DRB1.2digit.csv")

size = "1M"
size = "500K"

theme = ""
theme = "Pruning."

for (digit in c("2","4")) {
  df <- read.csv(paste0("patent/split_result/",size,"/compare/",theme,"compare.IMPvsNGS.splitimp.A.B.DRB1.",digit,"digit.csv"))
  out <- df[df$IID %in% ref$V1,]
  write.csv(out,paste0("patent/split_result/254sample/",size,".",theme,"compare.IMPvsNGS.splitimp.A.B.DRB1.",digit,"digit.csv"),row.names = F,quote = F)
}


#================================================================================
setwd("~/Desktop/KCDC/HLAimputation/patent/split_result/254sample/")

size = "1M"
size = "500K"

td <- read.csv(paste0(size,".compare.IMPvsNGS.splitimp.A.B.DRB1.2digit.csv"),header = T)
fd <- read.csv(paste0(size,".compare.IMPvsNGS.splitimp.A.B.DRB1.4digit.csv"),header = T)


td_pruning <- read.csv(paste0(size,".Pruning.compare.IMPvsNGS.splitimp.A.B.DRB1.2digit.csv"),header = T)
fd_pruning <- read.csv(paste0(size,".Pruning.compare.IMPvsNGS.splitimp.A.B.DRB1.4digit.csv"),header = T)
head(td)
head(fd)
head(fd_pruning)
type <- c("td_normal","fd_normal","td_pruning","fd_pruning")
A_match <- c(sum(td$A.match),sum(fd$A.match),sum(td_pruning$A.match),sum(fd_pruning$A.match))
B_match <- c(sum(td$B.match),sum(fd$B.match),sum(td_pruning$B.match),sum(fd_pruning$B.match))
DRB1_match <- c(sum(td$DRB1.match),sum(fd$DRB1.match),sum(td_pruning$DRB1.match),sum(fd_pruning$DRB1.match))

A_wrong <- c(sum(td$A.wrong),sum(fd$A.wrong),sum(td_pruning$A.wrong),sum(fd_pruning$A.wrong))
B_wrong <- c(sum(td$B.wrong),sum(fd$B.wrong),sum(td_pruning$B.wrong),sum(fd_pruning$B.wrong))
DRB1_wrong <- c(sum(td$DRB1.wrong),sum(fd$DRB1.wrong),sum(td_pruning$DRB1.wrong),sum(fd_pruning$DRB1.wrong))

A_empty <- c(sum(td$A.empty),sum(fd$A.empty),sum(td_pruning$A.empty),sum(fd_pruning$A.empty))
B_empty <- c(sum(td$B.empty),sum(fd$B.empty),sum(td_pruning$B.empty),sum(fd_pruning$B.empty))
DRB1_empty <- c(sum(td$DRB1.empty),sum(fd$DRB1.empty),sum(td_pruning$DRB1.empty),sum(fd_pruning$DRB1.empty))

accuracy <- c(sum(td$A.match)/(sum(td$A.match) + sum(td$A.wrong)),
              sum(fd$A.match)/(sum(fd$A.match) + sum(fd$A.wrong)),
              sum(td_pruning$A.match)/(sum(td_pruning$A.match) + sum(td_pruning$A.wrong)),
              sum(fd_pruning$A.match)/(sum(fd_pruning$A.match) + sum(fd_pruning$A.wrong)))

df <- data.frame(type,A_match,A_wrong,A_empty,B_match,B_wrong,B_empty,DRB1_match,DRB1_wrong,DRB1_empty)
df$A_accuracy <- df$A_match/(df$A_match + df$A_wrong)
df$B_accuracy <- df$B_match/(df$B_match + df$B_wrong)
df$DRB1_accuracy <- df$DRB1_match/(df$DRB1_match + df$DRB1_wrong)

df$overall <- (df$A_accuracy + df$B_accuracy +df$DRB1_accuracy)/3


write.table(df,paste0(size,".Accuracy.Result.txt"),col.names = T,row.names = F,quote = F,sep = '\t')

