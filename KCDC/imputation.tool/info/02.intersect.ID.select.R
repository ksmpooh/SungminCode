###Rsciprt
#Three args : [ref] [input] [output]
# Rscript --vanilla
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3){
	stop("At least three arguments : [ref] [input1] [output] ")
}

#setwd("/BDATA/smkim/JG/08.asso/OUTPUTs/plot/")

#library(qqman)

print("Read : ref")
ref <- read.table(args[1],header = F)
#ref <- ref[,c(2,5)]

print("Read : df")
df <- read.table(args[2],header = T)
#rownames(df) <- df$V1


print("Read :Merge")
out <- df[df$ID %in% ref$V1,]

colnames(out) <- c("ID","INFO","MAF")

write.table(out,args[3],col.names = T,row.names=F,quote =F,sep="\t")