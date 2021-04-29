#input
#10:12253597_G_A
#10:71091013_A_G
#10:71099913_T_C
#10:114752503_T_C
#10:114758349_C_T
#11:199256_G_A


#ref
#X:60014_T_C	rs370048753
#X:60014_T_G	rs370048753


#output

###Rsciprt
#Three args : [ref] [input] [output]
# Rscript --vanilla
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3){
	stop("At least three arguments : [ref] [input] [output] ")
}


print("Read : ref")
ref <- read.table(args[1],header = F)
#ref <- ref[,c(2,5)]

print("Read : df")
df <- read.table(args[2],header = F)
#rownames(df) <- df$V1


print("Read :Merge")
#out <- df[df$V1 %in% ref$V1,]
out <- merge(df,ref,by = "V1",all.x = T)

#colnames(out) <- c("ID","INFO","MAF")
out <- na.omit(out)
write.table(out,args[3],col.names = F,row.names=F,quote =F,sep="\t")
#write.table(out[,2],args[3],col.names = F,row.names=F,quote =F,sep="\t")