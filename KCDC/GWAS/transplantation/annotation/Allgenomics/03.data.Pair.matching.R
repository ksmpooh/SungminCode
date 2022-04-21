### score data frame
#ref
#     KBA_ID.x OriID.x    KBA_ID.y OriID.y
#1 NIH19KT5591 KR00033 NIH19KT5659 KD00033
#2 NIH19KT5595 KR00052 NIH19KT5660 KD00052
#3 NIH19KT5598 KR00056 NIH19KT5661 KD00056
#4 NIH19KT5602 KR00058 NIH19KT5664 KD00058


#Rscript
#Three args : [input.txt] [output]
#Rscript --vanilla 03.ploting.R [input.raw] [output.txt]
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2){
	stop("At least two arguments : [input.raw] [output.txt]")
}


df <- read.table(args[1],header=T)
df <- df[,c(1,7:ncol(df))]
df[1:5,1:10]
ref <-read.table("/BDATA/smkim/KR_allogenomics/Allogenomics_KR.KD.QCin.pairTable_ref.txt",header=T)
ref <- ref[,c(1,3)]


out <- merge(ref,df[df$FID %in% ref$KBA_ID.x,],by.x="KBA_ID.x",by.y="FID")
out <- merge(out,df[df$FID %in% ref$KBA_ID.y,],by.x="KBA_ID.y",by.y="FID")
write.table(out,args[2],col.names=T,row.names=F,quote=F,sep="\t")

