###Rsciprt
#Three args : [input.txt] [output]

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2){
	stop("At least two arguments : [input.txt] [output]")
}

setwd("/BDATA/myhwang/KBA_130K/11_UKB/RESULTs/assoMERGE_new/merge")

library(qqman)

df <- read.table(args[1],header = F)
colnames(df) <- c("ID","CHROM","POS","P")

png(filename = paste0(args[2],".png") ,width = 1207, height = 753)
manhattan(df, main=args[2], chr='CHROM', bp='POS', p='P', snp='ID', suggestiveline=T, genomewideline= -log10(5e-8))
#manhattan(df, main=args[2], chr='CHROM', bp='POS', p='P', snp='ID', suggestiveline=T, genomewideline= -log10(5e-8),ylim=c(0,-log10(1e-10)))

#manhattan(epacts_noNApval, main="epacts", chr='CHROM', bp='BEGIN', p='PVALUE', snp='MARKER_ID', suggestiveline=F, genomewideline= -log10(5e-8))
dev.off()