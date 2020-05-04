library(qqman)

epacts <- read.table("liver.organfailure.epacts.b.firth.result.txt")


###head False
colnames(epacts) <- c("CHROM","BEGIN", "END", "MARKER_ID", "NS", "AC", "CALLRATE", "MAF", "PVALUE", "BETA", "SEBETA", 'CHISQ', 'NS.CASE', 'NS.CTRL', 'AF.CASE', 'AF.CTRL')
epacts_pval_maf <- subset(epacts, select = c(CHROM, BEGIN, END, MARKER_ID, MAF, PVALUE))

#save(epacts_pval_maf, file = "./epact_pval.RData")

epacts_noNApval <- na.omit(epacts_pval_maf)
#save(epacts_noNApval, file="./epact_noNApvale.RData")

png(filename ="./liver.organfailure.epacts.b.firth.result.png" ,width = 1207, height = 753)
#manhattan(epacts_noNApval, main="epacts", chr='CHROM', bp='BEGIN', p='PVALUE', snp='MARKER_ID', suggestiveline=F, genomewideline= -log10(1e-6), ylim=c(0,-log10(1e-9)))
manhattan(epacts_noNApval, main="epacts", chr='CHROM', bp='BEGIN', p='PVALUE', snp='MARKER_ID', suggestiveline=F, genomewideline= -log10(5e-8),ylim=c(0,-log10(1e-9)) )
dev.off()

keep <- c(epacts_noNApval$MAF > 0.05)
table(keep)
epacts_result <- epacts_noNApval[keep,]

keep2 <- c(epacts_noNApval$MAF > 0.1)
table(keep2)
epacts_result2 <- epacts_noNApval[keep2,]

keep3 <- c(epacts_noNApval$MAF > 0.01)
table(keep3)
epacts_result3 <- epacts_noNApval[keep3,]

png(filename ="./liver.organfailure.epacts.b.firth.result.png" ,width = 1207, height = 753)
#manhattan(epacts_result, main="epacts", chr='CHROM', bp='BEGIN', p='PVALUE', snp='MARKER_ID', suggestiveline=F, genomewideline= -log10(1e-6), ylim=c(0,-log10(1e-9)))
manhattan(epacts_result, main="epacts", chr='CHROM', bp='BEGIN', p='PVALUE', snp='MARKER_ID', suggestiveline=F, genomewideline= -log10(5e-8))
#manhattan(epacts_result2, main="epacts", chr='CHROM', bp='BEGIN', p='PVALUE', snp='MARKER_ID', suggestiveline=F, genomewideline= -log10(1e-6))
#manhattan(epacts_result2, main="epacts", chr='CHROM', bp='BEGIN', p='PVALUE', snp='MARKER_ID', suggestiveline=F, genomewideline= -log10(1e-6), ylim=c(0,-log10(1e-8))) + scal
#e_y_continuous(breaks = seq(0,-log10(1e-8),-log10(1e-1)))
#manhattan(epacts_result3, main="epacts", chr='CHROM', bp='BEGIN', p='PVALUE', snp='MARKER_ID', suggestiveline=F, genomewideline= -log10(5e-5))
dev.off()


png(filename ="./Cohort.PDR.info04.Rplot0.1.png" ,width = 1207, height = 753)
#manhattan(epacts_result2, main="epacts", chr='CHROM', bp='BEGIN', p='PVALUE', snp='MARKER_ID', suggestiveline=F, genomewideline= -log10(1e-6), ylim=c(0,-log10(1e-9)))
manhattan(epacts_result2, main="epacts", chr='CHROM', bp='BEGIN', p='PVALUE', snp='MARKER_ID', suggestiveline=F, genomewideline= -log10(5e-8))
#manhattan(epacts_result2, main="epacts", chr='CHROM', bp='BEGIN', p='PVALUE', snp='MARKER_ID', suggestiveline=F, genomewideline= -log10(1e-6))
#manhattan(epacts_result2, main="epacts", chr='CHROM', bp='BEGIN', p='PVALUE', snp='MARKER_ID', suggestiveline=F, genomewideline= -log10(1e-6), ylim=c(0,-log10(1e-8))) + scal
#e_y_continuous(breaks = seq(0,-log10(1e-8),-log10(1e-1)))
#manhattan(epacts_result3, main="epacts", chr='CHROM', bp='BEGIN', p='PVALUE', snp='MARKER_ID', suggestiveline=F, genomewideline= -log10(5e-5))
dev.off()

png(filename ="./Cohort.PDR.info04.Rplot0.01.png" ,width = 1207, height = 753)
manhattan(epacts_result3, main="epacts", chr='CHROM', bp='BEGIN', p='PVALUE', snp='MARKER_ID', suggestiveline=F, genomewideline= -log10(5e-8))
dev.off()

