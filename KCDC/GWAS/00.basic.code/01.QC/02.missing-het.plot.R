###Rsciprt
#Three args : [ref] [input] [output]
#Three args : [missing.imiss] [het.het] [PDF.output]
# Rscript --vanilla input.R [missing.imiss] [het.het] [PDF.output]
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3){
	stop("At least three arguments : [missing.imiss] [het.het] [PDF.output] ")
}
#print("Read : ref")
#ref <- read.table(args[1],header = F)

print("Read : missing.imiss")
miss <-read.table(args[1],header = T)
print("Read : het.het")
het <- read.table(args[2], header = T)


miss <- cbind(miss, CR=((1 - miss$F_MISS)*100))
het <- cbind(het, HET=((het$N.NM. - het$O.HOM.)/het$N.NM.)*100)

lowSample <- merge(miss, het, by="FID")

pdf(args[3], height = 7, width = 10)
plot(lowSample$HET, lowSample$F_MISS, 
    #xlim=c(13,22), ylim=c(0,0.1), 
    xlab="heterozygosity rate",ylab="missing rate", 
    main="Missing vs heterozygosity", col=rgb(0,0,1,0.3), cex=1.5, pch=16)
abline(v=15.5, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(v=17, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(h=0.03, col=rgb(1,0,0,1), lty=3, lwd=2)
points(lowSample[lowSample$HET < 15.5 | 17 < lowSample$HET | 0.03 < lowSample$F_MISS,]$HET,
       lowSample[lowSample$HET < 15.5 | 17 < lowSample$HET | 0.03 < lowSample$F_MISS,]$F_MISS,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
dev.off()

rmList <- lowSample[0.03 < lowSample$F_MISS | lowSample$HET < 15.5 | 17 < lowSample$HET,]
#dim(rmList)

write.table(rmList[,c(1:2)], "rmLQSamples.txt", col.names= FALSE, row.names=FALSE, sep="\t", quote=FALSE)


