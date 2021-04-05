plink --bfile JG.2020.rep_snpolisher --missing --out MISS
plink --bfile JG.2020.rep_snpolisher --het --out HET



##missing-het plot

miss <-read.table("MISS.imiss",header = T)
het <- read.table("HET.het", header = T)
miss <- cbind(miss, CR=((1 - miss$F_MISS)*100))
het <- cbind(het, HET=((het$N.NM. - het$O.HOM.)/het$N.NM.)*100)

lowSample <- merge(miss, het, by="FID")

pdf("../PDF/JG.2nd.QC_SNPolisher_miss-het.pdf", height = 7, width = 10)
plot(lowSample$HET, lowSample$F_MISS, xlim=c(10,25), ylim=c(0,0.1), xlab="heterozygosity rate",
     ylab="missing rate", main="Missing vs heterozygosity", col=rgb(0,0,1,0.3), cex=1.5, pch=16)
abline(v=15, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(v=17, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(h=0.03, col=rgb(1,0,0,1), lty=3, lwd=2)
points(lowSample[lowSample$HET < 15 | 17 < lowSample$HET | 0.03 < lowSample$F_MISS,]$HET,
       lowSample[lowSample$HET < 15 | 17 < lowSample$HET | 0.03 < lowSample$F_MISS,]$F_MISS,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
dev.off()

rmList <- lowSample[0.03 < lowSample$F_MISS | lowSample$HET < 15 | 17 < lowSample$HET,]
#dim(rmList)

#write.table(rmList[,c(1:2)], "rmLQSamples.txt", col.names= FALSE, row.names=FALSE, sep="\t", quote=FALSE)