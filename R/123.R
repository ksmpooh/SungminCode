rm(list = ls(all = TRUE));
	
	## Settings ##
setwd("C:/Users/Administrator/Documents"); ## set your working directory ##
setwd("~/Downloads/eLD.v1.0/"); ## set your working directory ##
head <- "ExampleHLA.v1"; ## prefix of the input HLA allele data (in the SNP2HLA imputation reference format) ##
Thres <- 0.05; ## MAF threshold for combining the rare alleles ##
Iter <- 1000; ## No. iterations to estimate null distributin of epsilon ##

Data <-read.table(paste(head, ".bgl.phased", sep=""), header=F, sep=" ");prefix <- paste(head, ".bgl.phased.eps", sep="");
Fix <- (length(Data[1,])-2)/2;
HLA <- c("HLA_A_", "HLA_C_", "HLA_B_", "HLA_DRB1_", "HLA_DQA1_", "HLA_DQB1_", "HLA_DPA1_", "HLA_DPB1_");numHLA <- length(HLA);
	
PermOrd<- function(Fix) {
		flagR <- rep(F, Fix*2);
		Seq <- 1:Fix;
		ord <-  order(runif(Fix));
		flagR[Seq*2-1] <- ord*2-1;
		flagR[Seq*2] <- ord*2;
		return(flagR)
}
	
outfile <- paste(prefix, "_FreqThres", Thres, "_N", Fix, "_Iter", Iter, ".txt", sep="");
outfile1 <- paste(prefix, "_FreqThres", Thres, "_N", Fix, "_Iter", Iter, ".Matrix.txt", sep="");
outfile2 <- paste(prefix, "_FreqThres", Thres, "_N", Fix, "_Iter", Iter, "_Haplo.txt", sep="");

Geno <- Data[6:length(Data[,1]),3:length(Data[1,])]=="P";Geno <- Geno[, PermOrd(length(Geno[1,])/2)[1:(Fix*2)]];numSample <- length(Geno[1,]);
head(Geno)
Allele <- Data[6:length(Data[,1]),2];N <- length(Allele);
head(Allele)
Freq <- rowSums(Geno==T)/length(Geno[1,]);
head(Freq)

flagHLA <-matrix(rep(F, length(Allele)*numHLA), ncol=numHLA);
flagFreq <- Freq>=Thres;
for (i in 1:numHLA) {
	flagHLA[substr(Allele, 1, nchar(HLA[i]))==HLA[i],i] <- T;
}

NQC <- sum(flagFreq)+numHLA;
GenoQC <- matrix(rep(F, NQC*numSample), ncol=numSample);
AlleleQC <- rep(0, NQC);
flagHLAQC <-matrix(rep(F, length(AlleleQC)*numHLA), ncol=numHLA);

GenoQC[1:sum(flagFreq), ] <- Geno[flagFreq,]
AlleleQC[1:sum(flagFreq)] <- as.vector(Allele[flagFreq]);

counter <- sum(flagFreq);
for (i in 1:numHLA) {
	flagT <- flagHLA[,i] & !flagFreq;
	if (sum(flagT)>=2) {
		GenoQC[counter+i,] <- apply(Geno[flagT,], 2, any);
		AlleleQC[counter+i] <- paste(HLA[i], "rare", sep="");
	} else if (sum(flagT)==1) {
		GenoQC[counter+i,] <- Geno[flagT,];
		AlleleQC[counter+i] <- paste(HLA[i], "rare", sep="");
	} else {
		AlleleQC[counter+i] <- "-";
	}
	
	if (sum(GenoQC[counter+i,])==0) {
		AlleleQC[counter+i] <- "-";
	}
}

flagHLAQC <-matrix(rep(F, length(AlleleQC)*numHLA), ncol=numHLA);
FreqQC <- rowSums(GenoQC==T)/length(GenoQC[1,]);
for (i in 1:numHLA) {
	flagHLAQC[substr(AlleleQC, 1, nchar(HLA[i]))==HLA[i],i] <- T;
}

CalcEps<- function(Num1, Num2, Geno1, Geno2, numSample) {
	S_Exp <- 0;
	S_Obs <- 0;
	
	for (k in 1:Num1) {
		for (l in 1:Num2) {
			
			FreqExp <- sum(Geno1[k,])*sum(Geno2[l,])/numSample^2;
			FreqObs <- sum(Geno1[k,] & Geno2[l,])/numSample;
			
			S_Exp <- S_Exp + -1*FreqExp*log(FreqExp, 2);
			if (FreqObs>0) {
				S_Obs <- S_Obs + -1*FreqObs*log(FreqObs, 2);
			}
		}
	}
	eps <- 1-S_Obs/S_Exp;
	if (Num1==1 | Num2==1) {
		eps = 0;
	}
return(eps)
}

EPS <- matrix(numeric(choose(numHLA, 2)*5), ncol=5);
tmpEPS <- matrix(numeric(choose(numHLA,2)*Iter), ncol=choose(numHLA,2));
EPSmat <- matrix(numeric(numHLA^2), ncol=numHLA);

out <- paste("HLA1", "NumAllele", "HLA2", "NumAllele", "Epsilon", "Epsilon_mean_in_null", "Epsilon_sd_mean_in_null", "Z-score", sep="\t");
head(out)
write(out, file=outfile, append=F);

counter <- 1;
out <- "Iter";
for (i in 1:(numHLA-1)) {
	for (j in (i+1):numHLA) {
		EPS[counter,1] <- HLA[i];
		EPS[counter,2] <- sum(flagHLAQC[,i]);
		EPS[counter,3] <- HLA[j];
		EPS[counter,4] <- sum(flagHLAQC[,j]);
		
		Geno1 <- GenoQC[flagHLAQC[,i],];
		Geno2 <- GenoQC[flagHLAQC[,j],];
		Num1 <- sum(flagHLAQC[,i]);
		Num2 <- sum(flagHLAQC[,j]);
		
		if(Num1==1) {
			Geno1 <- matrix(as.vector(Geno1), ncol=length(Geno1));
		}
		if(Num2==1) {
			Geno2 <- matrix(as.vector(Geno2), ncol=length(Geno2));
		}
		
		tmpEPSval <- CalcEps(Num1, Num2, Geno1, Geno2, numSample);
		EPS[counter,5] <- tmpEPSval;
		EPSmat[i,j] <- as.numeric(tmpEPSval);
		EPSmat[j,i] <- as.numeric(tmpEPSval);
		
		for (k in 1:Iter) {
			tmpGeno2 <-  Geno2[, PermOrd(Fix)];
			if(Num2==1) {
				tmpGeno2 <- matrix(as.vector(tmpGeno2), ncol=length(tmpGeno2));
			}
			tmpEPS[k, counter] <- CalcEps(Num1, Num2, Geno1, tmpGeno2, numSample);
		}
		
		out <- paste(out, paste(HLA[i], "-", HLA[j], sep=""), sep="\t");
		
		counter <- counter+1;
	}
}
write.table(file=outfile, cbind(EPS, apply(tmpEPS, 2, mean), apply(tmpEPS, 2, sd), (as.numeric(EPS[,5])-apply(tmpEPS, 2, mean))/apply(tmpEPS, 2, sd)), quote=F, col.names=F, row.names=F, sep="\t", append=T);
write.table(EPSmat, file=outfile1, sep="\t", quote=F, row.names=F, col.names=HLA)

DataQC <- matrix(rep("-", Fix*(numHLA)), ncol=numHLA);
DataQC2 <- matrix(rep("-", Fix*(numHLA)*2), ncol=numHLA);
HapQC <- numeric(Fix*2);

for (i in 1:numHLA) {
	tmpAllele <- AlleleQC[flagHLAQC[,i]];
		
	for (j in 1:Fix) {
			A1 <- "NA";
			A2 <- "NA";
			tmpA1 <- tmpAllele[GenoQC[flagHLAQC[,i],j*2-1]];
			tmpA2 <- tmpAllele[GenoQC[flagHLAQC[,i],j*2]];
			if (length(tmpA1)==1) {
				A1 <- substr(tmpA1, nchar(tmpA1)-4, nchar(tmpA1));
			}
			if (length(tmpA2)==1) {가
				A2 <- substr(tmpA2, nchar(tmpA2)-4, nchar(tmpA2));
			}
			DataQC[j, i] <- paste(A1, A2, sep=" ");
			DataQC2[j*2-1, i] <- A1;
			DataQC2[j*2, i] <- A2;
			
	}
}


for (i in 1:(Fix*2)) {
	HapQC[i] <- paste(DataQC2[i,], collapse="-");
}
UniHapQC <- unique(HapQC);
UniHapQCFreq <- numeric(length(UniHapQC));
for (i in 1:length(UniHapQC)) {
	UniHapQCFreq[i] <- sum(HapQC==UniHapQC[i])/length(HapQC);
}
	
write.table(cbind(UniHapQC, UniHapQCFreq), file=outfile2, sep="\t", col.names=F, row.names=F, quote=F);