if(!requireNamespace("BiocManager",quitely = TRUE))
install.packages("BiocManager")  
BiocManager::install("affy")
BiocManager::install("DBI")
BiocManager::install("RSQLite")
BiocManager::install("hthgu133acdf")
library(affy)

setwd('C:/Users/user/git/SungminCode/snubi/R_MicroArray01/')
getwd()
d <- ReadAffy()
d
dl <- log2(exprs(d))
head(dl)
image(d)

###### Quantile normalization
q <- as.matrix(read.table("quantile.txt",header = TRUE, row.names = 1, sep = "\t"))
head(q)
qs <- matrix(ncol = ncol(q), nrow = nrow(q)) # sort된 유전자 발현량을 담을 벡터 만들기
qr <- matrix(ncol = ncol(q), nrow = nrow(q)) #sample 내에 유전자 발형량 순위를 담을 벡터 만들기

for (i in 1:ncol(q)){
  qs[,i] = sort(q[,i])
  qr[,i] = rank(q[,i])
}

qm <- apply(qs,1,mean)
qn <- matrix(ncol=ncol(q), nrow = nrow(q))
for (i in 1:length(qr)){
  r <- qr[i]
  qn[i] <- qm[r]
}

head(qn)

rownames(qn) = rownames(q)
colnames(qn) = colnames(q)

par(mfrow = c(1,2))
boxplot(q,main = "Before normalization")
boxplot(qn,main = "After normalization")


#######RMA

d_rma = rma(d)
head(d_rma)
dr = exprs(d_rma)
head(dr)

#####mas5
d_mas5 <- mas5(d)
dm <- (exprs(d_mas5))
head(dm)

#####expresso method
d_es <- expresso(d, bgcorrect.method = "none",normalize.method  = "quantiles",pmcorrect.method = "pmonly",summary.method = "medianpolish")
des = exprs(d_es)  #Expresso expression table 뽑기


#####plot


##density plot 
par(mfrow = c(1,2))
plotDensity(dl)
plotDensity(des)


## box plot
boxplot(dl)
boxplot(des)
boxplot(dl, col = rainbow(15), main = "before norm")
boxplot(des, col =rainbow(15), main = "after norm")

par(mfrow=c(1,1))
pairs(dl)		# scatter plot before normalization
pairs(des) 	# scatter plot after normalization

pairs(dl, lower.panel=panel.smooth,upper.panel = panel.smooth)
pairs(des, lower.panel=panel.smooth,upper.panel = panel.smooth)

###### MA plot ############
mva.pairs(dl)  # MA plot  before normalization
mva.pairs(des) # MA plot  after normalization