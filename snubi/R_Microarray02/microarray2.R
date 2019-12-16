if(!requireNamespace("BiocManager",quitely = TRUE))
  install.packages("BiocManager")  
BiocManager::install("affy")
#BiocManager::install("hgu113a.db")
BiocManager::install("hgu133plus2.db")

BiocManager::install("samr")
library(affy)
#library(hgu113a.db)
library(hgu133plus2.db)
library(samr)

setwd("c:/Users/user/git/SungminCode/snubi/R_Microarray02/")
no <- as.matrix(read.table("af_no.txt",header = TRUE, sep = '\t'))
mode(no)

head(no)

l = c(1,1,1,2,2,2)
g = factor(l)
t.test(no[1,]~g)

tp = vector()
tf = vector()

for(i in 1:nrow(no)){
  tmp = t.test(no[i,] ~ g, paired = F)
  tf[i] = log2(tmp$estimate[1] / tmp$estimate[2])
  tp[i] = tmp$p.value
}

id = which(tp<0.05)  ## index value in tp < 0.05
head(id)
length(id)
head(tp[id])

#다중 비교 문제 시 p-value를 그대로 쓰면 안되기에 보정한 값을 사용하게 된다.
##단순하게 p-value 값으로 그 유전자가 의미있는건지 확인할 수 없다... 다중 오류
pfdr = p.adjust(tp, method = 'fdr')
f.id = which(pfdr<0.05)
head(f.id)
head(pfdr[f.id])

plot(tf, -log10(tp))
points(tf[f.id],-log10(tp)[f.id], col = 'red')

f.no = cbind(no,tp,tf,pfdr)
dep = f.no[f.id,]

head(dep)
dim(dep)
write.table(dep,"ttest_fdr.txt",sep = "\t",row.names = T)

## probe를 gene symbol로 변환하기.
pid = row.names(no)[f.id]
gn = unlist(mget(pid[!is.na(pid)],hgu133plus2SYMBOL)) # annotation
head(gn)
deg = cbind(dep,gn)
head(deg)
degs = deg[order(deg[,9]),]
write.table(degs,"deg.txt",sep = '\t', row.names = T)


###samr

sg = c("1","1","1","2","2","2")
sg
sm = list(x = no, y =sg, logged2 = TRUE)
head(sm)
st = samr(sm,resp.type = "Two class unpaired",nperm = 100)
dt = samr.compute.delta.table(st)
head(dt)

d = 1.30
samr.plot(st,d)
st = samr.compute.siggenes.table(st,d,sm,dt)
names(st)
head(st$genes.up)
#write.table(st$genes.up,"sam.txt",sep = '\t',row.names = T)

