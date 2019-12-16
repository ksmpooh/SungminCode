BiocManager::install("MASS")
BiocManager::install("DAAG")

source("http://www.biostat.wics.edu/~kbroman/teaching/stat371/permfunc.R")
library(MASS)
library(DAAG)


data(birthwt)
attach(birthwt)
head(birthwt)

nor = birthwt[(birthwt[,1] == '0'),3]
head(nor)
low = birthwt[(birthwt[,1] == '1'),3]

tobs = t.test(nor,low)$statistic
tobs

##permutation test 

tperm = perm.test(nor, low, n.perm = 1000)
hist(tperm)
abline(v=abs(tobs), lty = 2, col = 2)

## empirical p-value

pvalue = mean(abs(tperm) >= abs(tobs))
pvalue



###cross Validation

head(birthwt)
bwt_kg = birthwt$bwt/1000
lwt_kg = birthwt$lwt * 0.45
smk = as.factor(birthwt$smoke)
mydata = data.frame(bwt_kg,lwt_kg,smk)
head(mydata)
dim(mydata)
train = mydata[-seq(1,180,by=20),]
dim(train)
attach(train)
test_in = mydata[seq(1,180,by= 20),c(2,3)]
dim(test_in)
head(test_in)
test_out = mydata[seq(1,180,by = 20),c(1)]


model = lm(bwt_kg~lwt_kg+smk,data = train)
summary(model)

# 10 -fold cross validation
cross = CVlm(train, m = 10, form.lm = model)

cross


p_out = predict(model, test_in)
p_out
test_out
plot(p_out,test_out)
abline(0,1,col = 'red',lwd = 0.5)
#########correction method
# bonferroni 
pbon = p.adjust(pvalue, method = "bonferroni")
pfdr = p.adjust(pvalue, method = 'fdr')

pbon
pfdr



## ²ÜÆÁ
install.packages("swirl")
library(swirl)
ls()
