#http://www.snubi.org/R_seminar/case_study.pdf

###############################
##########   CASE I  ##########
###############################

##### 1. CGDS-R package install

install.packages('cgdsr') # Package install
library(cgdsr) 

###############################

##### 2. Get list of cancer studies at server

mycgds = CGDS("http://www.cbioportal.org/")
test(mycgds)
cancerstudy <- getCancerStudies(mycgds)
head(cancerstudy)


cancerstudy$name
cancerstudy[165,1]
cancerstudy[165,]
cancerstudy[166,]

mycancerstudy = cancerstudy[165,1]

mycancerstudy = cancerstudy[178,1]

###############################




##### 3. Extract samples and features

getCaseLists(mycgds,mycancerstudy)[,1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[2,1]
mycaselist


mymutationprofile = getGeneticProfiles(mycgds,mycancerstudy)[6,1] 
head(mymutationprofile)
mymethylationprofile=getGeneticProfiles(mycgds,mycancerstudy)[5,1] 
###############################




##### 4. Get mutation profiles for BRCA1 and BRCA2

brca_mutation = getMutationData(mycgds,mycaselist,mymutationprofile,c('BRCA1','BRCA2') ) 
head(brca_mutation)
table(brca_mutation$gene_symbol) 
brca1_mutated_cases <- brca_mutation[which(brca_mutation$gene_symbol=='BRCA1'),3]
brca2_mutated_cases <- brca_mutation[which(brca_mutation$gene_symbol=='BRCA2'),3]
head(brca1_mutated_cases)
head(brca2_mutated_cases)
###############################




##### 4. Extract samples with BRCA1, BRCA2 methylation

brca_methylation = getProfileData(mycgds, c('BRCA1','BRCA2'),mymethylationprofile, mycaselist)
?getProfileData
head(brca_methylation)
brca1_methylation_cases = rownames(brca_methylation[which(brca_methylation$BRCA1>0.8),])
brca1_methylation_cases
brca2_methylation_cases = rownames(brca_methylation[which(brca_methylation$BRCA2>0.8),])
length(brca2_methylation_cases)
###############################


##### 5. Clinical data integration

myclinicaldata = getClinicalData(mycgds,mycaselist)
head(myclinicaldata)

myclinicaldata$OS_STATUS[myclinicaldata$OS_STATUS == ""] <- NA
head(myclinicaldata$OS_STATUS)
myclinicaldata$OS_MONTHS



total_sample <- rownames(myclinicaldata) 
type <- rep('Wild', length(total_sample))
names(type) <- total_sample
head(type)


brca1_mutated_cases = gsub("-",".",brca1_mutated_cases) 
brca2_mutated_cases = gsub("-",".",brca2_mutated_cases)
head(brca1_mutated_cases)

type[brca1_methylation_cases] <- "BRCA1_methylation"
type[brca2_methylation_cases] <- "BRCA2_methylation"
type[brca1_mutated_cases] <- "BRCA1_mutation"
type[brca2_mutated_cases] <- "BRCA2_mutation"

type <- type[names(type) %in% rownames(myclinicaldata)]
type <- factor(type, levels= c("Wild","BRCA1_mutation","BRCA2_mutation","BRCA2_methylation") )
###############################




##### 6. Survival Analysis


install.packages('survival')
library(survival)


out <- survfit(Surv(OS_MONTHS, OS_STATUS=="DECEASED")~type, data=myclinicaldata)

survdiff(Surv(OS_MONTHS, OS_STATUS=="DECEASED") ~ type,data=myclinicaldata)  #chisq °á°ú
coxph(Surv(OS_MONTHS, OS_STATUS=="DECEASED") ~ type, data=myclinicaldata)

color <- c("black","Skyblue","Blue", "Red") 
plot(out, col=color , main="Association of BRCA1/2 Mutations with Survival", xlab="Time,days", ylab="Proportion", lty=1:4, lwd=2)
legend("topright", levels(type), col=color, lty=1:4, lwd=3) 

###############################
# survminer ggplot 

install.packages("survminer")

#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/survminer")

library("survminer")

ggsurvplot_res <- ggsurvplot(out, data=myclinicaldata, pval=T, risk.table=T, conf.int=F, break.time.by = 30, legend.title = "Patient types", risk.table.fontsize = 2.5,surv.plot.height = 0.45, legend="right")

ggsurvplot_res





###############################
##########   CASE II  #########
###############################

##   expression data download link : 
##   http://www.snubi.org/R_seminar/case.zip

#####  1. Get the Expression data

setwd("c:/Users/user/Desktop/snubi/")
expression <- read.delim('expression.csv', sep=",", header=T, stringsAsFactors = F)

dim(expression)
length(unique(rownames(expression)))  


#####  2. Check variability of gene expression across patients, using MAD
myMad <- function(x){mad(x, na.rm=T)}
result <- apply(expression, 1, myMad)


#####  3. Get the highest variablity genes

result2 <- sort(result, decreasing=T)[1:100] 
exp_result <- expression[which(rownames(expression)%in%names(result2)),]



#####  4. Install NMF package

install.packages("NMF")
library(NMF)


#####  5. Make the values non-negative

nonNegativeFx <- function(x){
  a <- min(x)
  if(a<0){
    a <- a*(-1)
    b <- x + a}
  else{
    b <- x }
  return(b)
}



res2 <- apply(exp_result,1, nonNegativeFx) 
res2 <- as.data.frame(res2)



#####  6. Draw heatmap and NMF analysis

nmf_res <- nmf(res2, 2:6, nrun=10, seed=123456)

plot(nmf_res)

consensusmap(nmf_res, annCol = colnames(nmf_res), labCol = NA, labRow = NA)

