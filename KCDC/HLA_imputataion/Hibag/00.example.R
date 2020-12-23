if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("HIBAG")
library("HIBAG")

browseVignettes("HIBAG")

head(HLA_Type_Table)
dim(HLA_Type_Table)



hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
                 H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
                 H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
                 locus=hla.id, assembly="hg19")
hla
# divide HLA types randomly
set.seed(100)
hlatab <- hlaSplitAllele(hla, train.prop=0.5)
names(hlatab)
head(hlatab)
summary(hlatab)
??hlaSplitAllele
# "training" "validation"
summary(hlatab$training)
summary(hlatab$validation)
# SNP predictors within the flanking region on each side
summary(HapMap_CEU_Geno)
HapMap_CEU_Geno$snp.allele
region <- 500 # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
                        hla.id, region*1000, assembly="hg19")
length(snpid) # 275
# training and validation genotypes
train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
                            snp.sel=match(snpid, HapMap_CEU_Geno$snp.id),
                            samp.sel=match(hlatab$training$value$sample.id,
                                           HapMap_CEU_Geno$sample.id))
test.geno <- hlaGenoSubset(HapMap_CEU_Geno,
                           samp.sel=match(hlatab$validation$value$sample.id,
                                          HapMap_CEU_Geno$sample.id))
# train a HIBAG model
set.seed(100)
# please use "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hlatab$training, train.geno, nclassifier=4,
                        verbose.detail=TRUE)
summary(model)
# validation
pred <- predict(model, test.geno)
summary(pred)
# compare
(comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model,
                          call.threshold=0))
(comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model,
                          #HapMap_CEU_Geno 5
                          call.threshold=0.5))
# save the parameter file
mobj <- hlaModelToObj(model)
save(mobj, file="HIBAG_model.RData")
save(test.geno, file="testgeno.RData")
save(hlatab, file="HLASplit.RData")
# Clear Workspace
hlaClose(model) # release all resources of model
rm(list = ls())
######################################################################
# NOW, load a HIBAG model from the parameter file
mobj <- get(load("HIBAG_model.RData"))
model <- hlaModelFromObj(mobj)
# validation
test.geno <- get(load("testgeno.RData"))
hlatab <- get(load("HLASplit.RData"))
pred <- predict(model, test.geno)
# compare
(comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model,
                          call.threshold=0.5))
#########################################################################
# import a PLINK BED file
#
bed.fn <- system.file("extdata", "HapMap_CEU.bed", package="HIBAG")
fam.fn <- system.file("extdata", "HapMap_CEU.fam", package="HIBAG")
bim.fn <- system.file("extdata", "HapMap_CEU.bim", package="HIBAG")
hapmap.ceu <- hlaBED2Geno(bed.fn, fam.fn, bim.fn, assembly="hg19")
#########################################################################
# predict
#
hapmap.ceu
HapMap_CEU_Geno
pred <- predict(model, hapmap.ceu, type="response")
head(pred$value)
# sample.id allele1 allele2 prob
# 1 NA10859 01:01 03:01 0.9999992
# 2 NA11882 01:01 29:02 1.0000000
# ...
# delete the temporary files
unlink(c("HIBAG_model.RData", "testgeno.RData", "HLASplit.RData"), force=TRUE)
