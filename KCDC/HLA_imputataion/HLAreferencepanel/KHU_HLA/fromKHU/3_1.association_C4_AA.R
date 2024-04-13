library(dplyr);library(tidyr);library(optparse);library(bigreadr);
library(foreach);library(doParallel);library(lmtest);library(data.table)


####################################### arguement setting ##########################################
args_list <- list(
    make_option( "--disease" , type="character" , default= "SLE" , help="snpversion", metavar="character")
)

opt_parser <- OptionParser( option_list = args_list)
opt <- parse_args(opt_parser)

# 분석하시는 phenotype으로 해주시면 됩니다
# opt=list(); opt$disease="SLE"
###################################################################################################


setwd("~~/wholesample_imputation_result") )

# 제경우에는 아래 ped파일에 covariate 정보들이 있어서 불러왔습니다. 선생님 데이터로 고치시면 됩니다.!
ped <- read.table("/kimlab_wd/kah/0.data/1.4th_KCHIP/3.imputed/2.ped/2.4th.kchip.ra.vcfID.ver3.ped" ,  header=T, sep="\t")
# ped파일 형식은 아래와 같이 생겼습니다!
# > head(ped)
#             FID           IID FAT_ID MOT_ID SEX pheno RA RA_kin
# 1 BRG001_BRG001 BRG001_BRG001      0      0   2    RA  2      2
# 2 BRG003_BRG003 BRG003_BRG003      0      0   1    RA  2      2
# 3 BRG004_BRG004 BRG004_BRG004      0      0   2    RA  2      2
# 4 BRG005_BRG005 BRG005_BRG005      0      0   2    RA  2      2
# 5 BRG006_BRG006 BRG006_BRG006      0      0   2    RA  2      2
# 6 BRG007_BRG007 BRG007_BRG007      0      0   2    RA  2      2
#   RA_kin_with_omni RA_kin_w_omni_ver2 SLE_kin          PC1           PC2
# 1                2                  2      NA  0.007988699 -0.0090509930
# 2                2                  2      NA -0.004666920  0.0096160220
# 3                2                  2      NA -0.006730679  0.0105728000
# 4                2                  2      NA  0.011074450  0.0040453750
# 5                2                  2      NA  0.012846880  0.0009969028
# 6                2                  2      NA  0.005479033 -0.0008907822
#            PC3          PC4          PC5          PC6          PC7          PC8
# 1 -0.003603393 -0.015412200  0.004744021  0.009827589 -0.003971098 -0.001109956
# 2 -0.036126950 -0.004454542  0.020110690 -0.020720010 -0.003709694  0.001777732
# 3 -0.001532866 -0.025977620  0.017300390 -0.012624810 -0.002845906 -0.003472627
# 4  0.022344970 -0.003073701  0.012259560  0.005924254  0.012109840  0.008884985
# 5 -0.015044000  0.008804950  0.006273200 -0.011356020  0.012095590  0.007184493
# 6  0.005657466  0.015086760 -0.007748925  0.000788582  0.013418840  0.013645980
#            PC9          PC10
# 1  0.006983857  0.0132571600
# 2  0.003700590 -0.0149817100
# 3  0.007968257  0.0034289330
# 4  0.006174090 -0.0045441620
# 5 -0.005825741 -0.0006802852
# 6  0.012253170  0.0048812190


dosage_ori=fread2("./dosage_imputed_snps_only.txt", header=T,sep="\t") # 2.after_imputation.sh 에서 만든 파일입니다!
colnames(dosage_ori)[1] <- "[1]ID"
colnames(dosage_ori) = gsub("\\[\\d+\\]", "", colnames(dosage_ori)) %>% gsub(":DS", "", ., fixed=TRUE)
dosage_ori = select(dosage_ori, !c("CHROM", "REF", "ALT", "POS"))



########## AA Omnivus test

AA_glm <- function(){

AA_dosage=dosage_ori[grep("^AA", dosage_ori$ID),]
# AA_dosage=AA_dosage[1:100,]

rownames(AA_dosage)=AA_dosage$ID
AA_dosage=AA_dosage[-1]
AA_dosage=t(AA_dosage)
AA_dosage=as.data.frame(AA_dosage)
AA_dosage$ID=rownames(AA_dosage)
AA_dosage= relocate(AA_dosage, ID)
rownames(AA_dosage)=1:nrow(AA_dosage)

tt <- cbind( AA_dosage, ped[-c(1:4)])
tt[ , "SLE_kin"] = gsub( 1 , 0 , tt[, "SLE_kin"] ) %>% gsub(2, 1, .) %>% as.numeric()
tt[ , "RA_kin"] = gsub( 1 , 0 , tt[, "RA_kin"] ) %>% gsub(2, 1, .) %>% as.numeric()

AAs = unique( sapply( strsplit( grep("AA_", colnames(tt), value=T) , "_" ) , function(x) paste0( c( x[1], x[2], x[3] ), collapse="_") ) )

k = 0
result_for_AAomni = data.frame()
for(v in AAs){
    print(k)
    temp = grep(v,colnames(tt), value=T)
    temp = temp[nchar(temp) == min(nchar(temp))]

    temp = temp[colMeans(tt[, temp , drop=FALSE])/2  > 0.005 ]

    if(length(temp) >= 3){
    temp = temp[colMeans(tt[,temp])/2 != max(colMeans(tt[,temp])/2)]
    }


if( length(temp) !=0 ) {

    temp = paste0("`", temp ,"`")
    temp = paste0(temp, collapse = " + ")

#### 제가 SLE/RA 예시로 넣어놓았습니다 phenotype에 맞게 수정해주시면 됩니다!
    if(opt$disease=="SLE"){
        temp = paste0(" glm ( SLE_kin ~ ", temp , " + SEX + PC1 + PC2 + PC3 + PC4 + PC5 , family=binomial, data=tt ) ")
    }else if(opt$disease=="RA"){
        temp = paste0(" glm ( RA_kin ~ ", temp , " + SEX + PC1 + PC2 + PC3 + PC4 + PC5 , family=binomial, data=tt ) ")
    }
        full_glm=eval(parse(text = temp))

    if(opt$disease=="SLE"){
        nested_glm = glm( SLE_kin ~ SEX + PC1 + PC2 + PC3 + PC4 + PC5 ,  family=binomial, data=tt)
    }else if(opt$disease=="RA"){
        nested_glm = glm( RA_kin ~ SEX + PC1 + PC2 + PC3 + PC4 + PC5 ,  family=binomial, data=tt)
    }
        test = lrtest(nested_glm, full_glm)
        tmp_res = data.frame( ID = v, Df=test$Df[2] , chisq=test$Chisq[2], pval = test$Pr[2])

        result_for_AAomni = rbind(result_for_AAomni, tmp_res)
        k=k+1
}

}

return(result_for_AAomni)
}

result_for_AAomni1=AA_glm()
arrange(result_for_AAomni1, pval) %>% head (n=10)

write.table(result_for_AAomni1, paste0( "~~~/", opt$disease, "/result_for_AAomni1.txt" ), col.names=T, row.names=F, sep="\t", quote=F)


##### C4 analysis

setwd("~~~/wholesample_imputation_result") )

ped <- read.table("/kimlab_wd/kah/0.data/1.4th_KCHIP/3.imputed/2.ped/2.4th.kchip.ra.vcfID.ver3.ped" ,  header=T, sep="\t")

C4_dosage=fread2("./tableC4_diploid.txt", header=T,sep="\t")


C4_glm <- function(){
tt <- cbind( C4_dosage, ped[-c(1:4)])
tt[, "SLE_kin"] = gsub( 1, 0 , tt[, "SLE_kin"] )%>%  gsub(2, 1, .) %>% as.numeric()
tt[, "RA_kin"] = gsub( 1, 0 , tt[, "RA_kin"])%>%  gsub(2, 1, .) %>% as.numeric()

res_C4=data.frame();k=0

for( v in colnames(C4_dosage)[-1] ){
    print(k)


### 수정해주시면 됩니다!
if(opt$disease=="SLE"){
    temp=paste0( "glm( SLE_kin ~ `", v, "` + SEX + PC1 + PC2 + PC3 + PC4 + PC5 , family=binomial, data=tt) %>% summary() " )
}else if(opt$disease=="RA"){
    temp=paste0( "glm( RA_kin ~ `", v, "` + SEX + PC1 + PC2 + PC3 + PC4 + PC5 , family=binomial, data=tt) %>% summary() " )
}
    r=eval(parse(text=temp))

    if( length( which( names( table( rownames(r$coef)==paste0("`",v,"`") ) ) %in% "TRUE"  )) ==1 ){
    rest <- r$coef[rownames(r$coef)==paste0("`",v,"`"), ]
    }
    if( length( which( names( table( rownames(r$coef)==v ) ) %in% "TRUE"  )) ==1 ){
    rest <- r$coef[rownames(r$coef)==v, ]
    }

    rest = t( data.frame( rest ) )
    rownames(rest) = v
    res_C4= rbind(res_C4, rest)
    k=k+1
}

return(res_C4)
}

result_for_C41 = C4_glm()
colnames(result_for_C41) <- c("Est", "SE", "z", "pval")
result_for_C41$ID = rownames(result_for_C41 )
result_for_C41 = relocate(result_for_C41, ID)

write.table(result_for_C41, paste0( "~~~/", opt$disease, "/result_for_C41_dosage.txt" ), col.names=T, row.names=F, sep="\t", quote=F)


#################### combine Three

setwd(paste0( "~~~/wholesample_imputation_result") )

markers=read.table("./wholesample.imputed.markers", header=F)
markers=rename(markers, ID=V1)
imputed_snps <- read.table("./imputed_snps.txt", header=F)
imputed_markers= markers[which(markers$V2 %in% imputed_snps$V2),]
imputed_markers=rename(imputed_markers, POS=V2)
head(imputed_markers)

# merge three data
res_ref = fread2(paste0("zcat " , "~~~/", opt$disease, "/rvtest.", opt$disease,"1.MetaScore.assoc.gz" , "| grep -v '^##'"   ) )
res_ref=as.data.frame(res_ref)
res_ref_forsave=merge( imputed_markers[c("POS", "ID")], res_ref, by="POS") %>% arrange(PVALUE)
head(res_ref_forsave)
write.table(res_ref_forsave, paste0( "~~~/", opt$disease, "/result_for_rvtest1.txt" ), col.names = T, row.names = F, sep="\t", quote=F)
res_ref=merge( imputed_markers[,1:2], res_ref[c("POS", "PVALUE")], by="POS")
res_ref = rename(res_ref, pval=PVALUE)


res_AA = read.table(paste0( "~~~/"  , opt$disease, "/result_for_AAomni1.txt" ), header=T, sep="\t" )
res_c4_dosage = read.table( paste0( "~~~/"  , opt$disease, "/result_for_C41_dosage.txt" ) , header=T, sep="\t" )


head(res_ref); head(res_AA);

res_ref$label= "nonHLA_SNPs"
res_ref[grep("HLA", res_ref$ID), "label"] <- "HLA"
res_ref[grep("SNPS", res_ref$ID), "label"] <- "HLA_SNPS"
res_ref[grep("^AA", res_ref$ID), "label"] <- "AA"
res_ref[grep("copy", res_ref$ID), "label"] <- "C4_biallele"
res_AA$label="AA_omni"
res_c4_dosage$label="C4_dosage"


res_all = rbind( res_ref[c("ID", "pval", "label")] , res_AA[c("ID", "pval", "label")] ,  res_c4_dosage[c("ID", "pval", "label")])
res_all = arrange( res_all, pval)


# labeling
map_2f <- read.table("~~~/HLA_DICTIONARY_AA.hg38.imgt3320.orderchange.not2f.map", header=F)  # for HLA/AA
map_gg <- read.table("~~~/HLA_DICTIONARY_SNPS.hg38.Ggroup_4field.map", header=F)  # for SNPs
hlalocation = read.table("/kimlab_wd/yuo1996/C4_HLA_ref/try6/2.MakeReference_with_HLA/ver3/hla_location", header=F)
hlalocation = rbind( hlalocation, data.frame( V1= c("C4", "HERV") , V2=c(31980581, 31984684)) )


head(res_all)
nrow(res_all)

# imputed marker만!!
res = merge(res_all, imputed_markers[c("ID", "POS")], by="ID") %>% rename(pos=POS) %>% arrange(pos)
# AAomni + C4는 포함안되서 따로해줘야.
res= rbind(res, filter(res_all, label %in% c("AA_omni", "C4_binary", "C4_dosage")) %>% mutate(pos=NA) )
head(res)
tail(res)

filter(res, label=="HLA")
filter(res, label=="AA")
filter(res, label=="C4_binary")
filter(res, label=="C4_dosage")

for(i in 1:nrow(hlalocation)){
    h=hlalocation[i, "V1"]
    res[grep(h, res$ID), "pos"]  <- hlalocation[i, "V2"]
}

filter(res, label=="HLA")
filter(res, label=="AA")

head(map_2f)
for(i in 1:nrow(map_2f)){
    print(paste(i, "___", nrow(map_2f)) )
    h=map_2f[i, "V2"]
    res[grep(h, res$ID), "pos"]  <- map_2f[i, "V4"]
}


head(res)
# for AAomni
for(i in filter(res, label=="AA_omni")$ID  ){
    print(i)
    res[which(res$ID==i), "pos"]  <- map_2f[ grep( paste0(i, "_"), map_2f$V2), "V4"]
}


head(map_gg)
for(i in 1:nrow(map_gg)){
    print(paste(i, "___", nrow(map_gg)) )
    h=map_gg[i, "V2"]
    res[grep(h, res$ID), "pos"]  <- map_gg[i, "V4"]
}

filter(res, label=="C4_binary")
filter(res, label=="C4_dosage")
arrange(res, pval ) %>% head
arrange(res, pval ) %>% filter( label %in% c("HLA", "AA_omni" )) %>% head (n=20)
filter(res, label=="AA_omni")
table( is.na(res$pos) )


write.table(res, paste0("~~~/", opt$disease, "/result_all1_forplot.txt" ) , col.names=T, row.names=F, sep="\t", quote=F)
