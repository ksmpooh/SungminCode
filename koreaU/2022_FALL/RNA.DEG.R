# devtools
#system('sudo apt-get install -y libfribidi-dev', intern=TRUE)
#install.packages("devtools") #기본적인인 r library 설치

# tidyverse and DT libraries
#install.packages(c("tidyverse", "DT"))

# TCGAbiolinks to download RNA-seq data
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks", ask=F)

# DESEQ2 관련 라이브러리 설치
BiocManager::install("DESeq2", ask=F)

# WGCNA 관련 라이브러리 설치
BiocManager::install("WGCNA", ask=F)
BiocManager::install("limma", ask=F)
install.packages("cowplot")

# GSEA 관련 라이브러리 설치
install.packages("msigdbr")
install.packages(c("tidyverse", "DT"))

BiocManager::install("fgsea", ask=F)


# 필요한 라이브러리 불러오기
library(TCGAbiolinks)
library(tidyverse)
library(DT)
library(DESeq2)
library(WGCNA)


# 분석에 사용할 데이터 구성 설정
query <- GDCquery(
  project = "TCGA-LUSC",
  data.type = "Gene Expression Quantification", 
  data.category = "Transcriptome Profiling",
  workflow.type = "STAR - Counts"
)

# 다운로드
# 다운로드에 오랜 시간이 소요될수도 있습니다
'
GDCdownload(
  query = query, 
  method = "api"
)
'

setwd("~/Desktop/KU/2022_FALL/BDS/RNA.DEG.AN.pro/")
df <- GDCprepare(query = query)
head(df)
e = SummarizedExperiment::assay(df) %>% as_tibble()
g = SummarizedExperiment::rowData(df) %>% as_tibble()
s = SummarizedExperiment::colData(df) %>% as_tibble()

head(s)
write.csv(e,"e.csv")
write.csv(g,"g.csv")
write.csv(s,"s.csv")


s %>% dplyr::count(definition)
s %>% dplyr::count(shortLetterCode)
# 간단한 DESEQ2 분석을 위해, 소수의 샘플을 추출하겠습니다
s_tumor = s %>% dplyr::filter(shortLetterCode=="TP") %>% dplyr::slice(1:5) %>% pull(barcode)
s_normal = s %>% dplyr::filter(shortLetterCode=="NT") %>% dplyr::slice(1:5) %>% pull(barcode)
e1 = e %>% dplyr::select(s_tumor, s_normal)
dim(e1)
s1 = s %>% filter(barcode %in% c(s_tumor, s_normal))
dim(s1)

s1 = s1 %>% mutate(is_tumor = ifelse(shortLetterCode=="TP", 1, 0),
                   is_tumor = factor(is_tumor, levels=c(0, 1))) %>% column_to_rownames("barcode")

head(e1)
colnames(e1)
colnames(s1)
table(s1$is_tumor)
head(dds)

# DESEQ2 객체 만들기. DESEQ2 분석을 위해선, gene expression, sample information을 모두 하나의 객체로 취합하고 계산을 해야합니다.

dds <- DESeqDataSetFromMatrix(countData = e1,
                              colData = s1,
                              design= ~ is_tumor)

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

# 사용할 수 있는 intercept (비교대상)은 'is_tumor_1_vs_0'에 해당합니다. 이 컬럼 정보를 이용하여, 차이점을 살펴보겠습니다. 

res <- results(dds, name="is_tumor_1_vs_0")

# res 객체를 데이터프레임으로 변경합니다. 이때 gene information을 같이 추가합니다. 
res1 = data.frame(gene_id = g$gene_id, 
                  gene_name = g$gene_name, 
                  baseMean=res$baseMean, 
                  lfc=res$log2FoldChange, 
                  lfcSE=res$lfcSE, p=res$pvalue, 
                  padj=res$padj) %>% as_tibble()

# 어떤 유전자들이 크게 변화하였는지 살펴봅시다.
# lfc 는 log2 fold change로써 tumor와 normal 그룹에 따른 발현량 차이를 나타낸 컬럼입니다
# Deseq2는 wald test를 이용하여 p값을 산출하고, 이에 대한 다중비교 값을 padj 컬럼으로 제공합니다. 자세한 설명은 아래 링크를 참고하세요
# https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html

res1 %>% arrange(padj) %>% head()

# 얻은 결과를 시각화 하겠습니다
res1 %>% mutate(isSig = ifelse(padj < 0.05 & lfc > 0, "Up", "None"),
                isSig = ifelse(padj < 0.05 & lfc < 0, "Down", isSig)) %>%
  ggplot(aes(lfc, -log10(padj), fill=isSig)) + 
  geom_point(shape=21) + 
  scale_fill_manual(values=c("Up"="firebrick", "Down"="dodgerblue", "None"="grey"))



# msigdb download
library(msigdbr)
df_msigdb = msigdbr(species = "Homo sapiens")

# Hallmark pathway 만 써봅시다
df_msigdb_hallmarks = df_msigdb %>% filter(gs_cat=="H")
msigdbr_list = split(x = df_msigdb_hallmarks$gene_symbol, f = df_msigdb_hallmarks$gs_name)
# fgsea를 이용하여 GSEA 분석
library(fgsea)
input_ranks = res1 %>% filter(!is.na(lfc) & !is.na(padj)) %>% arrange(-abs(lfc)) %>% 
  filter(!duplicated(gene_name)) %>% pull(lfc, name=gene_name)

fgseaRes <- fgsea(msigdbr_list, input_ranks)
fgseaRes %>% as_tibble() %>% arrange(pval) %>% head(10)

plotEnrichment(msigdbr_list[["HALLMARK_E2F_TARGETS"]], input_ranks) + labs(title="HALLMARK_E2F_TARGETS")

'
Homework
위에 연습한 DESEQ2 방법을 바탕으로 흡연자와 비흡연자를 추출하여 DEG 분석을 시행하겠습니다. 아래 샘플을 사용하여 진행하면 됩니다.

다음에 대한 답을 제출하시면 됩니다. 답은 워드파일에 제출하시면 됩니다.

Q1. 몇개의 DEG (pdaj < 0.05 미만) 유전자들이 나옵니까?
  Q2. Hallmark GSEA를 시행하면, 흡연자 그룹에서 가장 up-regulation 되는 pathway 가 무엇인가요? GSEA plot을 같이 제출해주세요.
  '

table(s$paper_Smoking.Status)
# 여러분 질문을 위해 소수의 샘플을 추출하겠습니다
s_smoker = s %>% dplyr::filter(paper_Smoking.Status=="Current smoker") %>% 
  dplyr::slice(1:5) %>% 
  pull(barcode)
s_nonsmoker = s %>% dplyr::filter(paper_Smoking.Status=="Lifelong Non-smoker") %>% 
  dplyr::slice(1:5) %>% 
  pull(barcode)
head(s_smoker)
head(s_nonsmoker)


e1 = e %>% dplyr::select(s_smoker, s_nonsmoker)
dim(e1)
head(e1)
head(s1)
s1 = s %>% filter(barcode %in% c(s_smoker, s_nonsmoker))
dim(s1)
colnames(s1)
head(s1$shortLetterCode)
s1 = s1 %>% mutate(is_smoker = ifelse(paper_Smoking.Status=="Current smoker", 1, 0),
                   is_smoker = factor(is_smoker, levels=c(0, 1))) %>% column_to_rownames("barcode")



colnames(s1)
table(s1$is_smoker)
head(s1)
# DESEQ2 객체 만들기. DESEQ2 분석을 위해선, gene expression, sample information을 모두 하나의 객체로 취합하고 계산을 해야합니다.

dds <- DESeqDataSetFromMatrix(countData = e1,
                              colData = s1,
                              design= ~is_smoker)

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

# 사용할 수 있는 intercept (비교대상)은 'is_tumor_1_vs_0'에 해당합니다. 이 컬럼 정보를 이용하여, 차이점을 살펴보겠습니다. 

res <- results(dds, name="is_smoker_1_vs_0")

# res 객체를 데이터프레임으로 변경합니다. 이때 gene information을 같이 추가합니다. 
res1 = data.frame(gene_id = g$gene_id, 
                  gene_name = g$gene_name, 
                  baseMean=res$baseMean, 
                  lfc=res$log2FoldChange, 
                  lfcSE=res$lfcSE, p=res$pvalue, 
                  padj=res$padj) %>% as_tibble()

# 어떤 유전자들이 크게 변화하였는지 살펴봅시다.
# lfc 는 log2 fold change로써 tumor와 normal 그룹에 따른 발현량 차이를 나타낸 컬럼입니다
# Deseq2는 wald test를 이용하여 p값을 산출하고, 이에 대한 다중비교 값을 padj 컬럼으로 제공합니다. 자세한 설명은 아래 링크를 참고하세요
# https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html

res1 %>% arrange(padj) %>% head()
res1 %>% filter(padj<0.05) %>% nrow()
res11<-res1 %>% mutate(isSig = ifelse(padj < 0.05 & lfc > 0, "Up", "None"),
                isSig = ifelse(padj < 0.05 & lfc < 0, "Down", isSig)) 
table(res11$isSig)
# 얻은 결과를 시각화 하겠습니다
res1 %>% mutate(isSig = ifelse(padj < 0.05 & lfc > 0, "Up", "None"),
                isSig = ifelse(padj < 0.05 & lfc < 0, "Down", isSig)) %>%
  ggplot(aes(lfc, -log10(padj), fill=isSig)) + 
  geom_point(shape=21) + 
  scale_fill_manual(values=c("Up"="firebrick", "Down"="dodgerblue", "None"="grey"))

input_ranks = res1 %>% filter(!is.na(lfc) & !is.na(padj)) %>% arrange(-abs(lfc)) %>% 
  filter(!duplicated(gene_name)) %>% pull(lfc, name=gene_name)

fgseaRes <- fgsea(msigdbr_list, input_ranks)

fgseaRes %>% as_tibble() %>% arrange(-NES) %>% head(10)

# msigdb download
library(msigdbr)
df_msigdb = msigdbr(species = "Homo sapiens")

# Hallmark pathway 만 써봅시다
df_msigdb_hallmarks = df_msigdb %>% filter(gs_cat=="H")
msigdbr_list = split(x = df_msigdb_hallmarks$gene_symbol, f = df_msigdb_hallmarks$gs_name)
# fgsea를 이용하여 GSEA 분석
library(fgsea)
input_ranks = res1 %>% filter(!is.na(lfc) & !is.na(padj)) %>% arrange(-abs(lfc)) %>% 
  filter(!duplicated(gene_name)) %>% pull(lfc, name=gene_name)

fgseaRes <- fgsea(msigdbr_list, input_ranks)
fgseaRes %>% as_tibble() %>% arrange(pval) %>% head(10)
fgseaRes %>% as_tibble() %>% arrange(NES) %>% head(10)
fgseaRes %>% as_tibble() %>% arrange(-NES) %>% head(10)

plotEnrichment(msigdbr_list[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], input_ranks) + labs(title="HALLMARK_INTERFERON_GAMMA_RESPONSE")

##########



#############3
# 간단한 DESEQ2 분석을 위해, 소수의 샘플을 추출하겠습니다
s_tumor_wgcna = s %>% dplyr::filter(shortLetterCode=="TP") %>% dplyr::slice(1:50) %>% pull(barcode)

e2 = e %>% dplyr::select(any_of(s_tumor_wgcna)) %>%
  as.data.frame() %>% 
  mutate(rowMeans = rowMeans(.)) %>% 
  mutate(gene_name = g$gene_name) %>%
  arrange(-rowMeans) %>% 
  filter(!duplicated(gene_name)) %>%
  column_to_rownames('gene_name') %>% dplyr::select(-rowMeans)

e2 = log2(e2 + 1)
e2 = limma::normalizeQuantiles(e2)

## Check the genes and samples
gsg = goodSamplesGenes(t(e2), verbose = 5)
table(gsg$goodGenes)
table(gsg$goodSamples)

input_mat = e2[gsg$goodGenes==T, gsg$goodSamples == T]
input_mat = input_mat[rowMeans(input_mat) >= 10,]
dim(input_mat)

Power = seq(2,15,by=1)
sft <- pickSoftThreshold(t(input_mat), powerVector=Power, networkType="signed", verbose=5)

p1 <- ggplot(sft$fitIndices, aes(Power, SFT.R.sq)) + 
  geom_point() + theme_minimal(base_size = 7) + 
  geom_hline(yintercept = 0.8, color="red") + 
  labs(title = 'Scale independence',
       x='Soft Threshold (power)', 
       y='Scale Free Topology Model Fit,signed R^2') 

p2 <- ggplot(sft$fitIndices, aes(Power, mean.k.)) + 
  geom_point() + theme_minimal(base_size = 7) +
  geom_hline(yintercept = 100, color="red") + 
  labs(title = 'Mean connectivity',
       x='Soft Threshold (power)', 
       y='Mean Connectivity')

cowplot::plot_grid(p1, p2, labels = c('A', 'B'), ncol=2)

sft_cut = min(sft$fitIndices[sft$fitIndices$SFT.R.sq >= 0.8 & sft$fitIndices$mean.k. <= 100,]$Power)
print (paste('The sft cutoff for this analysis is', sft_cut, sep=' '))

type <- 'signed'
minModuleSize <- 100
deepSplit <- 2

TOMsim <- TOMsimilarityFromExpr(t(input_mat), 
                                corType = "pearson", 
                                networkType = type, 
                                power = sft_cut, 
                                TOMType = "signed", 
                                nThreads = 2)

dissTOMA1 <- 1 - TOMsim
geneTree = hclust(as.dist(dissTOMA1), method="average")

# 상관계수를 토대로 상관관계가 높은 모듈을 선별
tree = cutreeHybrid(dendro = geneTree,
                    minClusterSize= minModuleSize,
                    pamStage=FALSE,
                    cutHeight = 0.999,
                    deepSplit= deepSplit,
                    distM= dissTOMA1 )

# 분류된 모듈 중에 상관관계가 비슷한 것들 통합
merged = mergeCloseModules(exprData= t(input_mat),
                           colors = tree$labels,
                           cutHeight=0.05,
                           verbose = 3)


res = setNames(as.data.frame(cbind(merged$colors)), c('Module_number'))
res$gene_name <- rownames(input_mat)

res$Merged_color = labels2colors(res$Module_number)


res %>% head()
table(res$Module_number)[2:16]
max(table(res$Module_number)[2:16])

res %>%select(-Merged_color) %>% pivot_wider(names_from = Module_number,values_from = gene_name)
writexl::write_xlsx(res,"co-expression.gene.xlsx")
#response to other organism
#detection of stimulus
'Homework
우리는 총 20개의 co-expression module을 구했습니다. 여러분은 이 정보를 바탕으로 다음의 질문에 답을 구하고 과제 파일에 적어서 제출하면 됩니다.

Q3. Module 중 가장 많은 유전자 갯수를 갖는 모듈의 번호는 무엇인가요?

Q4. 각 모듈은 주로 어떤 유전자를 갖고 있나요? Amigo 사이트에 들어가서, 유전자 리스트를 넣으면, 가장 p값이 낮은 pathway를 찾아서 과제 파일에 작성하면 됩니다.'