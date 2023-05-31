# devtools
#system('sudo apt-get install -y libfribidi-dev', intern=TRUE)
#install.packages("devtools") #�⺻������ r library ��ġ

# tidyverse and DT libraries
#install.packages(c("tidyverse", "DT"))

# TCGAbiolinks to download RNA-seq data
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks", ask=F)

# DESEQ2 ���� ���̺귯�� ��ġ
BiocManager::install("DESeq2", ask=F)

# WGCNA ���� ���̺귯�� ��ġ
BiocManager::install("WGCNA", ask=F)
BiocManager::install("limma", ask=F)
install.packages("cowplot")

# GSEA ���� ���̺귯�� ��ġ
install.packages("msigdbr")
install.packages(c("tidyverse", "DT"))

BiocManager::install("fgsea", ask=F)


# �ʿ��� ���̺귯�� �ҷ�����
library(TCGAbiolinks)
library(tidyverse)
library(DT)
library(DESeq2)
library(WGCNA)


# �м��� ����� ������ ���� ����
query <- GDCquery(
  project = "TCGA-LUSC",
  data.type = "Gene Expression Quantification", 
  data.category = "Transcriptome Profiling",
  workflow.type = "STAR - Counts"
)

# �ٿ�ε�
# �ٿ�ε忡 ���� �ð��� �ҿ�ɼ��� �ֽ��ϴ�
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
# ������ DESEQ2 �м��� ����, �Ҽ��� ������ �����ϰڽ��ϴ�
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

# DESEQ2 ��ü �����. DESEQ2 �м��� ���ؼ�, gene expression, sample information�� ��� �ϳ��� ��ü�� �����ϰ� ����� �ؾ��մϴ�.

dds <- DESeqDataSetFromMatrix(countData = e1,
                              colData = s1,
                              design= ~ is_tumor)

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

# ����� �� �ִ� intercept (�񱳴��)�� 'is_tumor_1_vs_0'�� �ش��մϴ�. �� �÷� ������ �̿��Ͽ�, �������� ���캸�ڽ��ϴ�. 

res <- results(dds, name="is_tumor_1_vs_0")

# res ��ü�� ���������������� �����մϴ�. �̶� gene information�� ���� �߰��մϴ�. 
res1 = data.frame(gene_id = g$gene_id, 
                  gene_name = g$gene_name, 
                  baseMean=res$baseMean, 
                  lfc=res$log2FoldChange, 
                  lfcSE=res$lfcSE, p=res$pvalue, 
                  padj=res$padj) %>% as_tibble()

# � �����ڵ��� ũ�� ��ȭ�Ͽ����� ���캾�ô�.
# lfc �� log2 fold change�ν� tumor�� normal �׷쿡 ���� ������ ���̸� ��Ÿ�� �÷��Դϴ�
# Deseq2�� wald test�� �̿��Ͽ� p���� �����ϰ�, �̿� ���� ���ߺ� ���� padj �÷����� �����մϴ�. �ڼ��� ������ �Ʒ� ��ũ�� �����ϼ���
# https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html

res1 %>% arrange(padj) %>% head()

# ���� ����� �ð�ȭ �ϰڽ��ϴ�
res1 %>% mutate(isSig = ifelse(padj < 0.05 & lfc > 0, "Up", "None"),
                isSig = ifelse(padj < 0.05 & lfc < 0, "Down", isSig)) %>%
  ggplot(aes(lfc, -log10(padj), fill=isSig)) + 
  geom_point(shape=21) + 
  scale_fill_manual(values=c("Up"="firebrick", "Down"="dodgerblue", "None"="grey"))



# msigdb download
library(msigdbr)
df_msigdb = msigdbr(species = "Homo sapiens")

# Hallmark pathway �� �Ế�ô�
df_msigdb_hallmarks = df_msigdb %>% filter(gs_cat=="H")
msigdbr_list = split(x = df_msigdb_hallmarks$gene_symbol, f = df_msigdb_hallmarks$gs_name)
# fgsea�� �̿��Ͽ� GSEA �м�
library(fgsea)
input_ranks = res1 %>% filter(!is.na(lfc) & !is.na(padj)) %>% arrange(-abs(lfc)) %>% 
  filter(!duplicated(gene_name)) %>% pull(lfc, name=gene_name)

fgseaRes <- fgsea(msigdbr_list, input_ranks)
fgseaRes %>% as_tibble() %>% arrange(pval) %>% head(10)

plotEnrichment(msigdbr_list[["HALLMARK_E2F_TARGETS"]], input_ranks) + labs(title="HALLMARK_E2F_TARGETS")

'
Homework
���� ������ DESEQ2 ����� �������� �����ڿ� �������ڸ� �����Ͽ� DEG �м��� �����ϰڽ��ϴ�. �Ʒ� ������ ����Ͽ� �����ϸ� �˴ϴ�.

������ ���� ���� �����Ͻø� �˴ϴ�. ���� �������Ͽ� �����Ͻø� �˴ϴ�.

Q1. ��� DEG (pdaj < 0.05 �̸�) �����ڵ��� ���ɴϱ�?
  Q2. Hallmark GSEA�� �����ϸ�, ������ �׷쿡�� ���� up-regulation �Ǵ� pathway �� �����ΰ���? GSEA plot�� ���� �������ּ���.
  '

table(s$paper_Smoking.Status)
# ������ ������ ���� �Ҽ��� ������ �����ϰڽ��ϴ�
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
# DESEQ2 ��ü �����. DESEQ2 �м��� ���ؼ�, gene expression, sample information�� ��� �ϳ��� ��ü�� �����ϰ� ����� �ؾ��մϴ�.

dds <- DESeqDataSetFromMatrix(countData = e1,
                              colData = s1,
                              design= ~is_smoker)

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

# ����� �� �ִ� intercept (�񱳴��)�� 'is_tumor_1_vs_0'�� �ش��մϴ�. �� �÷� ������ �̿��Ͽ�, �������� ���캸�ڽ��ϴ�. 

res <- results(dds, name="is_smoker_1_vs_0")

# res ��ü�� ���������������� �����մϴ�. �̶� gene information�� ���� �߰��մϴ�. 
res1 = data.frame(gene_id = g$gene_id, 
                  gene_name = g$gene_name, 
                  baseMean=res$baseMean, 
                  lfc=res$log2FoldChange, 
                  lfcSE=res$lfcSE, p=res$pvalue, 
                  padj=res$padj) %>% as_tibble()

# � �����ڵ��� ũ�� ��ȭ�Ͽ����� ���캾�ô�.
# lfc �� log2 fold change�ν� tumor�� normal �׷쿡 ���� ������ ���̸� ��Ÿ�� �÷��Դϴ�
# Deseq2�� wald test�� �̿��Ͽ� p���� �����ϰ�, �̿� ���� ���ߺ� ���� padj �÷����� �����մϴ�. �ڼ��� ������ �Ʒ� ��ũ�� �����ϼ���
# https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html

res1 %>% arrange(padj) %>% head()
res1 %>% filter(padj<0.05) %>% nrow()
res11<-res1 %>% mutate(isSig = ifelse(padj < 0.05 & lfc > 0, "Up", "None"),
                isSig = ifelse(padj < 0.05 & lfc < 0, "Down", isSig)) 
table(res11$isSig)
# ���� ����� �ð�ȭ �ϰڽ��ϴ�
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

# Hallmark pathway �� �Ế�ô�
df_msigdb_hallmarks = df_msigdb %>% filter(gs_cat=="H")
msigdbr_list = split(x = df_msigdb_hallmarks$gene_symbol, f = df_msigdb_hallmarks$gs_name)
# fgsea�� �̿��Ͽ� GSEA �м�
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
# ������ DESEQ2 �м��� ����, �Ҽ��� ������ �����ϰڽ��ϴ�
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

# �������� ���� ������谡 ���� ����� ����
tree = cutreeHybrid(dendro = geneTree,
                    minClusterSize= minModuleSize,
                    pamStage=FALSE,
                    cutHeight = 0.999,
                    deepSplit= deepSplit,
                    distM= dissTOMA1 )

# �з��� ��� �߿� ������谡 ����� �͵� ����
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
�츮�� �� 20���� co-expression module�� ���߽��ϴ�. �������� �� ������ �������� ������ ������ ���� ���ϰ� ���� ���Ͽ� ��� �����ϸ� �˴ϴ�.

Q3. Module �� ���� ���� ������ ������ ���� ����� ��ȣ�� �����ΰ���?

Q4. �� ����� �ַ� � �����ڸ� ���� �ֳ���? Amigo ����Ʈ�� ����, ������ ����Ʈ�� ������, ���� p���� ���� pathway�� ã�Ƽ� ���� ���Ͽ� �ۼ��ϸ� �˴ϴ�.'