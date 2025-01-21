## Supplementary table 3
rm(list = ls())
library(openxlsx)
library(tidyverse)

setwd("~/Dropbox/SMC_AD_WGS_paper/Data/")

abbr = readxl::read_excel('~/Dropbox/SMC_AD_WGS_paper/Data/SMC_cwas_results_20240416/abbreviation_label.xlsx')
abbr = abbr[9:nrow(abbr),]

# Function to capitalize the first letter of each string
capitalize_first_letter <- function(x, low_other = F) {
  if(low_other){
    paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
  }else{
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  }
}

# all category
all_cate = read_tsv("SMC_cwas_results_20240416/burden_test/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.category_info.txt")
# FS after
select_cate = read_tsv("SMC_cwas_results_20240416/burden_test/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.category_info.txt")

enrich_cate = read_tsv("SMC_cwas_results_20240416/burden_test/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.burden_test.txt")
burden_cate = read_tsv("SMC_cwas_results_20240416/burden_shift/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.burdenshift_p0.05_cutoff10.txt") %>%
  filter(Category_set %in% c('is_coding', 'is_PTV', 'is_missense', "is_noncoding","is_promoter",   "is_intron",
                             "is_3primeUTR", "is_5primeUTR", "is_splice",
                             "is_intergenic"))
burden_cate$Category_set = gsub(x = burden_cate$Category_set,
                                pattern = 'is_',
                                replacement = '', fixed = T)
burden_cate$Category_set = capitalize_first_letter(x = burden_cate$Category_set, low_other = F)
burden_cate$Category_set = gsub(x = burden_cate$Category_set,
                                pattern = '_wo_',
                                replacement = ' without ',
                                fixed = T)
colnames(burden_cate)[1] = 'Genomic region'



###########
d0 = read.delim('SMC_cwas_results_20240416/burden_test/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.category_info.txt')

# cate1 = ("SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_noncoding_0.7_4_S107_CS.categorization_result.tsv")
d1 = read_csv('SMC_cwas_results_20240416/dawn/Kmeans_result/DAWN.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.intergenic.C10_R0.22_S2_L2_seed42.cluster_risk_pvalue_table.csv') %>%
  dplyr::select(-FDR)
d2 = read_csv("SMC_cwas_results_20240416/dawn/Kmeans_result/DAWN.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.intergenic.C10_R0.22_S2_L2_seed42.cluster_dawn_fdr_table.csv") %>% dplyr::rename(Cluster.idx = Name) %>%
  dplyr::select(-p.value)


# cluster_annotation.csv
d3 = read_csv("SMC_cwas_results_20240416/dawn/Kmeans_result/DAWN.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.intergenic.C10_R0.22_S2_L2_seed42.cluster_annotation.csv")

c1 = merge(d1,
           d2,
           by = 'Cluster.idx')
fdr = c1 %>%
  arrange(Cluster.idx) %>%
  dplyr::rename(cluster = Cluster.idx)

get_df = function(iteration, x){
  t0 = data.frame(cluster = iteration,
                  categories = x[,iteration])
  colnames(t0) = c('cluster', 'categories')
  t0 = t0 %>% filter(!is.na(categories)) %>% as.data.frame()
  return(t0)
}

res = lapply(unique(colnames(d3)), get_df, x = d3)
df = as.data.frame(do.call(rbind, res))

df_split <- df %>%
  separate(categories, into = c("variant_type", "gene_set", "functional_score", "gencode", "functional_annotation"), sep = "_",
           remove = F)

df_split2 = merge(df_split,
                  fdr %>%
                    dplyr::select(-indicator),
                  by = 'cluster')

write.table(df_split2,
            'SMC_cwas_results_20240416/DAWN_L2.cluster_list.txt',
            quote = F, row.names = F, col.names = T, sep = '\t')

# Grouping by 'cluster' and summarizing each component with unique comma-separated values
df_split$gene_set = ifelse(df_split$gene_set %in% abbr$Abbr,
                              paste('Path', df_split$gene_set, sep = '.'),
                              df_split$gene_set)
df_split = df_split %>%
  mutate(categories = paste(variant_type, gene_set, functional_score, gencode, functional_annotation,
                          sep = '_'))
df_new <- df_split %>%
  mutate(cluster = as.numeric(cluster)) %>% # Convert cluster to numeric
  group_by(cluster) %>%
  summarise(
    variant_type = toString(sort(unique(variant_type))),
    gene_set = toString(sort(unique(gene_set))),
    functional_score = toString(sort(unique(functional_score))),
    gencode = toString(sort(unique(gencode))),
    functional_annotation = toString(sort(unique(functional_annotation))),
    .groups = 'drop' # Optionally drop the grouping structure after summarising
  )


dawn_cluster = merge(df_new, fdr,
                     by = 'cluster') %>%
  select(-c(Pvalue, indicator)) %>%
  dplyr::rename(DAWN_RR = Risk , DAWN_FDR = FDR) %>%
  arrange(DAWN_FDR)


#cate = read_tsv(cate1) %>% dplyr::rename(sample = ...1)

TableS3 = createWorkbook()

# Table S4. Association Test for Rare Noncoding Variants with AD
# Key
openxlsx::addWorksheet(TableS3, "Table S3 Overview")
openxlsx::addWorksheet(TableS3, "Table S3A")
openxlsx::addWorksheet(TableS3, "Table S3B")
openxlsx::addWorksheet(TableS3, "Table S3C")
openxlsx::addWorksheet(TableS3, "Table S3D")
openxlsx::addWorksheet(TableS3, "Table S3E")
openxlsx::addWorksheet(TableS3, "Table S3F")
openxlsx::addWorksheet(TableS3, "Table S3G")
openxlsx::addWorksheet(TableS3, "Table S3H")

############  Table S3A	Annotation terms used in association test ########### 
S3A = rbind(data.frame(gene_set = unique(all_cate$gene_set)) %>%
              pivot_longer(gene_set),
            data.frame(functional_score = unique(all_cate$functional_score)) %>%
              pivot_longer(functional_score),
            data.frame(functional_annotation = unique(all_cate$functional_annotation)) %>%
              pivot_longer(functional_annotation)) %>%
  dplyr::rename(annotation_type = name, annotation_name = value)

S3A$annotation_name = ifelse(S3A$annotation_name %in% abbr$Abbr, paste('Path', S3A$annotation_name,
                                                                       sep = '.'),
                             S3A$annotation_name)

openxlsx::writeData(wb = TableS3, sheet = "Table S3A", x = S3A, startCol = 1, startRow = 1, rowNames = F)


############  Table S3B	All categories from rare variants ########### 

all_cate = read_tsv('SMC_cwas_results_20240416/burden_test/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.burden_test.txt') %>% select(-"Relative_Risk"   , -"P"  , -"P_1side",-"Z_1side"       )

openxlsx::writeData(wb = TableS3, sheet = "Table S3B", x = all_cate, startCol = 1, startRow = 1, rowNames = F)

############  Table S3C	Risk score from selected categories ########### 

res_cs = read_tsv("SMC_cwas_results_20240416/risk_score_after_fs/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.lasso_results_tf0.8_thres_3.txt") %>% filter(Domain != "all") %>% filter(Seed == 'average') %>% select(-Seed)


openxlsx::writeData(wb = TableS3, sheet = "Table S3C", x = res_cs, startCol = 1, startRow = 1, rowNames = F)


########### Table S3D	Single category enrichment test result ########### 
enrich_cate$gene_set = ifelse(enrich_cate$gene_set %in% abbr$Abbr,
                              paste('Path', enrich_cate$gene_set, sep = '.'),
                              enrich_cate$gene_set)
enrich_cate = enrich_cate %>%
  mutate(Category = paste(variant_type, gene_set, functional_score, gencode, functional_annotation,
                          sep = '_'))
openxlsx::writeData(wb = TableS3, sheet = "Table S3D", x = enrich_cate, startCol = 1, startRow = 1, rowNames = F)

###########  Table S3E	Comparison results of the significant count for each category by genomic region ########### 
openxlsx::writeData(wb = TableS3, sheet = "Table S3E", x = burden_cate, startCol = 1, startRow = 1, rowNames = F)

############  Table S3F	Clusters from sub-network analysis ########### 

openxlsx::writeData(wb = TableS3, sheet = "Table S3F", x = dawn_cluster, startCol = 1, startRow = 1, rowNames = F)


##### Table S3G	Variants from risk clusters for all samples ################

## DAWN
pheno = read_tsv("pheno_SMC_batch1-7_RCL40_visual_1559.tsv") %>% 
  select(-FAMILY)

dawn_rnv_list = read_csv("SMC_cwas_results_20240416/DAWN_L2.RNV_list.csv") %>%
  separate(Variant, into = c("chr", "locus", "ref", "alt"), sep = ":",
           remove = F)

dawn_rnv_list = dawn_rnv_list %>%
  select(cluster, Variant, GENE, `Gene distance`= gene_dist, Variant_locus, SAMPLE, PHENOTYPE) %>%
  unique()

colnames(dawn_rnv_list) = gsub(x = colnames(dawn_rnv_list),
                               pattern = '_',
                               replacement = ' ',
                               fixed = T)
# Apply the function to the vector
colnames(dawn_rnv_list) <- sapply(colnames(dawn_rnv_list), capitalize_first_letter, low_other=T)
dawn_rnv_list = dawn_rnv_list %>% arrange(Cluster)
dawn_rnv_list_final = dawn_rnv_list %>% distinct()
dawn_rnv_list_final$Phenotype = ifelse(dawn_rnv_list_final$Phenotype=='case',
                                       'Case',
                                       'Control')

openxlsx::writeData(wb = TableS3, sheet = "Table S3G", x = dawn_rnv_list_final, startCol = 1, startRow = 1, rowNames = F)



##### Table S3H	Unique variants from risk clusters ################
DAWN_RNV_annot = read_csv("SMC_cwas_results_20240416/DAWN_L2.RNV_annot.csv") %>% arrange(cluster)
DAWN_RNV_annot = merge(DAWN_RNV_annot,
                       dawn_rnv_list_final %>%
                         dplyr::select(Variant, `Variant locus`,
                                       `Gene distance`),
                       by = 'Variant',
                       all.x = T)
colnames(DAWN_RNV_annot) = gsub(x = colnames(DAWN_RNV_annot),
                                pattern = '2023.',
                                replacement = '',
                                fixed = T)
colnames(DAWN_RNV_annot) = ifelse(colnames(DAWN_RNV_annot)=='Xiong.CRE.Ex',
                                  'Xiong.CRE.Exc',
                                  colnames(DAWN_RNV_annot))
colnames(DAWN_RNV_annot) = ifelse(colnames(DAWN_RNV_annot)=='Xiong.CRE.In',
                                  'Xiong.CRE.Inh',
                                  colnames(DAWN_RNV_annot))
colnames(DAWN_RNV_annot) = gsub(x = colnames(DAWN_RNV_annot),
                                pattern = '.',
                                replacement = ' ',
                                fixed = T)
colnames(DAWN_RNV_annot) = gsub(x = colnames(DAWN_RNV_annot),
                                pattern = '_',
                                replacement = ' ',
                                fixed = T)
DAWN_RNV_annot = DAWN_RNV_annot %>%
  dplyr::select(-'Xiong CRE Vas')
colnames(DAWN_RNV_annot) <- sapply(colnames(DAWN_RNV_annot), capitalize_first_letter,
                                   low_other = F)
DAWN_RNV_annot = DAWN_RNV_annot %>%
  distinct() %>%
  dplyr::rename('Gene' = 'GENE')


## Sample count by variant & phenotype
tt = dawn_rnv_list_final %>%
  group_by(Cluster, Variant,
           Phenotype) %>%
  dplyr::count()
dawn_rnv_count = tt %>%
  spread(key = Phenotype, value = n)
dawn_rnv_count$Case = ifelse(is.na(dawn_rnv_count$Case), 0, dawn_rnv_count$Case)
dawn_rnv_count$Control = ifelse(is.na(dawn_rnv_count$Control), 0, dawn_rnv_count$Control)
colnames(dawn_rnv_count) = c('Cluster', 'Variant', 'case_count', 'control_count')


Table_S3H = DAWN_RNV_annot %>% left_join(dawn_rnv_count) %>%
  dplyr::select('Cluster', 'Variant',
                'Xiong CRE Ast', 'Xiong CRE Exc', 'Xiong CRE Inh', 'Xiong CRE Micro', 'Xiong CRE Oligo', 'Xiong CRE OPC', 'Sun CRE Micro',
                'Gene', 'Gene distance', 'Variant locus', 'Cell type marker', 'GO term',
                'case_count', 'control_count')

openxlsx::writeData(wb = TableS3, sheet = "Table S3H", x = Table_S3H, startCol = 1, startRow = 1, rowNames = F)

saveWorkbook(wb = TableS3, file = "../Tables/Supplementary Table 3 v4.1.xlsx", overwrite=T)



