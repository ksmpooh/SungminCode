library(tidyverse)

version = 'v2.1'

extrafont::font_import(pattern = "Arial", prompt = F)
extrafont::loadfonts()

setwd('~/Dropbox/SMC_AD_WGS_paper/Data/SMC_cwas_results_20240416/')
# Set colors

# Function definition
theme_step1 <- function(base_size = 11, base_family = "",
                        base_line_size = base_size / 22,
                        base_rect_size = base_size / 22) {
  theme(title = element_text(family = 'Arial', size = 18, color = 'black'), text = element_text(family = 'Arial', size = 16, color = 'black'),
        axis.title = element_text(family = 'Arial', size = 18, color = 'black'), axis.text = element_text(family = 'Arial', size = 16, color = 'black'), 
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = NA), axis.line = element_line(colour = "black", size = rel(1)),
        legend.background = element_rect(color = 'black'), legend.title = element_text(family = 'Arial', size = 16),
        legend.text = element_text(family = 'Arial', size = 14),
        legend.direction = "vertical", 
        legend.box = c("horizontal", "vertical"),
        legend.spacing.x = unit(0.1, 'cm'),
        plot.margin = unit(c(0.25, 1, 1, 0.5), 'cm'),
        axis.title.y = element_text(margin = margin(r = 10, unit = "pt"))) }

fn_DEG <- "~/Desktop/AD_WGS/Annotation/Matyhs_2023_DecontX_DEG.txt.gz"
df_DEG_orig <- data.table::fread(fn_DEG)
df_DEG = df_DEG_orig

df_DEG = df_DEG %>% filter(!cluster_id %in% c('SMC', 'Fib', 'CAMs', 'T cells',
                                              'Exc' ,'Inh',  "Inh LAMP5", "Inh PAX6" , "Inh PVALB" ,"Inh SST"  , "Inh VIP", "Exc IT"  ))

create_phenotype_column <- function(df) {
  
  # 항목 이름과 매칭되는 괄호 내용을 포함한 리스트
  matching_list <- list(
    "amyloid" = "Overall amyloid level",
    "Apoe_e4" = "ApoeE genotype",
    "arteriol_scler" = "Arteriolosclerosis",
    "bradysc_lv" = "Bradykinesia score",
    "caa_4gp" = "Cerebral amyloid angiopathy",
    "cancer_bl" = "Cancer at baseline",
    "chd_cogact_freq" = "Cognitive activity - childhood",
    "ci_num2_gct" = "Gross chronic infarcts",
    "ci_num2_gtt" = "Gross infarcts",
    "ci_num2_mct" = "Chronic microinfarcts",
    "ci_num2_mtt" = "Microinfarcts",
    "cogdx" = "Final cognitive consensus diagnosis",
    "cogdx_stroke" = "Clinical diagnosis of stroke",
    "cogn_ep_lv" = "Episodic memory",
    "cogn_global_lv" = "Global cognitive function",
    "cogn_po_lv" = "Perceptual orientation",
    "cogn_ps_lv" = "Perceptual speed",
    "cogn_se_lv" = "Semantic memory",
    "cogn_wo_lv" = "Working memory",
    "cognep_random_slope" = "Rate of change of episodic memory",
    "cogng_random_slope" = "Rate of change of global cognitive function",
    "cognpo_random_slope" = "Rate of change of perceptual orientation",
    "cognps_random_slope" = "Rate of change of perceptual speed",
    "cognse_random_slope" = "Rate of change of semantic memory",
    "cognwo_random_slope" = "Rate of change of working memory",
    "cvda_4gp2" = "Cerebral atherosclerosis",
    "diabetes_sr_rx_bl" = "Diabetes",
    "dlbdx" = "Lewy Body disease",
    "dxpark" = "Clinical diagnosis of Parkinson’s disease",
    "gaitsc_lv" = "Gait score",
    "gpath" = "Global AD pathology",
    "gpath_3neocort" = "Global measure of neocortical AD pathology",
    "gpath_CDR_score" = "Global AD pathology cognitive resilience score",
    "headinjrloc_bl" = "Head injury",
    "heart_bl" = "Heart",
    "hypertension_bl" = "Hypertension",
    "lifetime_cogact_freq_bl" = "Cognitive activity - lifetime (total)",
    "ma_adult_cogact_freq" = "Cognitive activity - middel age",
    "msex" = "Sex",
    "nft" = "Neurofibrillary tangle burden",
    "nft_CDR_score" = "Neurofibrillary tangle burden CDR score",
    "nft_CR_score" = "Neurofibrillary tangle burden CR score",
    "nft_mf" = "Neurofibrillary tangle burden in the midfrontal cortex",
    "parksc_lv" = "Parkinsonian summary score",
    "phys5itemsum_bl" = "Physical activity - baseline",
    "phys5itemsum_lv" = "Physical activity - last valid",
    "plaq_d" = "Diffuse plaque burden",
    "plaq_d_mf" = "Diffuse plaque burden in the midfrontal cortex",
    "plaq_n" = "Neuritic plaque burden",
    "plaq_n_CDR_score" = "Neuritic plaque burden CDR score",
    "plaq_n_CR_score" = "Neuritic plaque burden CR score",
    "plaq_n_mf" = "Neuritic plaque burden in the midfrontal cortex",
    "soc_net_bl" = "Social network - baseline",
    "social_isolation_avg" = "Social isolation score - average",
    "social_isolation_lv" = "Social isolation score - last valid",
    "stroke_bl" = "Stroke",
    "tangles" = "Tangle density",
    "tangles_CDR_score" = "Tangle density CDR score",
    "tangles_CR_score" = "Tangle density CR score",
    "tdp_st4" = "TDP-43 stage",
    "ya_adult_cogact_freq" = "Cognitive activity - young adult "
  )
  
  # 매칭되는 항목을 찾아 phenotype 값 설정
  matched_indices <- match(df$coef, names(matching_list))
  
  # 매칭되는 항목이 없는 경우에는 coef 값을 그대로 사용하도록 설정
  df$phenotype <- ifelse(is.na(matched_indices), df$coef, unlist(matching_list)[matched_indices])
  
  return(df)
}

# 새로운 phenotype 열 초기화
df_DEG$phenotype <- ""
# 함수 호출하여 phenotype 열 생성
df_DEG <- create_phenotype_column(df_DEG)

# 주어진 Phenotype 그룹에 따라 그룹명을 설정하는 함수
get_group_name <- function(phenotype) {
  # AD pathology
  ad_pathology <- c("nft", "nft_mf", "tangles", "niareagansc", "gpath_3neocort", 
                    "cogdx", "apoe_genotype", "gpath", "plaq_d_mf", "braaksc", 
                    "ceradsc", 
                    # "msex",
                    "amyloid", "plaq_d", "plaq_n", "plaq_n_mf", 
                    "gpath_CR_score","Apoe_e4",
                    'nft_CDR_score', 'nft_CR_score', "plaq_n_CDR_score",  "plaq_n_CR_score"  ,"tangles_CDR_score", "tangles_CR_score", 'gpath_CDR_score')
  # Lewy body / TDP
  lewy_tdp <- c("tdp_st4", "dlbdx")
  # Vascular pathology
  vascular_pathology <- c("arteriol_scler", "caa_4gp", "cvda_4gp2", "ci_num2_gct", 
                          "ci_num2_mct", "ci_num2_gtt", "ci_num2_mtt")
  # Medical conditions
  medical_conditions <- c("hypertension_bl", "cancer_bl", "diabetes_sr_rx_bl", 
                          "headinjrloc_bl", "heart_bl", "stroke_bl", "cogdx_stroke")
  # Motor and gait
  motor_gait <- c("bradysc_lv", "gaitsc_lv", "parksc_lv", "dxpark")
  # Lifestyle
  lifestyle <- c("cogn_global_lv", "cogn_ep_lv", "cogn_po_lv", "cogn_ps_lv", 
                 "cogn_se_lv", "cogn_wo_lv", "cognep_random_slope", 
                 "cogng_random_slope", "cognpo_random_slope", "cognps_random_slope", 
                 "cognse_random_slope", "cognwo_random_slope", "age_death_CR_score",'chd_cogact_freq', 'lifetime_cogact_freq_bl','ma_adult_cogact_freq', "phys5itemsum_bl","phys5itemsum_lv" ,'soc_net_bl','social_isolation_avg',  "soc_net_bl","social_isolation_avg"  ,"social_isolation_lv"  ,   "ya_adult_cogact_freq"   )
  
  # 그룹에 속하는지 확인 후 그룹명 반환
  if (phenotype %in% ad_pathology) {
    return("AD pathology")
  } else if (phenotype %in% lewy_tdp) {
    return("Lewy / TDP")
  } else if (phenotype %in% vascular_pathology) {
    return("Vascular pathology")
  } else if (phenotype %in% medical_conditions) {
    return("Medical conditions")
  } else if (phenotype %in% motor_gait) {
    return("Motor and gait")
  } else if (phenotype %in% lifestyle) {
    return("Lifestyle")
  } else {
    return("Others")
  }
}

# phenotype_group 열을 추가
df_DEG$phenotype_group <- sapply(df_DEG$coef, get_group_name)

# 결과 확인
head(df_DEG)

df_DEG %>% filter(phenotype_group == 'Others') %>% pull(coef) %>% unique()


df_DEG$phenotype %>% unique() %>% length
df_DEG$coef %>% unique() %>% length


gene_list = c('NTS', 'PMAIP1', 'GRP', 'DTD2', 'SNORA89', 'HECTD1', 'COCH', 'ARHGAP5', 'AP4S1', 'HECTD1', 'COCH', 'ARHGAP5', 'AP4S1', 'PACRG', 'NOL4', 'LOC107985157', 'LOC105372058', 'PMAIP1', 'GRP', 'VAPB', 'TUBB1', 'STX16', 'ATP5F1E', 'TBCCD1', 'RTP4', 'NMRAL2P', 'MIR1248', 'RFC4', 'SNORA4', 'SNORA63', 'SNORA63B', 'SNORA81', 'LPP', 'IGF2BP2-AS1', 'IGF2BP2', 'EPHB3', 'EIF4A2', 'SNORD2', 'CRYGS', 'SLC1A3', 'RANBP3L', 'NUP155', 'GDNF-AS1', 'CAPSL', 'SAMD5', 'RAET1K', 'C6orf118', 'QKI', 'CAHM', 'PACRG-AS1', 'PACRG', 'C6orf118', 'PDE10A', 'PACRG-AS1', 'SLC1A1', 'PUM3', 'PLGRKT', 'INSL6')

gene_list = c('APOE','APP', 'VGF')

# gene_list = c('MIR3123', 'MBL1P', 'NTS', 'LOC105369877', 'LOC107985156', 'DTD2', 'LOC105370438', 'LOC105378100', 'LOC105369171','C6orf118', 'VAPB', 'LOC107987284', 'ATP5F1E', 'MIR1248', 'RFC4', 'SNORA4', 'SNORA63', 'SNORA63B', 'SNORA81', 'IGF2BP2', 'EIF4A2', 'SNORD2', 'SLC1A3', 'LOC105378043', 'QKI', 'LOC105378115', 'LOC105378104', 'CAHM', 'PACRG-AS1', 'LOC72968', 'LOC105378116')


df_DEG_t = df_DEG %>% filter(gene %in% gene_list)

df_DEG_t$significance_label <- ifelse(df_DEG_t$p_adj.loc < 0.001 & abs(df_DEG_t$logFC) > 0.02, "***",
                                      ifelse(df_DEG_t$p_adj.loc < 0.01 & abs(df_DEG_t$logFC) > 0.02, "**",
                                             ifelse(df_DEG_t$p_adj.loc < 0.05 & abs(df_DEG_t$logFC) > 0.02, "*", "")))

# df_DEG_t$significance_label <- ifelse(df_DEG_t$p_adj.loc < 0.001 & abs(df_DEG_t$logFC) > 0.02, "*", "")


unique_values <- c(
  'Exc L2-3 CBLN2 LINC02306', 'Exc L3-4 RORB CUX2', 'Exc L3-5 RORB PLCH1', 'Exc L4-5 RORB IL1RAPL2', 'Exc L4-5 RORB GABRG1', 'Exc L5-6 RORB LINC02196', 'Exc L6 THEMIS NFIA', "Exc L5_6 IT Car3", "Exc L5 ET", "Exc L5_6 NP", "Exc L6 CT", "Exc L6b", "Exc NRGN","Exc RELN CHD7", 'Inh L1-6 LAMP5 CA13', 'Inh LAMP5 NRG1 (Rosehip)', 'Inh LAMP5 RELN', 'Inh L1 PAX6 CA4', 'Inh L1-2 PAX6 SCGN', 'Inh GPC5 RIT2', 'Inh L5-6 PVALB STON2', 'Inh PVALB CA8 (Chandelier)', 'Inh PVALB HTR4', 'Inh PVALB SULF1', 'Inh CUX2 MSR1', 'Inh ENOX2 SPHKAP', 'Inh FBN2 EPB41L4A', 'Inh L3-5 SST MAFB', 'Inh L5-6 SST TH', 'Inh L6 SST NPY', 'Inh ALCAM TRPM3', 'Inh PTPRK FAM19A1', 'Inh RYR3 TSHZ2', 'Inh SGCD PDE3A', 'Inh SORCS1 TTN', 'Inh VIP ABI3BP', 'Inh VIP CLSTN2', 'Inh VIP THSD7B', 'Inh VIP TSHZ2',
  "Ast","Oli","OPC","Mic","End","Per"     )

# 요인(factor)으로 변환
df_DEG_t$cluster_id <- factor(df_DEG_t$cluster_id, levels = unique_values)

df_DEG_t <- df_DEG_t %>% mutate(inh_subcalss = case_when(
  cluster_id %in% c('Exc L2-3 CBLN2 LINC02306', 'Exc L3-4 RORB CUX2', 'Exc L3-5 RORB PLCH1', 'Exc L4-5 RORB IL1RAPL2', 'Exc L4-5 RORB GABRG1', 'Exc L5-6 RORB LINC02196', 'Exc L6 THEMIS NFIA') ~ 'IT',
  cluster_id %in% c("Exc L5_6 IT Car3") ~ 'L5 IT Car3',
  cluster_id %in% c("Exc L5 ET") ~ 'L5 ET',
  cluster_id %in% c("Exc L5_6 NP") ~ 'Exc L5/6 NP',
  cluster_id %in% c("Exc L6 CT") ~ 'L6 CT',
  cluster_id %in% c("Exc L6b") ~ 'L6b',
  cluster_id %in% c("Exc NRGN","Exc RELN CHD7") ~ 'L4 IT',
  cluster_id %in% c('Inh L1-6 LAMP5 CA13', 'Inh LAMP5 NRG1 (Rosehip)', 'Inh LAMP5 RELN') ~ 'LAMP5',
  cluster_id %in% c('Inh L1 PAX6 CA4',
                    'Inh L1-2 PAX6 SCGN') ~ 'PAX6',
  cluster_id %in% c('Inh GPC5 RIT2', 'Inh L5-6 PVALB STON2', 'Inh PVALB CA8 (Chandelier)', 'Inh PVALB HTR4', 'Inh PVALB SULF1') ~ 'PVALB',
  cluster_id %in% c('Inh CUX2 MSR1', 'Inh ENOX2 SPHKAP', 'Inh FBN2 EPB41L4A', 'Inh L3-5 SST MAFB', 'Inh L5-6 SST TH', 'Inh L6 SST NPY') ~ 'SST',
  cluster_id %in% c('Inh ALCAM TRPM3', 'Inh PTPRK FAM19A1', 'Inh RYR3 TSHZ2', 'Inh SGCD PDE3A', 'Inh SORCS1 TTN', 'Inh VIP ABI3BP', 'Inh VIP CLSTN2', 'Inh VIP THSD7B', 'Inh VIP TSHZ2') ~ 'VIP',
  
  T ~ 'Non-neuronal'
))

# df_DEG %>% filter(inh_subtype == 'others') %>% pull(cluster_id) %>% unique()

sub_class <- c('IT', 'L5 IT Car3', 'L5 ET', 'Exc L5/6 NP', 'L6 CT', 'L6b', 'L4 IT', 'LAMP5', 'PAX6', 'PVALB', 'SST', 'VIP', 'Non-neuronal'  )

# 요인(factor)으로 변환
df_DEG_t$inh_subcalss <- factor(df_DEG_t$inh_subcalss, levels = sub_class)


ex_fig1_a = df_DEG_t %>% 
  filter(inh_subcalss %in% c('IT', 'L5 IT Car3', 'L5 ET', 'Exc L5/6 NP', 'L6 CT', 'L6b', 'L4 IT')) %>%
  mutate(logFC = ifelse(logFC > 1, 1,
                        ifelse(logFC < -1, -1 , logFC) )) %>%
  ggplot(aes(x = gene, 
             y = phenotype,
             # y = gene, 
             fill = logFC)) +
  geom_tile(color = "grey", lwd = 1.5, linetype = 1) +
  geom_text(aes(label = significance_label, fontface = "bold", angle = 90), color = "black", size = 5, hjust = 0.5, vjust =0.8,family = 'Arial') +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_gradient2(low = "#009E73", mid = "white", high = "#D55E00"# , limits = c(-3,3)
  ) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 5)) +
  facet_nested(
    phenotype_group  ~  inh_subcalss  + cluster_id,  #
    scales = "free",
    space = "free",
    switch = 'both'#, nest_line = element_line(linetype = 2)
  ) +
  theme_step1() +
  xlab("") +
  ylab("") + 
  labs(title = "") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour = "black", fill = "white", 
                                        size = 1.5, linetype = "solid"))

ex_fig1_b = df_DEG_t %>% 
  filter(inh_subcalss %in% c('LAMP5', 'PAX6', 'PVALB', 'SST', 'VIP')) %>%
  mutate(logFC = ifelse(logFC > 1, 1,
                        ifelse(logFC < -1, -1 , logFC) )) %>%
  ggplot(aes(x = gene, 
             y = phenotype,
             # y = gene, 
             fill = logFC)) +
  geom_tile(color = "grey", lwd = 1.5, linetype = 1) +
  geom_text(aes(label = significance_label, fontface = "bold", angle = 90), color = "black", size = 5, hjust = 0.5, vjust =0.8,family = 'Arial') +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_gradient2(low = "#009E73", mid = "white", high = "#D55E00"# , limits = c(-3,3)
  ) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 5)) +
  facet_nested(
    phenotype_group  ~  inh_subcalss  +cluster_id, 
    scales = "free",
    space = "free",
    switch = 'both'#, nest_line = element_line(linetype = 2)
  ) +
  theme_step1() +
  xlab("") +
  ylab("") + 
  labs(title = "") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour = "black", fill = "white", 
                                        size = 1.5, linetype = "solid"))

ex_fig1_c = df_DEG_t %>% 
  filter(inh_subcalss %in% c('Non-neuronal')) %>%
  mutate(logFC = ifelse(logFC > 1, 1,
                        ifelse(logFC < -1, -1 , logFC) )) %>%
  ggplot(aes(x = gene, 
             y = phenotype,
             # y = gene, 
             fill = logFC)) +
  geom_tile(color = "grey", lwd = 1.5, linetype = 1) +
  geom_text(aes(label = significance_label, fontface = "bold", angle = 90), color = "black", size = 5, hjust = 0.5, vjust =0.8,family = 'Arial') +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_gradient2(low = "#009E73", mid = "white", high = "#D55E00"# , limits = c(-3,3)
  ) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 5)) +
  facet_nested(
    phenotype_group  ~  inh_subcalss  + cluster_id, 
    scales = "free",
    space = "free",
    switch = 'both'#, nest_line = element_line(linetype = 2)
  ) +
  theme_step1() +
  xlab("") +
  ylab("") + 
  labs(title = "") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour = "black", fill = "white", 
                                        size = 1.5, linetype = "solid"))



ex_fig1 = ggarrange(ex_fig1_a, ex_fig1_b,ex_fig1_c, nrow = 1, labels = c('a','b','c'), font.label = list(size = 28), label.y = 1.01,
          widths = c(1.2,1.4,0.5))

ggsave(paste0('~/Dropbox/SMC_AD_WGS_paper/Figures/Extended Data Figure/Extended_figure3_a_',version,'.pdf'), ex_fig1_a, width = 120,height = 30, limitsize = FALSE)
ggsave(paste0('~/Dropbox/SMC_AD_WGS_paper/Figures/Extended Data Figure/Extended_figure3_b_',version,'.pdf'), ex_fig1_b, width = 120,height = 30, limitsize = FALSE)
ggsave(paste0('~/Dropbox/SMC_AD_WGS_paper/Figures/Extended Data Figure/Extended_figure3_c_',version,'.pdf'), ex_fig1_c, width = 60,height = 30, limitsize = FALSE)
ggsave(paste0('~/Dropbox/SMC_AD_WGS_paper/Figures/Extended Data Figure/Extended_figure3_AD_gene_',version,'.pdf'), ex_fig1, width = 100,height = 30, limitsize = FALSE)

# ggsave("~/Dropbox/ADWGS/CWAS/Figures/C73_Nearest_gene_v10.1.pdf", p, width = 100,height = 25, limitsize = FALSE)



{
  dawn_c = "Cluster73"
  
  list = read_tsv("table.DAWN_L2.RNV_carrier_list.txt")
  sig_cluster = list %>% select(-1,-2,-3) %>% colnames()
  
  print(paste0("",dawn_c))
  
  pheno_cwas = list
  
  pheno_age = readxl::read_xlsx("~/Dropbox/ADWGS/sample_information/WGS_1824_phenotype_update.231218.xlsx",sheet = 6) %>% select(WGS_SerialNo., age)
  
  pheno = read.csv("~/Dropbox/ADWGS/sample_information/WGS_label_240206_MR_Cog_add.csv") %>% mutate(sex = ifelse(sex == "M", 1,0)) %>% filter(RCL40_visual %in% c(1,0) & DX %in% c("DAT","MCI","CU") & WGS_SerialNo. %in% pheno_cwas$SAMPLE) %>% left_join(pheno_age)
  
  pheno$sex <- as.numeric(pheno$sex)
  pheno$KlunkCL <- as.numeric(pheno$KlunkCL)
  pheno$RCL_global <- as.numeric(pheno$RCL_global)
  
  # case 만 추출
  case = pheno %>% filter(RCL40_visual == 1)
  ctrl = pheno %>% filter(RCL40_visual == 0)
  
  carr = case
  all = pheno
  
  cluster_c = list %>% filter(.data[[dawn_c]]  == 1 & PHENOTYPE == 'case') %>% pull(SAMPLE)
  
  
  # carrier setting
  b = carr %>% mutate(carrier = ifelse(WGS_SerialNo. %in% cluster_c, 'RNV carriers', 'Non-carriers')) %>% 
    mutate(is_carrier = ifelse(WGS_SerialNo. %in% cluster_c, 1, 0))
  
  # Verbal memory
  model = aov( Verbal_memory ~ is_carrier + education + age , data=b)
  
  d1 <- summary(model)[[1]]$`Pr(>F)`[1]
  b$carrier %>% table()
  
  # Visual_memory
  model = aov( Visual_memory ~ is_carrier + education + age , data=b)
  d2 <- summary(model)[[1]]$`Pr(>F)`[1]
  b$carrier %>% table()
  
  
  # tmp  = b %>% pivot_longer(cols = c(Verbal_memory, Visual_memory), names_to = 'phenotype_domain', values_to = 'phenotype_score')
  
  library(ggrepel)
  
  # 새로운 데이터 프레임 생성
  new_data <- data.frame(carrier = rep(1, 5), Verbal_memory = rep(-2, 5), label = rep("EIF4A2", 5))
  gene_label = gene_label %>% mutate(carrier = 1)
  # ggplot에 적용
  library(ggrepel)
  
  ggplot(b, aes(x = carrier, y = Verbal_memory)) +   
    geom_violin(aes(fill = carrier)) + 
    geom_boxplot(width = 0.25) + 
    geom_point(data = gene_label, aes(x = 2.3, y = Verbal_memory), shape = 1, size = 3) + 
    geom_signif(tip_length = 0, xmin = 1, xmax = 2, annotations = c(paste0('P = ', signif(d1, 2))), y_position = max(b$Verbal_memory + 0.1, na.rm = TRUE), textsize = 5) +
    geom_text_repel(data = gene_label, aes(label = risk_gene), fill = "white", color = "black", vjust = 1.5, hjust = 0, box.padding = unit(0.35, "lines")) +
    scale_fill_manual(values = c("green", "red")) +  
    theme(axis.title.x = element_blank(), 
          legend.position = 'None') +
    labs(y = 'Verbal memory', x = "") +
    theme_step1()
  
  
  
  ggplot(b, aes(x = carrier, y = Visual_memory)) +   
    geom_violin(aes(fill = carrier)) + 
    geom_boxplot(width = 0.25) + 
    geom_point(data = gene_label, aes(x = 2.3, y = Visual_memory), shape = 1, size = 3) + 
    geom_signif(tip_length = 0, xmin = 1, xmax = 2, annotations = c(paste0('P = ', signif(d1, 2))), y_position = max(b$Visual_memory + 0.1, na.rm = TRUE), textsize = 5) +
    geom_text_repel(data = gene_label, aes(label = risk_gene), fill = "white", color = "black", vjust = 1.5, hjust = 0, box.padding = unit(0.35, "lines")) +
    scale_fill_manual(values = c("green", "red")) +  
    theme(axis.title.x = element_blank(), 
          legend.position = 'None') +
    labs(y = 'Verbal memory', x = "") +
    theme_step1()
  
  
  # 첫 번째 그래프: Verbal memory
  p3_f_1 = ggplot(b, aes(x=carrier, y=Verbal_memory)) +   #theme_step1() + 
    geom_violin(aes(fill=carrier)) + geom_boxplot(width=0.25) + #geom_jitter(alpha=.1) +
    geom_signif(tip_length = 0, xmin=1, xmax=2, annotations= c(paste0('P = ', signif(d1, 2))), y_position = max(b$Verbal_memory +0.1,na.rm = T), textsize = 5) +
    scale_fill_manual(values = c(green, red))  +
    theme(axis.title.x = element_blank(),
          legend.position = 'None') +
    labs(#title = "Cluster 192", 
      y = 'Verbal memory',x = "") +
    theme_step1()
  
  p3_f_1
  
  # 두 번째 그래프: Visual memory
  p3_f_1 <- ggplot(b, aes(x=carrier, y=Visual_memory)) +   #theme_step1() + 
    geom_violin(aes(fill=carrier)) + geom_boxplot(width=0.25) + #geom_jitter(alpha=.1) +
    geom_signif(tip_length = 0, xmin=1, xmax=2, annotations= c(paste0('P = ', signif(d2, 2))), y_position = max(b$Visual_memory +0.2,na.rm = T), textsize = 5) +
    scale_fill_manual(values = c(green, red))  +
    theme(axis.title.x = element_blank(),
          legend.position = 'None') +
    labs(#title = "Cluster 192", 
      y = 'Visual memory',x = "") +
    theme_step1()
  
  # 그래프 배열
  ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(1, 1), align = "h")
  
  
}




## I. pheno check --------------------------------------------------------------
{
  
  
  cluster_target = c('WGS_1420', 'WGS_1564', 
                     'WGS_1063', 
                     'WGS_0721')
  c73 = pheno %>% filter(WGS_SerialNo. %in% cluster_c)

  gene_c = c73 %>% mutate(gene_carrier = ifelse(WGS_SerialNo. %in% cluster_target, 1, 0))
  
  gene_label = gene_c %>% filter(gene_carrier == 1)
  gene_label = gene_label %>% mutate(risk_gene = case_when(
    WGS_SerialNo. == 'WGS_1420' ~ 'VAPB',
    WGS_SerialNo. == 'WGS_1564' ~ 'EIF4A2, TBCCD1, LPP',
    WGS_SerialNo. == 'WGS_1063' ~ 'SAMD5',
    WGS_SerialNo. == 'WGS_0721' ~ 'SLC1A1, PLGRKT',
  ) )


  # Visual_memory
  model = aov( Visual_memory ~ gene_carrier + education + age , data=gene_c)
  
  wilcox.test(Verbal_memory ~ gene_carrier , data=gene_c)
  
  summary(model)[[1]]$`Pr(>F)`[1]
  
  #wilco
  
}

########### STR coverage ########### 
{
  filter = dplyr::filter
  
  red = '#D55E00'
    orange = '#E69F00'
      green = '#009E73'
  
  pca_rm = read_tsv("~/Dropbox/SMC_AD_WGS_paper/Data/STR/eh.v5_coverage_240419_1558_PCA.txt")
  pca_rm_samples = pca_rm %>% filter(pc_remove == 1) %>% pull(sample)
  
  eh.coverage.mean = data.table::fread("~/Dropbox/SMC_AD_WGS_paper/Data/STR/eh.v5.sample_coverage_240419.txt.gz")
  eh.coverage.mean = eh.coverage.mean %>% filter(!sample %in%pca_rm_samples)
  
  eh.coverage.mean %>% dim()
  pheno = read_tsv("~/Dropbox/SMC_AD_WGS_paper/Data/WGS_1824_phenotype_update.240416.txt") %>% dplyr::rename(sample = IID)
  eh.coverage.mean = eh.coverage.mean %>% left_join(pheno)
  
  # shapiro.test(pull(subset(eh.coverage.mean,RCL40_visual == 1), eh.coverage.mean))
  
  eh.coverage.mean$RCL40_visual <- as.character(eh.coverage.mean$RCL40_visual)
  eh.coverage.mean$RCL40_visual <- as.factor(eh.coverage.mean$RCL40_visual)
  
  eh.coverage.mean$eh.coverage.mean = as.numeric(eh.coverage.mean$eh.coverage.mean)
  
  eh.coverage.mean$batch = as.character(eh.coverage.mean$batch)
  eh.coverage.mean$batch = as.factor(eh.coverage.mean$batch)
  
  # glm
  d = glm(eh.coverage.mean ~ RCL40_visual + batch, data=eh.coverage.mean)
  d = summary(d)
  d = d$coefficients[2,4]
  
  p1 = eh.coverage.mean %>% ggplot(aes(RCL40_visual,eh.coverage.mean)) +   
    geom_violin(aes(fill=RCL40_visual),width=0.6,position=position_dodge(width=1),adjust = 0.6) + 
    geom_boxplot(aes(),width=0.1,position=position_dodge(width=1))  +
    geom_signif(tip_length = 0, xmin=1, xmax=2, annotations= c(paste0('P = ', signif(d, 2))), y_position = max(eh.coverage.mean$eh.coverage.mean +1,na.rm = T), textsize = 5)  +
    scale_fill_manual(values = c(green, red))+
    theme_step1() +
    xlab("Amyloid beta") +
    ylab("STR mean coverage by sample")+ 
    scale_x_discrete(labels = c('1' = "Aβ-positive", '0' = "Aβ-negative"))  +
    # ggtitle("All samples (n=1,515)") +
    theme(legend.position = 'none')
  
  
  eh.coverage.mean= eh.coverage.mean %>% filter(DX %in% c('DAT','CU'))
 
  
  d = glm(eh.coverage.mean ~ DX + batch, data=eh.coverage.mean)
  d = summary(d)
  d = d$coefficients[2,4]
  
  p2 =eh.coverage.mean %>% ggplot(aes(DX,eh.coverage.mean)) +   
    geom_violin(aes(fill=DX),width=0.6,position=position_dodge(width=1),adjust = 0.6) + 
    geom_boxplot(aes(),width=0.1,position=position_dodge(width=1))  +
    geom_signif(tip_length = 0, xmin=1, xmax=2, annotations= c(paste0('P = ', signif(d, 2))), y_position = max(eh.coverage.mean$eh.coverage.mean +1,na.rm = T), textsize = 5)  +
    scale_fill_manual(values = c(green, red))+
    theme_step1() +
    xlab("Diagnosis") +
    ylab("STR mean coverage by sample")+ 
    # scale_x_discrete(labels = c('1' = "Aβ-positive", '0' = "Aβ-negative"))  +
    # ggtitle("All samples (n=1,515)") +
    theme(legend.position = 'none')
  p2
  
  p = ggarrange( p1, p2, ncol = 2 ,font.label = list(size = 20), label.y = 1.01, labels = c('a','b'))
  
  ggsave( '~/Dropbox/SMC_AD_WGS_paper/Figures/Extended Data Figure/old/STR_coverage.pdf',p, width = 10, height = 5, device = 'pdf')
  
  
}

########### STR PCA ########### 
{
  pca = read_tsv('STR/eh.v5_coverage_240419_1558_PCA.txt')
  
  pca$batch = as.character(pca$batch)
  
  sd_pc1 <- sd(pca$V1)
  sd_pc2 <- sd(pca$V2)
  
  p = ggplot(pca, aes(V1, V2, color = batch)) +
    geom_point(size=3,alpha = 0.5) +
    xlab(paste0("PC1: ",'64',"% variance")) +
    ylab(paste0("PC2: ",'4',"% variance")) + 
    # coord_fixed()+
    # xlim(c(-2.5,2.5))+
    # ylim(c(-1.8, 1.8))+
    ggtitle('PCA from EH calling coverage') +
    # scale_color_manual(values = c(hc_col,  cd_col, tn_col)) +
    scale_shape_manual(values = c(17))+
    # theme(axis.text.x = element_text(size = 12),
    #       axis.text.y = element_text(size = 12),
    #       axis.title = element_text(size = 14),
    #       legend.title = element_text(size = 14),
    #       legend.text = element_text(size = 14),
    #       #panel.border = element_rect(size = 0.8),
    #       plot.title = element_text(size = 16)) +
    theme_step1() +
    geom_vline(xintercept = 4*sd_pc1, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 4*sd_pc2, color = "red", linetype = "dashed") +
    geom_vline(xintercept = -4*sd_pc1, color = "red", linetype = "dashed") +
    geom_hline(yintercept = -4*sd_pc2, color = "red", linetype = "dashed") +
    geom_hline(yintercept = -220, color = "grey", linetype = "dashed")
  
  ggsave( '~/Dropbox/SMC_AD_WGS_paper/Figures/Extended Data Figure/old/STR_PCA.pdf',p, width = 8, height = 7, device = 'pdf')
  
}

################## original vs inverse transform ################
{
  library(tidyverse)
  setwd("~/Dropbox/SMC_AD_WGS_working/Data/STR/")
  data = read_tsv("eh_DX_single_str_test_inv_noPC_240419_pca_rm.txt")
  
  # 두 변수 간의 Pearson 상관 계수와 p 값 계산
  cor_test <- cor.test(data$allele_count_inv_norm_logistic_p, data$allele_count_logistic_p)
  
  cat("Pearson 상관 계수 =", cor_test$estimate, "\n")
  cat("p 값 =", cor_test$p.value, "\n")
  
  
  
  # 선형 회귀 모델 적합
  linear_model <- lm(allele_count_inv_norm_logistic_p ~ allele_count_logistic_p, data = data)
  
  # 회귀분석 결과 요약
  summary(linear_model)
  
  # R^2 값
  r_squared <- summary(linear_model)$r.squared
  cat("R^2 =", r_squared, "\n")
  
  # p-value
  p_value <- summary(linear_model)$coefficients[2, 4]
  cat("p-value =", p_value, "\n")
  if (p_value == 0) {
    p_value <- "< 2.2e-16"
  }
  
  data <- data %>% mutate(
    `-log10 (untransformed p-value)` = -log10(allele_count_logistic_p),
    `-log10 (rank-based inverse normal transformed p-value)` = -log10(allele_count_logistic_p)
  )
  
  # 산점도와 회귀선 그리기
  p <- ggplot(data, aes(x =  `-log10 (untransformed p-value)`, y = `-log10 (rank-based inverse normal transformed p-value)`)) +
    geom_point() +  # 산점도 추가
    geom_smooth(method = "lm", se = FALSE, color = "red") +  
    labs(title = "Scatterplot with Linear Regression Line",
         x = "allele_count_logistic_p", y = "allele_count_inv_norm_logistic_p") + 
    geom_text(x = Inf, y = -Inf, label = paste("R^2 =", round(r_squared, 3), "\np-value ", p_value), 
                hjust = 1, vjust = 0, size = 4, parse = TRUE) +
    theme_step1()
  
  # ggsave( '~/Dropbox/SMC_AD_WGS_paper/Figures/Extended Data Figure/old/STR_inverse_liner.pdf',p, width = 8, height = 7, device = 'pdf')
  
}


################## STR DX model 1 ################
{
  dat.test = read_tsv("~/Dropbox/SMC_AD_WGS_working/Data/STR/eh_DX_single_str_test_inv_noPC_240419_pca_rm.txt")
  
  
  { 
    p = 'allele_count_inv_norm'
    # allele_count
    input = 'allele_count'
    
    ## allele_count
    df_cond <- dat.test %>% 
      mutate(standard_mean_difference = ((.data[[paste0(input, "_mean_ad")]] - .data[[paste0(input, "_mean_control")]]) / .data[[paste0(input, "_sd")]])) %>%
      mutate(p_val_adj = .data[[paste0(p, "_logistic_p")]]) %>%
      mutate(p_val_adj = ifelse(p_val_adj == 0, 1.0e-300, p_val_adj)) %>%
      mutate(enrichment = ifelse(standard_mean_difference > 0 , 'case', ifelse(standard_mean_difference < 0 , 'control', 'nan')))
    
    x_range = max(abs(df_cond$standard_mean_difference), na.rm = T) * 1.05
    y_range = max(-log10(df_cond$p_val_adj),na.rm = T) * 1.1
    # gene_up = df_cond %>% filter(standard_mean_difference > 0) %>% slice_min(., p_val_adj,n = 5) %>% pull(gene)
    # gene_down = df_cond %>% filter(standard_mean_difference < 0) %>% slice_min(., p_val_adj, n = 5) %>% pull(gene)
    # res.gene = c(gene_up, gene_down)
    
    library(ggplot2)
    
    # 특정 index 선택
    specific_indices <- c(293559, 293569)
    
    df_cond <- df_cond %>%
      mutate(color_group = ifelse(index %in% specific_indices, "STRs identified in \namyloid beta positivity", as.character(cut(p_val_adj, breaks = c(-Inf, 0.001, 0.005, 0.01,0.05, Inf), 
                                                                                                    labels = c("p < 0.001", "0.001 ≤ p < 0.005","0.005 ≤ p < 0.01","0.01 ≤ p < 0.05", "p ≥ 0.05")))))
    
    # 그래프 생성
  p = ggplot(df_cond, aes(standard_mean_difference, -log10(p_val_adj), color = color_group)) + 
      coord_cartesian(ylim = c(0, 15)) +
      geom_point(size = 0.8) + 
      scale_color_manual(values = c(
        "p < 0.001" = "#9D0B0B", 
        "0.001 ≤ p < 0.005" = "#DA2D2D", 
        "0.005 ≤ p < 0.01" = "#EB8242", 
        "0.01 ≤ p < 0.05" = "#F6DA63",
        "p ≥ 0.05" = "grey50",
        "Specific Indices" = "#640D6B"
      ), guide = FALSE) +
      theme_minimal(base_size = 10) +
      geom_hline(yintercept = -log10(0.05/293751), linetype="dashed", color = 'red') +
      geom_vline(xintercept = 0, linetype="dashed",color = 'grey') +
      guides(color = guide_legend(title = "p-value"),
             y = guide_axis(title = "-log10(p_val)"))  +
      theme(
        # legend.position = 'none'
      ) +
      geom_point(data = df_cond %>% filter(index %in% specific_indices), aes(x = standard_mean_difference, y = -log10(p_val_adj)), color = "#640D6B")+
    xlim(c(-max(df_cond$standard_mean_difference), max(df_cond$standard_mean_difference))) +
    theme_step1()
    
    
    ggsave(
      filename = paste('~/Dropbox/SMC_AD_WGS_paper/Figures/Extended Data Figure/old/eh.DX_model1.png'),
      plot = p,
      width = 8,
      height = 5
    )
    
    
    
  }
  
}
