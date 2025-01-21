{ 
  # Last updated: 
# project : Figure 3. 


# Setting environment
rm(list=ls())
sessionInfo()
setwd('~/Dropbox/SMC_AD_WGS_paper/Data/SMC_cwas_results_20240416/')

options(stringsAsFactors = F)

library(ggh4x)
library(ggsignif)
library(grid)
library(tidyverse)
library(readxl)
library(gtable)
library(cowplot)
library(ggpubr)
library(extrafont)
library(ggcorrplot)
library(corrplot) 
library(ggplotify)
library(patchwork)
library(gridExtra)
library(igraph)

rename = dplyr::rename
select = dplyr::select

red = '#D55E00'
  orange = '#E69F00'
    green = '#009E73'
      
    case = 925
    ctrl = 634
      
}
# font_import()
# y
# loadfonts()
if(T){
  extrafont::font_import(pattern = "Arial", prompt = F)
  extrafont::loadfonts()
}


`%notin%` <- Negate(`%in%`)

version = 'v9.5'


# Data prepare
{
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
        axis.title.y = element_text(margin = margin(r = 10, unit = "pt")))
}

resize_heights <- function(g, heights = rep(1, length(idpanels))){
  idpanels <- unique(g$layout[grepl("panel",g$layout$name), "t"])
  g$heights <- unit.c(g$heights)
  g$heights[idpanels] <- unit.c(do.call(unit, list(heights, 'null')))
  g
}



# Load data
## Phenotype data
### korean AD dataset



## AD CWAS
### Risk score
#### noncoding
res_cs = read_tsv("risk_score_after_fs/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.lasso_results_tf0.8_thres_3.txt") %>% filter(Domain != "all")


### Burden shift

# burdenShift = data.table::fread("~/Dropbox/ADWGS/CWAS/outputs/CWAS_SMC_batch1-7_RCL40_visual_AF_DB_MAF0.1_annot5.1_1559/CS_noncoding_tf0.7_thr4_S107/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_noncoding_0.7_4_S107_CS.binom_pvals.txt.gz")
# burdenShift = data.table::fread("burden_shift/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.binom_pvals.txt.gz")
bs_res = data.table::fread("burden_shift/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.burdenshift_p0.05_cutoff10.txt")

### DAWN
# dawn_annot = read_tsv("SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220.DAWN_cluster_filtered_intergenic_k253_l5.25_corr0.22_size2.cluster_list.tsv") 
# result = readxl::read_xlsx('../Tables/Supplementary Table 3.xlsx',sheet = 8)

# risk_cate = read_tsv('SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_noncoding_0.7_4_S107_CS.intergenic_k253_l5.25_corr0.22_size2.tsv')

# dawn_corr = read_csv("SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_noncoding_0.7_4_S107_CS.intergenic_k253_l5.25_corr0.22_size2.category_correlation.csv")
# 
# dawn_corr_matrix = read.delim('~/Dropbox/SMC_AD_WGS_paper/Data/SMC_cwas_results_20240416/')
dawn_corr = read_rds("DAWN_L2.correlation_matrix.rds") #%>% select(-1)
rownames(dawn_corr) <- colnames(dawn_corr)
dawn_corr_matrix <- as.matrix(dawn_corr)

}

# Figure 3
## B. Risk score with feature selection --------------------------------------------------------------------------------
{

conversion_rules <- c(
  "coding" = "Coding",
  "PTV" = "PTV",
  "missense" = "Missense",
  "noncoding" = "Noncoding",
  "promoter" = "Promoter",
  "intron" = "Intron",
  "5primeUTR" = "5`UTR",
  "3primeUTR" = "3`UTR",
  "splice" = "Splice site",
  "intergenic" = "Intergenic"
)


data = res_cs %>% filter(Seed == "average") %>%
  mutate(value = R2*100) %>%
  mutate(Formatted_Perm_P = sprintf("%.2e", Perm_P))

data$Domain <- conversion_rules[data$Domain]
data$Domain <- factor(data$Domain, levels = c( "Coding", "PTV", "Missense", "Noncoding", "Promoter", "Intron", "5`UTR", "3`UTR","Splice Region", "Intergenic"))
data <- data %>%  rename( category = Domain,p_value = Perm_P)
# 데이터프레임 내 category 변수를 factor로 변환하고 역순으로 순서를 지정
data$category <- factor(data$category, levels = rev(levels(data$category)))

# ggplot 그래프 생성 및 y 축 순서를 역순으로 설정
p3_b = ggplot(data, aes(x = value, y = category)) +
  geom_bar(stat = "identity", aes(fill = ifelse(p_value >= 0.05, "Gray", "Blue"))) +
  labs(x = "R2 (%)", y = "Category") +
  coord_cartesian(xlim = c(ifelse(min(data$value) > 0, 0, min(data$value) - 0), max(data$value) + 1)) +  # x 축 범위
  geom_text(aes(label = sprintf("P = %.2e", p_value),
                fontface = ifelse(p_value < 0.05, "bold", "plain"),
                color = text),
            hjust = -0.05,
            color = "black",
            size = 5,
            data = subset(data, value > 0 & p_value < 0.05)) +  
  geom_text(aes(label = sprintf("P = %.2e", p_value),
                fontface = ifelse(p_value < 0.05, "bold", "plain"),
                color = text),
            hjust = -0.05,
            color = "grey",
            size = 5,
            data = subset(data, value > 0 & p_value >= 0.05)) +
  geom_text(aes(label = sprintf("P = %.2e", p_value),
                fontface = ifelse(p_value < 0.05, "bold", "plain"),
                color = text),
            hjust = 0,
            x = 0.2,
            color = "black", 
            size = 5,
            data = subset(data, value < 0 & p_value < 0.05)) +
  geom_text(aes(label = sprintf("P = %.2e", p_value),
                fontface = ifelse(p_value < 0.05, "bold", "plain"),
                color = text),
            hjust = 0,
            x = 0.2,
            color = "grey", 
            size = 5,
            data = subset(data, value < 0 & p_value >= 0.05)) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("Gray" = "gray", "Blue" = red)) +
  theme_step1()
p3_b

# ggsave("~/Dropbox/Paper_Korean_AD_WGS/Figures/Figure3_1.pdf",p1,width = 10 ,height = 9)
# ggsave("~/Dropbox/Paper_Korean_AD_WGS/Figures/Figure3_1.png",p1,width = 10 ,height = 9)
}

## C. Burden test ------------------------------------------------------------------------------------------
{
  options(stringsAsFactors = F)
  library(tidyverse)
  library(readr)
  library(cowplot)
  library(ggsignif)
  library(ggtext)
  library(ggplot2)
  
  if(F){
    extrafont::font_import(pattern = "Arial", prompt = FALSE)
    extrafont::loadfonts()
  }
  
  
  source('~/Dropbox/ADWGS/CWAS/Scripts/func_pval_to_power10.R')
  source('~/Dropbox/ADWGS/CWAS/Scripts/func_allocate_var.R')
  
  date = '20240221'
  file_path = ''
  ## B. Volcano plot
  if(T){
    dt1 = data.table::fread(paste('burden_test/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.burden_test.txt', sep = '/'), header = T, data.table = F)
    
    # dt1 %>% filter(P < 0.05 & Relative_Risk > 1) %>% pull(gencode) %>% table() %>% view()
    
    
    category_count = data.table::fread(paste0(file_path, 'burden_test/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.category_counts.txt'), data.table = F)
    categories = read_tsv('burden_test/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.category_info.txt')
    
    categories = merge(categories %>% dplyr::select(-c(variant_type, gene_set, functional_score, gencode, functional_annotation)), category_count, by = 'Category')
    
    dt1 = merge(dt1, categories, by = 'Category')
    
    # dt1 = dt1 %>% filter(is_noncoding==1)
    # 
    # dt1 = dt1 %>% filter(Raw_counts>=7) # eight..?
    
    ## Relative risk visualization using binom results
    head(dt1) # gencode = genomic functional_annotation / functional_annotation = regulatory track
    dt1$Direction = ifelse(dt1$Relative_Risk>1, 'case', 'control') # No. 1.
    
    dt1 = dt1 %>% mutate(log2_RR = log2(Relative_Risk))
    table(dt1$Direction)
    # x axis intervals & modify Inf to numerical numbers
    minrr = dt1 %>% filter(is.finite(log2_RR)) %>% pull(log2_RR) %>% min()
    maxrr = dt1 %>% filter(is.finite(log2_RR)) %>% pull(log2_RR) %>% max()
    real_max = max(abs(minrr), abs(maxrr))
    max_lim = real_max + 0.5
    dt1$log2_RR = ifelse(dt1$log2_RR == Inf, max_lim,
                         ifelse(dt1$log2_RR == -Inf, -max_lim,
                                dt1$log2_RR))
    
    eff_num = 1105
    pcutoff = 0.05/eff_num
    eff_num_chr = as.character(eff_num)
    eff_num_chr = strsplit(eff_num_chr, split = '', fixed = T)[[1]]
    eff_num_chr = paste(eff_num_chr[1:(length(eff_num_chr)-3)],
                        paste(eff_num_chr[(length(eff_num_chr)-2):length(eff_num_chr)], collapse = ''),
                        sep = ',')
    
    my_exp = pval_to_power10(pcutoff)
    
    # dt1 <- dt1 %>%
    #   mutate(Risk_cluster = case_when(
    #     Category %in% risk_cate$Category ~ T,
    #     TRUE ~ F
    #   ))
    dt1 = dt1 %>% mutate(is_intergenic = ifelse(gencode == 'IntergenicRegion', 1, 0))
    
    
    
    p3_c = dt1 %>%
      mutate(is_inter = case_when(
        is_intergenic == 1 & P < 0.05 & Relative_Risk > 1 ~ T,
        is_intergenic == 1 & P > 0.05 & Relative_Risk > 1 ~ F,
        is_intergenic == 0 & P < 0.05 & Relative_Risk > 1 ~ F,
        is_intergenic == 0 & P > 0.05 & Relative_Risk > 1 ~ F,
        T ~ F
      )) %>%
      #mutate(Risk_cluster = ifelse(Category %in% risk_cate, 1, 0)) %>%
      ggplot(aes(x = log2_RR, -log10(P))) +
      geom_vline(xintercept = log2(1), color = 'grey80') +
      geom_point(aes(
        # fill = Risk_cluster 
        fill = is_inter 
                     ),pch=21, size=3, show.legend = T, color = 'grey50', stroke = 0) +
      guides(colour = "none") + # Turning off color legend
      theme_bw() +
      ylim(c(0, -log10(pcutoff)+0.25)) +
      xlim(c(-4,4)) +
      theme(#text=element_text(family="Arial"),
        plot.title = element_blank(),
        # legend.text = element_text(size = 8),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.margin = ggplot2::margin(0.2, 0.2, 0.2, 0.2, "cm"),
        legend.position = "none",
        # legend.position = c(0.3, 1),
        legend.justification = c(1, 2.2),
        legend.background = element_rect(color = 'black'),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.1, "cm"),
        legend.title = element_blank()
        ) +
      labs(x = expression('Relative Risk' ~ (log[2])),
           y = expression('P' ~ (-log[10])),
           fill = 'Risk category') +
      scale_x_continuous(breaks = c(-max_lim, -2, 0, 2, max_lim),
                         labels = c('-Inf', '-2', '0', '2', 'Inf')) +
      geom_hline(yintercept = -log10(0.05/eff_num), color = red, linetype = "dashed") +
      geom_richtext(data = data.frame(x = -Inf,
                                      y = -log10(0.05/eff_num)+0.25,
                                      label = paste(' P=', my_exp, sep = '')),
                    aes(x = x, y = y, label = label),
                    label.color = NA, fill = NA, color = red,
                    size = 5,
                    hjust = 'left',
                    family="Arial",
                    inherit.aes = FALSE) +
      geom_text(data = data.frame(x = -Inf,
                                  y=-log10(0.05/eff_num)+0.1,
                                  label = paste('', eff_num_chr, 'effective tests', sep = ' ')),
                aes(x = x, y = y, label = label),
                color = red,
                size = 5,
                family="Arial",
                hjust = 'left') +
      geom_hline(yintercept = -log10(0.05), color = 'black', linetype = "dashed") + # P=0.05
      geom_richtext(data = data.frame(x = -Inf,
                                      y = (-log10(0.05))+0.12,
                                      label = paste(' P=', pval_to_power10(0.05), sep = '')),
                    aes(x = x, y = y, label = label),
                    label.color = NA,
                    fill = NA,
                    color = 'black',
                    family="Arial",
                    size = 5,
                    hjust = 'left',
                    inherit.aes = FALSE)+
      # geom_point(data = subset(dt1,Risk_cluster == TRUE), aes(x = log2_RR, -log10(P)), color = '#D83F31')+
      scale_fill_manual(values = c("TRUE" = "#D83F31", "FALSE" = "lightgrey")) + # Setting colors for "T" and "F"
      theme_step1() 
    p3_c
    
    # ggsave('~/Dropbox/SMC_AD_WGS_paper/Figures/olds/volcano_plot.pdf',width = 6,height = 6)
    
  }
}



## D. Burden shfit --------------------------------------------------------------

if(T){
  maxP = 0.05
  count_cutoff = 10
  
  # Function for counting categories
  countCats <- function(pvals, pvalThresh) {
    nCase <- length(pvals[pvals>0 & pvals<=pvalThresh])
    nControl <- length(pvals[pvals<0 & abs(pvals)<=pvalThresh])
    return(c("case"=nCase,"control"=nControl))
  }

  # burden test
  dt1 = read.delim('burden_test/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.burden_test.txt')
  # category_counts.txt
  dt2 = read.delim('burden_test/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.category_counts.txt')

  dt3 = merge(dt1, dt2, by = 'Category')

  # gencode 리스트
  # gencode_list <- c(  "3primeUTR" , "5primeUTR" , "coding"   ,  "intergenic" ,"intron","missense" , "noncoding" , "promoter" ,  "PTV"   ,     "splice" )

  gencode_list <- c("3PrimeUTRsRegion"   ,   "5PrimeUTRsRegion"  ,      "IntergenicRegion"     ,  "IntronRegion"    , "PromoterRegion"    ,     "SpliceSiteRegion"    )
  
  # # 결과를 저장할 빈 데이터 프레임 생성
  # combined_ggpermCounts <- data.frame()
  # combined_obs <- data.frame()
  # 
  # # 각 gencode에 대해 행위 수행
  # for (gen in gencode_list) {
  #   # 특정 gencode에 대한 데이터 필터링
  #   burdenResTrim_gencode <- dt3 %>%
  #     filter(Raw_counts >= count_cutoff) %>%
  #     filter(gencode == gen)
  # 
  #   # nObsCase 및 nObsControl 계산
  #   nObsCase <- nrow(burdenResTrim_gencode[burdenResTrim_gencode$P <= maxP & burdenResTrim_gencode$Relative_Risk > 1, ])
  #   nObsControl <- nrow(burdenResTrim_gencode[burdenResTrim_gencode$P <= maxP & burdenResTrim_gencode$Relative_Risk < 1, ])
  # 
  #   # burdenShift에서 필요한 열 선택
  #   c_idx <- colnames(burdenShift) %in% burdenResTrim_gencode$Category
  #   burdenShiftTrim <- burdenShift[, ..c_idx]
  # 
  #   # permCounts 계산
  #   permCounts <- data.frame(t(apply(burdenShiftTrim, 1, countCats, pvalThresh = maxP)))
  # 
  #   # nPermCase 및 nPermControl 계산
  #   nPermCase <- nrow(permCounts[permCounts$case >= nObsCase, ])
  #   pCase <- nPermCase / nrow(permCounts)
  #   nPermControl <- nrow(permCounts[permCounts$control >= nObsControl, ])
  #   pControl <- nPermControl / nrow(permCounts)
  # 
  #   # ggpermCounts 데이터 프레임 생성
  #   ggpermCounts <- data.frame("N_signif_tests" = c(permCounts$case, permCounts$control))
  #   myfontsize <- 8
  #   larger <- max(c(nObsCase, nObsControl))
  #   ggpermCounts$scaled_cnts <- scale(ggpermCounts$N_signif_tests)
  # 
  #   # obs 데이터 프레임 생성
  #   mean_val <- mean(ggpermCounts$N_signif_tests, na.rm = TRUE)
  #   std_dev <- sd(ggpermCounts$N_signif_tests, na.rm = TRUE)
  #   obs <- data.frame(
  #     nObsCase = nObsCase,
  #     nObsControl= nObsControl,
  #     norm_nObsCase = (nObsCase - mean_val) / std_dev,
  #                     pCase = pCase,
  #                     norm_nObsControl = (nObsControl - mean_val) / std_dev,
  #                     pControl = pControl)
  # 
  #   # gencode 정보 추가
  #   ggpermCounts$gen <- gen
  #   obs$gen <- gen
  # 
  #   # 결과를 합침
  #   combined_ggpermCounts <- bind_rows(combined_ggpermCounts, ggpermCounts)
  #   combined_obs <- bind_rows(combined_obs, obs)
  # 
  # }
  # 
  # ggpermCounts = combined_ggpermCounts
  # obs = combined_obs
  # obs$y_position = 1.5
  # 
  # ggpermCounts$gen <- factor(ggpermCounts$gen, levels =  gencode_list <- c("3PrimeUTRsRegion"   ,   "5PrimeUTRsRegion"  ,      "IntergenicRegion"     ,  "IntronRegion"    , "PromoterRegion"    ,     "SpliceSiteRegion"    ))
  # obs$gen <- factor(obs$gen, levels =  gencode_list <- c("3PrimeUTRsRegion"   ,   "5PrimeUTRsRegion"  ,      "IntergenicRegion"     ,  "IntronRegion"    , "PromoterRegion"    ,     "SpliceSiteRegion"    ))
  # 
  # write_rds(ggpermCounts,"240416_Burden_shift_results.RDS")
  # write_rds(obs,"240416_Burden_shift_p_val.RDS")
  larger = 10
  ggpermCounts =  readRDS("240416_Burden_shift_results.RDS")
  obs =  readRDS("240416_Burden_shift_p_val.RDS")
  pCase =0
  pControl = 0
  
  # p3_d = ggplot(ggpermCounts, aes(x = scaled_cnts, y = gen)) +
  #   geom_violin(width = 0.6, adjust = 2, fill = 'grey') + 
  #   geom_boxplot(width = 0.2, outlier.shape = NA) +
  #   facet_wrap(~ gen, 
  #              # scales = "free_y", 
  #              scales = "free_y", 
  #              ncol = 1, strip.position = "left") +
  #   geom_vline(data = obs, aes(xintercept = nObsCase), size = 1, colour = red, linetype = 1) +
  #   geom_vline(data = obs, aes(xintercept = nObsControl), size = 1, colour = green, linetype = 1) +
  #   geom_text(data = obs, aes(x = nObsCase+ 0.2,y = y_position +0.1, label = paste("p = ", round(pCase, 4))),
  #             vjust = 1.5, hjust = 0, color = red, size = 6) +
  #   theme(strip.background = element_blank(),
  #         strip.text = element_blank(), # strip 제목 제거
  #         axis.text.x = element_text(hjust = 0, angle = 0),
  #         axis.title.y = element_blank()) +
  #   theme_step1() 
  # p3_d
  
  # 'PromoterRegion', 'UTRsRegion', 'IntronRegion', 'IntergenicRegion'
  
  bs_PromoterRegion = ggpermCounts %>%
    filter(gen == 'PromoterRegion') %>% 
    ggplot(aes(x = N_signif_tests)) +
    geom_density(size = 0.25, fill = "gray92", bw = 0.7) +
    scale_x_continuous(breaks = seq(round(max(ggpermCounts$N_signif_tests) / 4), round(max(ggpermCounts$N_signif_tests) / 4) * 4, round(max(ggpermCounts$N_signif_tests) / 4))) +
    theme_step1() +
    geom_vline(mapping = aes(xintercept = pull(subset(obs, gen == 'PromoterRegion'), nObsCase)), size = 0.7, colour = red, linetype = 2) +
    geom_vline(mapping = aes(xintercept = pull(subset(obs, gen == 'PromoterRegion'), nObsControl)), size = 0.7, colour = green, linetype = 2) +
    scale_y_continuous(n.breaks = 2.5) +
    labs(y = 'Density', x = 'Number of significant tests',
         title = "Promoter") +
    geom_richtext(data = data.frame(x = pull(subset(obs, gen == 'PromoterRegion'), nObsCase),
                                    y = Inf,
                                    label = paste('p=', sprintf("%.2g", pull(subset(obs, gen == 'PromoterRegion'), pCase)), sep = '')),
                  aes(x = x, y = Inf, label = label),
                  vjust = 1, hjust = 0,
                  label.color = NA, fill = NA, color = red,
                  size = 5 ,
                  family = "Arial",
                  inherit.aes = FALSE) +
    geom_richtext(data = data.frame(x = pull(subset(obs, gen == 'PromoterRegion'), nObsControl) + (larger * 0.02),
                                    y = Inf,
                                    label = paste('p=', ifelse(pControl == 1, 1, sprintf("%.2g", pull(subset(obs, gen == 'PromoterRegion'), pControl))), sep = '')),
                  aes(x = x, y = Inf, label = label),
                  vjust = 0, hjust = 0.5,
                  family = "Arial",
                  label.color = NA, fill = NA, color = green,
                  size = 5 ,
                  inherit.aes = FALSE) +
    theme(plot.title = element_text(hjust = 0, margin = margin(b = 20))) +
    coord_cartesian(clip = 'off')
  
  bs_5UTRsRegion = ggpermCounts %>%
    filter(gen == '5PrimeUTRsRegion') %>% 
    ggplot(aes(x = N_signif_tests)) +
    geom_density(size = 0.25, fill = "gray92", bw = 0.7) +
    scale_x_continuous(breaks = seq(round(max(ggpermCounts$N_signif_tests) / 4), round(max(ggpermCounts$N_signif_tests) / 4) * 4, round(max(ggpermCounts$N_signif_tests) / 4))) +
    theme_step1() +
    geom_vline(mapping = aes(xintercept = pull(subset(obs, gen == '5PrimeUTRsRegion'), nObsCase)), size = 0.7, colour = red, linetype = 2) +
    geom_vline(mapping = aes(xintercept = pull(subset(obs, gen == '5PrimeUTRsRegion'), nObsControl)), size = 0.7, colour = green, linetype = 2) +
    scale_y_continuous(n.breaks = 2.5) +
    labs(y = 'Density', x = 'Number of significant tests',
         title = "5` UTR") +
    geom_richtext(data = data.frame(x = pull(subset(obs, gen == '5PrimeUTRsRegion'), nObsCase),
                                    y = Inf,
                                    label = paste('p=', sprintf("%.2g", pull(subset(obs, gen == '5PrimeUTRsRegion'), pCase)), sep = '')),
                  aes(x = x + 5, y = Inf, label = label),
                  vjust = 1, hjust = 0.3,
                  label.color = NA, fill = NA, color = red,
                  size = 5 ,
                  family = "Arial",
                  inherit.aes = FALSE) +
    geom_richtext(data = data.frame(x = pull(subset(obs, gen == '5PrimeUTRsRegion'), nObsControl) + (larger * 0.02),
                                    y = Inf,
                                    label = paste('p=', ifelse(pControl == 1, 1, sprintf("%.2g", pull(subset(obs, gen == '5PrimeUTRsRegion'), pControl))), sep = '')),
                  aes(x = x, y = Inf, label = label),
                  vjust = 0, hjust = 0.5,
                  family = "Arial",
                  label.color = NA, fill = NA, color = green,
                  size = 5 ,
                  inherit.aes = FALSE) +
    theme(plot.title = element_text(hjust = 0, margin = margin(b = 20))) +
    coord_cartesian(clip = 'off')
  
  bs_3UTRsRegion = ggpermCounts %>%
    filter(gen == '3PrimeUTRsRegion') %>% 
    ggplot(aes(x = N_signif_tests)) +
    geom_density(size = 0.25, fill = "gray92", bw = 0.7) +
    scale_x_continuous(breaks = seq(round(max(ggpermCounts$N_signif_tests) / 4), round(max(ggpermCounts$N_signif_tests) / 4) * 4, round(max(ggpermCounts$N_signif_tests) / 4))) +
    theme_step1() +
    geom_vline(mapping = aes(xintercept = pull(subset(obs, gen == '3PrimeUTRsRegion'), nObsCase)), size = 0.7, colour = red, linetype = 2) +
    geom_vline(mapping = aes(xintercept = pull(subset(obs, gen == '3PrimeUTRsRegion'), nObsControl)), size = 0.7, colour = green, linetype = 2) +
    scale_y_continuous(n.breaks = 2.5) +
    labs(y = 'Density', x = 'Number of significant tests',
         title = "3` UTR") +
    geom_richtext(data = data.frame(x = pull(subset(obs, gen == '3PrimeUTRsRegion'), nObsCase),
                                    y = Inf,
                                    label = paste('p=', sprintf("%.2g", pull(subset(obs, gen == '3PrimeUTRsRegion'), pCase)), sep = '')),
                  aes(x = x + 5, y = Inf, label = label),
                  vjust = 1, hjust = 0.5,
                  label.color = NA, fill = NA, color = red,
                  size = 5 ,
                  family = "Arial",
                  inherit.aes = FALSE) +
    geom_richtext(data = data.frame(x = pull(subset(obs, gen == '3PrimeUTRsRegion'), nObsControl) + (larger * 0.02),
                                    y = Inf,
                                    label = paste('p=', ifelse(pControl == 1, 1, sprintf("%.2g", pull(subset(obs, gen == '3PrimeUTRsRegion'), pControl))), sep = '')),
                  aes(x = x, y = Inf, label = label),
                  vjust = 0, hjust = 0.5,
                  family = "Arial",
                  label.color = NA, fill = NA, color = green,
                  size = 5 ,
                  inherit.aes = FALSE) +
    theme(plot.title = element_text(hjust = 0, margin = margin(b = 20))) +
    coord_cartesian(clip = 'off')
  bs_3UTRsRegion
  
  
  bs_IntronRegion = ggpermCounts %>%
    filter(gen == 'IntronRegion') %>% 
    ggplot(aes(x = N_signif_tests)) +
    geom_density(size = 0.25, fill = "gray92", bw = 0.7) +
    scale_x_continuous(breaks = seq(round(max(ggpermCounts$N_signif_tests) / 4), round(max(ggpermCounts$N_signif_tests) / 4) * 4, round(max(ggpermCounts$N_signif_tests) / 4))) +
    theme_step1() +
    geom_vline(mapping = aes(xintercept = pull(subset(obs, gen == 'IntronRegion'), nObsCase)), size = 0.7, colour = red, linetype = 2) +
    geom_vline(mapping = aes(xintercept = pull(subset(obs, gen == 'IntronRegion'), nObsControl)), size = 0.7, colour = green, linetype = 2) +
    scale_y_continuous(n.breaks = 2.5) +
    labs(y = 'Density', x = 'Number of significant tests',
         title = "Intron") +
    geom_richtext(data = data.frame(x = pull(subset(obs, gen == 'IntronRegion'), nObsCase),
                                    y = Inf,
                                    label = paste('p=', sprintf("%.2g", pull(subset(obs, gen == 'IntronRegion'), pCase)), sep = '')),
                  aes(x = x + 3.5, y = Inf, label = label),
                  vjust = 0, hjust = 0.5, 
                  label.color = NA, fill = NA, color = red,
                  size = 5 ,
                  family = "Arial",
                  inherit.aes = FALSE) +
    geom_richtext(data = data.frame(x = pull(subset(obs, gen == 'IntronRegion'), nObsControl) + (larger * 0.02),
                                    y = Inf,
                                    label = paste('p=', ifelse(pControl == 1, 1, sprintf("%.2g", pull(subset(obs, gen == 'IntronRegion'), pControl))), sep = '')),
                  aes(x = x, y = Inf, label = label),
                  vjust = 1, hjust = 0,
                  family = "Arial",
                  label.color = NA, fill = NA, color = green,
                  size = 5 ,
                  inherit.aes = FALSE) +
    theme(plot.title = element_text(hjust = 0, margin = margin(b = 20))) +
    coord_cartesian(clip = 'off')
  
  
  bs_SpliceSiteRegion = ggpermCounts %>%
    filter(gen == 'SpliceSiteRegion') %>% 
    ggplot(aes(x = N_signif_tests)) +
    geom_density(size = 0.25, fill = "gray92", bw = 0.7) +
    scale_x_continuous(breaks = seq(round(max(ggpermCounts$N_signif_tests) / 4), round(max(ggpermCounts$N_signif_tests) / 4) * 4, round(max(ggpermCounts$N_signif_tests) / 4))) +
    theme_step1() +
    geom_vline(mapping = aes(xintercept = pull(subset(obs, gen == 'SpliceSiteRegion'), nObsCase)), size = 0.7, colour = red, linetype = 2) +
    geom_vline(mapping = aes(xintercept = pull(subset(obs, gen == 'SpliceSiteRegion'), nObsControl)), size = 0.7, colour = green, linetype = 2) +
    scale_y_continuous(n.breaks = 2.5) +
    labs(y = 'Density', x = 'Number of significant tests',
         title = "Splice Site") +
    geom_richtext(data = data.frame(x = pull(subset(obs, gen == 'SpliceSiteRegion'), nObsCase),
                                    y = Inf,
                                    label = paste('p=', sprintf("%.2g", pull(subset(obs, gen == 'SpliceSiteRegion'), pCase)), sep = '')),
                  aes(x = x, y = Inf, label = label),
                  vjust = 0, hjust = 0.5,
                  label.color = NA, fill = NA, color = red,
                  size = 5 ,
                  family = "Arial",
                  inherit.aes = FALSE) +
    geom_richtext(data = data.frame(x = pull(subset(obs, gen == 'SpliceSiteRegion'), nObsControl) + (larger * 0.02),
                                    y = Inf,
                                    label = paste('p=', ifelse(pControl == 1, 1, sprintf("%.2g", pull(subset(obs, gen == 'SpliceSiteRegion'), pControl))), sep = '')),
                  aes(x = x, y = Inf, label = label),
                  vjust = 0, hjust = 0.5,
                  family = "Arial",
                  label.color = NA, fill = NA, color = green,
                  size = 5 ,
                  inherit.aes = FALSE) +
    theme(plot.title = element_text(hjust = 0, margin = margin(b = 20))) +
    coord_cartesian(clip = 'off')
  
  
  bs_IntergenicRegion = ggpermCounts %>%
    filter(gen == 'IntergenicRegion') %>% 
    ggplot(aes(x = N_signif_tests)) +
    geom_density(size = 0.25, fill = "gray92", bw = 0.7) +
    scale_x_continuous(breaks = seq(round(max(ggpermCounts$N_signif_tests) / 4), round(max(ggpermCounts$N_signif_tests) / 4) * 4, round(max(ggpermCounts$N_signif_tests) / 4))) +
    theme_step1() +
    geom_vline(mapping = aes(xintercept = pull(subset(obs, gen == 'IntergenicRegion'), nObsCase)), size = 0.7, colour = red, linetype = 2) +
    geom_vline(mapping = aes(xintercept = pull(subset(obs, gen == 'IntergenicRegion'), nObsControl)), size = 0.7, colour = green, linetype = 2) +
    scale_y_continuous(n.breaks = 2.5) +
    labs(y = 'Density', x = 'Number of significant tests',
         title = "Intergenic") +
    geom_richtext(data = data.frame(x = pull(subset(obs, gen == 'IntergenicRegion'), nObsCase),
                                    y = Inf,
                                    label = paste('p=', sprintf("%.2g", pull(subset(obs, gen == 'IntergenicRegion'), pCase)), sep = '')),
                  aes(x = x, y = Inf, label = label),
                  vjust = 0, hjust = 0.5,
                  label.color = NA, fill = NA, color = red,
                  size = 5 ,
                  family = "Arial",
                  inherit.aes = FALSE) +
    geom_richtext(data = data.frame(x = pull(subset(obs, gen == 'IntergenicRegion'), nObsControl) + (larger * 0.02),
                                    y = Inf,
                                    label = paste('p=', ifelse(pControl == 1, 1, sprintf("%.2g", pull(subset(obs, gen == 'IntergenicRegion'), pControl))), sep = '')),
                  aes(x = x, y = Inf, label = label),
                  vjust = 0, hjust = 0.5,
                  family = "Arial",
                  label.color = NA, fill = NA, color = green,
                  size = 5 ,
                  inherit.aes = FALSE) +
    theme(plot.title = element_text(hjust = 0, margin = margin(b = 20))) +
    coord_cartesian(clip = 'off')
  p3_d_1 = ggarrange(bs_IntergenicRegion,bs_IntronRegion,  bs_5UTRsRegion, ncol = 1)
  p3_d_2 = ggarrange(bs_PromoterRegion,bs_SpliceSiteRegion,  bs_3UTRsRegion, ncol = 1)
  p3_d = ggarrange(p3_d_1,NULL,p3_d_2, nrow = 1, widths = c(1,0.00001,1))
  
  # p3_d = ggarrange(bs_IntergenicRegion,bs_PromoterRegion,bs_IntronRegion,bs_SpliceSiteRegion,bs_5UTRsRegion,bs_3UTRsRegion,ncol = 2, nrow = 3)
  p3_d
  

}



## 

## E. Intergenic network ------------------------------------------------------- 
if(T){
  library(igraph)
  library(ggraph)
  if(T){
    # layout
    a1 = read_csv('dawn/Kmeans_result/train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.intergenic.C10_R0.22_S2_L2_seed42.graph_layout.csv') %>%
      dplyr::select(Cluster = `Cluster.index`,
                    X.pos, Y.pos, zvalue)
    
    # cluster_risk_pvalue_table
    a2 = read_csv('dawn/Kmeans_result/DAWN.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.intergenic.C10_R0.22_S2_L2_seed42.cluster_risk_pvalue_table.csv') %>%
      dplyr::select(Cluster = `Cluster.idx`,
                    Risk,
                    Pvalue)
    # cluster_dawn_fdr_table
    a3 = read_csv('dawn/Kmeans_result/DAWN.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.intergenic.C10_R0.22_S2_L2_seed42.cluster_dawn_fdr_table.csv') %>%
      dplyr::select(Cluster = Name,
                    FDR)
    
    m0 = merge(a2, a3, by = 'Cluster')
    m1 = merge(a1, m0, by = 'Cluster')
    
    if(T){
      adjacency = read.delim(file = 'dawn/Kmeans_result/DAWN.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.intergenic.C10_R0.22_S2_L2_seed42.adjacency_matrix.csv',
                             sep = ',', header = F)
      adjacency = adjacency[2:nrow(adjacency),]
      cols = adjacency[,1]
      adjacency = as.matrix(adjacency[,2:ncol(adjacency)])
      colnames(adjacency) = cols
      rownames(adjacency) = cols
      g = graph_from_adjacency_matrix(adjacency, mode = "undirected")
      igraph::V(g)$name = cols
    }
    
    zvalues = m1$zvalue
    normalized_zval <- (zvalues - min(zvalues)) / (max(zvalues) - min(zvalues))
    igraph::V(g)$norm_z = normalized_zval
    
    sig_clusters = m1 %>% filter(m1$FDR<0.05) %>%
      filter(Risk>1) %>% pull(Cluster)
    keep_clusters = sig_clusters
    
    cl_list = list()
    # Get the neighbors of the selected nodes (first-degree neighbors)
    for(k in as.character(keep_clusters)){
      # get neighbors
      idx = adjacency[,paste(k)]
      idx2 = which(idx==1)
      cl_list[[k]] =  colnames(adjacency)[idx2] %>% unique() %>% sort()
    }
    
    # g1 = intersect(cl_list[['145']], cl_list[['186']])
    # g2 = intersect(g1, cl_list[['188']])
    
    # neigh = m1 %>%
    #   filter(Risk>1) %>%
    #   filter(Cluster %in% g2) %>%
    #   arrange(FDR) %>% head(8) %>%
    #   pull(Cluster)
    
    
    # neigh = cl_list[[k]]
    
    neigh = unlist(cl_list)
    
    #as.numeric(neigh) %in% sig_clusters
    
    # filter cluster
    if(T){
      nodes_to_keep <- c(keep_clusters, as.numeric(neigh)) %>% unique() %>% sort() %>% as.character()
      #nodes_to_keep = nodes_to_keep[nodes_to_keep %in% sig_clusters]
    }
    
    # Create a subgraph with the selected nodes and their first-degree neighbors
    if(T){
      name_list = igraph::V(g)$name
      set.seed(1234)
      random_10 = sample(name_list, 15)
      nodes_to_keep = c(nodes_to_keep,random_10)
      # non_match = name_list[!(name_list %in% m1[m1$Pvalue<0.05,]$Cluster)]
      non_match = name_list[!(name_list %in% nodes_to_keep)]
    }
    # name_list = igraph::V(g)$name
    # non_match = name_list[!(name_list %in% nodes_to_keep)] 
    # 
    # delete manually
    subgraph = delete_vertices(g, v = as.character(non_match))
    plot(subgraph)
    
    df = igraph::as_data_frame(subgraph)
    
    
    layout = create_layout(subgraph, layout = 'kk')
    
    ### draw
    if(T){
    set.seed(1234)
    if(T){
      layout2 = merge(layout, m1 %>%
                        mutate(name = Cluster),
                      by = 'name')
    }else{
      layout2 = merge(layout, m1 %>%
                        mutate(name = Cluster),
                      by = 'name')
      layout2$x = layout2$X.pos
      layout2$y = layout2$Y.pos
    }
    
    
    # cluster size
    if(T){
      cluster = read_csv('dawn/Kmeans_result/DAWN.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.intergenic.C10_R0.22_S2_L2_seed42.cluster_annotation.csv')
      size_df = data.frame(stringsAsFactors = F)
      for(k in colnames(cluster)){
        cats = cluster[,as.character(k)]
        cats = cats[!is.na(cats)]
        size_df = rbind(size_df, data.frame(name = k,
                                            size = length(cats)))
      }
    }
    
    layout2 = merge(layout2, size_df,
                    by = 'name')
    
    p3_e = ggraph(subgraph, layout = layout2) +
      geom_edge_link(color = 'grey', alpha = 1, width=0.1) +
      geom_node_point(aes(fill = norm_z, size = size), shape = 21, color = "black", show.legend = T, stroke = 0.2) +
      
      geom_point(aes(x = x, y = y,
                     fill = norm_z, size = size),
                 data = layout2 %>% filter(FDR<0.05),
                 shape = 21, color = "black", show.legend = F, stroke = 1) +
      # geom_node_label(aes(label = name)) +
      theme_classic() +
      scale_size(range = c(2,6),
                 breaks = c(2, 6, 8)) +
      scale_edge_width(guide = 'none', range = c(0.2, 0.8)) + 
      scale_fill_gradient(
        low = 'white',
        high = '#FF0000',
        space = "Lab",
        limits = c(0,1)) +
      labs(fill = 'Z-score', size = 'Number of\ncategories',
           title = "") + #Intergenic category subnetwork
      theme_step1() +
      #geom_text(aes(x = x, y = y,
      #              label = name)) +
      theme(axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),axis.title.y = element_blank()
              ) 
    p3_e
    }
    
    c_idx = colnames(cluster) %in% layout2$name
    a0 = cluster[,c_idx]
  }
}


## F. correlation ------------------------------------------------------------------------------------------

pdf(file = "~/Dropbox/SMC_AD_WGS_paper/Figures/old/Figure3f_v2.pdf", width = 13, height = 10)
# 
# # dawn_corr_matrix = read.delim('~/Dropbox/SMC_AD_WGS_paper/Data/SMC_cwas_results_20240416/DAWN_L2.correlation_matrix.txt')
# 
p <- corrplot(dawn_corr_matrix,
              # method="color",
               order = 'hclust',
               outline = T,
              # addrect = 6,
               cl.lim = c(0, 1),is.corr=FALSE,
               col = colorRampPalette(c("white", "#CE5A67"))(500),
               tl.col = "black", tl.srt = 45,
               tl.cex= 1.2,
               hclust.method = "complete", #complete average
) %>% corrRect(c(1,3)) %>% corrRect(c(4,7)) %>% corrRect(c(7,11)) %>%corrRect(c(12,15))# %>% corrRect(c(22,24)) %>% corrRect(c(27,28)) %>% corrRect(c(30,31))
p
dev.off()

# "ward.D", "ward.D2", "single", "complete", "average", "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
{
#    p3_f_1 <-   ggcorrplot(t(dawn_corr_matrix),
#                      method = 'circle',
#                 hc.order = TRUE,
#                 outline.col = "black",
#                 ggtheme = theme_gray(),
#                 colors = c("#6D9EC1", "white", "#CE5A67"),
#                 
#                 show.legend = TRUE,
#                 hc.method = "complete", # mcquitty, ward.D
#                 pch.col = "black",
#                 pch.cex = 0) + #mcquitty
#   theme_step1() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.ticks.y=element_blank()) +
#      theme(plot.margin = margin(2, 1, 2, 1, "cm"))
#    p3_f_1
# 
#   # Add a black rectangle around the entries for variables 3 to 5
# # Adjust the coordinates if necessary after checking the plot
# p3_f = p3_f_1 +
#   geom_rect(aes(xmin = 0.5, xmax = 3.5, ymin = 0.5, ymax = 3.5),
#               color = "black", fill = NA, size = 1)  +
#   geom_rect(aes(xmin = 3.5, xmax = 7.5, ymin = 3.5, ymax = 7.5),
#             color = "black", fill = NA, size = 1) +
#   geom_rect(aes(xmin = 6.5, xmax = 11.5, ymin = 6.5, ymax = 11.5),
#             color = "black", fill = NA, size = 1) +
#   # geom_rect(aes(xmin = 8.5, xmax = 11.5, ymin = 8.5, ymax = 11.5),
#   #           color = "black", fill = NA, size = 1) +
#   geom_rect(aes(xmin = 11.5, xmax = 15.5, ymin = 11.5, ymax = 15.5),
#             color = "black", fill = NA, size = 1)
# p3_f

}



## H. cluster 73 -------------------------------------------------------------------

# C73
{
  dawn_c = "Cluster73"

list = read_tsv("table.DAWN_L2.RNV_carrier_list.txt")
sig_cluster = list %>% dplyr::select(-1,-2,-3) %>% colnames()

print(paste0("",dawn_c))

pheno_cwas = list

# pheno_age = readxl::read_xlsx("~/Dropbox/ADWGS/sample_information/WGS_1824_phenotype_update.231218.xlsx",sheet = 6) %>% select(WGS_SerialNo., age)

pheno = read_tsv("~/Dropbox/SMC_AD_WGS_paper/Data/WGS_1824_phenotype_update.240416.txt") %>% 
  # mutate(sex = ifelse(sex == "M", 1,0)) %>% 
  filter(RCL40_visual %in% c(1,0) & DX %in% c("DAT","MCI","CU") & WGS_SerialNo. %in% pheno_cwas$SAMPLE) #%>% left_join(pheno_age)

pheno$sex <- as.factor(pheno$sex)
pheno$batch <- as.factor(pheno$batch)
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

b$is_carrier = as.factor(b$is_carrier)
b$sex = as.factor(b$sex)

# Verbal memory
model = lm( Verbal_memory ~ is_carrier + education + age + sex , data=b)

d1 <- summary(model)$coefficients["is_carrier1", ][4]

model = lm( Visual_memory ~ is_carrier + education + age + sex , data=b)

d2 <- summary(model)$coefficients["is_carrier1", ][4]

b$carrier %>% table()


# tmp  = b %>% pivot_longer(cols = c(Verbal_memory, Visual_memory), names_to = 'phenotype_domain', values_to = 'phenotype_score')

p3_h = b %>%
  mutate(carrier = case_when(
    carrier == 'Non-carriers' ~ 'Non-carriers\n(n=876)',
    carrier == 'RNV carriers' ~ 'C73 carrier\n(n=49)'
  )) %>% 
  mutate(carrier = factor(carrier, levels = c('Non-carriers\n(n=876)', 'C73 carrier\n(n=49)'))) %>%
  ggplot(aes(x=carrier, y=Verbal_memory)) +   #theme_step1() + 
  geom_violin(aes(fill=carrier),width=0.7) + geom_boxplot(width=0.2) + #geom_jitter(alpha=.1) +
  geom_signif(tip_length = 0, xmin=1, xmax=2, annotations= c(paste0('P = ', signif(d1, 2))), y_position = max(b$Verbal_memory +0.1,na.rm = T), textsize = 5, family = 'Arial') +
  scale_fill_manual(values = c(green, red))  +
  theme(axis.title.x = element_blank(),
        legend.position = 'None') +
  labs(#title = "Cluster 192", 
    y = 'Verbal memory',x = "")  +
  theme_step1()+
  theme(plot.margin = margin(0, 2, 0, 2, "cm"))
p3_h

# p3_h_2 = b %>%
#   mutate(carrier = case_when(
#     carrier == 'Non-carriers' ~ 'Non-carriers\n(n=876)',
#     carrier == 'RNV carriers' ~ 'C73 carrier\n(n=49)'
#   )) %>% 
#   mutate(carrier = factor(carrier, levels = c('Non-carriers\n(n=876)', 'C73 carrier\n(n=49)'))) %>%
#   ggplot(aes(x=carrier, y=Visual_memory)) +   #theme_step1() + 
#   geom_violin(aes(fill=carrier),width=0.7) + geom_boxplot(width=0.2) + #geom_jitter(alpha=.1) +
#   geom_signif(tip_length = 0, xmin=1, xmax=2, annotations= c(paste0('P = ', signif(d2, 2))), y_position = max(b$Visual_memory +0.1,na.rm = T), textsize = 5, family = 'Arial') +
#   scale_fill_manual(values = c(green, red))  +
#   theme(axis.title.x = element_blank(),
#         legend.position = 'None') +
#   labs(#title = "Cluster 192", 
#     y = 'Visual memory',x = "")  +
#   # theme(plot.margin = margin(0, 2, 0, 2, "cm"))+
#   theme_step1()
# p3_h_2

}


### J. DEG plot -----------------------------------------------------------------
{
if(F){
  fn_DEG <- "~/Desktop/AD_WGS/Annotation/Matyhs_2023_DecontX_DEG.txt.gz"
df_DEG_orig <- data.table::fread(fn_DEG)
df_DEG = df_DEG_orig

# df_DEG = df_DEG %>% filter(!cluster_id %in% c('SMC', 'Fib', 'CAMs', 'T cells',
#                                               'Exc' ,'Inh',  "Inh LAMP5", "Inh PAX6" , "Inh PVALB" ,"Inh SST"  , "Inh VIP", "Exc IT"  ))

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


gene_list = c( 'SAMD5', 'SLC1A1', 'PLGRKT', 'TBCCD1', 'LPP', 'VAPB', 'EIF4A2')
# gene_list = c('MIR3123', 'MBL1P', 'NTS', 'LOC105369877', 'LOC107985156', 'DTD2', 'LOC105370438', 'LOC105378100', 'LOC105369171','C6orf118', 'VAPB', 'LOC107987284', 'ATP5F1E', 'MIR1248', 'RFC4', 'SNORA4', 'SNORA63', 'SNORA63B', 'SNORA81', 'IGF2BP2', 'EIF4A2', 'SNORD2', 'SLC1A3', 'LOC105378043', 'QKI', 'LOC105378115', 'LOC105378104', 'CAHM', 'PACRG-AS1', 'LOC72968', 'LOC105378116')


df_DEG_t = df_DEG %>% filter(gene %in% gene_list)

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

df_DEG_t <- df_DEG_t %>% mutate(cluster_id = case_when(
  cluster_id == 'Inh LAMP5 NRG1 (Rosehip)' ~ 'Inh LAMP5 NRG1 \n(Rosehip)',
  cluster_id == 'Inh PVALB CA8 (Chandelier)' ~ 'Inh PVALB CA8 \n(Chandelier)',
  cluster_id == 'Inh ENOX2 SPHKAP' ~ 'Inh ENOX2 \nSPHKAP',
  T ~ cluster_id
))


unique_values <- c('Inh L1-6 LAMP5 CA13', 
                   'Inh LAMP5 NRG1 \n(Rosehip)',
                   'Inh PVALB HTR4',
                   'Inh PVALB SULF1',
                   'Inh PVALB CA8 \n(Chandelier)' , 
                   'Inh VIP CLSTN2', 
                   'Inh VIP ABI3BP',
                   'Inh ENOX2 \nSPHKAP')

# unique_values <- c('Inh L1-6 LAMP5 CA13', 
#                    'Inh LAMP5 NRG1 (Rosehip)',
#                    'Inh PVALB HTR4',
#                    'Inh PVALB SULF1',
#                    'Inh PVALB CA8 (Chandelier)' , 
#                    'Inh VIP CLSTN2', 
#                    'Inh VIP ABI3BP',
#                    'Inh ENOX2 SPHKAP')

# 요인(factor)으로 변환
df_DEG_t$cluster_id <- factor(df_DEG_t$cluster_id, levels = unique_values)

phenotype_target = c('Overall amyloid level', 
                           'Neurofibrillary tangle burden', 
                           'Neuritic plaque burden', 
                           'Global measure of neocortical AD pathology', 
                           'Global AD pathology', 
                           'Final cognitive consensus diagnosis', 
                           'Working memory', 
                           'Semantic memory', 
                           'Rate of change of working memory', 
                           'Rate of change of semantic memory', 
                           'Rate of change of perceptual speed', 
                           'Rate of change of global cognitive function', 
                           'Rate of change of episodic memory', 
                           'Perceptual speed', 
                           'Global cognitive function', 
                           'Episodic memory')

df_DEG_t$significance_label <- ifelse(df_DEG_t$p_adj.loc < 0.001 & abs(df_DEG_t$logFC) > 0.02, "***",
                                      ifelse(df_DEG_t$p_adj.loc < 0.01 & abs(df_DEG_t$logFC) > 0.02, "**",
                                             ifelse(df_DEG_t$p_adj.loc < 0.05 & abs(df_DEG_t$logFC) > 0.02, "*", "")))

write_csv(df_DEG_t,'../Cluster73_DEG.csv')

}

df_DEG_t = read_csv('../Cluster73_DEG.csv')

phenotype_target = c(
                     'Global AD pathology', #
                     'Overall amyloid level', # 
                     'Neuritic plaque burden', #
                     'Neurofibrillary tangle burden', #
                     'Global measure of neocortical AD pathology', #
                     'Final cognitive consensus diagnosis', #
                     'Global cognitive function', 
                     'Episodic memory',
                     'Semantic memory', 
                     'Perceptual speed', 
                     'Working memory', 
                     'Rate of change of global cognitive function', 
                     'Rate of change of episodic memory' ,
                     'Rate of change of semantic memory', 
                     'Rate of change of perceptual speed', 
                     'Rate of change of working memory' )


df_DEG_t <- df_DEG_t %>% filter(phenotype %in% phenotype_target) 

df_DEG_t$phenotype <- factor(df_DEG_t$phenotype , levels = rev(phenotype_target))

id_target = c('Inh L1-6 LAMP5 CA13', 
              'Inh LAMP5 NRG1 \n(Rosehip)',
              'Inh PVALB CA8 \n(Chandelier)' , 
              'Inh PVALB HTR4',
              'Inh PVALB SULF1',
              'Inh VIP ABI3BP', 
              'Inh ENOX2 \nSPHKAP', 
              'Inh VIP CLSTN2' ) 



# df_DEG_t$phenotype <- factor(df_DEG_t$phenotype , levels = rev(id_target))


p3_j = df_DEG_t %>% 
  filter(
  ( gene %in% c('VAPB') & cluster_id %in% c('Inh VIP ABI3BP' ) ) |
    ( gene %in% c('EIF4A2') & cluster_id %in% c('Inh PVALB CA8 \n(Chandelier)', 'Inh L1-6 LAMP5 CA13', 'Inh PVALB SULF1', 'Inh PVALB HTR4' ) ) |
    ( gene %in% c('TBCCD1') & cluster_id %in% c('Inh PVALB SULF1', 'Inh PVALB HTR4' ) ) |
    ( gene %in% c('LPP') & cluster_id %in% c('Inh PVALB SULF1', 'Inh PVALB HTR4') ) |
    ( gene %in% c('SAMD5') & cluster_id %in% c('Inh ENOX2 \nSPHKAP') ) |
    ( gene %in% c('SLC1A1') & cluster_id %in% c('Inh L1-6 LAMP5 CA13', 'Inh LAMP5 NRG1 \n(Rosehip)') )|
    ( gene %in% c('PLGRKT') & cluster_id %in% c('Inh VIP CLSTN2' , 'Inh LAMP5 NRG1 \n(Rosehip)') )
) %>% 
  filter(cluster_id %in% c('Inh VIP ABI3BP', 
                           'Inh PVALB CA8 \n(Chandelier)' , 'Inh L1-6 LAMP5 CA13', 'Inh PVALB HTR4',
                           'Inh ENOX2 \nSPHKAP', 'Inh LAMP5 NRG1 \n(Rosehip)',
                           'Inh VIP CLSTN2', 'Inh PVALB HTR4', 'Inh PVALB SULF1') ) %>%
  mutate(cluster_id = factor(cluster_id, levels = c('Inh L1-6 LAMP5 CA13','Inh LAMP5 NRG1 \n(Rosehip)',
                                                    'Inh PVALB HTR4',  'Inh PVALB SULF1','Inh PVALB CA8 \n(Chandelier)' , 
                                                    'Inh ENOX2 \nSPHKAP', 
                                                    'Inh VIP CLSTN2',  
                                                    'Inh VIP ABI3BP'))) %>% 
  mutate(logFC = ifelse(logFC > 1, 1,
                        ifelse(logFC < -1, -1 , logFC) )) %>%
  ggplot(aes(x = gene, 
             y = phenotype,
             # y = gene, 
             fill = logFC)) +
  geom_tile(color = "grey", lwd = 1.5, linetype = 1) +
  geom_text(aes(label = significance_label, fontface = "bold", angle = 90), color = "black", size = 5, hjust = 0.5, vjust =0.8,family = 'Arial') +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_gradient2(low = green, mid = "white", high = red , limits = c(-1.1,1.1)
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
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))+
  theme(strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour = "black", fill = "white", 
                                        size = 1.5, linetype = "solid"))
p3_j
}
# ggsave('~/Dropbox/SMC_AD_WGS_paper/Figures/old/DEG.pdf', p3_i, width = 25, height = 10)

### Plotting  ------------------------------------------------------------------

{
  # overview
p3_a = NULL
# p3_e = NULL
p3_f = NULL
  
p3_g = NULL
# p3_h_3 = NULL
p3_i = NULL
# p3_i = NULL

# heights = c(1.2,1.2,1.1,1)
# widths = c(1)

p3_bc = ggarrange(p3_b,p3_c, ncol = 1, nrow = 2, labels = c('b','c'), font.label = list(size = 28), label.y = 1.01,heights = c(1,1.3) )

p3_1 = ggarrange(p3_a, p3_bc, p3_d,  ncol = 3, nrow = 1, labels = c('a','','d'), font.label = list(size = 28), label.y = 1.01,
                 widths = c(1,1,1.5))

p3_gh =  ggarrange(p3_g,p3_h,  ncol = 1, nrow = 2, labels = c('g','h'), font.label = list(size = 28), label.y = 1.01,heights = c(1.1,1)  )

p3_2 = ggarrange(p3_e, p3_f ,p3_gh, ncol = 3, nrow = 1, labels = c('e','f',''), font.label = list(size = 28), label.y = 1.01,
                 widths = c(0.8,1.5,0.8))

# p3_h =  ggarrange(p3_h_1, p3_h_2, ncol = 1, nrow = 2, labels = c('h',''), font.label = list(size = 28), label.y = 1.01 )

p3_3 = ggarrange(p3_i, p3_j, nrow = 1, ncol = 2, labels = c('i','j'), font.label = list(size = 28), label.y = 1.01, widths = c(1,1.2))

# p3_g = ggarrange(p3_g_1, p3_g_2, p3_g_3, nrow = 1,widths = c(1,1,2) )
# p3_h = ggarrange(p3_h_1, p3_h_2, p3_h_3, nrow = 1,widths = c(1,1,2))
# 
# p3_gh = ggarrange(p3_g,p3_h, ncol = 1, labels = c('g','h'), font.label = list(size = 28), label.y = 1.01)
# 
# p3_2 = ggarrange(p3_f,p3_gh, nrow = 1, labels = c('f',''), font.label = list(size = 28), label.y = 1.01,widths = c(1,2))


# p3 = ggarrange(p3_1,p3_2,p3_3, nrow = 3, ncol = 1, height = c(1.2,1,1))
#   
# ggsave(
#   filename = paste('../../Figures/Figure3', version, 'pdf', sep = '.'),
#   plot = p3,
#   width = 30,
#   height = 30,
#   device = "pdf"
# )

p3 <- plot_grid(p3_1, p3_2, p3_3, nrow = 3, ncol = 1, align = 'h', hjust = 1, vjust = 1,
                rel_heights = c(1.2,1,1))

# 플롯 저장
ggsave(
  filename = paste('../../Figures/Figure3', version, 'pdf', sep = '.'),
  plot = p3,
  width = 23,
  height = 29,
  device = "pdf"
)

}

