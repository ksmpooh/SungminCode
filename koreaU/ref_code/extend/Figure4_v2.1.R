{ 
  # Last updated: 
  # project : Figure 3. 
  
  
  # Setting environment
  rm(list=ls())
  sessionInfo()
  setwd('~/Dropbox/SMC_AD_WGS_working/Data/')
  
  options(stringsAsFactors = F)
  
  # install.packages("~/Downloads/Kmisc_0.5.0.tar.gz", repos = NULL, type = "source")
  
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
  library(broom)
  library(RColorBrewer)   
  library(MetBrewer)
  library(ggpie)
  library(ggbreak)
  
  mutate = dplyr::mutate
  case_when = dplyr::case_when

}


version <- 'v2.2.2_tmp'
  
if(T){
  extrafont::font_import(pattern = "Arial", prompt = F)
  extrafont::loadfonts()
}



# Data prepare
{
  # Set colors
  red = '#D55E00'
  orange = '#E69F00'
  green = '#009E73'
  
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
  
  
  pheno <- read_tsv("WGS_1824_phenotype_update.240416.txt") %>% dplyr::rename(sample = IID)
  pheno$sex <- as.factor(pheno$sex)
  pheno$batch <- as.factor(pheno$batch)
  
  
  
}

## Fig. 4d, e, f, g distribution --------------------------------------------------------
{
  bed = read_csv('STR/EH.1515.293751.bed.csv')
  bed %>% dim
  
  # 데이터프레임에 비율 계산 추가
  bed_summary <- bed %>%
    group_by(`Genomic annotation`) %>% 
    summarise(count = n()) %>%
    mutate(percent = count / sum(count) * 100,
           label = paste(`Genomic annotation`, "(", round(percent, 1), "%)"))
  
  # 원형 차트 그리기
  p4_d <- ggplot(bed_summary, aes(x = "", y = count, fill = `Genomic annotation`, label = label)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    theme_void() +
    theme(legend.position = "right") +
    scale_fill_brewer(palette = "Set3") +
    geom_text(position = position_stack(vjust = 0.5)) +
    labs(fill = "Genomic annotation", y = NULL) +
    guides(fill = guide_legend(title = "Genomic annotation")) +
    theme_step1()+
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank())
  p4_d
  
  fig_4c_1 <- bed %>% 
    mutate(motif_length = nchar(motif)) %>%
    filter(motif_length <= 15) %>%
    group_by(motif_length) %>% 
    summarise(count = n()) %>%
    ggplot(aes(motif_length,count)) +
    theme_step1()+
    geom_bar(stat = 'identity', fill = '#9AC8CD', color = 'black') +
    scale_x_continuous(breaks = 2:15) 
  fig_4c_1
    
  
  fig_4c_2 =  bed %>% 
    mutate(motif_length = nchar(motif)) %>%
    filter(motif_length >= 7 & motif_length <= 18) %>%
    group_by(motif_length) %>% 
    summarise(count = n()) %>%
    ggplot(aes(motif_length,count)) +
    theme_step1()+
    geom_bar(stat = 'identity', fill = '#9AC8CD', color = 'black') +
    scale_x_continuous(breaks = 7:18) 
  fig_4c_2
  
  
  p4_e = fig_4c_1 + annotation_custom(
    ggplotGrob(fig_4c_2), 
    xmin = 5, xmax = 17, ymin = 1000, ymax = 120000
  )
  p4_e
  
  #p4_f
  p4_f = bed %>% 
    filter(ref_length <= 20) %>%
    group_by(ref_length) %>% 
    summarise(count = n()) %>%
    ggplot(aes(ref_length,count)) +
    theme_step1()+
    geom_bar(stat = 'identity', fill = '#9AC8CD', color = 'black') +
    scale_x_continuous(breaks = 2:25) 
  p4_f
  
  median_ref_length = read_tsv("STR/eh.median_ref_length.1515.tsv")
  
  p4_g <- median_ref_length %>% 
    # filter(!id %in% c('chr19-42663234-AT', 'chr11-122018396-TA')) %>%
    filter(median_allele_count < 20 &median_allele_count > -20 ) %>%
    ggplot(aes(median_allele_count)) +
    geom_histogram(bins = 40, fill = '#9AC8CD', color = 'black')+
    # coord_cartesian(ylim = c(0, 200000)) +
    theme_step1() 
  p4_g
  
  
}

## Fig. 4h --------------------------------------------------------------------
{
  {
    # dat.test <- read_tsv('Dropbox/ADWGS/STR/Tables/strling_single_STR_test_QC1.2_DX_noMCI_allele_count_240227.tsv')
    dat.test <- read_tsv('STR/eh_RCL40_visual_single_str_test_inv_noPC_240419_PC_rm.txt')
    
    # dat.test = dat.test %>% mutate(allele_count_inv_norm_FDR__logistic_p = p.adjust(allele_count_inv_norm_logistic_p, method = "fdr"))
    
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
      
      p4_h = df_cond %>% 
        ggplot(aes(standard_mean_difference, -log10(p_val_adj))) + 
        # scale_y_break(c(7.5, 33)) +
        coord_cartesian(ylim = c(0, 15)) +
      geom_point(aes(color = cut(p_val_adj, breaks = c(-Inf, 0.001, 0.005, 0.01,0.05, Inf), 
                                 labels = c("p < 0.001", "0.001 ≤ p < 0.005","0.005 ≤ p < 0.01","0.01 ≤ p < 0.05", "p ≥ 0.05"))),
                 size = 0.8) + 
        scale_color_manual(values = c(
          "p < 0.001" = "#9D0B0B", 
          "0.001 ≤ p < 0.005" = "#DA2D2D", 
          "0.005 ≤ p < 0.01" = "#EB8242", 
          "0.01 ≤ p < 0.05" = "#F6DA63",
          "p ≥ 0.05" = "grey50"), guide = FALSE) + # 색상 설정
        # ggtitle(paste0('By amyloid beta, STR tract lengths \nunder single STR association test')) +
        theme_minimal(base_size = 10) +
        geom_hline(yintercept = -log10(0.05/293751), linetype="dashed", color = 'red') +
        geom_vline(xintercept = 0, linetype="dashed",color = 'grey') +
        # theme(plot.title = element_text(size = 14),
        #       axis.title = element_text(size = 12),
        #       axis.text = element_text(size = 10, color = 'black')) +
        guides(color = guide_legend(title = "p-value"),
               y = guide_axis(title = "-log10(p_val)")) +  # 범례 설정
        theme_step1() +
        # xlim(c(-1,1)) +
        theme(
          # legend.position = c(0.9,0.9)
          legend.position = 'none'
          # legend.position = c(0.3, 0.7)
        ) +
        geom_point(data = df_cond, aes(x = 0.70685109, y = 15)) 
      p4_h
      # p
      
      ggsave(
        filename = paste('../Figures/old/Figure4h', version, 'png', sep = '.'),
        plot = p4_h,
        width = 5,
        height = 5
      )
      
      
      
    }
      {
      
      result_01 <- df_cond %>% filter(p_val_adj < 0.05 & standard_mean_difference != 0) %>% group_by(sign = if_else(standard_mean_difference > 0, "case_enrichment", "control_enrichment")) %>%
        summarise(count = n()) %>% mutate(p_range = "p < 0.05" )
      
      result_005 <- df_cond %>% filter(p_val_adj < 0.01 & standard_mean_difference != 0) %>% group_by(sign = if_else(standard_mean_difference > 0, "case_enrichment", "control_enrichment")) %>%
        summarise(count = n())  %>% mutate(p_range = "p < 0.01" )
      
      result_001 <- df_cond %>% filter(p_val_adj < 0.005 & standard_mean_difference != 0) %>% group_by(sign = if_else(standard_mean_difference > 0, "case_enrichment", "control_enrichment")) %>%
        summarise(count = n())  %>% mutate(p_range = "p < 0.005" )
      
      result_0005 <- df_cond %>% filter(p_val_adj < 0.001 & standard_mean_difference != 0) %>% group_by(sign = if_else(standard_mean_difference > 0, "case_enrichment", "control_enrichment")) %>%
        summarise(count = n())  %>% mutate(p_range = "p < 0.001" )
      # 
      # result_0001 <- df_cond %>% filter(p_val_adj < 0.001) %>% group_by(sign = if_else(standard_mean_difference >= 0, "case_enrichment", "control_enrichment")) %>%
      #   summarise(count = n())  %>% mutate(p_range = "p < 0.001" )
      # 
      # 결과들을 리스트로 만들고 bind_rows로 합침
      result_combined <- bind_rows(result_01, result_005, result_001, result_0005)
      
      
      
      # df_num = df_cond %>%
      #   mutate(level = case_when(
      #     p_val_adj < 0.001 ~ "p < 0.001",
      #     p_val_adj < 0.005 ~ "p < 0.005",
      #     p_val_adj < 0.01 ~ "p < 0.01",
      #     p_val_adj < 0.05 ~ "p < 0.05",
      #     p_val_adj >= 0.05 ~ "p ≥ 0.05" 
      #   )) %>%
      #   group_by(level, sign = if_else(standard_mean_difference >= 0, "case_enrichment", "control_enrichment")) %>%
      #   summarise(count = n()) %>%
      #   ungroup()
      
      # df_num %>% view
    }
    
    #df_num %>% mutate()
    
    df_num <- result_combined %>%
      group_by(p_range) %>%
      mutate(ratio = if_else(sign == "case_enrichment",
                             count / lead(count),
                             count / lag(count))) %>%
      ungroup()
    
    # p ≤ 0.001에 대한 ratio 값을 설정
    # df_num$ratio[which(df_num$level == "p ≤ 0.001")] <- -2.5
    # 
    # new_data <- tribble(
    #   ~level,    ~sign,             ~count,  ~ratio,
    #   "p ≤ 0.001", "case_enrichment", 5,      0
    # )
    
    # 새로운 데이터 추가
    #df_num <- bind_rows(df_num, new_data)
    df_num <- df_num %>% filter(!is.na(sign))
    
    desired_order <- c( "p < 0.05", "p < 0.01", "p < 0.005", "p < 0.001")
    # desired_order <- rev(desired_order)
    
    df_num %>%
      filter(sign == 'case_enrichment') %>%
      mutate(level = factor(p_range, levels = desired_order)) %>%
      ggplot(aes(ratio - 1, level, fill = level)) +  # fill aesthetic을 추가하여 p_range에 따라 색상 지정
      geom_bar(stat = "identity") +
      scale_x_continuous(labels = function(x) ifelse(x == -1, "Inf", x + 1)) +
      #scale_x_continuous(labels = function(x) (x + 1)) +
      scale_fill_manual(values = c("p < 0.001" = "#9D0B0B", 
                                   "p < 0.005" = "#DA2D2D", 
                                   "p < 0.01" = "#EB8242",
                                   "p < 0.05" = "#F6DA63"), guide = FALSE) +
      ggtitle("Ratio of STR length \n(cases enrichment)") +
      guides(color = guide_legend(title = "p-value"),
             x = guide_axis(title = "Association p-value threshold"),
             y = guide_axis(title = "ratio of # of STRs longer in \ncases versus longer in controls"))  +
      theme_step1() +
      coord_flip()
    
  }
}

## Fig. 4i --------------------------------------------------------------------
{
  result_df = read_tsv( "STR/eh_RCL40_visual_single_str_thr_freq_ref_test_PC_rm.txt")
  
  
  result_df$threshold = as.character(result_df$threshold)
  
  p4_i <- result_df %>%
    filter(threshold != 1) %>%
    mutate(threshold = factor(threshold, levels = unique(threshold))) %>%
    ggplot(aes(x = threshold, y = fold_difference - 1, fill = as.factor(frequency))) +
    geom_bar(stat = "identity", position = "dodge",color = "black") +
    scale_fill_manual(values = c("1" = "#F0F5FC", "5" = "#BED3E0", "10" = "#6BA5C8", "100" = "#307BB0", "Inf" = "#124D97")) +
    labs(x = "STR tract length threshold", y = "Fold difference (case/control)", fill = "Max number of observations") +
    scale_y_continuous(labels = function(x)  x + 1, expand = expansion(add = c(0.01, 0.02))) +
    # scale_y_continuous(labels = function(x)  x + 0.9, expand = expansion(add = c(-0.04, 0.05))) +
    theme_step1() +
    # geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
    theme(
      legend.position = "bottom" 
      #axis.line.x = element_blank()
    ) 
  p4_i
}

## Fig. 4j --------------------------------------------------------------------
{
  expan.test.sample <- read_tsv('STR/eh_RCL40_visual_expansion_count_sample_240419_PC_rm_figure.tsv')
  
  
  expan.test.sample$batch <- as.factor(expan.test.sample$batch)
  expan.test.sample$sex <- as.factor(expan.test.sample$sex)
  
  
  d2 = lm(is_case ~ outlier_count + is_female + batch + age + s_avg_depth +
            PC1 + PC2 + PC3, data = expan.test.sample )
  
  d2 = coef(summary(d2))[2,4]
  
  
  p4_j <- expan.test.sample %>%
    filter(outlier_count < 75) %>%
    mutate(Amyloid = ifelse(
      is_case == 1 , "Case", "Control"
    )) %>% 
    ggplot(aes(factor(Amyloid, levels = c("Control", "Case")), outlier_count, group = Amyloid)) +
    geom_violin(aes(fill = Amyloid)) + 
    geom_boxplot(width = 0.25) + 
    geom_signif(comparisons = list(c("Control", "Case")),
                map_signif_level = TRUE, textsize = 5,
                annotations = c(paste0('P = ', signif(d2, 2))),
                y_position = max(75, na.rm = TRUE),
                family = 'Arial') +
    theme(axis.title.x = element_blank(), legend.position = 'None') +
    labs(y = 'STR outlier count per samples', x = "Amyloid beta") +
    scale_fill_manual(values = c("Control" = green, "Case" = red)) +
    ggtitle('STR outlier per samples')+
    theme_step1()
  p4_j
}


## Fig. 4k --------------------------------------------------------------------
{
  n_Case = 895
  n_Control = 620
  
  
  expan.test.sample<- read_tsv("STR/eh_RCL40_visual_expansion_count_sample_240419_PC_rm.txt")
  expan.test.sample <- expan.test.sample %>% left_join(pheno %>% select(sample, RCL_global, KlunkCL))
  
  
  threshold_list <- c(10,20,30,40)
  
  result_df <- data.frame(threshold = numeric(),
                          ratio = numeric(), stringsAsFactors = FALSE)
  
  
  for(i in threshold_list){
    if(i == 10){
      expan.test.sample.case = expan.test.sample %>% filter(outlier_count < i & is_case == 1) %>% pull(sample) %>% unique %>% length()
      expan.test.sample.ctrl = expan.test.sample %>% filter(outlier_count < i & is_case == 0) %>% pull(sample) %>% unique %>% length()
      
      ratio <- (expan.test.sample.case / expan.test.sample.ctrl) / ((n_Case - expan.test.sample.case) / (n_Control - expan.test.sample.ctrl))
      result_df <- rbind(result_df, data.frame(threshold = paste0("<",i), ratio = ratio))
      
      # expan.test.sample.case = expan.test.sample %>% filter(outlier_count >= i & is_case == 1) %>% pull(sample) %>% unique %>% length()
      # expan.test.sample.ctrl = expan.test.sample %>% filter(outlier_count >= i & is_case == 0) %>% pull(sample) %>% unique %>% length()
      #  
      # ratio <- (expan.test.sample.case / expan.test.sample.ctrl) / ((n_Case - expan.test.sample.case) / (n_Control - expan.test.sample.ctrl))
      # result_df <- rbind(result_df, data.frame(threshold = paste0('≥',i), ratio = ratio))
      next
      
    }
    expan.test.sample.case = expan.test.sample %>% filter(outlier_count >= i & is_case == 1) %>% pull(sample) %>% unique %>% length()
    expan.test.sample.ctrl = expan.test.sample %>% filter(outlier_count >= i & is_case == 0) %>% pull(sample) %>% unique %>% length()
    
    ratio <- (expan.test.sample.case / expan.test.sample.ctrl) / ((n_Case - expan.test.sample.case) / (n_Control - expan.test.sample.ctrl))
    result_df <- rbind(result_df, data.frame(threshold = paste0('≥',i), ratio = ratio))
  }
  result_df
  
  thr_lvl = c("<10" ,"≥20","≥30","≥40")
  
  p4_k = result_df %>% 
    mutate(threshold = factor(threshold, levels = thr_lvl)) %>%
    ggplot() +
    geom_bar(aes(threshold ,ratio - 1, fill = threshold),stat = "identity", show.legend = FALSE) +
    geom_text(aes(threshold,label = round(ratio, 2), y = ratio - 1), vjust = -0.5, color = "black", size = 6, family = 'Arial') +
    scale_fill_manual(values = c("<10" = "#F6DA63", 
                                 #"≥10" = "#F6DA63", 
                                 "≥20" = "#EB8242", "≥30" = "#DA2D2D", "≥40" = "#9D0B0B")) +
    scale_y_continuous(labels = function(x)  x + 1, expand = expansion(add = c(0.05, 0.05))) +
    theme_step1() +
    theme(
      legend.position=
        #axis.line.x = element_blank()
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(x = "Burden of STR expansions", y = "Odd ratio for AD", fill = "Max number of observations") 
  p4_k
  
}


## Fig. 4l --------------------------------------------------------------------
{
  expan.test.sample <- expan.test.sample %>%
    # filter(is_case == 1) %>%
    mutate(
      over40 = ifelse(outlier_count >= 40, 'Over_40_STR_outliers',
                      ifelse(outlier_count < 10, 'Under_10_STR_outliers',NA))
    )  #%>%filter(!is.na(over40))
  expan.test.sample %>% dim
  
  expan.test.sample.target = expan.test.sample %>% filter(!is.na(over40))
  expan.test.sample.target$batch <- factor(expan.test.sample.target$batch)
  expan.test.sample.target$sex <- factor(expan.test.sample.target$sex)
  
  ###### RCL_global KlunkCL
  p = summary(lm(  RCL_global ~ over40 + is_female + batch + age + s_avg_depth +
                     PC1 + PC2 + PC3 ,expan.test.sample.target))
  
  p1 = p$coefficients[2,4]
  p1
  
  ###### RCL_global KlunkCL
  p = summary(lm(  RCL_global ~ over40 + is_female + batch + age + s_avg_depth +
                     PC1 + PC2 + PC3 ,expan.test.sample.target %>% filter(DX == 'CU')))
  
  p2 = p$coefficients[2,4]
  p2
  
  
  
  p4_l_1 = expan.test.sample %>%
    filter(!is.na(over40)) %>% 
    ggplot(aes(factor(over40, levels = c("Under_10_STR_outliers", "Over_40_STR_outliers")),RCL_global)) +
    geom_violin(aes(fill = over40), width = 0.7) + 
    geom_boxplot(width = 0.1) + 
    geom_signif(comparisons = list(c("Over_40_STR_outliers", "Under_10_STR_outliers")), 
                map_signif_level = TRUE, textsize = 5, 
                annotations = c(paste0('P = ', signif(p1, 2))), 
                y_position = max(expan.test.sample$RCL_global -10 , na.rm = TRUE), 
                family = 'Arial') +
    theme(axis.title.x = element_blank(), legend.position = 'None') +
    labs(y = 'Amyloid beta level (RCL global)', x = "Amyloid beta") +
    scale_fill_manual(values = c("Over_40_STR_outliers" = red, "Under_10_STR_outliers" = green)) +
    # ggtitle('Amyloid beta difference')+ 
    scale_x_discrete(labels = c("Under_10_STR_outliers" = "Outlier count < 10 \n(n=181)", "Over_40_STR_outliers" = "Outlier count ≥ 40 \n(n=82)"))+
    theme_step1()
  p4_l_1
  
  
  p4_l_2 = expan.test.sample %>%
    filter(!is.na(over40) & DX == 'CU') %>% 
    ggplot(aes(factor(over40, levels = c("Under_10_STR_outliers", "Over_40_STR_outliers")),RCL_global)) +
    geom_violin(aes(fill = over40), width = 0.7) + 
    geom_boxplot(width = 0.1) + 
    geom_signif(comparisons = list(c("Over_40_STR_outliers", "Under_10_STR_outliers")), 
                map_signif_level = TRUE, textsize = 5, 
                annotations = c(paste0('P = ', signif(p2, 2))), 
                y_position = max(expan.test.sample$RCL_global -50, na.rm = TRUE), 
                family = 'Arial') +
    theme(axis.title.x = element_blank(), legend.position = 'None') +
    labs(y = 'Amyloid beta level (KlunkCL)', x = "Amyloid beta") +
    scale_fill_manual(values = c("Over_40_STR_outliers" = red , "Under_10_STR_outliers" =  green)) +
    # ggtitle('Amyloid beta difference') + 
    scale_x_discrete(labels = c("Under_10_STR_outliers" = "<10", "Over_40_STR_outliers" = "≥40")) +
    theme_step1()
  p4_l_2
  
  
  
  
}

# Fig.4m --------------------------------------------------------------------
{
  GO_df <- read_tsv('STR/eh_STR_over40_goterm.tsv')
  
  GO_df$Description <- factor(GO_df$Description, levels = (unique(GO_df$Description)))
  
  p4_m <- GO_df %>%
    tail(10) %>%
    ggplot(aes(-log10(p.adjust),Description)) +
    geom_point(colour = "#124076", size = 7) +
    theme_step1() +
    ylab("") +
    ggtitle("STRs over 40 expansions")
  p4_m
  
  
}

# manhateen
# {
#   
#   res<-data.table::fread("logistic_regression_result")
#   res$SV_start<-as.numeric(res$SV_start)
#   
#   fig3 <- ggmanh::manhattan_plot(x = res, pval.colname = "P", chr.colname = "SV_chrom", pos.colname = "SV_start", ylim = 10, signif = c(2.4e-06, 1.6e-04), chr.col=c("black", "lightgrey"), label.font.size =6) + theme(plot.title = element_text(hjust = 0.5))  + theme_step1()
#   
#   fig3 <- fig3 + theme(plot.title = element_text(hjust = 0.5))  + theme_step1()
#   
#   p_sv =  fig3
#   ggsave( '~/Dropbox/SMC_AD_WGS_paper/Figures/old/Figure4.manh.pdf',p_sv)
#   
# }


## Plotting --------------------------------------------------------------------
{
  # overview
p4_a = NULL
# CNV manh
p4_b = NULL
p4_c = NULL
p4_h = NULL

p4_1 = ggarrange(p4_a,p4_b,p4_c, ncol = 3 ,labels = c('a', 'b','c'), font.label = list(size = 28), label.y = 1.01, widths = c(1.7,1,1))

p4_defg = ggarrange( p4_d,p4_e,p4_g,  p4_f, ncol = 2, nrow = 2, labels = c('d', 'e','f','g'), font.label = list(size = 28), label.y = 1.01, widths = c(1,1.8))
p4_defg

p4_2 = ggarrange(p4_defg, p4_h, p4_i ,ncol = 3 ,labels = c('','h','i'), font.label = list(size = 28), label.y = 1.01, widths = c(1.5,0.8,0.7))

p4_l =  ggarrange( p4_l_1, p4_l_2, ncol = 1 ,font.label = list(size = 28), label.y = 1.01)

p4_3 = ggarrange(p4_j, p4_k, p4_l, p4_m ,ncol = 4 ,labels = c('j', 'k','l','m'), font.label = list(size = 28), label.y = 1.01, widths = c(1,1,0.8,1.7))

# p4_4 = ggarrange(p4_k, ncol = 1 ,labels = c('k'), font.label = list(size = 28), label.y = 1.01, widths = c(1))


p4 = ggarrange(p4_1,p4_2,p4_3,heights = c(1.2,1,1),ncol = 1)

ggsave(
  filename = paste('../Figures/Figure4', version, 'pdf', sep = '.'),
  plot = p4,
  width = 24,
  height = 20,
  device = "pdf"
)


}


