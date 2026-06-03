library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(ggplot2)

setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/bandage/final_pattern")

df <- read_tsv(
  "region4.haplotype_representative.canonical.compare.add.tsv",
  show_col_types = FALSE
)

# --------------------------------------------------
# 1) Final_path 기준으로 합치기
# --------------------------------------------------
df_sum <- df %>%
  mutate(
    Final_path = str_trim(Final_path)
  ) %>%
  group_by(region_id, region_name, Final_path) %>%
  summarise(
    KPP_count  = sum(KPPD_count, na.rm = TRUE),
    HPRC_count = sum(HPRC_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    KPP_freq  = KPP_count / sum(KPP_count),
    HPRC_freq = HPRC_count / sum(HPRC_count)
  ) %>%
  arrange(
    KPP_count == 0,
    desc(KPP_freq),
    desc(KPP_count),
    desc(HPRC_freq),
    desc(HPRC_count)
  )

# --------------------------------------------------
# 2) 빈도 표기 함수
#    - 0이면 "0"
#    - 아니면 소수 셋째자리까지
# --------------------------------------------------
fmt_freq3 <- function(x) {
  ifelse(is.na(x), NA_character_,
         ifelse(x == 0, "0", sprintf("%.3f", x)))
}

# --------------------------------------------------
# 3) parser
#    M4 규칙:
#    h  1개 -> 가운데 1개
#    h2/hh -> 1/3, 2/3 지점 2개
# --------------------------------------------------
parse_m4_suffix <- function(s) {
  s <- str_trim(s)
  
  has_a <- str_detect(s, "a")
  has_b <- str_detect(s, "b")
  
  h_count <- case_when(
    str_detect(s, "h2") ~ 2L,
    str_detect(s, "hh") ~ 2L,
    str_detect(s, "h")  ~ 1L,
    TRUE ~ 0L
  )
  
  list(
    has_a = has_a,
    has_b = has_b,
    h_count = h_count
  )
}

# --------------------------------------------------
# 4) gene body 안에 H segment 삽입
# --------------------------------------------------
make_gene_with_inner_h <- function(x_start, width, y_val,
                                   gene_fill = "C4A/B",
                                   h_count = 0L) {
  out <- list()
  
  if (h_count == 0) {
    out[[1]] <- tibble(
      y = y_val,
      xmin = x_start,
      xmax = x_start + width,
      fill_group = gene_fill
    )
    return(bind_rows(out))
  }
  
  if (h_count == 1) {
    h_w <- width * 0.16
    left_w <- (width - h_w) / 2
    
    x1 <- x_start
    x2 <- x1 + left_w
    x3 <- x2 + h_w
    x4 <- x_start + width
    
    out[[1]] <- tibble(
      y = y_val,
      xmin = x1,
      xmax = x2,
      fill_group = gene_fill
    )
    out[[2]] <- tibble(
      y = y_val,
      xmin = x2,
      xmax = x3,
      fill_group = "H"
    )
    out[[3]] <- tibble(
      y = y_val,
      xmin = x3,
      xmax = x4,
      fill_group = gene_fill
    )
    return(bind_rows(out))
  }
  
  if (h_count == 2) {
    h_w <- width * 0.11
    
    h1_center <- x_start + width / 3
    h2_center <- x_start + 2 * width / 3
    
    h1_xmin <- h1_center - h_w / 2
    h1_xmax <- h1_center + h_w / 2
    h2_xmin <- h2_center - h_w / 2
    h2_xmax <- h2_center + h_w / 2
    
    out[[1]] <- tibble(
      y = y_val,
      xmin = x_start,
      xmax = h1_xmin,
      fill_group = gene_fill
    )
    out[[2]] <- tibble(
      y = y_val,
      xmin = h1_xmin,
      xmax = h1_xmax,
      fill_group = "H"
    )
    out[[3]] <- tibble(
      y = y_val,
      xmin = h1_xmax,
      xmax = h2_xmin,
      fill_group = gene_fill
    )
    out[[4]] <- tibble(
      y = y_val,
      xmin = h2_xmin,
      xmax = h2_xmax,
      fill_group = "H"
    )
    out[[5]] <- tibble(
      y = y_val,
      xmin = h2_xmax,
      xmax = x_start + width,
      fill_group = gene_fill
    )
    return(bind_rows(out))
  }
  
  bind_rows(out)
}

# --------------------------------------------------
# 5) Final_path -> block table
# --------------------------------------------------
make_blocks <- function(path_string, y_val) {
  toks <- strsplit(path_string, "->", fixed = TRUE)[[1]]
  
  gene_gap <- 0.14
  
  widths <- c(
    M1 = 1.05,
    M2 = 1.05,
    M3 = 0.95,
    M4_single = 1.05,
    M4_double = 1.55
  )
  
  cursor <- 7.45
  out <- list()
  idx <- 1
  
  for (tok in toks) {
    
    if (tok == "M1") {
      out[[idx]] <- make_gene_with_inner_h(
        x_start = cursor, width = widths["M1"], y_val = y_val,
        gene_fill = "C4A/B", h_count = 0
      )
      cursor <- cursor + widths["M1"] + gene_gap
      idx <- idx + 1
      next
    }
    
    if (tok == "M1H") {
      out[[idx]] <- make_gene_with_inner_h(
        x_start = cursor, width = widths["M1"], y_val = y_val,
        gene_fill = "C4A/B", h_count = 1
      )
      cursor <- cursor + widths["M1"] + gene_gap
      idx <- idx + 1
      next
    }
    
    if (tok == "M2") {
      out[[idx]] <- make_gene_with_inner_h(
        x_start = cursor, width = widths["M2"], y_val = y_val,
        gene_fill = "C4A/B", h_count = 0
      )
      cursor <- cursor + widths["M2"] + gene_gap
      idx <- idx + 1
      next
    }
    
    if (tok == "M2H") {
      out[[idx]] <- make_gene_with_inner_h(
        x_start = cursor, width = widths["M2"], y_val = y_val,
        gene_fill = "C4A/B", h_count = 1
      )
      cursor <- cursor + widths["M2"] + gene_gap
      idx <- idx + 1
      next
    }
    
    if (tok == "M3") {
      out[[idx]] <- make_gene_with_inner_h(
        x_start = cursor, width = widths["M3"], y_val = y_val,
        gene_fill = "C4A/B", h_count = 0
      )
      cursor <- cursor + widths["M3"] + gene_gap
      idx <- idx + 1
      next
    }
    
    if (tok == "M3H") {
      out[[idx]] <- make_gene_with_inner_h(
        x_start = cursor, width = widths["M3"], y_val = y_val,
        gene_fill = "C4A/B", h_count = 1
      )
      cursor <- cursor + widths["M3"] + gene_gap
      idx <- idx + 1
      next
    }
    
    if (startsWith(tok, "M4_")) {
      suffix <- sub("^M4_", "", tok)
      info <- parse_m4_suffix(suffix)
      
      m4_w <- if (info$has_a && info$has_b) widths["M4_double"] else widths["M4_single"]
      
      out[[idx]] <- make_gene_with_inner_h(
        x_start = cursor, width = m4_w, y_val = y_val,
        gene_fill = "C4A/B", h_count = info$h_count
      )
      
      cursor <- cursor + m4_w + gene_gap
      idx <- idx + 1
      next
    }
    
    if (tok == "M4") {
      out[[idx]] <- make_gene_with_inner_h(
        x_start = cursor, width = widths["M4_single"], y_val = y_val,
        gene_fill = "C4A/B", h_count = 0
      )
      cursor <- cursor + widths["M4_single"] + gene_gap
      idx <- idx + 1
      next
    }
  }
  
  bind_rows(out)
}

# --------------------------------------------------
# 6) plotting 함수
# --------------------------------------------------
plot_region4 <- function(dat, title_text = NULL,
                         outfile_prefix = "region4_plot",
                         width = 13.5, height = 7) {
  
  dat2 <- dat %>%
    mutate(
      haplotype_label = str_replace_all(Final_path, "->", "-"),
      kpp_zero = KPP_count == 0
    ) %>%
    arrange(
      kpp_zero,
      desc(KPP_freq),
      desc(KPP_count),
      desc(HPRC_freq),
      desc(HPRC_count)
    ) %>%
    mutate(
      y_id = row_number(),
      y = factor(
        paste0(haplotype_label, "__", y_id),
        levels = rev(paste0(haplotype_label, "__", y_id))
      )
    )
  
  block_df <- purrr::map2_dfr(
    dat2$Final_path,
    dat2$y,
    ~ make_blocks(.x, .y)
  )
  
  max_end <- max(block_df$xmax, na.rm = TRUE)
  
  baseline_df <- dat2 %>%
    mutate(
      x = 7.30,
      xend = max_end + 0.05
    )
  
  anno_df <- dat2 %>%
    mutate(
      path_lab = haplotype_label,
      kpp_freq_lab  = fmt_freq3(KPP_freq),
      hprc_freq_lab = fmt_freq3(HPRC_freq),
      x_path      = 0.40,
      x_kpp_count = 3.50,
      x_kpp_freq  = 4.65,
      x_hprc_count = 5.80,
      x_hprc_freq  = 6.95
    )
  
  header_y <- nlevels(dat2$y) + 0.75
  
  header_df <- tibble(
    x = c(0.40, 3.50, 4.65, 5.80, 6.95),
    lab = c("Haplotype", "KPP\nCount", "KPP\nFreq", "HPRC\nCount", "HPRC\nFreq")
  )
  
  group_cols <- c(
    "C4A/B" = "#8DB33E",
    "H"     = "#C63D2F"
  )
  
  legend_df <- tibble(
    fill_group = factor(c("C4A/B", "H"), levels = c("C4A/B", "H")),
    x = 0, y = 0
  )
  
  # 두께 조금 더 굵게
  block_rect_df <- block_df %>%
    mutate(
      y_num = as.numeric(y),
      ymin = y_num - 0.22,
      ymax = y_num + 0.22
    )
  
  p <- ggplot() +
    geom_text(
      data = anno_df,
      aes(x = x_path, y = y, label = path_lab),
      hjust = 0,
      size = 4
    ) +
    geom_text(
      data = anno_df,
      aes(x = x_kpp_count, y = y, label = KPP_count),
      hjust = 1,
      size = 4
    ) +
    geom_text(
      data = anno_df,
      aes(x = x_kpp_freq, y = y, label = kpp_freq_lab),
      hjust = 1,
      size = 4
    ) +
    geom_text(
      data = anno_df,
      aes(x = x_hprc_count, y = y, label = HPRC_count),
      hjust = 1,
      size = 4
    ) +
    geom_text(
      data = anno_df,
      aes(x = x_hprc_freq, y = y, label = hprc_freq_lab),
      hjust = 1,
      size = 4
    ) +
    geom_text(
      data = header_df,
      aes(x = x, y = header_y, label = lab),
      hjust = c(0, 1, 1, 1, 1),
      fontface = "bold",
      size = 4,
      lineheight = 0.9
    ) +
    geom_segment(
      data = baseline_df,
      aes(x = x, xend = xend, y = y, yend = y),
      linewidth = 3.4,
      color = "grey80",
      lineend = "butt"
    ) +
    geom_rect(
      data = block_rect_df,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = fill_group
      ),
      color = "black",
      linewidth = 0.28
    ) +
    geom_point(
      data = legend_df,
      aes(x = x, y = y, fill = fill_group),
      shape = 22,
      size = 4,
      alpha = 0
    ) +
    scale_fill_manual(
      values = group_cols,
      breaks = c("C4A/B", "H"),
      labels = c("C4A/B", "H"),
      name = NULL
    ) +
    labs(x = NULL, y = NULL, title = title_text) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = c(0.88, 0.98),
      legend.justification = c(1, 1),
      legend.direction = "horizontal",
      legend.background = element_blank(),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.margin = margin(20, 10, 10, 10)
    ) +
    coord_cartesian(
      xlim = c(0, max_end + 0.08),
      ylim = c(0.5, header_y + 0.35),
      clip = "off"
    )
  
  print(p)
  
  ggsave(
    paste0(outfile_prefix, ".pdf"),
    p,
    width = width,
    height = height
  )
  ggsave(
    paste0(outfile_prefix, ".png"),
    p,
    width = width,
    height = height,
    dpi = 300
  )
  
  invisible(p)
}

# --------------------------------------------------
# 7) 데이터 분리
# --------------------------------------------------
df_top <- df_sum %>%
  filter(KPP_count > 2)

df_rest <- df_sum %>%
  filter(KPP_count <= 2)

# --------------------------------------------------
# 8) 그림 생성
# --------------------------------------------------
plot_region4(
  dat = df_sum,
  title_text = "All haplotypes",
  outfile_prefix = "region4_all_rect_hcount_final",
  width = 10,
  height = 6
)

plot_region4(
  dat = df_top,
  title_text = "Top-frequency haplotypes (KPP count > 2)",
  outfile_prefix = "region4_top_rect_hcount_final",
  width = 10,
  height = 6
)

plot_region4(
  dat = df_rest,
  title_text = "Remaining haplotypes (KPP count <= 2)",
  outfile_prefix = "region4_rest_rect_hcount_final",
  width = 10,
  height = 6
)

