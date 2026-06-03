library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(tidyr)
library(ggplot2)
library(gggenes)
library(grid)

package_version(tidyvesse)
packageVersion("rnaturalearthdata")
library(sf)
library(rnaturalearth)
library(ggplot2)
library(ggspatial)

install.packages(c("sf","rnaturalearth","rnaturalearthdata","ggspatial"))



setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/bandage/final_pattern")

# ==================================================
# 0) input
# ==================================================
infile <- "region2.final_compare.expanded.tsv"

df <- read_tsv(
  infile,
  show_col_types = FALSE
)

# ==================================================
# 1) style
# ==================================================
token_widths <- c(
  "RHD"     = 1.45,
  "TMEM50A" = 0.42,
  "RHCE"    = 1.55,
  "RHD2"    = 0.90
)

gene_cols <- c(
  "RHD"     = "#B7BEA3",
  "TMEM50A" = "#A9D6E5",
  "RHCE"    = "#EFE9C9"
)

gap_touch  <- 0.02
gap_half   <- 0.16
gap_one    <- 0.34
contig_gap <- 0.36
x_start0   <- 7.45

# ==================================================
# 2) helpers
# ==================================================
fmt_freq3 <- function(x) {
  ifelse(is.na(x), NA_character_,
         ifelse(x == 0, "0", sprintf("%.3f", x)))
}

is_multicontig <- function(x) {
  str_detect(x, fixed(" + "))
}

# ==================================================
# 3) haplotype label 정리
#    - R/L 제거
#    - RHD2 -> RHD
#    - 특정 패턴은 이름 치환
#    - -> 대신 -
# ==================================================
clean_haplotype_label <- function(x) {
  groups <- str_split(x, fixed(" + "))
  
  map_chr(groups, function(grps) {
    cleaned_groups <- map_chr(grps, function(g) {
      toks <- str_split(g, fixed("->"))[[1]] %>% str_trim()
      key <- paste(toks, collapse = "->")
      
      if (key == "RHD->TMEM50A->L->R->RHCE") {
        return("RHD-TMEM50A(inv)-RHCE")
      }
      
      if (key == "RHD->R->RHD2->L->TMEM50A->RHCE") {
        return("RHD-Sub-TMEM50A-RHCE")
      }
      
      toks <- toks[!toks %in% c("R", "L")]
      toks <- ifelse(toks == "RHD2", "RHD", toks)
      
      paste(toks, collapse = "-")
    })
    
    paste(cleaned_groups, collapse = " + ")
  })
}

# ==================================================
# 4) sort
# ==================================================
df2 <- df %>%
  mutate(
    final_class_expanded = str_trim(final_class_expanded),
    haplotype_clean = clean_haplotype_label(final_class_expanded),
    kpp_zero = KPPD_count == 0,
    multicontig_flag = is_multicontig(final_class_expanded)
  ) %>%
  arrange(
    kpp_zero,
    desc(KPPD_freq),
    desc(KPPD_count),
    desc(HPRC_freq),
    desc(HPRC_count)
  ) %>%
  mutate(
    row_id = row_number(),
    y_label = haplotype_clean,
    y = factor(
      paste0(haplotype_clean, "__", row_id),
      levels = rev(paste0(haplotype_clean, "__", row_id))
    )
  )

# ==================================================
# 5) parser
# ==================================================
parse_path_groups <- function(path_string) {
  groups <- str_split(path_string, fixed(" + "))[[1]]
  groups <- str_trim(groups)
  
  tibble(
    group_id = seq_along(groups),
    group_string = groups
  ) %>%
    mutate(
      tokens = map(group_string, ~ str_split(.x, fixed("->"))[[1]] %>% str_trim())
    )
}

# ==================================================
# 6) token pattern -> visible genes + gaps
#    R/L은 spacing rule만 담당
#    유전자 위치는 이전 버전처럼 유지
# ==================================================
convert_group_to_layout <- function(tokens) {
  key <- paste(tokens, collapse = "->")
  
  gene <- function(name, forward = TRUE, label = name) {
    tibble(
      gene_raw = name,
      gene_draw = label,
      forward = forward
    )
  }
  
  make_layout <- function(visible, gaps) {
    list(visible = visible, gaps = gaps)
  }
  
  if (key == "RHD->R->L->TMEM50A->RHCE") {
    visible <- bind_rows(
      gene("RHD", TRUE, "RHD"),
      gene("TMEM50A", TRUE, "TMEM50A"),
      gene("RHCE", FALSE, "RHCE")
    )
    return(make_layout(visible, c(gap_one, gap_touch)))
  }
  
  if (key == "RHD->TMEM50A->L->R->RHCE") {
    visible <- bind_rows(
      gene("RHD", TRUE, "RHD"),
      gene("TMEM50A", FALSE, "TMEM50A"),
      gene("RHCE", FALSE, "RHCE")
    )
    return(make_layout(visible, c(gap_touch, gap_one)))
  }
  
  if (key == "TMEM50A->RHCE") {
    visible <- bind_rows(
      gene("TMEM50A", TRUE, "TMEM50A"),
      gene("RHCE", FALSE, "RHCE")
    )
    return(make_layout(visible, c(gap_touch)))
  }
  
  if (key == "RHCE") {
    visible <- bind_rows(
      gene("RHCE", FALSE, "RHCE")
    )
    return(make_layout(visible, numeric(0)))
  }
  
  if (key == "RHD") {
    visible <- bind_rows(
      gene("RHD", TRUE, "RHD")
    )
    return(make_layout(visible, numeric(0)))
  }
  
  if (key == "RHD->TMEM50A->RHCE") {
    visible <- bind_rows(
      gene("RHD", TRUE, "RHD"),
      gene("TMEM50A", TRUE, "TMEM50A"),
      gene("RHCE", FALSE, "RHCE")
    )
    return(make_layout(visible, c(gap_touch, gap_touch)))
  }
  
  if (key == "R->L->TMEM50A") {
    visible <- bind_rows(
      gene("TMEM50A", TRUE, "TMEM50A")
    )
    return(make_layout(visible, numeric(0)))
  }
  
  if (key == "RHD->R->RHD2->L->TMEM50A->RHCE") {
    visible <- bind_rows(
      gene("RHD", TRUE, "RHD"),
      gene("RHD2", TRUE, "RHD"),
      gene("TMEM50A", TRUE, "TMEM50A"),
      gene("RHCE", FALSE, "RHCE")
    )
    return(make_layout(visible, c(gap_half, gap_half, gap_touch)))
  }
  
  visible_tokens <- tokens[tokens %in% c("RHD", "TMEM50A", "RHCE", "RHD2")]
  
  if (length(visible_tokens) == 0) {
    return(make_layout(tibble(), numeric(0)))
  }
  
  visible <- map_dfr(visible_tokens, function(tok) {
    if (tok == "RHCE") {
      tibble(gene_raw = tok, gene_draw = "RHCE", forward = FALSE)
    } else if (tok == "RHD2") {
      tibble(gene_raw = tok, gene_draw = "RHD", forward = TRUE)
    } else {
      tibble(gene_raw = tok, gene_draw = tok, forward = TRUE)
    }
  })
  
  gaps <- rep(gap_touch, max(0, nrow(visible) - 1))
  make_layout(visible, gaps)
}

# ==================================================
# 7) one full row -> blocks
# ==================================================
make_blocks <- function(path_string, y_val,
                        token_widths,
                        x_start = 7.45,
                        contig_gap = 0.36) {
  
  groups_df <- parse_path_groups(path_string)
  
  cursor <- x_start
  out <- list()
  idx <- 1
  
  for (g in seq_len(nrow(groups_df))) {
    toks <- groups_df$tokens[[g]]
    layout <- convert_group_to_layout(toks)
    
    vis <- layout$visible
    gaps <- layout$gaps
    
    if (nrow(vis) > 0) {
      for (i in seq_len(nrow(vis))) {
        raw_name  <- vis$gene_raw[i]
        draw_name <- vis$gene_draw[i]
        fwd       <- vis$forward[i]
        w         <- unname(token_widths[[raw_name]])
        
        out[[idx]] <- tibble(
          y = y_val,
          gene = draw_name,
          gene_raw = raw_name,
          xmin = cursor,
          xmax = cursor + w,
          forward = fwd,
          group_id = g
        )
        
        cursor <- cursor + w
        idx <- idx + 1
        
        if (i < nrow(vis)) {
          cursor <- cursor + gaps[i]
        }
      }
    }
    
    if (g < nrow(groups_df)) {
      cursor <- cursor + contig_gap
    }
  }
  
  bind_rows(out)
}

gene_df <- purrr::map2_dfr(
  df2$final_class_expanded,
  df2$y,
  ~ make_blocks(
    path_string  = .x,
    y_val        = .y,
    token_widths = token_widths,
    x_start      = x_start0,
    contig_gap   = contig_gap
  )
)

# ==================================================
# 8) multicontig gap marks
# ==================================================
make_gap_marks <- function(path_string, y_val,
                           token_widths,
                           x_start = 7.45,
                           contig_gap = 0.36) {
  
  groups_df <- parse_path_groups(path_string)
  
  if (nrow(groups_df) <= 1) {
    return(tibble())
  }
  
  cursor <- x_start
  out <- list()
  idx <- 1
  
  for (g in seq_len(nrow(groups_df))) {
    toks <- groups_df$tokens[[g]]
    layout <- convert_group_to_layout(toks)
    
    vis <- layout$visible
    gaps <- layout$gaps
    
    if (nrow(vis) > 0) {
      for (i in seq_len(nrow(vis))) {
        raw_name <- vis$gene_raw[i]
        w <- unname(token_widths[[raw_name]])
        cursor <- cursor + w
        if (i < nrow(vis)) {
          cursor <- cursor + gaps[i]
        }
      }
    }
    
    if (g < nrow(groups_df)) {
      out[[idx]] <- tibble(
        y = y_val,
        x = cursor + contig_gap / 2
      )
      cursor <- cursor + contig_gap
      idx <- idx + 1
    }
  }
  
  bind_rows(out)
}

gap_df <- purrr::map2_dfr(
  df2$final_class_expanded,
  df2$y,
  ~ make_gap_marks(
    path_string  = .x,
    y_val        = .y,
    token_widths = token_widths,
    x_start      = x_start0,
    contig_gap   = contig_gap
  )
)

# ==================================================
# 9) baseline / annotation
# ==================================================
max_end <- max(gene_df$xmax, na.rm = TRUE)

baseline_df <- df2 %>%
  mutate(
    x = 7.25,
    xend = max_end + 0.04
  )

anno_df <- df2 %>%
  mutate(
    path_lab = y_label,
    kpp_count_lab = ifelse(multicontig_flag & KPPD_count == 1,
                           paste0(KPPD_count, "*"),
                           as.character(KPPD_count)),
    kpp_freq_lab  = fmt_freq3(KPPD_freq),
    hprc_freq_lab = fmt_freq3(HPRC_freq),
    x_path        = 0.40,
    x_kpp_count   = 3.60,
    x_kpp_freq    = 4.75,
    x_hprc_count  = 5.90,
    x_hprc_freq   = 7.05
  )

header_y <- nlevels(df2$y) + 0.75

header_df <- tibble(
  x = c(0.40, 3.60, 4.75, 5.90, 7.05),
  lab = c("Haplotype", "KPP\nCount", "KPP\nFreq", "HPRC\nCount", "HPRC\nFreq")
)

legend_df <- tibble(
  gene = factor(
    c("RHD", "TMEM50A", "RHCE"),
    levels = c("RHD", "TMEM50A", "RHCE")
  ),
  x = 0,
  y = 0
)

# ==================================================
# 10) plot
# ==================================================
p <- ggplot() +
  geom_text(
    data = anno_df,
    aes(x = x_path, y = y, label = path_lab),
    hjust = 0,
    size = 4
  ) +
  geom_text(
    data = anno_df,
    aes(x = x_kpp_count, y = y, label = kpp_count_lab),
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
    linewidth = 3.0,
    color = "grey82",
    lineend = "butt"
  ) +
  geom_gene_arrow(
    data = gene_df,
    aes(
      xmin = xmin,
      xmax = xmax,
      y = y,
      fill = gene,
      forward = forward
    ),
    color = "black",
    linewidth = 0.3,
    arrowhead_height = unit(5.0, "mm"),
    arrow_body_height = unit(5.0, "mm"),
    arrowhead_width = unit(1.8, "mm"),
    show.legend = TRUE
  ) +
  geom_segment(
    data = gap_df,
    aes(x = x, xend = x, y = y - 0.24, yend = y + 0.24),
    linewidth = 0.7,
    linetype = "dashed",
    color = "grey35"
  ) +
  geom_point(
    data = legend_df,
    aes(x = x, y = y, fill = gene),
    shape = 22,
    size = 4,
    alpha = 0
  ) +
  scale_fill_manual(
    values = gene_cols,
    breaks = c("RHD", "TMEM50A", "RHCE"),
    labels = c("RHD", "TMEM50A", "RHCE"),
    name = NULL
  ) +
  labs(x = NULL, y = NULL) +
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
    plot.margin = margin(35, 10, 10, 10)
  ) +
  coord_cartesian(
    xlim = c(0, max_end + 0.08),
    ylim = c(0.5, header_y + 0.35),
    clip = "off"
  )

p

ggsave("region2_finalpath_ggenes_spacingrule_final_v3.pdf", p, width = 15, height = 6.5)
ggsave("region2_finalpath_ggenes_spacingrule_final_v3.png", p, width = 10, height = 6, dpi = 300)
