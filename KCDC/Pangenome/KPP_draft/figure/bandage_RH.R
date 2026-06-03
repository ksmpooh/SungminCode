library(tidyverse)
library(gggenes)
library(grid)

theme_step1 <- function(base_size = 11, base_family = "",
                        base_line_size = base_size / 22,
                        base_rect_size = base_size / 22) {
  theme(
    title = element_text(family = "Arial", size = 18, color = "black"),
    text = element_text(family = "Arial", size = 16, color = "black"),
    axis.title = element_text(family = "Arial", size = 18, color = "black"),
    axis.text = element_text(family = "Arial", size = 16, color = "black"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(family = "Arial", size = 14),
    legend.spacing.x = unit(0.1, "cm"),
    plot.margin = unit(c(0.4, 0.8, 0.6, 0.4), "cm"),
    axis.title.y = element_text(margin = margin(r = 10, unit = "pt"))
  )
}

df <- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/bandage/plot/RH_KPPD_HPRC_merged_path.final.compare.tsv")

# --------------------------------------------------
# 1) merged_path -> display label
# --------------------------------------------------
rename_path <- function(x) {
  case_when(
    x == "RHD->middle->TMEM50A->RHCE"   ~ "RHD;TMEM50A-RHCE",
    x == "TMEM50A->RHCE"                ~ "TMEM50A-RHCE",
    x == "RHD->TMEM50A->middle->RHCE"   ~ "RHD-TMEM50A;RHCE",
    x == "RHD->TMEM50A->middle->RHD"    ~ "RHD-TMEM50A;RHD",
    x == "RHD->TMEM50A"                 ~ "RHD-TMEM50A",
    x == "RHCE->middle->TMEM50A->RHCE"  ~ "RHCE;TMEM50A-RHCE",
    x == "RHD->RHCE"                    ~ "RHD-RHCE",
    TRUE                                ~ x
  )
}

# --------------------------------------------------
# 2) 정렬
#    - KPPD freq 내림차순
#    - KPPD=0 아래
#    - unknown 맨 마지막
# --------------------------------------------------
df2 <- df %>%
  mutate(
    merged_path = str_trim(merged_path),
    display_path = rename_path(merged_path),
    is_unknown = merged_path == "unknown",
    kppd_zero = KPPD_count == 0
  ) %>%
  arrange(
    is_unknown,
    kppd_zero,
    desc(KPPD_freq),
    desc(KPPD_count),
    desc(HPRC_freq),
    desc(HPRC_count)
  ) %>%
  mutate(
    y = factor(display_path, levels = rev(display_path))
  )

# --------------------------------------------------
# 3) row 전체 baseline
# --------------------------------------------------
baseline_df <- df2 %>%
  mutate(
    x = 7.20,
    xend = 12.55
  )

# --------------------------------------------------
# 4) path -> gene block
#    - baseline 위에 gene block만 그림
#    - middle은 baseline만 보이게 공간만 확보
#    - TMEM50A는 짧게
# --------------------------------------------------
make_blocks <- function(path_string, y_val) {
  genes <- strsplit(path_string, "->", fixed = TRUE)[[1]]
  
  if (length(genes) == 1 && genes[1] == "unknown") {
    return(tibble())
  }
  
  add_gene <- function(gene, xmin, xmax, forward = TRUE) {
    tibble(
      y = y_val,
      gene = gene,
      xmin = xmin,
      xmax = xmax,
      forward = forward
    )
  }
  
  widths <- c(
    RHD = 1.55,
    middle = 0.95,
    TMEM50A = 0.32,
    RHCE = 1.65
  )
  
  gap <- 0.03
  cursor <- 7.45
  
  out <- list()
  idx <- 1
  
  for (g in genes) {
    if (g == "RHD") {
      out[[idx]] <- add_gene("RHD", cursor, cursor + widths["RHD"], TRUE)
      cursor <- cursor + widths["RHD"] + gap
      idx <- idx + 1
      
    } else if (g == "RHCE") {
      out[[idx]] <- add_gene("RHCE", cursor, cursor + widths["RHCE"], FALSE)
      cursor <- cursor + widths["RHCE"] + gap
      idx <- idx + 1
      
    } else if (g == "TMEM50A") {
      out[[idx]] <- add_gene("TMEM50A", cursor, cursor + widths["TMEM50A"], TRUE)
      cursor <- cursor + widths["TMEM50A"] + gap
      idx <- idx + 1
      
    } else if (g == "middle") {
      cursor <- cursor + widths["middle"] + gap
    }
  }
  
  bind_rows(out)
}

gene_df <- purrr::map2_dfr(
  df2$merged_path,
  df2$y,
  ~ make_blocks(.x, .y)
)

# --------------------------------------------------
# 5) unknown 위치 조정
# --------------------------------------------------
unknown_df <- df2 %>%
  filter(merged_path == "unknown") %>%
  mutate(
    x = (baseline_df$x[1] + baseline_df$xend[1]) / 2
  )

# --------------------------------------------------
# 6) 왼쪽 표
# --------------------------------------------------
anno_df <- df2 %>%
  mutate(
    path_lab = display_path,
    kppd_freq_lab = sprintf("%.2f", KPPD_freq),
    hprc_freq_lab = sprintf("%.2f", HPRC_freq),
    
    x_path       = 0.40,
    x_kppd_count = 3.40,
    x_kppd_freq  = 4.40,
    x_hprc_count = 5.50,
    x_hprc_freq  = 6.60
  )

header_y <- nlevels(df2$y) + 0.75

header_df <- tibble(
  x = c(0.40, 3.40, 4.40, 5.50, 6.60),
  lab = c("Haplotype", "KPP\nCount", "KPP\nFreq", "HPRC\nCount", "HPRC\nFreq")
)

# --------------------------------------------------
# 7) 색상
# --------------------------------------------------
gene_cols <- c(
  "RHD" = "#b7bea3",
  "RHCE" = "#efe9c9",
  "TMEM50A" = "#a9d6e5"
)

# legend용 dummy
legend_df <- tibble(
  gene = factor(c("RHD", "RHCE", "TMEM50A"), levels = c("RHD", "RHCE", "TMEM50A")),
  x = c(0, 0, 0),
  y = c(0, 0, 0)
)

# --------------------------------------------------
# 8) plot
# --------------------------------------------------
p <- ggplot() +
  # 왼쪽 텍스트 표
  geom_text(
    data = anno_df,
    aes(x = x_path, y = y, label = path_lab),
    hjust = 0,
    size = 4
  ) +
  geom_text(
    data = anno_df,
    aes(x = x_kppd_count, y = y, label = KPPD_count),
    hjust = 1,
    size = 4
  ) +
  geom_text(
    data = anno_df,
    aes(x = x_kppd_freq, y = y, label = kppd_freq_lab),
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
  
  # 헤더
  geom_text(
    data = header_df,
    aes(x = x, y = header_y, label = lab),
    hjust = c(0, 1, 1, 1, 1),
    fontface = "bold",
    size = 4,
    lineheight = 0.9
  ) +
  
  # row 전체 baseline
  geom_segment(
    data = baseline_df,
    aes(x = x, xend = xend, y = y, yend = y),
    linewidth = 3.1,
    color = "grey78",
    lineend = "butt"
  ) +
  
  # 유전자 블록
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
    arrowhead_height = unit(3, "mm"),
    arrowhead_width = unit(1.5, "mm"),
    show.legend = TRUE
  ) +
  
  # unknown 텍스트
  geom_text(
    data = unknown_df,
    aes(x = x, y = y, label = "unknown"),
    size = 4
  ) +
  
  # legend 유지용 dummy
  geom_point(
    data = legend_df,
    aes(x = x, y = y, fill = gene),
    shape = 22,
    size = 4,
    alpha = 0
  ) +
  
  scale_fill_manual(
    values = gene_cols,
    breaks = c("RHD", "RHCE", "TMEM50A"),
    labels = c("RHD", "RHCE", "TMEM50A"),
    name = NULL
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 13) +
  theme_step1() +
  theme(
    legend.position = c(0.75, 0.97),
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
    xlim = c(0, 15.0),
    ylim = c(0.5, header_y + 0.35),
    clip = "off"
  )

p

ggsave(
  "/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/bandage/plot/RH_KPPD_HPRC_merged_path.final.onepanel_lefttable.pdf",
  p,
  width = 15,
  height = 6
)