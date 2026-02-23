library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)


theme_step1 <- function(base_size = 11, base_family = "",
                        base_line_size = base_size / 22,
                        base_rect_size = base_size / 22) {
  theme(title = element_text(family = 'Arial', size = 18, color = 'black'), text = element_text(family = 'Arial', size = 16, color = 'black'),
        axis.title = element_text(family = 'Arial', size = 18, color = 'black'), axis.text = element_text(family = 'Arial', size = 16, color = 'black'), 
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = NA), axis.line = element_line(colour = "black", size = rel(1)),
        #legend.background = element_rect(color = 'black'), 
        legend.title = element_text(family = 'Arial', size = 16),
        legend.text = element_text(family = 'Arial', size = 14),
        legend.direction = "vertical", 
        legend.box = c("horizontal", "vertical"),
        legend.spacing.x = unit(0.1, 'cm'),
        plot.margin = unit(c(0.25, 1, 1, 0.5), 'cm'),
        axis.title.y = element_text(margin = margin(r = 10, unit = "pt"))) }


setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/PCA")
df <- read_table("PCA_ethnic.txt")
df <- read_table("PCA_ethnic.txt")
head(df)
sample_1gkp<- read_tsv("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/PCA/1kgp.ethnic.info.tsv")
head(sample_1gkp)
dim(sample_1gkp)
sample_1gkp %>% select(1,7) %>%
  mutate(
    superpop5 = case_when(
      str_detect(`Superpopulation name`, "African")    ~ "AFR",
      str_detect(`Superpopulation name`, "American")   ~ "AMR",
      str_detect(`Superpopulation name`, "East Asia")  ~ "EAS",
      str_detect(`Superpopulation name`, "European")   ~ "EUR",
      str_detect(`Superpopulation name`, "South Asia") ~ "SAS",
      TRUE                                             ~ NA_character_
    )
  ) -> sample_1gkp
colnames(sample_1gkp) <- c("ID","rawg","GROUP") 

dup_iids <- df %>%
  count(IID) %>%
  filter(n > 1) %>%
  pull(IID)

dup_iids

#sample_1gkp
dim(sample_1gkp)
df %>% 
  filter(IID != "CHM13v2") %>%
  mutate(cohort = ifelse(str_detect(IID,"KPPD"),"KPPD","other")) %>%
  ggplot(aes(x=PC1,y=PC2,color=cohort)) + 
  geom_point()

library(dplyr)
library(stringr)
library(ggplot2)

df %>%
  filter(IID != "CHM13v2") %>%
  mutate(
    rownum = row_number(),
    cohort = case_when(
      rownum <= 43            ~ "HPRC",
      str_detect(IID, "HIFI") ~ "CPC",
      str_detect(IID, "RY")   ~ "CPC",
      str_detect(IID, "KPPD") ~ "KPPD",
      rownum <= 226           ~ "HPRC",
      TRUE                    ~ "1KGP"
    ),
    cohort = factor(cohort, levels = c("1KGP", "CPC", "HPRC", "KPPD"))
  ) %>%
  ggplot(aes(x = PC1, y = PC2, color = cohort)) +
  
  ## 1KGP: 가장 뒤 (배경)
  geom_point(
    data = ~ filter(.x, cohort == "1KGP"),
    size = 0.7,
    alpha = 0.6
  ) +
  
  ## CPC: 그 위
  geom_point(
    data = ~ filter(.x, cohort == "CPC"),
    size = 1.6,
    alpha = 0.85
  ) +
  
  ## HPRC: 그 위
  geom_point(
    data = ~ filter(.x, cohort == "HPRC"),
    size = 1.8,
    alpha = 0.9
  ) +
  
  ## KPPD: 최상단
  geom_point(
    data = ~ filter(.x, cohort == "KPPD"),
    size = 2.2,
    alpha = 1
  ) +
  
  scale_color_manual(
    values = c(
      "1KGP" = "grey60",
      "CPC"  = "#2E7D32",
      "HPRC" = "#C00000",
      "KPPD" = "#1F4E9E"
    )
  ) + 
  theme_step1() +
  theme(
    legend.title = element_blank()
  )


duplicated(df$IID) %>% table
df %>%
  filter(IID != "CHM13v2") %>%
  mutate(rownum = row_number(),cohort = case_when(
      rownum <= 43            ~ "HPRC",
      str_detect(IID, "HIFI") ~ "CPC",
      str_detect(IID, "RY")   ~ "CPC",
      str_detect(IID, "KPPD") ~ "KPPD",
      rownum <= 226           ~ "HPRC",
      TRUE                    ~ "1KGP"
    )) %>% duplicated(IID)




###
df<- read_table("PCA_ethnic.kppdip.txt")

head(sample_1gkp)
head(df)
df %>%
  filter(IID != "CHM13v2") %>%
  mutate(
    rownum = row_number(),
    cohort = case_when(
      str_detect(IID, "KPPD") ~ "KPPD",
      IID %in% dup_iids ~ "HPRC",
      TRUE                    ~ "1KGP"
    )
  ) %>% 
  left_join(sample_1gkp %>% rename(IID = ID)) %>%
  mutate(GROUP = ifelse(cohort %in% c("KPPD","HPRC"),cohort,GROUP)) %>% 
  #count(GROUP)
  ggplot(aes(x = PC2, y = PC1, color = GROUP)) +
  geom_point() + 
 
  theme_step1() +
  theme(
    legend.title = element_blank())

library(dplyr)
library(ggplot2)

pop5 <- c("AFR","AMR","EAS","EUR","SAS")

df_pca_plot <- df %>%
  filter(IID != "CHM13v2") %>%
  mutate(
    cohort = case_when(
      str_detect(IID, "KPPD") ~ "KPPD",
      IID %in% dup_iids       ~ "HPRC",
      TRUE                    ~ "1KGP"
    )
  ) %>%
  left_join(sample_1gkp %>% rename(IID = ID), by = "IID") %>%
  mutate(
    GROUP = ifelse(cohort %in% c("KPPD","HPRC"), cohort, GROUP),
    GROUP = factor(GROUP, levels = c(pop5, "HPRC", "KPPD"))
  )

ggplot(df_pca_plot, aes(PC1, PC2)) +
  
  ## 5개 superpop: 배경(옅게)
  geom_point(
    data = ~ filter(.x, GROUP %in% pop5),
    aes(color = GROUP),
    size = 0.7,
    alpha = 0.25
  ) +
  
  ## HPRC: 삼각형(세모)
  geom_point(
    data = ~ filter(.x, GROUP == "HPRC"),
    aes(color = GROUP, shape = GROUP),
    size = 2.0,
    alpha = 0.95
  ) +
  
  ## KPPD: 네모
  geom_point(
    data = ~ filter(.x, GROUP == "KPPD"),
    aes(color = GROUP, shape = GROUP),
    size = 2.2,
    alpha = 1
  ) +
  
  scale_color_manual(
    values = c(
      "AFR"  = "grey70",
      "AMR"  = "grey70",
      "EAS"  = "grey70",
      "EUR"  = "grey70",
      "SAS"  = "grey70",
      "HPRC" = "#C00000",
      "KPPD" = "#1F4E9E"
    )
  ) +
  scale_shape_manual(
    values = c(
      "HPRC" = 17,  # triangle
      "KPPD" = 15   # square
    )
  ) +
  theme_step1() +
  theme(
    legend.title = element_blank()
  )
  


pop5 <- c("AFR","AMR","EAS","EUR","SAS")

ggplot(df_pca_plot, aes(PC1, PC2)) +
  
  ## 1) AFR/AMR/EAS/EUR/SAS: 배경 (연하게)
  geom_point(
    data = ~ filter(.x, GROUP %in% pop5),
    aes(color = GROUP),
    size = 0.7,
    alpha = 0.35
  ) +
  
  ## 2) HPRC: 삼각형
  geom_point(
    data = ~ filter(.x, GROUP == "HPRC"),
    aes(color = GROUP, shape = GROUP),
    size = 2.0,
    alpha = 0.95,
    show.legend = c(color = FALSE, shape = TRUE)
  ) +
  
  ## 3) KPPD: 네모
  geom_point(
    data = ~ filter(.x, GROUP == "KPPD"),
    aes(color = GROUP, shape = GROUP),
    size = 2.2,
    alpha = 1,
    show.legend = c(color = FALSE, shape = TRUE)
  ) +
  
  ## color scale
  scale_color_manual(
    values = c(
      "AFR"  = "#E69F00",
      "AMR"  = "#009E73",
      "EAS"  = "#56B4E9",
      "EUR"  = "#CC79A7",
      "SAS"  = "#F0E442",
      "HPRC" = "#C00000",
      "KPPD" = "#1F4E9E"
    ),
    breaks = pop5   # legend에는 pop5만
  ) +
  
  
  
  ## shape scale (HPRC/KPPD만)
  scale_shape_manual(
    values = c(
      "HPRC" = 17,
      "KPPD" = 15
    ),
    breaks = c("HPRC","KPPD")
  ) +
  
  theme_step1() +
  theme(
    legend.title = element_blank()
  )

ggplot(df_pca_plot, aes(PC1, PC2)) +
  
  ## 배경: 5 superpop (연하게 유지)
  geom_point(
    data = ~ filter(.x, GROUP %in% pop5),
    aes(color = GROUP),
    size = 0.7,
    alpha = 0.35
  ) +
  
  ## HPRC
  geom_point(
    data = ~ filter(.x, GROUP == "HPRC"),
    aes(color = GROUP, shape = GROUP),
    size = 2.0,
    alpha = 0.95,
    show.legend = c(color = FALSE, shape = TRUE)
  ) +
  
  ## KPPD
  geom_point(
    data = ~ filter(.x, GROUP == "KPPD"),
    aes(color = GROUP, shape = GROUP),
    size = 2.2,
    alpha = 1,
    show.legend = c(color = FALSE, shape = TRUE)
  ) +
  
  scale_color_manual(
    values = c(
      "AFR"  = "#E69F00",
      "AMR"  = "#009E73",
      "EAS"  = "#56B4E9",
      "EUR"  = "#CC79A7",
      "SAS"  = "#F0E442",
      "HPRC" = "#C00000",
      "KPPD" = "#1F4E9E"
    ),
    breaks = pop5
  ) +
  
  scale_shape_manual(
    values = c(
      "HPRC" = 17,
      "KPPD" = 15
    ),
    breaks = c("HPRC","KPPD")
  ) +
  
  ## 🔴 여기서 legend만 선명하게 설정
  guides(
    color = guide_legend(
      override.aes = list(alpha = 1, size = 3)
    ),
    shape = guide_legend(
      override.aes = list(alpha = 1, size = 3)
    )
  ) +
  
  theme_step1() +
  theme(
    legend.title = element_blank()
  )






##### last all


df<- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/PCA/PCA_ethnic.kppdip_hprcdip.txt")

head(sample_1gkp)
head(df)
df %>%filter(IID == "CHM13v2")  

library(dplyr)
library(ggplot2)

pop5 <- c("AFR","AMR","EAS","EUR","SAS")

df_pca_plot <- df %>%
  mutate(
    cohort = case_when(
      str_detect(IID, "KPPD") ~ "KPPD",
      str_detect(IID, "HPRC") ~ "HPRC",
      TRUE                    ~ "1KGP"
    )
  ) %>%
  left_join(sample_1gkp %>% rename(IID = ID), by = "IID") %>%
  mutate(
    GROUP = ifelse(cohort %in% c("KPPD","HPRC"), cohort, GROUP),
    GROUP = factor(GROUP, levels = c(pop5, "HPRC", "KPPD"))
  )


pop5 <- c("AFR","AMR","EAS","EUR","SAS")


ggplot(df_pca_plot, aes(PC1, PC2)) +
  
  ## 배경: 5 superpop (연하게 유지)
  geom_point(data = ~ filter(.x, GROUP %in% pop5),aes(color = GROUP),size = 0.7,alpha = 0.35) +
  
  ## HPRC
  geom_point(data = ~ filter(.x, GROUP == "HPRC"),aes(color = GROUP, shape = GROUP),size = 2.0,alpha = 0.95,show.legend = c(color = FALSE, shape = TRUE)) +
  
  ## KPPD
  geom_point(
    data = ~ filter(.x, GROUP == "KPPD"),
    aes(color = GROUP, shape = GROUP),
    size = 2.2,
    alpha = 1,
    show.legend = c(color = FALSE, shape = TRUE)
  ) +
  
  scale_color_manual(
    values = c(
      "AFR"  = "#E69F00",
      "AMR"  = "#009E73",
      "EAS"  = "#56B4E9",
      "EUR"  = "#CC79A7",
      "SAS"  = "#F0E442",
      "HPRC" = "#C00000",
      "KPPD" = "#1F4E9E"
    ),
    breaks = pop5
  ) +
  
  scale_shape_manual(
    values = c(
      "HPRC" = 17,
      "KPPD" = 15
    ),
    breaks = c("HPRC","KPPD")
  ) +
  
  ## 🔴 여기서 legend만 선명하게 설정
  guides(
    color = guide_legend(
      override.aes = list(alpha = 1, size = 3)
    ),
    shape = guide_legend(
      override.aes = list(alpha = 1, size = 3)
    )
  ) +
  
  theme_step1() +
  theme(legend.title = element_blank())



ggplot(df_pca_plot, aes(PC1, PC2)) +
  
  ## 배경: 5 superpop
  geom_point(
    data = ~ filter(.x, GROUP %in% pop5),
    aes(color = GROUP),
    size = 0.7,
    alpha = 0.35
  ) +
  
  ## HPRC: 삼각형 (fill 색 + darkgrey 테두리)
  geom_point(
    data = ~ filter(.x, GROUP == "HPRC"),
    aes(fill = GROUP, shape = GROUP),
    size = 2.1,
    alpha = 1,
    colour = "lightgrey",
    stroke = 0.3,
    show.legend = c(fill = FALSE, shape = TRUE)
  ) +
  
  ## KPPD: 네모 (fill 색 + darkgrey 테두리)
  geom_point(
    data = ~ filter(.x, GROUP == "KPPD"),
    aes(fill = GROUP, shape = GROUP),
    size = 2.4,
    alpha = 1,
    colour = "lightgrey",
    stroke = 0.4,
    show.legend = c(fill = FALSE, shape = TRUE)
  ) +
  
  ## 색상 정의 (background + fill 공용)
  scale_color_manual(
    values = c(
      "AFR"  = "#E69F00",
      "AMR"  = "#009E73",
      "EAS"  = "#56B4E9",
      "EUR"  = "#CC79A7",
      "SAS"  = "#F0E442"
    ),
    breaks = pop5
  ) +
  
  scale_fill_manual(
    values = c(
      "HPRC" = "#C00000",
      "KPPD" = "#1F4E9E"
    )
  ) +
  
  ## shape: 테두리 가능한 도형으로 변경
  scale_shape_manual(
    values = c(
      "HPRC" = 24,  # filled triangle
      "KPPD" = 22   # filled square
    ),
    breaks = c("HPRC","KPPD")
  ) +
  
  ## legend는 선명하게
  guides(
    color = guide_legend(
      override.aes = list(alpha = 1, size = 3)
    ),
    shape = guide_legend(
      override.aes = list(alpha = 1, size = 3)
    )
  ) +
  
  theme_step1() +
  theme(
    legend.title = element_blank()
  )
