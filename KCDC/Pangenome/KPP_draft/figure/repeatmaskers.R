### repeatmasker pattern
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/repeat_masker")
df <- read_table("merge.repeat_size.txt")
hprc <- read_table("merge.repeat_size.hprc.txt")
head(df)
head(hprc)

df$cohort <- "KPPD"
hprc$cohort <- "HPRC"

df %>% 
  ggplot(aes(x=Category,y=bp,color=Category)) + 
  geom_violin() + 
  facet_grid(~Category,scales = "free")


df %>% 
  ggplot(aes(x=Category,y=bp,color=Category)) + 
  geom_boxplot() +
  theme_step1() + 
  theme(legend.position = 'none'
  ) +
  facet_wrap(~Category,scales = "free",nrow = 1)

table(df$Category)
df %>% rbind(hprc) %>% 
  mutate(
    Category = factor(
      Category,
      levels = c("LINE","SINE","LTR","Satellite","DNA_transposon","Simple_repeat","Low_complexity","Other"))) %>%
  ggplot(aes(x = cohort, y = bp, color = cohort,fill=cohort)) +
  ## boxplot
  geom_violin(width = 0.9,alpha = 0.6 ) +
  #geom_boxplot(width = 0.15,outlier.shape = NA,linewidth = 0.7) +
  
  ## quasirandom points
  geom_quasirandom(
    width = 0.25,
    alpha = 0.6,
    size = 1
  ) +
  
  facet_wrap(
    ~ Category,
    scales = "free",
    nrow = 2
  ) +
  scale_y_continuous(
    labels = function(x) x / 1e6
  ) + 
  theme_step1() +
  theme(
    legend.position = "rigth",
    
    ## facet title 제거
    #strip.text = element_blank()
    
    ## x축 텍스트 정리 (facet이므로 크게 의미는 없음)
    axis.title.x = element_blank()
  ) + labs(y = "Length (Mb)")

"HPRC" = "red", "KPPD" = "blue"


df %>% rbind(hprc) %>%
  ggplot(aes(x=Category,y=bp,color=Category,fill=cohort)) + 
  geom_boxplot() +
  theme_step1() + 
  theme(legend.position = 'none'
  ) +
  facet_wrap(~Category,scales = "free",nrow = 1)


#####
df %>%
  rbind(hprc) %>%
  mutate(
    cohort = factor(cohort, levels = c("KPPD", "HPRC")),
    Category = factor(
      Category,
      levels = c(
        "LINE","SINE","LTR","Satellite",
        "DNA_transposon","Simple_repeat",
        "Low_complexity","Other"
      )
    )
  ) %>%
  ggplot(aes(x = cohort, y = bp, color = cohort, fill = cohort)) +
  
  geom_violin(
    width = 0.9,
    alpha = 0.6,
    linewidth = 0
  ) +
  
  geom_quasirandom(
    width = 0.25,
    alpha = 0.6,
    size = 1
  ) +
  
  facet_wrap(
    ~ Category,
    scales = "free",
    nrow = 2
  ) +
  
  scale_y_continuous(
    labels = function(x) x / 1e6
  ) +
  
  scale_color_manual(
    values = c("HPRC" = "red", "KPPD" = "blue")
  ) +
  scale_fill_manual(
    values = c("HPRC" = "red", "KPPD" = "blue")
  ) +
  
  labs(y = "Length (Mb)") +
  
  theme_step1() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )
