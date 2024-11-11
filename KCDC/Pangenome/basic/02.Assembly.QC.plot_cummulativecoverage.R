setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/construction/03.hifiasm/")

library(tidyverse)

library(ggplot2)
library(plyr)
library(dplyr)
head(df)
df <- read.csv("length.csv", header = TRUE)
hprc <- read.csv("merge.HPRC.phaseI.NGx.length.csv", header = F)
grch <- read.csv("length_grch38.csv", header = F)
grch <- read.csv("lnegth_grch38_plamary_scaf.csv", header = F)
head(hprc)
head(grch)
colnames(hprc) <- colnames(df)
colnames(grch) <- colnames(df)
grch %>% mutate(Tyep = "GRCh38")

hprc %>% mutate(Type = case_when(
  Sample == "CHM13Y.fa" ~ "CHM13",
  str_detect(Sample,"GR") ~ "GRCh38",
  TRUE ~ "HPRC")) %>% filter(Type != "GRCh38") %>%
  rbind(grch %>% mutate(Type = "GRCh38")) %>% #head()
  rbind(df %>% mutate(Type = "KoGES")) %>% #head
  mutate(Type = as.factor(Type)) %>%# count(Type)
  mutate(Coverage = Coverage*3) %>%
  ggplot(aes(x=Coverage, y=Length/1000000, color=Type, group=Sample)) +
  geom_vline(xintercept = 1.5, linetype="dotted", linewidth=0.5) +
  geom_step(alpha=0.8) +
  #scale_color_manual(values = colors) +
  #scale_x_continuous(limits = c(0, 1.25), breaks = seq(0, 1.25, by = 0.25)) +
  labs(x = "Cumulative coverage (Gb)", y = "Length (Mb)") +
  theme_bw()



hprc %>%
  mutate(Type = case_when(
    Sample == "CHM13Y.fa" ~ "CHM13",
    str_detect(Sample, "GR") ~ "GRCh38",
    TRUE ~ "HPRC"
  )) %>%
  filter(Type != "GRCh38") %>%
  rbind(grch %>% mutate(Type = "GRCh38")) %>%
  rbind(df %>% mutate(Type = "KoGES")) %>%
  mutate(Type = as.factor(Type)) %>%
  mutate(Coverage = Coverage * 3) %>%
  ggplot(aes(x = Coverage, y = Length / 1000000, color = Type, group = Sample)) +
  geom_vline(xintercept = 1.5, linetype = "dotted", linewidth = 0.5) +
  geom_step(size = 0.5,alpha=0.5) +
  #scale_color_manual(values = c("GRCh38" = "black", "CHM13" = "black", "HPRC" = "red", "KoGES" = "blue")) +
  scale_alpha_manual(values = c("GRCh38" = 1, "CHM13" = 1, "HPRC" = 0.3, "KoGES" = 0.3)) +
  labs(x = "Cumulative coverage (Gb)", y = "Length (Mb)") +
  theme_bw() +
  guides(alpha = "none")

library(dplyr)
library(ggplot2)

# Data preparation
processed_data <- hprc %>%
  mutate(Type = case_when(
    Sample == "CHM13Y.fa" ~ "CHM13",
    str_detect(Sample, "GR") ~ "GRCh38",
    TRUE ~ "HPRC"
  )) %>%
  filter(Type != "GRCh38") %>%
  rbind(grch %>% mutate(Type = "GRCh38")) %>%
  rbind(df %>% mutate(Type = "KoGES")) %>%
  mutate(Type = as.factor(Type)) %>%
  mutate(Coverage = Coverage * 3)

# Plot
ggplot() +
  # Plot HPRC and KoGES first with lower alpha
  geom_step(
    data = processed_data %>% filter(Type %in% c("HPRC", "KoGES")),
    aes(x = Coverage, y = Length / 1000000, color = Type, group = Sample),
    size = 0.5, alpha = 0.2
  ) +
  # Plot GRCh38 and CHM13 with higher alpha and black color, no legend
  geom_step(
    data = processed_data %>% filter(Type %in% c("GRCh38", "CHM13")),
    aes(x = Coverage, y = Length / 1000000, group = Sample),
    color = "black", size = 0.5, alpha = 1, show.legend = FALSE
  ) +
  geom_vline(xintercept = 1.5, linetype = "dotted", linewidth = 0.5) +
  scale_color_manual(values = c("HPRC" = "red", "KoGES" = "blue")) +
  labs(x = "Cumulative coverage (Gb)", y = "Length (Mb)") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  # Annotations for CHM13 and GRCh38
  annotate("text", x = 0.8, y = 200, label = "CHM13", color = "black", size = 5, hjust = 0) +
  annotate("text", x = 0.8, y = 100, label = "GRCh38", color = "black", size = 5, hjust = 0)





ggsave("coverage.png", plot=coveragePlot)
pdf("coverage.pdf",width=8,height=4,paper='special')
print(coveragePlot)
dev.off()
head(df)
df %>% count(Length < LengthSum)
df %>% mutate(Coverage2 = Coverage/10*3) %>% head()
df %>% mutate(Coverage2 = LengthSum/Length*100*300000000) %>% head()

df %>% mutate(Coverage = Coverage/10*3) %>%  #head()
  ggplot(aes(x=Coverage, y=Length/1000000, color=Type, group=Sample)) +
  geom_vline(xintercept = 1.5, linetype="dotted", linewidth=0.5) +
  geom_step(alpha=0.5) +
  scale_color_manual(values = colors) +
  scale_x_continuous(limits = c(0, 3.0), breaks = seq(0, 3, by = 1)) +
  labs(x = "Cumulative coverage", y = "Length (Mb)") +
  theme_bw()
