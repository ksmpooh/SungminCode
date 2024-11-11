#!/home/junkim/00_education/240731_PacBioWorkshop/00_mamba/mambaforge/envs/Rplot/bin/Rscript
# conda activate Rplot
library(ggplot2)
library(plyr)
library(dplyr)

df <- read.csv("ReadLengthQual.csv")

df$sample <- as.factor(df$sample)
df$sample <- factor(df$sample,level = c("GMPB010_1", "GMPB010_4", "GMPB011_1", "GMPB011_4", "GMPB034_1", "GMPB012_1", "GMPB012_4", "GMPB013_1", "GMPB013_4", "GMPB015_1", "GMPB020_1", "GMPB024_1", "GMPB024_2", "GMPB024_3", "GMPB025_1", "GMPB025_2", "GMPB025_3", "GMPB035_1", "GMPB035_4", "GMPB036_1", "GMPB036_4"))

length_violin <- ggplot(df, aes(x=sample, y=length/1000)) +
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1, outlier.size = NA, outlier.shape = NA) +
  labs(x = "Sample", y = "Read length (kb)") +
  ylim(1,40) +
  theme_bw()

qual_violin <- ggplot(df, aes(x=sample, y=mean_quality)) +
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1, outlier.size = NA, outlier.shape = NA) +
  labs(x = "Sample", y = "Mean read quality") +
  ylim(1,50) +
  theme_bw()

#plot <- plot_grid(length_violin, qual_violin, ncol = 1)

ggsave("readLength.png", plot=length_violin)
ggsave("readQual.png", plot=qual_violin)

pdf("readLength.pdf",width=20,height=2,paper='special')
print(length_violin)
dev.off()

pdf("readQual.pdf",width=20,height=2,paper='special')
print(qual_violin)
dev.off()
