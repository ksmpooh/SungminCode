#!/home/junkim/00_education/240731_PacBioWorkshop/00_mamba/mambaforge/envs/Rplot/bin/Rscript
# conda activate Rplot
library(ggplot2)
library(plyr)
library(dplyr)

df <- read.csv("length.csv", header = TRUE)
df$Type <- factor(df$Type, levels=c("Assembly"))

colors <- c("Assembly" = "black")

coveragePlot <- ggplot(data=df, aes(x=Coverage, y=Length/1000000, color=Type, group=Sample)) +
    geom_vline(xintercept = 0.5, linetype="dotted", linewidth=0.5) +
    geom_step(alpha=0.5) +
    scale_color_manual(values = colors) +
    scale_x_continuous(limits = c(0, 1.25), breaks = seq(0, 1.25, by = 0.25)) +
    labs(x = "Cumulative coverage", y = "Length (Mb)") +
    theme_bw()

ggsave("coverage.png", plot=coveragePlot)
pdf("coverage.pdf",width=8,height=4,paper='special')
    print(coveragePlot)
dev.off()
