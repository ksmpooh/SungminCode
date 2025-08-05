library(tidyverse)
setwd("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/KBA_QC/")
pca <- read.table("PCA_ethnic.txt",header = T)

samplegnomad<- read.table("~/Desktop/KCDC/transplantation/ethical/1000GP_Phase3.sample",header = T)
head(samplegnomad)
case <- pca %>% select("FID") %>% filter(grepl("NIH",FID))
#case <- pca %>% select("FID")

#control <- tera %>% select("FID")
case$FID <- as.factor(case$FID)
#control$FID <- as.factor(control$FID)

gnomad <- subset(samplegnomad,select = c("ID","GROUP"))
colnames(gnomad) <- c("FID","GROUP")
case$GROUP <- "KOTRY"

table(df$GROUP)
df <- rbind(gnomad,case)

head(df)
df <- merge(pca,df,by = "FID")

plot(df$PC1,df$PC2,col = rgb(1,1,1,0.1),xlab = "PC1",ylab = "PC2",main="KOTRY Ethnic PCA",
     #     xlim = c(-0.2,0.2),ylim = c(-0.2,0.2),
     #cex.main = 1,
     cex = 1,pch = 16
)
points(df[df$GROUP == "AFR",]$PC1,df[df$GROUP == "AFR",]$PC2,col = rgb(0,0,1,0.5), cex = 0.7 , pch = 16)
points(df[df$GROUP == "AMR",]$PC1,df[df$GROUP == "AMR",]$PC2,col = rgb(1,1,0,0.5), cex = 0.7 , pch = 16)
points(df[df$GROUP == "EUR",]$PC1,df[df$GROUP == "EUR",]$PC2,col = rgb(1,0,1,0.5), cex = 0.7 , pch = 16)

points(df[df$GROUP == "SAS",]$PC1,df[df$GROUP == "SAS",]$PC2,col = rgb(0,1,0,0.5), cex = 0.7 , pch = 16)
points(df[df$GROUP == "EAS",]$PC1,df[df$GROUP == "EAS",]$PC2,col = rgb(0,1,1,0.5), cex = 0.7 , pch = 16)
#points(df[df$GROUP == "CONTROL",]$PC1,df[df$GROUP == "CONTROL",]$PC2,col = rgb(0,0,0,1), cex = 0.7 , pch = 16)
points(df[df$GROUP == "KOTRY",]$PC1,df[df$GROUP == "KOTRY",]$PC2,col = rgb(1,0,0,0.5), cex = 0.7 , pch = 16)
#abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)

color <- c(
  #  rgb(0,0,0,1),
  rgb(1,0,0,1),
  rgb(0,1,0,1),
  rgb(0,0,1,1),
  rgb(1,1,0,1),
  rgb(0,1,1,1),
  rgb(1,0,1,1))

list <- c("KOTRY","SAS","AFR","AMR","EAS","EUR")

legend("topright",list,col = color,cex = 0.8,pch = 16,pt.cex = 0.5)
#legend("bottomright",list,col = color,cex = 1,pch = 16)
#legend(y = -0.05,x=-0.15,list,col = color,cex = 1,pch = 16)

legend("bottom",c("# of Markers : 5,226"),box.lwd = 0,box.col = "white",bg = "white")
dev.off()

#####
library(tidyverse)
library(patchwork)

# 작업 디렉토리 및 데이터 불러오기
setwd("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/KBA_QC/")
pca <- read.table("PCA_ethnic.txt", header = TRUE)
samplegnomad <- read.table("~/Desktop/KCDC/transplantation/ethical/1000GP_Phase3.sample", header = TRUE)

# case: NIH 샘플
case <- pca %>%
  filter(str_detect(FID, "NIH")) %>%
  select(FID) %>%
  mutate(GROUP = "KOTRY")

# control (gnomAD) 샘플
gnomad <- samplegnomad %>%
  select(FID = ID, GROUP)

# 병합
df <- bind_rows(gnomad, case) %>%
  left_join(pca, by = "FID")

# GROUP 순서 지정
df$GROUP <- factor(df$GROUP, levels = c("KOTRY", "SAS", "AFR", "AMR", "EAS", "EUR"))

# 색상 정의
group_colors <- c(
  "KOTRY" = rgb(1, 0, 0, 0.5),
  "SAS"   = rgb(0, 1, 0, 0.5),
  "AFR"   = rgb(0, 0, 1, 0.5),
  "AMR"   = rgb(1, 1, 0, 0.5),
  "EAS"   = rgb(0, 1, 1, 0.5),
  "EUR"   = rgb(1, 0, 1, 0.5)
)

# PC1 vs PC2
p1 <- ggplot(df, aes(x = PC1, y = PC2, color = GROUP)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = group_colors) +
  labs(title = "PC1 vs PC2", x = "PC1", y = "PC2", color = "Group") +
  theme_bw()

# PC2 vs PC3
p2 <- ggplot(df, aes(x = PC2, y = PC3, color = GROUP)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = group_colors) +
  labs(title = "PC2 vs PC3", x = "PC2", y = "PC3", color = "Group") +
  theme_bw()

# 두 플롯 나란히 출력
p1 + p2
