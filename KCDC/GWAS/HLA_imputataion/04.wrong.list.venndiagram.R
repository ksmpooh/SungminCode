### HLA imputation result


setwd("c:/Users/user/Desktop/KCDC/HLAimputation/Final/")

pan <- read.csv("HLA.imputation.2digit.result.compare.ngs.and.sm_usingPan.Kor.csv")
pan <- pan[1:255,]
han <- read.csv("HLA.imputation.2digit.result.compare.ngs.and.sm_usingHan.csv")
han <- han[1:255,]

head(pan)
head(han)

pan.wrong <- pan[pan$A.wrong > 0 | pan$B.wrong > 0 | pan$DRB.wrong > 0  ,]
han.wrong <- han[han$A.wrong > 0 | han$B.wrong > 0 | han$DRB.wrong > 0  ,]
#pan.wrong[pan.wrong$A.wrong > 0 & pan.wrong$B.wrong > 0 & pan.wrong$DRB.wrong > 0  ,]
pan[pan$A.wrong > 0 | pan$B.wrong > 0 | pan$DRB.wrong > 0  ,]$ID

library(VennDiagram)
library(dplyr)
library(magrittr) # %>%
library(RColorBrewer)
library(coefficientalpha)

myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(
    pan[pan$A.wrong > 0 ,]$ID %>% unlist(),
    pan[pan$B.wrong > 0 ,]$ID %>% unlist(),
    pan[pan$DRB.wrong> 0 ,]$ID %>% unlist()
  ),
  category.names = c("HLA.A","HLA.B","HLA.DRB1"),
  filename = 'pan.wrong.list.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 700 , 
  width = 700 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 20, 135),
  cat.dist = c(0.055, 0.055, 0.055),
  cat.fontfamily = "sans",
  rotation = 1
  
)

venn.diagram(
  x = list(
    han[han$A.wrong > 0 ,]$ID %>% unlist(),
    han[han$B.wrong > 0 ,]$ID %>% unlist(),
    han[han$DRB.wrong> 0 ,]$ID %>% unlist()
  ),
  category.names = c("HLA.A","HLA.B","HLA.DRB1"),
  filename = 'han.wrong.list.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 700 , 
  width = 700 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 20, 135),
  cat.dist = c(0.055, 0.055, 0.055),
  cat.fontfamily = "sans",
  rotation = 1
)


myCol
venn.diagram(
  x = list(
    han[han$A.wrong > 0 ,]$ID %>% unlist(),
    pan[pan$A.wrong > 0 ,]$ID %>% unlist()
  ),
  category.names = c("HAN.A","PAN + KOR.A"),
  filename = 'han.pan.A.list.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 300 , 
  width = 300 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#B3E2CD","#FDCDAC"),
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  
  # Set names
  cat.cex = 0.4,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = 1,
  #cat.dist = 2,
  cat.fontfamily = "sans",
  #rotation = 1
  
)

venn.diagram(
  x = list(
    han[han$B.wrong > 0 ,]$ID %>% unlist(),
    pan[pan$B.wrong > 0 ,]$ID %>% unlist()
  ),
  category.names = c("HAN.B","PAN + KOR.B"),
  filename = 'han.pan.B.list.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 400 , 
  width = 400 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#B3E2CD","#FDCDAC"),
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
    
  # Set names
  cat.cex = 0.4,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = 5,
  cat.dist = -0.0001,
  cat.fontfamily = "sans",
  #rotation = 1
)

venn.diagram(
  x = list(
    han[han$DRB.wrong> 0 ,]$ID %>% unlist(),
    pan[pan$DRB.wrong> 0 ,]$ID %>% unlist()
  ),
  category.names = c("HAN.DRB1","PAN + KOR.DRB1"),
  filename = 'han.pan.DRB1,wrong.list.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 400 , 
  width = 400 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#B3E2CD","#FDCDAC"),
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.4,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = 0,
  #cat.dist = 2,
  cat.fontfamily = "sans"
  #rotation = 
  
)
