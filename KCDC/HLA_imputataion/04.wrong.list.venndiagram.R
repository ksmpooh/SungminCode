### HLA imputation result

setwd("c:/Users/user/Desktop/KCDC/HLAimputation/255sample/")
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/Final/")

pan <- read.csv("01.pan/HLA.imputation.4digit.result.compare.ngs.and.sm.csv")
pan <- pan[1:255,]
han <- read.csv("02.han/HLA.imputation.4digit.result.compare.ngs.and.sm.csv")
han <- han[1:255,]

head(pan)
head(han)

pan.wrong <- pan[pan$A.wrong > 0 | pan$B.wrong > 0 | pan$DRB1.wrong > 0  ,]
han.wrong <- han[han$A.wrong > 0 | han$B.wrong > 0 | han$DRB1.wrong > 0  ,]
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
    pan[pan$DRB1.wrong> 0 ,]$ID %>% unlist()
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
    han[han$DRB1.wrong> 0 ,]$ID %>% unlist()
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
  height = 500 , 
  width = 500 , 
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
  height = 500 , 
  width = 500 , 
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
    han[han$DRB1.wrong> 0 ,]$ID %>% unlist(),
    pan[pan$DRB1.wrong> 0 ,]$ID %>% unlist()
  ),
  category.names = c("HAN.DRB1","PAN + KOR.DRB1"),
  filename = 'han.pan.DRB1,wrong.list.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 500 , 
  width = 500 , 
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
