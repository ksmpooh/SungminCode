### HLA imputation result

#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/255sample/")
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/Final/")
setwd("c:/Users/user/Desktop/KCDC/HLAimputation/20201026/")

pan <- read.csv("impute4/pan/compare.IMPvsNGS.A.B.DRB1.4digit.csv")
han <- read.csv("impute4/han/compare.IMPvsNGS.A.B.DRB1.4digit.csv")
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
    pan[pan$A.wrong > 0 ,]$IID %>% unlist(),
    pan[pan$B.wrong > 0 ,]$IID %>% unlist(),
    pan[pan$DRB1.wrong> 0 ,]$IID %>% unlist()
  ),
  category.names = c("A","B","DRB1"),
  filename = 'pan.wrong.list.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 400 , 
  width = 400 , 
  resolution = 200,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .7,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.7,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(-27, 20, 135),
  cat.pos = c(0, 10, 0),
  #cat.dist = c(0.055, 0.055, 0.055),
  cat.dist = c(-0.005, 0.005, -0.015),
  cat.fontfamily = "sans",
  rotation = 1
  
)

venn.diagram(
  x = list(
    han[han$A.wrong > 0 ,]$IID %>% unlist(),
    han[han$B.wrong > 0 ,]$IID %>% unlist(),
    han[han$DRB1.wrong> 0 ,]$IID %>% unlist()
  ),
  category.names = c("A","B","DRB1"),
  filename = 'han.wrong.list.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 400 , 
  width = 400 , 
  resolution = 200,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .7,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.7,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0, 0),
  cat.dist = c(0.01, 0.01, 0.01),
  cat.fontfamily = "sans",
  rotation = 1
)


myCol
venn.diagram(
  x = list(
    han[han$A.wrong > 0 ,]$IID %>% unlist(),
    pan[pan$A.wrong > 0 ,]$IID %>% unlist()
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
  cat.cex = 0.6,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = 1,
  #cat.dist = 2,
  cat.fontfamily = "sans",
  #rotation = 1
  
)

venn.diagram(
  x = list(
    han[han$B.wrong > 0 ,]$IID %>% unlist(),
    pan[pan$B.wrong > 0 ,]$IID %>% unlist()
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
  cat.cex = 0.6,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = 5,
  cat.dist = -0.0001,
  cat.fontfamily = "sans",
  #rotation = 1
)

venn.diagram(
  x = list(
    han[han$DRB1.wrong> 0 ,]$IID %>% unlist(),
    pan[pan$DRB1.wrong> 0 ,]$IID %>% unlist()
  ),
  category.names = c("HAN.DRB1","PAN + KOR.DRB1"),
  filename = 'han.pan.DRB1.wrong.list.png',
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
  cat.cex = 0.6,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = -5,
  cat.dist = c(-0.3,-0.4),
  cat.fontfamily = "sans"
  #rotation = 
)



####################################################
### 2 type 2 tool 2 ref


setwd("c:/Users/user/Desktop/KCDC/HLAimputation/20201026/")

impute4.pan <- read.csv("impute4/pan/compare.IMPvsNGS.A.B.DRB1.4digit.csv")
impute4.han <- read.csv("impute4/han/compare.IMPvsNGS.A.B.DRB1.4digit.csv")
cookHLA.pan <- read.csv("cookHLA/pan/compare.IMPvsNGS.A.B.DRB1.4digit.csv")
cookHLA.han <- read.csv("cookHLA/han/compare.IMPvsNGS.A.B.DRB1.4digit.csv")


cookHLA.pan.wrong <- cookHLA.pan[cookHLA.pan$A.wrong > 0 | cookHLA.pan$B.wrong > 0 | cookHLA.pan$DRB1.wrong > 0  ,]
cookHLA.han.wrong <- cookHLA.han[cookHLA.han$A.wrong > 0 | cookHLA.han$B.wrong > 0 | cookHLA.han$DRB1.wrong > 0  ,]
impute4.pan.wrong <- impute4.pan[impute4.pan$A.wrong > 0 | impute4.pan$B.wrong > 0 | impute4.pan$DRB1.wrong > 0  ,]
impute4.han.wrong <- impute4.han[impute4.han$A.wrong > 0 | impute4.han$B.wrong > 0 | impute4.han$DRB1.wrong > 0  ,]

head(cookHLA.han)
head(cookHLA.pan)


library(VennDiagram)
library(dplyr)
library(magrittr) # %>%
library(RColorBrewer)
library(coefficientalpha)

myCol <- brewer.pal(4, "Pastel2")


venn.diagram(
  x = list(
    impute4.pan[impute4.pan$A.wrong > 0 ,]$IID %>% unlist(),
    impute4.han[impute4.han$A.wrong > 0 ,]$IID %>% unlist(),
    cookHLA.pan[cookHLA.pan$A.wrong > 0 ,]$IID %>% unlist(),
    cookHLA.han[cookHLA.han$A.wrong > 0 ,]$IID %>% unlist()
  ),
  #category.names = c("A","B","DRB1"),
  category.names = c("impute4.pan.A","impute4.han.A","cookHLA.pan.A","cookHLA.han.A"),
  filename = 'tool.A.wrong.list.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 200,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .7,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.7,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0, 0, 0),
  cat.dist = c(-0.1, -0.1, 0.06, 0.06),
  cat.fontfamily = "sans",
  #rotation = 1
  
)




myCol <- brewer.pal(4, "Pastel2")
venn.diagram(
  x = list(
    impute4.pan[impute4.pan$B.wrong > 0 ,]$IID %>% unlist(),
    impute4.han[impute4.han$B.wrong > 0 ,]$IID %>% unlist(),
    cookHLA.pan[cookHLA.pan$B.wrong > 0 ,]$IID %>% unlist(),
    cookHLA.han[cookHLA.han$B.wrong > 0 ,]$IID %>% unlist()
  ),
  #category.names = c("A","B","DRB1"),
  category.names = c("impute4.pan.B","impute4.han.B","cookHLA.pan.B","cookHLA.han.B"),
  filename = 'tool.B.wrong.list.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 200,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .7,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.7,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0, 0, 0),
  cat.dist = c(-0.1, -0.1, 0.06, 0.06),
  cat.fontfamily = "sans",
  #rotation = 1
  
)



myCol <- brewer.pal(4, "Pastel2")
venn.diagram(
  x = list(
    impute4.pan[impute4.pan$DRB1.wrong > 0 ,]$IID %>% unlist(),
    impute4.han[impute4.han$DRB1.wrong > 0 ,]$IID %>% unlist(),
    cookHLA.pan[cookHLA.pan$DRB1.wrong > 0 ,]$IID %>% unlist(),
    cookHLA.han[cookHLA.han$DRB1.wrong > 0 ,]$IID %>% unlist()
  ),
  #category.names = c("A","B","DRB1"),
  category.names = c("impute4.pan.DRB1","impute4.han.DRB1","cookHLA.pan.DRB1","cookHLA.han.DRB1"),
  filename = 'tool.DRB1.wrong.list.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 200,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .7,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0, 0, 0),
  cat.dist = c(-0.1, -0.1, 0.06, 0.06),
  cat.fontfamily = "sans",
  #rotation = 1
  
)




#######################chcek list all 20201214
a <- intersect(impute4.pan[impute4.pan$A.wrong > 0 ,]$IID,
          impute4.han[impute4.han$A.wrong > 0 ,]$IID)
b <- intersect(cookHLA.pan[cookHLA.pan$A.wrong > 0 ,]$IID,
          cookHLA.han[cookHLA.han$A.wrong > 0 ,]$IID)
c1 <- intersect(a,b)

a <- intersect(impute4.pan[impute4.pan$B.wrong > 0 ,]$IID,
               impute4.han[impute4.han$B.wrong > 0 ,]$IID)
b <- intersect(cookHLA.pan[cookHLA.pan$B.wrong > 0 ,]$IID,
               cookHLA.han[cookHLA.han$B.wrong > 0 ,]$IID)
c2 <- intersect(a,b)

a <- intersect(impute4.pan[impute4.pan$DRB1.wrong > 0 ,]$IID,
               impute4.han[impute4.han$DRB1.wrong > 0 ,]$IID)
b <- intersect(cookHLA.pan[cookHLA.pan$DRB1.wrong > 0 ,]$IID,
               cookHLA.han[cookHLA.han$DRB1.wrong > 0 ,]$IID)
c3 <- intersect(a,b)

intersect(c1,c2)
intersect(c2,c3)
intersect(c3,c1)
