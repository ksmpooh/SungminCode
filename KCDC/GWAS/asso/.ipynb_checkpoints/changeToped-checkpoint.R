setwd("C:/Users/user/Desktop/KCDC/")
df <- read.table("association/KCHIP130K_transformation_hematological_20190903.txt",
                 header = T,sep = " ")

colnames(df)
df$FAM_ID <- df$id
df$IND_ID <- df$id
df$FAT_ID <- 0
df$MOT_ID <- 0

df <- subset(df,select = c(30,31,32,33,2:29))

write.table(df,"association/KCHIP130K_transformation_hematological_20190903.ped",col.names = T,row.names = F,quote = F)
