##### /home/epigenome/바탕화면/smkim


## R
library(tidyverse)
library(readxl)

out <- NULL
df <-read_excel("revio&ont_md5sum.xlsx",sheet = 1,skip = 3)
df$orderNumber <- "HN00209115"
df$platform <- "Revio"
head(df)
out <- rbind(out,df)
head(out)

df <-read_excel("revio&ont_md5sum.xlsx",sheet = 2,skip = 3)
df$orderNumber <- "HN00205800"
df$platform <- "Revio"

head(df)
out <- rbind(out,df)
head(out)


df <-read_excel("revio&ont_md5sum.xlsx",sheet = 3,skip = 3)
df$orderNumber <- "HN00202145"
df$platform <- "Revio"

head(df)
out <- rbind(out,df)
head(out)

write.table(out[,c(2,1,3,4)],"Revio_sample_md5sum.txt",col.names = T,row.names = F,quote = F,sep = "\t")

library(tidyverse)
library(readxl)

out <- NULL
df <-read_excel("revio&ont_md5sum.xlsx",sheet = 4,skip = 3)
df$orderNumber <- "HN00205849"
df$platform <- "Nanopore"
head(df)
out <- rbind(out,df)
head(out)

df <-read_excel("revio&ont_md5sum.xlsx",sheet = 5,skip = 3)
df$orderNumber <- "HN00207351"
df$platform <- "Nanopore"

head(df)
out <- rbind(out,df)
head(out)


write.table(out[,c(2,1,3,4)],"Nanopore_sample_md5sum.txt",col.names = T,row.names = F,quote = F,sep = "\t")







###


grep HN00202145 Revio_sample_md5sum.txt | grep HDD1 | awk '{print $1,$2}' | sed 's/HDD1/./g' > HN00202145_HDD1.md5sum

cp HN00202145_HDD1.md5sum /media/epigenome/One\ Touch/HN00202145_HDD1/

