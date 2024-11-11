##### /home/epigenome/ë°???????ë©´¯smkim


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


#write.table(out[,c(2,1,3,4)],"Nanopore_sample_md5sum.txt",col.names = T,row.names = F,quote = F,sep = "\t")







###


#grep HN00202145 Revio_sample_md5sum.txt | grep HDD1 | awk '{print $1,$2}' | sed 's/HDD1/./g' > HN00202145_HDD1.md5sum

#cp HN00202145_HDD1.md5sum /media/epigenome/One\ Touch/HN00202145_HDD1/


## data check 

nano <- readxl::read_xlsx("~/Desktop/KCDC/pangenome/00.datacheck/2023_pro_long.xlsx",sheet = 1)
revio48 <- readxl::read_xlsx("~/Desktop/KCDC/pangenome/00.datacheck/2023_pro_long.xlsx",sheet = 3)
revio18 <- readxl::read_xlsx("~/Desktop/KCDC/pangenome/00.datacheck/2023_pro_long.xlsx",sheet = 2,skip = 1)

head(nano)
head(revio48)
head(revio18)
nano %>% select(1,2,3,4,5) -> nano
revio48 %>% select(1,2,3,7,4) -> revio48
#revio18 %>% select(1,2,3,7,4) -> revio18
colnames(revio48) <- colnames(nano)
colnames(revio18) <- colnames(nano)

nano$platform <- "Nanopore"
revio48$platform <- "Revio_b1"
revio18$platform <- "Revio_b2"

nano %>% rbind(revio48) %>% rbind(revio18)-> df
#df$`Read N50 (bp)` <- as.numeric(df$`Read N50 (bp)`)

#df$`Average Read Length (bp)` <- as.numeric(df$`Average Read Length (bp)`)
head(df)
ggplot(df,aes(y=`Total Bases (bp)`,x=platform,fill=platform)) +
  geom_boxplot() + 
  theme(legend.position = 'none') -> a
a
ggplot(df,aes(y=`Total Reads`,x=platform,fill=platform)) +
  geom_boxplot() + 
  theme(legend.position = 'none') -> b
b
ggplot(df,aes(y=`Read N50 (bp)`,x=platform,fill=platform)) +
  geom_boxplot() + 
  theme(legend.position = 'none') -> c
c
ggplot(df,aes(y=`Average Read Length (bp)`,x=platform,fill=platform)) +
  geom_boxplot() + 
  theme(legend.position = 'none')-> d
d
ggarrange(b,d,c,a, nrow = 2,ncol = 2)
c
head(df) 


head(df)
ref <- read.table("~/Desktop/KCDC/pangenome/00.datacheck/KBA.Long_Revio_Nanopore_short.IDmatchinagtable.txt",header = T) %>% na.omit()
head(ref)
ref %>% select(Nanopore,Revio) %>% pivot_longer(1:2) %>% rename(ID = value) -> ref
df %>% mutate(ID = str_split_fixed(`Sample ID`,"_",2)[,1]) %>% filter(ID %in% ref$ID) %>% group_by(platform) %>%
  summarise(mean(`Average Read Length (bp)`),mean(`Read N50 (bp)`)) -> a

