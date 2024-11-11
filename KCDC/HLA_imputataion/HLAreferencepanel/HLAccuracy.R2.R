#### R2
setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/aggR/RS/")

flist = grep(list.files("./"),pattern = "Multi", value=TRUE,invert = T)
flist
a <- NULL
#df <- read.table("g4.KMHC.DR2.txt")
head(a)
for (i in flist) {
  tmp <- read.table(i)
  tmp$g <- i
  a <- rbind(a,tmp)
}
flist = grep(list.files("./"),pattern = "Multi", value=TRUE)
multi <- NULL
for (i in flist) {
  tmp <- read.table(i)
  tmp$g <- i
  multi <- rbind(multi,tmp)
}
head(multi)
head(a)
colnames(a) <- c("ID","AF","DR2","g")
colnames(multi) <- c("ID","AF",'DR2',"ER2","g")

multi %>% select(-"ER2") %>% rbind(a) %>%
  mutate(fold = str_split_fixed(g,"\\.",4)[,1]) %>% 
  mutate(Ref = str_split_fixed(g,"\\.",4)[,2]) %>% 
  mutate(HLAtype = str_split_fixed(ID,"HLA_",2)[,2]) %>% 
  mutate(Gene = str_split_fixed(HLAtype,"\\*",2)[,1]) %>%
  filter(AF != 0,Ref != "1KGP") %>% select(-ID,-g) %>%
  mutate(Ref = factor(Ref,c("KMHC","Multi","Han","PanKor"),labels = c("KMHC","Multi-ehthnic","Han Chinese","Pan-Kor"))) -> df

head(df)
freq_ref <- readxl::read_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/KMHC.HLAtype.freq.xlsx") %>% select(-n,-prop)
head(freq_ref)
colnames(freq_ref)[2]<- "HLAtype"

df %>% left_join(freq_ref) -> df



df %>%
  group_by(Ref,fold,Gene,freq) %>%
  summarise(R2 = mean(DR2)) -> df1

df %>% group_by(Ref,fold,freq) %>% #head()
  summarise(R2 = mean(DR2)) %>% #dim()
  mutate(Gene = "Overall") -> df2
head(df1)
head(df2)


df1 %>% ungroup() %>% count(Ref)

df1 %>% rbind(df2) %>%
  ggplot(aes(x=freq,y=R2,fill=Ref)) +
  geom_boxplot() +
  facet_wrap(~Gene,nrow = 3)

