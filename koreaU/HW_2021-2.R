### 고려대학원 2021-2 : ㄷ데이터시각화 과제
library(tidyverse)
library(readxl)
library(grid)
library(gridExtra)
setwd("~/Desktop/KU/2021_Fall/biovisual/hw/ScienceDirect_files_19Oct2021_06-39-57.618/")


df <- read_excel("1-s2.0-S0092867421008576-mmc1.xlsx",sheet = 2)
head(df)
colnames(df)
df <- as.data.frame(df)



df %>% select(Type,Ethnicity,Gender,Smoking.History_modified,Stage) %>%
  group_by(Type,Stage) %>%
  summarise()
    










df %>% filter(Type == "Tumor") %>% # stage 간 NAT는 NA
  select(Ethnicity,Gender,Smoking.History_modified,Stage) %>%
  #count(Ethnicity) 
  #count(Gender)
  #count(Smoking.History_modified)
  #count(Stage)
  mutate_at('Ethnicity',str_replace_all,c(han = 'asian',tssdidnotcollectthisinformation = 'NA')) %>%
  mutate_at('Ethnicity',str_replace,"white\\(caucasian\\)","caucasian") %>%
  #count(Ethnicity)
  mutate_at("Stage",str_replace_all,c(IA = 'I',IB = 'I',IIA = 'II',IIB = 'II',IIIA = 'III',IIIB = 'III')) %>%
  #count(Stage)
  gather() %>%
  filter(value != "NA") %>% 
  group_by(key,value) %>%
  summarise(n = n()) %>%
  mutate(Frequency = n/sum(n)) %>% #head()
  arrange(key,-Frequency) %>% 
  #arrange(c(1:4,6,5,9,7,8,11,10,12,13)) %>%
  arrange(c(1:4,5,6,7,8,9,11,10,12,13)) %>%
  as.data.frame() %>%
  mutate(Color = c("darkred","deepskyblue3","aquamarine4","black",
                   "darkred","deepskyblue3",
                   "darkred","deepskyblue3","aquamarine4",
                   "darkred","deepskyblue3","aquamarine4","darkslateblue")) %>%
  ggplot(aes(x=key,y=Frequency,fill=Color)) + 
  geom_bar(position="fill",stat= 'identity') +
  theme(legend.position ="none") 
  


df %>% filter(Type == "Tumor") %>% # stage 간 NAT는 NA
  select(Ethnicity,Gender,Smoking.History_modified,Stage) %>%
  rename("Smoking\nHistory" = Smoking.History_modified) %>%
  mutate_at('Ethnicity',str_replace_all,c(han = 'asian',tssdidnotcollectthisinformation = 'NA')) %>%
  mutate_at('Ethnicity',str_replace,"white\\(caucasian\\)","caucasian") %>%
  mutate(Ethnicity = str_to_title(Ethnicity),Gender = str_to_title(Gender),)
  mutate_at("Stage",str_replace_all,c(IA = 'I',IB = 'I',IIA = 'II',IIB = 'II',IIIA = 'III',IIIB = 'III')) %>%
  gather() %>%
  filter(value != "NA") %>% 
  group_by(key,value) %>%
  summarise(n = n()) %>%
  mutate(Frequency = n/sum(n)) %>% #head()
  arrange(key,-Frequency) %>% 
  as.data.frame() %>%
  mutate(Color = c(1,2,3,4,2,1,2,3,1,3,4,2,1)) %>%
  arrange(desc(Color)) %>%
  #mutate(Color = c("red","blue","green","black","blue","red","blue","green","red","green","black","red","blue")) %>%
  #arrange(c(1:4,6,5,7:13)) %>%
  #arrange(c(2,3,4,1,6,5,7,8,9,11,10,13,12)) %>%
  as.data.frame() %>% #head()
  ggplot(aes(x=key,y=Frequency,fill=Color)) + 
  #geom_bar(position="fill",stat= 'identity') +
  geom_bar(position="fill",stat= 'identity') +
  scale_y_continuous(labels = label_percent(accuracy =  1)) + 
  #scale_fill_manual(values = "red","green","blue","black") + 
  labs(y = "Frequency",x="") + 
  xlab(c("Ethinicity","Gender","Smoking\nhistory","Stage")) + 
  theme(legend.position ="none") 

library(scales)





a <- df %>% filter(Type == "Tumor") %>% # stage 간 NAT는 NA
  select(Ethnicity,Gender,Smoking.History_modified,Stage) %>%
  rename("Smoking\nHistory" = Smoking.History_modified) %>%
  mutate_at('Ethnicity',str_replace_all,c(han = 'asian',tssdidnotcollectthisinformation = 'NA')) %>%
  mutate_at('Ethnicity',str_replace,"white\\(caucasian\\)","caucasian") %>%
  mutate(Ethnicity = str_to_title(Ethnicity),Gender = str_to_title(Gender)) %>%
  mutate_at("Stage",str_replace_all,c(IA = 'I',IB = 'I',IIA = 'II',IIB = 'II',IIIA = 'III',IIIB = 'III')) %>%
  gather() %>%
  filter(value != "Na") %>%
  filter(value != "NA") %>%
  group_by(key,value) %>%
  summarise(n = n()) %>%
  mutate(Frequency = n/sum(n)) %>% #head()
  arrange(key,-Frequency) %>% 
  as.data.frame() %>% #head
  #mutate(Color = c(1,2,3,4,2,1,2,3,1,3,4,2,1)) %>%
#  arrange(desc(Color)) %>%
  mutate(Color = c("red","blue","green","black","blue","red","blue","green","red","green","black","red","blue")) %>%
  #arrange(c(1:4,6,5,7:13)) %>%
  #arrange(c(2,3,4,1,6,5,7,8,9,11,10,13,12)) %>%
  as.data.frame()#head()


b <- df %>% filter(Type == "Tumor") %>% # stage 간 NAT는 NA
  select(Ethnicity,Gender,Smoking.History_modified,Stage) %>%
  rename("Smoking\nHistory" = Smoking.History_modified) %>%
  mutate_at('Ethnicity',str_replace_all,c(han = 'asian',tssdidnotcollectthisinformation = 'NA')) %>%
  mutate_at('Ethnicity',str_replace,"white\\(caucasian\\)","caucasian") %>%
  mutate(Ethnicity = str_to_title(Ethnicity),Gender = str_to_title(Gender)) %>%
  mutate_at("Stage",str_replace_all,c(IA = 'I',IB = 'I',IIA = 'II',IIB = 'II',IIIA = 'III',IIIB = 'III')) %>%
  gather() %>%
  filter(value != "Na") %>%
  filter(value != "NA") %>%
  group_by(key,value) %>%
  summarise(n = n()) %>%
  mutate(Frequency = n/sum(n)) %>% #head()
  arrange(key,-Frequency) %>% 
  as.data.frame() #%>% #head
  #mutate(Color = c(1,2,3,4,2,1,2,3,1,3,4,2,1)) %>%
  #  arrange(desc(Color)) %>%
  #mutate(Color = c("red","blue","green","black","blue","red","blue","green","red","green","black","red","blue")) %>%
  #arrange(c(1:4,6,5,7:13)) %>%
  #arrange(c(2,3,4,1,6,5,7,8,9,11,10,13,12)) %>%
 # as.data.frame()#head()

a

head(a)
table(a$value)
a$Color <- 'Red'
a1 <- a %>% filter(key == 'Ethnicity')
a2 <- a %>% filter(key == 'Stage')

b1 <- b %>% filter(key == 'Ethnicity')
b2 <- b %>% filter(key == 'Stage')

#ggplot(a,mapping = aes(x=key,fill=interaction(key,value)))+

ggplot(a,mapping = aes(x=key,fill=interaction(key,value)))+
    geom_bar(position = 'fill') +
  guides(fill = guide_legend(title = 'Category and subcategory'))
#https://stackoverflow.com/questions/33376750/specifying-the-number-of-factors-in-legend-columns


head(b)
b
b <- b %>% mutate(Color= c("red","blue","green","black","blue","red","blue","green","red","green","black","red","blue")) 
ggplot(b,aes(x=key,y=Frequency,fill=value)) +
  geom_bar(stat= 'identity') + 
  theme(legend.position = "bottom",legend.justification = "left") +
  guides(fill=guide_legend(title = 'test'))
  #geom_bar(position="fill",stat= 'identity')+
  #theme(legend.position ="bottom",legend.justification = "left")
  #guides(fill = guide_legend(title = 'test'))#
  #guides(fill=guide_legend(override.aes = list(fill=c(1,2,3,4,5,6,7,8,9,10,11,12,13)))) +
  #guides(fill=guide_legend(override.aes = list(fill=c(5,3,2,1,5,6,7,8,9,10,11,12,13))))

  


ggplot(a,aes(x=key,fill=factor(value,levels = c('Caucasian','Asian','Slavic','Black','Female','Male','Non-Smoker','Reformed Smoker','Smoker',"IV","III","II","I")))) + 
  geom_bar(position = 'fill')+
  #theme(legend.position ="none") + 
  theme(legend.position ="bottom") + 
  scale_fill_manual(values = c("darkred","deepskyblue3","darkseagreen4","black","darkred","deepskyblue3","darkred","deepskyblue3","darkseagreen4","darkred","deepskyblue3","darkseagreen4","darkslateblue")) +
  labs(y = "Frequency",x="") +
  guides(fill=guide_legend(title = "test"))+
  scale_fill_manual(values = c(hcl))  


ggplot() + 
  geom_bar(a1,mapping = aes(x=key,fill = value)) +
  geom_bar(position = 'fill')
ggplot(a,aes(x=key,fill=value)) +
  geom_bar(position = 'fill') +
  theme(legend.position = "left",legend.margin = margin(1,1,1,1))
  
    guides(fill=guide_legend(override.aes = list(fill=c(1,2,3,4))))
  
ggplot() + 
  geom_bar(a1,mapping = aes(x=key,size = 'Ethnicity')) + 
  geom_bar(a2,mapping = aes(x=key,colour = 'Stage')) + 
  geom_bar(position = 'fill')

#ggplot(a,aes(x=key,fill=value)) + 
ggplot(a,aes(x=key,fill=factor(value,levels = c('Caucasian','Asian','Slavic','Black','Female','Male','Non-Smoker','Reformed Smoker','Smoker',"IV","III","II","I")))) + 
  geom_bar(position = 'fill')+
  #theme(legend.position ="none") + 
  theme(legend.position ="bottom") + 
  scale_fill_manual(values = c("darkred","deepskyblue3","darkseagreen4","black","darkred","deepskyblue3","darkred","deepskyblue3","darkseagreen4","darkred","deepskyblue3","darkseagreen4","darkslateblue")) +
  labs(y = "Frequency",x="") + 
  guides(fill=guide_legend(override.aes = list(fill = c(1,2,3,4)),title="test"))
  #guides(fill = guide_legend(byrow = T,ncol = key))

p <- ggplot(a,aes(x=key,fill=factor(value,levels = c('Caucasian','Asian','Slavic','Black','Female','Male','Non-Smoker','Reformed Smoker','Smoker',"IV","III","II","I")))) + 
  geom_bar(position = 'fill')+
  theme(legend.position ="none") + 
  #theme(legend.position ="bottom") + 
  scale_fill_manual(values = c("darkred","deepskyblue3","darkseagreen4","black","darkred","deepskyblue3","darkred","deepskyblue3","darkseagreen4","darkred","deepskyblue3","darkseagreen4","darkslateblue")) +
  labs(y = "Frequency",x="")
  #guides(fill = guide_legend(title = c("test","2"),nrow = 4))
  #scale_fill_hue("test",guide=guide_legend(order = 1))
  #scale_color_manual(guide_legend())
  #guides(fill=guide_legend(title="New1"))

p1 <-ggplot(a %>% filter(key == "Ethnicity"), aes(x = key,fill=factor(value,levels = c('Caucasian','Asian','Slavic','Black')))) + 
  geom_bar() +
  scale_fill_manual(values = c("darkred","deepskyblue3","darkseagreen4","black")) +
  theme(legend.position ="bottom",legend.justification = 'left',legend.box.margin = margin(0,0,0,0)) +
  guides(fill=guide_legend(title="Ethnicity"))


p1 <- ggplot(a %>% filter(key == "Ethnicity"), aes(x = key,fill=factor(value,levels = c('Caucasian','Asian','Slavic','Black')))) + 
  geom_bar() +
  scale_fill_manual(values = c("darkred","deepskyblue3","darkseagreen4","black")) +
  theme(legend.position ="bottom") +
  guides(fill=guide_legend(title="Ethnicity"))
p2 <- ggplot(a %>% filter(key == "Gender"), aes(x = key,fill=factor(value,levels = c("Female","Male")))) + 
  geom_bar() +
  scale_fill_manual(values = c("darkred","deepskyblue3")) +
  theme(legend.position ="bottom") +
  guides(fill=guide_legend(title="Gender"))

p2 <- ggplot(a %>% filter(key == "Gender"), aes(x = key,fill=factor(value,levels = c("Female","Male")))) + 
  geom_bar() +
  scale_fill_manual(values = c("darkred","deepskyblue3")) +
  theme(legend.position ="bottom",legend.box.margin = margin(0,0,0,0)) +
  guides(fill=guide_legend(title="Gender"))

p3 <- ggplot(a %>% filter(key == "Smoking\nHistory"), aes(x = key,fill=factor(value,levels = c('Non-Smoker','Reformed Smoker','Smoker')))) + 
  geom_bar() +
  scale_fill_manual(values = c("darkred","deepskyblue3","darkseagreen4")) +
  theme(legend.position ="bottom") +
  guides(fill=guide_legend(title="Smoking\nHistory"))


p4 <- ggplot(a %>% filter(key == "Stage"), aes(x = key,fill=factor(value,levels = c("I","II","III","IV")))) + 
  geom_bar() +
  scale_fill_manual(values = c("darkslateblue","darkseagreen4","deepskyblue3","darkred")) +
  theme(legend.position ="bottom") +
  guides(fill=guide_legend(title="Stage"))

dev.off()
cowplot::get_legend(p1)
p1 <- cowplot::get_legend(p1)
p2 <- cowplot::get_legend(p2)
p3 <- cowplot::get_legend(p3)
p4 <- cowplot::get_legend(p4)

plot()
ggplot() + geom_blank()
  
ggplot() + theme_void()
legend("")
p2 <- ggplot(a %>% filter(key == "Gender"), aes(x = key,fill=factor(value,levels = c("Female","Male")))) + 
  geom_bar() +
  scale_fill_manual(values = c("darkred","deepskyblue3")) +
  theme(legend.position ="bottom",legend.box.margin = margin(10,0,0,0)) +
  guides(fill=guide_legend(title="Gender"))
p1 <- cowplot::get_legend(p1)
p2 <- cowplot::get_legend(p2)

plot_grid(p4,p2,p3,p4,ncol = 1,align = 'hv',rel_heights = c(1,1,1,1))

plot_grid(p,p1,p2,p3,p4,ncol = 1,rel_heights = c(1,0.1,0.1,0.1,0.1))


plot.new()
legend("bottomleft",title="test",c('Caucasian','Asian','Slavic','Black'),box.lwd = 0,box.col = "white",bg = "white",horiz = T,
       col = c("darkred","deepskyblue3","darkseagreen4","black"),pch = 15)
p1 <- ggplot(a %>% filter(key == "Ethnicity"), aes(x = key,fill=factor(value,levels = c('Caucasian','Asian','Slavic','Black')))) + 
  geom_bar() +
  scale_fill_manual(values = c("darkred","deepskyblue3","darkseagreen4","black")) +
  theme(legend.position ="bottom") +
  guides(fill=guide_legend(title="Ethnicity"))

legend("topright",list,col = color,cex = 1,pch = 16)


### part 1

a <- df %>% filter(Type == "Tumor") %>% # stage 간 NAT는 NA
  select(Ethnicity,Gender,Smoking.History_modified,Stage) %>%
  rename("Smoking\nHistory" = Smoking.History_modified) %>%
  mutate_at('Ethnicity',str_replace_all,c(han = 'asian',tssdidnotcollectthisinformation = 'NA')) %>%
  mutate_at('Ethnicity',str_replace,"white\\(caucasian\\)","caucasian") %>%
  mutate(Ethnicity = str_to_title(Ethnicity),Gender = str_to_title(Gender)) %>%
  mutate_at("Stage",str_replace_all,c(IA = 'I',IB = 'I',IIA = 'II',IIB = 'II',IIIA = 'III',IIIB = 'III')) %>%
  gather() %>%
  filter(value != "Na") %>%
  filter(value != "NA") %>%
  group_by(key,value) %>%
  as.data.frame()

p <- ggplot(a,aes(x=key,fill=factor(value,levels = c('Caucasian','Asian','Slavic','Black','Female','Male','Non-Smoker','Reformed Smoker','Smoker',"IV","III","II","I")))) + 
  geom_bar(position = 'fill')+
  theme(legend.position ="none") + 
  #theme(legend.position ="bottom") + 
  scale_fill_manual(values = c("darkred","deepskyblue3","darkseagreen4","black","darkred","deepskyblue3","darkred","deepskyblue3","darkseagreen4","darkred","deepskyblue3","darkseagreen4","darkslateblue")) +
  labs(y = "Frequency",x="")


p1 <- ggplot(a %>% filter(key == "Ethnicity"), aes(x = key,fill=factor(value,levels = c('Caucasian','Asian','Slavic','Black')))) + 
  geom_bar() +
  scale_fill_manual(values = c("darkred","deepskyblue3","darkseagreen4","black")) +
  theme(legend.position ="bottom",legend.justification = "left") +
  guides(fill=guide_legend(title="Ethnicity  "))
p2 <- ggplot(a %>% filter(key == "Gender"), aes(x = key,fill=factor(value,levels = c("Female","Male")))) + 
  geom_bar() +
  scale_fill_manual(values = c("darkred","deepskyblue3")) +
  theme(legend.position ="bottom",legend.justification = "left") +
  guides(fill=guide_legend(title="Gender   "))
p3 <- ggplot(a %>% filter(key == "Smoking\nHistory"), aes(x = key,fill=factor(value,levels = c('Non-Smoker','Reformed Smoker','Smoker')))) + 
  geom_bar() +
  scale_fill_manual(values = c("darkred","deepskyblue3","darkseagreen4")) +
  theme(legend.position ="bottom",legend.justification = "left") +
  guides(fill=guide_legend(title="Smoking \nHistory"))
p4 <- ggplot(a %>% filter(key == "Stage"), aes(x = key,fill=factor(value,levels = c("I","II","III","IV")))) + 
  geom_bar() +
  scale_fill_manual(values = c("darkslateblue","darkseagreen4","deepskyblue3","darkred")) +
  theme(legend.position ="bottom",legend.justification = "left") +
  guides(fill=guide_legend(title="Stage      "))

p1 <- cowplot::get_legend(p1)
p2 <- cowplot::get_legend(p2)
p3 <- cowplot::get_legend(p3)
p4 <- cowplot::get_legend(p4)

g1 <- plot_grid(p,p1,p2,p3,p4,ncol = 1,rel_heights = c(1,0.1,0.1,0.1,0.1))
g1

plot_grid(g1,ncol = 1)





##test : 
library(patchwork)
library(cowplot)
a <- df %>% filter(Type == "Tumor") %>% # stage 간 NAT는 NA
  select(Ethnicity,Gender,Smoking.History_modified,Stage) %>%
  rename("Smoking\nHistory" = Smoking.History_modified) %>%
  mutate_at('Ethnicity',str_replace_all,c(han = 'asian',tssdidnotcollectthisinformation = 'NA')) %>%
  mutate_at('Ethnicity',str_replace,"white\\(caucasian\\)","caucasian") %>%
  mutate(Ethnicity = str_to_title(Ethnicity),Gender = str_to_title(Gender)) %>%
  mutate_at("Stage",str_replace_all,c(IA = 'I',IB = 'I',IIA = 'II',IIB = 'II',IIIA = 'III',IIIB = 'III')) %>%
  gather() %>%
  filter(value != "Na") %>%
  filter(value != "NA") %>%
  group_by(key,value) %>%
  as.data.frame()

#level1  <- c('Caucasian','Asian','Slavic','Black','Female','Male','Non-Smoker','Reformed Smoker','Smoker',"IV","III","II","I")
level1  <- c('Ethnicity','Caucasian','Asian','Slavic','Black',' ','Gender','Female','Male','  ','Smoking History','Non-Smoker','Reformed Smoker','Smoker','   ','Stage',"IV","III","II","I")
a$value <- factor(a$value,levels = level1)
ggplot(a,aes(x=key,fill = factor(value, levels = level1)))+
  #scale_fill_manual(values = c("darkred","deepskyblue3","darkseagreen4","black","darkred","deepskyblue3","darkred","deepskyblue3","darkseagreen4","darkred","deepskyblue3","darkseagreen4","darkslateblue")) +
  geom_bar(position ="fill") + 
  scale_fill_manual(values = c("white","darkred","deepskyblue3","darkseagreen4","black","white","white","darkred","deepskyblue3","white","white","darkred","deepskyblue3","darkseagreen4","white","white","darkred","deepskyblue3","darkseagreen4","darkslateblue")) +
  theme(legend.position = "bottom",
        legend.justification = "left",
        legend.key = element_rect(fill=NA),
        legend.title = element_blank()) +
  guides(fill=guide_legend(order = 1)) +  ##
  labs(y = "Frequency",x="")

ggplot(a,aes(x=key)) + 
  geom_bar(aes(fill = level1),position= "fill") +
  scale_fill_manual(aesthetics = "fill",values = level1,
                    breaks = c("darkred","deepskyblue3","darkseagreen4","black"),
                    b)
library(ggplot2)
library(ggnewscale)
library(tidyverse)
library(relayer)
ggplot(a,aes(x=key))+
  geom_bar(aes(fill=factor(value,level1)),position="fill") +
  scale_fill_manual(aesthetics = "fill",values = c("darkred","deepskyblue3","darkseagreen4","black","darkred","deepskyblue3","darkred","deepskyblue3","darkseagreen4","darkred","deepskyblue3","darkseagreen4","darkslateblue"),
                    breaks = c('Caucasian','Asian','Slavic','Black'),name = "1st") +
  new_scale_fill() + 
  geom_bar(aes(fill=factor(value,level1)),position="fill") %>% rename_geom_aes(new_aes = c(fill = "fill")) +
  scale_fill_manual(aesthetics = "fill",values = c("darkred","deepskyblue3","darkseagreen4","black","darkred","deepskyblue3","darkred","deepskyblue3","darkseagreen4","darkred","deepskyblue3","darkseagreen4","darkslateblue"),
                    breaks = c("Female","Male"),name = "2nd") +
  guides(fill= guide_legend(ncol = 2))+
#        fill2 = guide_legend(ncol=1)) +
  theme(legend.position = "bottom")


rename_ge

#geom_bar(aes(fill2=value),position="fill") %>% rename_geom +

ggplot(diamonds, aes(color)) +
  geom_bar(aes(fill = cut)) + 
  scale_fill_manual(aesthetics = "fill", values = cut.values,
                    breaks = cut.levs[1:2], name = "First Grouop:") +
  new_scale_fill() +
  geom_bar(aes(fill2 = cut)) %>% rename_geom_aes(new_aes = c(fill = "fill2")) +
  scale_fill_manual(aesthetics = "fill2", values = cut.values,
                    breaks = cut.levs[-(1:2)], name = "Second Group:") +
  guides(fill=guide_legend(ncol = 2)) +
  theme(legend.position="bottom")
  

  
  scale_fill_manual(aesthetics = "fill",values = c("darkred","deepskyblue3","darkseagreen4","black","darkred","deepskyblue3","darkred","deepskyblue3","darkseagreen4","darkred","deepskyblue3","darkseagreen4","darkslateblue"),
                    breaks = level1[1:4],name="1st") +  
  scale_fill_manual(aesthetics = "fill",values = c("darkred","deepskyblue3","darkseagreen4","black","darkred","deepskyblue3","darkred","deepskyblue3","darkseagreen4","darkred","deepskyblue3","darkseagreen4","darkslateblue"),
                  breaks = level1[5:6],name="2nd")  



  scale_fill_manual(aesthetics = "fill", values = cut.values,
                    breaks = cut.levs[1:2], name = "First Grouop:") +
  scale_fill_manual(aesthetics = "fill2", values = cut.values,
                    breaks = cut.levs[-(1:2)], name = "Second Group:")

  

geom_bar(aes(x = key,fill=factor(value,levels = level1))) + 
  geom_bar(position = 'fill')

             
             ggplot(a,aes(x=key,fill=factor(value,levels = c('Caucasian','Asian','Slavic','Black','Female','Male','Non-Smoker','Reformed Smoker','Smoker',"IV","III","II","I")))) + 
  geom_bar(position = 'fill')+
  theme(legend.position ="none") + 
  #theme(legend.position ="bottom") + 
  scale_fill_manual(values = c("darkred","deepskyblue3","darkseagreen4","black","darkred","deepskyblue3","darkred","deepskyblue3","darkseagreen4","darkred","deepskyblue3","darkseagreen4","darkslateblue")) +
  labs(y = "Frequency",x="")
