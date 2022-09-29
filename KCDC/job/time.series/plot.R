## plot

library(readxl)
library(tidyverse)
library(zoo)
library(stringr)
setwd("~/Desktop/KCDC/job/time_series_data/")


#df <- read_excel("../all.cities.xlsx")
library(lubridate)
library(timetk)
library(extrafont)
library(plotly)
library(ggplot2)
library(scales)



df <- read_excel("./all.cities_v2_20220921.xlsx")
font_import()
interactive <- FALSE
head(df)

head(df)
a <- df %>% select(date,city,generation,value,type) %>%
  na.omit()
cities = df$city %>% unique()
types = df$type %>% unique()
dates = df$date %>% unique()
ages = df$generation %>% unique()
par(mfrow=c(6,3))

for (t in types) {
  df %>% select(date,city,generation,value,type) %>%
    na.omit() %>%
    filter(type==t) %>% 
    ggplot(aes(x=date,y=value)) +
    geom_line(aes(color=generation)) + 
    scale_x_datetime(breaks = date_breaks("1 year"),labels = date_format("%Y")) +
    facet_wrap(~city,ncol = 3,scales = 'free') + 
    geom_vline(xintercept=as.POSIXct("2020-01-01"), linetype="dashed", color = "red") +
    labs(x = "Frequency", y = "Date",title = t,color="")
  
  ggsave(paste0("./plot/",t,".png"),width = 10,height = 12)
}
p
dev.off()

df %>% select(date,city,generation,value,type) %>%
  na.omit() %>%
  filter(type=="수두") %>% 
  ggplot(aes(x=date,y=value)) +
  geom_line(aes(color=generation)) + 
  scale_x_datetime(breaks = date_breaks("1 year"),labels = date_format("%Y")) +
  facet_wrap(~city,ncol = 3,scales = 'free') + 
  geom_vline(xintercept=as.POSIXct("2020-01-01"), linetype="dashed", color = "red") +
  labs(x = "Frequency", y = "Date",title = "하이",color="")

ggsave("수두.png",width = 10,height = 15)



head(df)
df %>% select(date,year,city,generation,value,type) %>% 
  mutate(date = as.Date(date)) %>%
  pivot_wider(names_from = generation) %>%
  na.omit() %>%
  #filter(type=="수두",total >= 100,city=="경기") %>%
  filter(type=="수두") %>%  
  pivot_longer(cols = ages,names_to = "generation",values_to = "value") %>%
  filter(generation != "total") %>%
  ggplot(aes(x=date,y=value,fill = generation)) + 
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~city,ncol = 3,scales = 'free') + 
  geom_vline(xintercept=as.Date("2020-01-01"), linetype="dashed", color = "red") +
  labs(x = "", y = "",title = "하이",color="")

  
for (t in types) {
  df %>% select(date,year,city,generation,value,type) %>% 
    mutate(date = as.Date(date)) %>%
    pivot_wider(names_from = generation) %>%
    na.omit() %>%
    #filter(type=="수두",total >= 100,city=="경기") %>%
    filter(type==t) %>%  
    pivot_longer(cols = ages,names_to = "generation",values_to = "value") %>%
    filter(generation != "total") %>%
    ggplot(aes(x=date,y=value,fill = generation)) + 
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~city,ncol = 3,scales = 'free') + 
    geom_vline(xintercept=as.Date("2020-01-01"), linetype="dashed", color = "red") +
    labs(x = "", y = "",title = t,color="")
  
  ggsave(paste0("./plot/",t,"_hist.png"),width = 10,height = 12)
}

###### 
library(tidyverse)
library(tidymodels)
library(modeltime)
library(timetk)
library(lubridate)
library(readxl)
library(extrafont)
library(plotly)
library(zoo)
library(stringr)
font_import()
interactive <- FALSE
par(family = 'D2Coding')
par(family = 'NanumGothic')
par(family = "AppleGothic")

setwd("~/Desktop/KCDC/job/time_series_data/")
df <-read_excel("all.cities_v2_20220921.xlsx")
head(df)
df %>% filter(type == "total",type == "수두") %>%
  group_by(city) %>%
  plot_time_series(date,value, #.color_var = generation,
                   .smooth = F,
                   .title = "수두",
                   .facet_ncol = 3, .facet_scales = "free", .interactive = interactive)



df %>% select(date,city,type,value,generation) %>% 
  filter(type == '수두') %>%
  group_by(city) %>%
  plot_time_series(date,value, .color_var = generation,
                   .smooth = F,
                   .title = "수두",
                   .facet_ncol = 3, .facet_scales = "free", .interactive = interactive) + 
  theme_set(theme_gray(base_family='NanumGothic')) 
  
