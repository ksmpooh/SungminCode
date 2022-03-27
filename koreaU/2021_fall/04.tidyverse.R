# tidyverse
library(dslabs)
data(murders)
library(tidyverse)

murders


data(murders)
head(murders)

murders <-mutate(murders,rate=total/population * 100000)
filter(murders,rate <= 0.71)
new_table <- select(murders,state,region,rate)
head(new_table)

library(dplyr)
library(dslabs)
#data(murders)
murders <- mutate(murders, population_in_millions = population / 10^6)

select(murders, state, population) %>% head()

filter(murders, state == "New York")
no_florida <- filter(murders, state != "Florida")
filter(murders, state %in% c("New York", "Texas"))
filter(murders, population < 5000000 & region == "Northeast")

murders %>% select(state, region, rate) %>% filter(rate <= 0.71)
16 %>% sqrt()
16 %>% sqrt() %>% log2()
16 %>% sqrt() %>% log(base = 2)
murders %>% select(state, region, rate) %>% filter(rate <= 0.71)
murders <- mutate(murders, rate =  total / population * 100000, 
                  rank = rank(-rate))
my_states <- filter(murders, region %in% c("Northeast", "West") & 
                      rate < 1)

select(my_states, state, rate, rank)
head(murders)
pull(murders)


#######수업

library(tidyverse)
library(dslabs)

data(murders)

filter(murders, region == "South")
murders %>% filter(region == "South")

murders$state
murders[,'state',drop = F]
murders %>% select(state)
murders %>% pull(state)


murders %>% select(state,abb)
murders %>% select(abb,state)
murders %>% select(abb,State = state )


head(murders)
murders <- murders %>% mutate(rate2 = total/population* 1000000)

paste("Today","Is","Monday")
head(murders)

#murders$paste_cal <-paste(murders$state,murders$abb)
murders <- murders %>% mutate(paste_col = paste(state,abb))

murders %>% select(rate) %>% mutate(rate2 = )
                              
murders %>% filter(region == "West")%>% mutate(rate = total/population * 1000000)



murders %>% 
  filter(region == "West")%>% 
  mutate(rate = total/population * 1000000)


sum(1,2,3,4)
head(murders)

murders %>% 
  filter(region == "Northeast") %>% 
  select(total) %>% 
  sum()

murders %>% 
  filter(region == "Northeast") %>% 
  pull(total) %>% 
  sum()

murders %>% 
  filter(region == "Northeast") %>% 
  select(total,population) %>% 
  colSums()

data("heights")
head(heights)
heights %>%
  filter(sex == "Female") %>%
  summarize(median = median(height), minimum = min(height), maximum = max(height))

heights %>%
  group_by(sex) %>%
  summarize(median = median(height), minimum = min(height), maximum = max(height))
# tibble 이라는 데이터 형태로 변경됨

murders %>%
  group_by(region) %>%
  summarise(median = median(total),minimum = min(rate),maximum = max(rate),mean = mean(rate))

murders %>% arrange(desc(population))
murders %>% arrange(population,decreasing = T)

murders %>%
  top_n(3,population)


df <- data.frame(a = seq(2,10),b=seq(3,11),c=seq(14,22))
colMeans(df)

apply(df, 2, mean)


library(readxl)
read_excel("",sheet = 2)

readxl::excel_sheets() # 각 sheet 이름이 나옴
