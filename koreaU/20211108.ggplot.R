### 수업
library(tidyverse)
library(cowplot)

ggplot()


data("population", package = 'tidyr')

tidy1 = population %>% filter(country %in% c("Republic of Korea","China","Canada"))


ggplot(tidy1) + aes(x=year,y=population,color=country) + geom_point()


#geom_bar(stat = 'count') #

#facet_grid(scales = )

population %>% filter(country %in% c("Republic of Korea","China","Canada")) %>%
  ggplot() + 
  aes(x = year, y=population, color = country) +
  geom_point() +
  geom_line(aes(group=country)) + 
  scale_color_brewer(palette = 'Set1') + 
  #theme_bw() +
  theme_bw(base_size = 7) +
  #theme(axis.line = element_line(size = 0.5))+
  theme(axis.line.x = element_line(size = 0.5))
  
  

p1 <- population %>% 
  filter(country %in% c('Republic of Korea', 'China', 'Canada')) %>%
  ggplot(.) + 
  aes(x = year, y = population, color=country) + 
  geom_point() + 
  geom_line(aes(group = country)) + 
  scale_color_brewer(palette='Set1') + 
  theme_bw(base_size = 7) + 
  theme(axis.line.x = element_line(size = 1))

p2 <- population %>% 
  filter(country %in% c('Chile','Brazil','Mexico')) %>%
  ggplot(.) + 
  aes(x = year, y = population, color=country) + 
  geom_point() + 
  geom_line(aes(group = country)) + 
  scale_color_brewer(palette='Set1') + 
  theme_bw(base_size = 7) + 
  theme(axis.line.x = element_line(size = 1))
#scale_y_discrete(value=c('Mother','Father'))
#scale_y_log10()

cowplot::plot_grid(p1,p2, ncol = 2,labels = c("A","B"))
cowplot::plot_grid(p1,p2,p1,p2, ncol = 2,labels = c("A","B","C","D"))
cowplot::plot_grid(plot_grid(p1,p2, ncol = 2,labels=c("A","B")),
                   plot_grid(p1,p2,p1, ncol = 3,labels=c("C","D","E"))
                   ,ncol=1)

cowplot::plot_grid(plot_grid(p1,p2, ncol = 2,labels=c("A","B")),
                   plot_grid(p1,p2,p1, ncol = 3,labels=c("C","D","E"))
                   ,ncol = 1,rel_heights = c(1,2))


cowplot::plot_grid(plot_grid(p1,NULL, ncol = 2,labels=c("A","B")),
                   plot_grid(p1,p2,p1, ncol = 3,labels=c("C","D","E"))
                   ,ncol = 1,rel_heights = c(1,2))





