setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/construction/03.hifiasm/")

theme_step1 <- function(base_size = 11, base_family = "",
                        base_line_size = base_size / 22,
                        base_rect_size = base_size / 22) {
  theme(title = element_text(family = 'Arial', size = 18, color = 'black'), text = element_text(family = 'Arial', size = 16, color = 'black'),
        axis.title = element_text(family = 'Arial', size = 18, color = 'black'), axis.text = element_text(family = 'Arial', size = 16, color = 'black'), 
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = NA), axis.line = element_line(colour = "black", size = rel(1)),
        legend.background = element_rect(color = 'black'), legend.title = element_text(family = 'Arial', size = 16),
        legend.text = element_text(family = 'Arial', size = 14),
        legend.direction = "vertical", 
        legend.box = c("horizontal", "vertical"),
        legend.spacing.x = unit(0.1, 'cm'),
        plot.margin = unit(c(0.25, 1, 1, 0.5), 'cm'),
        axis.title.y = element_text(margin = margin(r = 10, unit = "pt"))) }


library(tidyverse)

library(ggplot2)
library(plyr)
library(dplyr)
head(df)
df <- read.csv("length.csv", header = TRUE)
hprc <- read.csv("merge.HPRC.phaseI.NGx.length.csv", header = F)
grch <- read.csv("length_grch38.csv", header = F)
grch <- read.csv("lnegth_grch38_plamary_scaf.csv", header = F)
head(hprc)
head(grch)
colnames(hprc) <- colnames(df)
colnames(grch) <- colnames(df)
grch %>% mutate(Tyep = "GRCh38")

hprc %>% mutate(Type = case_when(
  Sample == "CHM13Y.fa" ~ "CHM13",
  str_detect(Sample,"GR") ~ "GRCh38",
  TRUE ~ "HPRC")) %>% filter(Type != "GRCh38") %>%
  rbind(grch %>% mutate(Type = "GRCh38")) %>% #head()
  rbind(df %>% mutate(Type = "KoGES")) %>% #head
  mutate(Type = as.factor(Type)) %>%# count(Type)
  mutate(Coverage = Coverage*3) %>%
  ggplot(aes(x=Coverage, y=Length/1000000, color=Type, group=Sample)) +
  geom_vline(xintercept = 1.5, linetype="dotted", linewidth=0.5) +
  geom_step(alpha=0.8) +
  #scale_color_manual(values = colors) +
  #scale_x_continuous(limits = c(0, 1.25), breaks = seq(0, 1.25, by = 0.25)) +
  labs(x = "Cumulative coverage (Gb)", y = "Length (Mb)") +
  theme_bw()



hprc %>%
  mutate(Type = case_when(
    Sample == "CHM13Y.fa" ~ "CHM13",
    str_detect(Sample, "GR") ~ "GRCh38",
    TRUE ~ "HPRC"
  )) %>%
  filter(Type != "GRCh38") %>%
  rbind(grch %>% mutate(Type = "GRCh38")) %>%
  rbind(df %>% mutate(Type = "KoGES")) %>%
  mutate(Type = as.factor(Type)) %>%
  mutate(Coverage = Coverage * 3) %>%
  ggplot(aes(x = Coverage, y = Length / 1000000, color = Type, group = Sample)) +
  geom_vline(xintercept = 1.5, linetype = "dotted", linewidth = 0.5) +
  geom_step(size = 0.5,alpha=0.5) +
  #scale_color_manual(values = c("GRCh38" = "black", "CHM13" = "black", "HPRC" = "red", "KoGES" = "blue")) +
  scale_alpha_manual(values = c("GRCh38" = 1, "CHM13" = 1, "HPRC" = 0.3, "KoGES" = 0.3)) +
  labs(x = "Cumulative coverage (Gb)", y = "Length (Mb)") +
  theme_bw() +
  guides(alpha = "none")

library(dplyr)
library(ggplot2)

# Data preparation
processed_data <- hprc %>%
  mutate(Type = case_when(
    Sample == "CHM13Y.fa" ~ "CHM13",
    str_detect(Sample, "GR") ~ "GRCh38",
    TRUE ~ "HPRC"
  )) %>%
  filter(Type != "GRCh38") %>%
  rbind(grch %>% mutate(Type = "GRCh38")) %>%
  rbind(df %>% mutate(Type = "KoGES")) %>%
  mutate(Type = as.factor(Type)) %>%
  mutate(Coverage = Coverage * 3)

# Plot
ggplot() +
  # Plot HPRC and KoGES first with lower alpha
  geom_step(
    data = processed_data %>% filter(Type %in% c("HPRC", "KoGES")),
    aes(x = Coverage, y = Length / 1000000, color = Type, group = Sample),
    size = 0.5, alpha = 0.2
  ) +
  # Plot GRCh38 and CHM13 with higher alpha and black color, no legend
  geom_step(
    data = processed_data %>% filter(Type %in% c("GRCh38", "CHM13")),
    aes(x = Coverage, y = Length / 1000000, group = Sample),
    color = "black", size = 0.5, alpha = 1, show.legend = FALSE
  ) +
  geom_vline(xintercept = 1.5, linetype = "dotted", linewidth = 0.5) +
  scale_color_manual(values = c("HPRC" = "red", "KoGES" = "blue")) +
  labs(x = "Cumulative coverage (Gb)", y = "Length (Mb)") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  # Annotations for CHM13 and GRCh38
  annotate("text", x = 0.8, y = 200, label = "CHM13", color = "black", size = 5, hjust = 0) +
  annotate("text", x = 0.8, y = 100, label = "GRCh38", color = "black", size = 5, hjust = 0)





ggsave("coverage.png", plot=coveragePlot)
pdf("coverage.pdf",width=8,height=4,paper='special')
print(coveragePlot)
dev.off()
head(df)
df %>% count(Length < LengthSum)
df %>% mutate(Coverage2 = Coverage/10*3) %>% head()
df %>% mutate(Coverage2 = LengthSum/Length*100*300000000) %>% head()

df %>% mutate(Coverage = Coverage/10*3) %>%  #head()
  ggplot(aes(x=Coverage, y=Length/1000000, color=Type, group=Sample)) +
  geom_vline(xintercept = 1.5, linetype="dotted", linewidth=0.5) +
  geom_step(alpha=0.5) +
  scale_color_manual(values = colors) +
  scale_x_continuous(limits = c(0, 3.0), breaks = seq(0, 3, by = 1)) +
  labs(x = "Cumulative coverage", y = "Length (Mb)") +
  theme_bw()

#### Busco visual

library(tidyverse)
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/construction/03.hifiasm/busco/prep/")
flist = grep(list.files("./"),pattern = "txt", value=TRUE)
flist
busco <- NULL
for (i in flist) {
  tmp <- read_table(i,col_names = T)
  tmp$filename = i
  #  tmp %>% filter(X6 == "PASS") %>% select(X3,X4,X5) -> tmp
  busco <- rbind(busco,tmp)
}
head(busco)
colnames(busco) <- c("Complete","Complete_single_copy","Complete_duplicated","Fragmented","Missing","total","filename")
busco_percet <- busco
busco_percet[1:5] <- lapply(busco_percet[1:5], function(x) x / busco_percet$total * 100)

head(busco_percet)


t2t<- data.frame(
  Complete = 95.7,
  Complete_single_copy = 94.1,
  Complete_duplicated = 1.6,
  Fragmented = 1.1,
  Missing = 3.2,
  total = 13870,
  filename = "T2T-CHM13"
) 
head(t2t)

head(busco_percet)
busco_percet %>% rbind(t2t) %>% mutate(filename = ifelse(str_detect(filename,"T2T"),filename,"KPP")) %>% 
  group_by(filename) %>% #count(filename)
  summarise(
    Complete = sprintf("%.2f ± %.2f", mean(Complete), sd(Complete)),
    Complete_single_copy = sprintf("%.2f ± %.2f", mean(Complete_single_copy), sd(Complete_single_copy)),
    Complete_duplicated = sprintf("%.2f ± %.2f", mean(Complete_duplicated), sd(Complete_duplicated)),
    Fragmented = sprintf("%.2f ± %.2f", mean(Fragmented), sd(Fragmented)),
    Missing = sprintf("%.2f ± %.2f", mean(Missing), sd(Missing))
  ) %>% writexl::write_xlsx("/Users/ksmpooh/Desktop/KCDC/pangenome/construction/03.hifiasm/busco.stats.mean.xlsx")

writexl::write_xlsx(busco,"/Users/ksmpooh/Desktop/KCDC/pangenome/construction/03.hifiasm/busco.stats.xlsx")

library(grid)
library(ggforce)
busco_percet %>% rbind(t2t) %>% mutate(facet = ifelse(str_detect(filename,"T2T"),filename,"KPP")) %>%  select(-Complete,-total) %>%
  pivot_longer(Complete_single_copy:Missing) %>% mutate(
    name = factor(name, levels = c("Missing", "Fragmented", "Complete_duplicated", "Complete_single_copy")), # y축 순서 지정
#    filename = ifelse(str_detect(filename, "T2T"), filename, "") # T2T만 x축 텍스트 표시
  ) %>%
  ggplot(aes(x=filename,y=value,fill=name)) +  
  geom_bar(stat = 'identity',position = 'fill') +
  labs(y="Busco (%)") + 
  theme_step1() + 
  facet_row(~ facet,scales = "free", space = "free") + 
  theme(strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")
  
    

#### Quast

quast <- read.table("/Users/ksmpooh/Desktop/KCDC/pangenome/construction/03.hifiasm/quast/03.Quast/NIH23F1144740.h2.quast/report.tsv",sep="\t",comment.char = "")
head(quast)

quast %>% filter(V1 %in% c("Assembly","# contigs","Largest contig","Total length","GC (%)","N50","N90")) -> quast
head(quast)
quast %>% pivot_wider(names_from = V1,values_from = V2) %>%
  rename(`Number contigs`=`# contigs`)

setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/construction/03.hifiasm/quast/03.Quast/")
flist = grep(list.files("./"),pattern = "NIH", value=TRUE)
flist
quast <- NULL
for (i in flist) {
  tmp <- read.table(paste0(i,"/report.tsv"),sep="\t",comment.char = "")
  tmp %>% filter(V1 %in% c("Assembly","# contigs","Largest contig","Total length","GC (%)","N50","N90")) %>%
    pivot_wider(names_from = V1,values_from = V2) %>%
    rename(`Number contigs`=`# contigs`) -> tmp
  #  tmp %>% filter(X6 == "PASS") %>% select(X3,X4,X5) -> tmp
  quast <- rbind(quast,tmp)
}

head(quast)

t2t_quast<- data.frame(
  Assembly = "T2T-CHM13",
  `Number contigs` = '24',
  "Largest contig" = 248387328,
  "Total length" = 3055677272,
  "GC (%)" = 41.0,
  N50 = 154259566,
  N90 = 57882016
) 
colnames(t2t_quast) <- colnames(quast)
head(quast)
quast %>% rbind(t2t_quast) %>% mutate(across(-Assembly, as.numeric)) %>%
  mutate(Assembly = ifelse(str_detect(Assembly,"T2T"),Assembly,"KPP")) %>% 
  group_by(Assembly) %>% #count(filename)
  summarise(
    `Number contigs` = sprintf("%.2f ± %.2f", mean(`Number contigs`), sd(`Number contigs`)),
    `Largest contig` = sprintf("%.2f ± %.2f", mean(`Largest contig`), sd(`Largest contig`)),
    `Total length` = sprintf("%.2f ± %.2f", mean(`Total length`), sd(`Total length`)),
    `GC (%)` = sprintf("%.2f ± %.2f", mean(`GC (%)`), sd(`GC (%)`)),
    N50 = sprintf("%.2f ± %.2f", mean(N50), sd(N50)),
    N90 = sprintf("%.2f ± %.2f", mean(N90), sd(N90))
  ) %>% writexl::write_xlsx("/Users/ksmpooh/Desktop/KCDC/pangenome/construction/03.hifiasm/quast.stats.mean.xlsx")

writexl::write_xlsx(quast,"/Users/ksmpooh/Desktop/KCDC/pangenome/construction/03.hifiasm/quast.stats.xlsx")


head(quast)

quast %>%
  ggplot(aes(x="",y= `Number contigs`)) +
  geom_violin()

quast %>% mutate(across(-Assembly, as.numeric)) %>%
  pivot_longer(`Number contigs`:N90) %>% #head()
  mutate(
    # 패싯 그룹 설정: 위, 아래 그룹으로 나눔
    facet_group = case_when(
      name %in% c("Number contigs", "Largest contig", "Total length") ~ "Group 1: Assembly Stats",
      TRUE ~ "Group 2: Other Stats"
    ),
    # facet 순서 지정
    name = factor(name, levels = c("Number contigs", "Largest contig", "Total length", 
                                   "GC (%)", "N50", "N90"))) %>%
  ggplot(aes(x=name,y=value,fill = name)) +
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.5,alpha=0.8)  + 
  #labs(y="Busco (%)") + 
  theme_step1() + 
  #facet_grid(facet_group ~ name, scales = "free", space = "free") + # facet을 행(row)으로 분리
    
  facet_row(facet_group~ name,scales = "free", space = "free") + 
  #scale_y_continuous(labels = scales::comma) + 
  scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) + # 천 단위 쉼표 추가
  theme(strip.text = element_blank(),
        #axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal")

