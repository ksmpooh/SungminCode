## Figure 3
#setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/figure/")
theme_step1 <- function(base_size = 11, base_family = "",
                        base_line_size = base_size / 22,
                        base_rect_size = base_size / 22) {
  theme(title = element_text(family = 'Arial', size = 18, color = 'black'), text = element_text(family = 'Arial', size = 16, color = 'black'),
        axis.title = element_text(family = 'Arial', size = 18, color = 'black'), axis.text = element_text(family = 'Arial', size = 16, color = 'black'), 
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = NA), axis.line = element_line(colour = "black", size = rel(1)),
        #legend.background = element_rect(color = 'black'), 
        legend.title = element_text(family = 'Arial', size = 16),
        legend.text = element_text(family = 'Arial', size = 14),
        legend.direction = "vertical", 
        legend.box = c("horizontal", "vertical"),
        legend.spacing.x = unit(0.1, 'cm'),
        plot.margin = unit(c(0.25, 1, 1, 0.5), 'cm'),
        axis.title.y = element_text(margin = margin(r = 10, unit = "pt"))) }

library(tidyverse)
library(ggplot2)
library(googledrive)
library(googlesheets4)


##### variant
url <- "https://docs.google.com/spreadsheets/d/1Y8cqzCizaQ9vvAzoAEhqqGW7GaplxRRnmPsQ0dK6o30/edit?usp=sharing"
drive_auth()

props <- sheet_properties(url)
props$name
#13,"1.3.asm_correcteness1_inspector"
#14,"1.3.asm_correcteness1_merqury"   
#11 "1.2.assembly_contiguity" 

graph_variant_count_GRCh38 <- read_sheet(url,sheet = "2.1.variant_landscape",skip = 1)[1:24,]
tail(graph_variant_count_GRCh38)
graph_variant$Scale_SVtype
graph_variant$Tag

graph_variant_length_GRCh38 <- read_sheet(url,sheet = "2.1.variant_landscape",skip = 28)[1:24,]
#graph_variant_length_GRCh38


graph_variant_count_CHM13 <- read_sheet(url,sheet = "2.1.variant_landscape",skip = 55)[1:24,]
graph_variant_count_CHM13

variant_colors <- c(
  "DEL" = "#D55E00",  # 붉은 계열 (deletion → loss)
  "INS" = "#009E73",  # 초록 계열 (insertion → gain)
  "INV" = "#CC79A7",  # 보라 계열 (structural rearrangement)
  "SNV" = "#56B4E9"   # 파랑 계열 (base-level variation)
)

#graph_variant_length_CHM13 <- read_sheet(url,sheet = "2.1.variant_landscape",skip = 82)[1:24,]
#graph_variant_length_CHM13

graph_variant_count_GRCh38 %>% as.data.frame() %>% #head()
  pivot_longer(3:ncol(graph_variant_count_GRCh38)) %>% 
  filter(value != "NA") %>% mutate(value = as.numeric(value)) %>%
  mutate(size = str_split_fixed(Scale_SVtype,"_",2)[,1]) %>%
  mutate(type = str_split_fixed(Scale_SVtype,"_",2)[,2]) %>%
  filter(size == "Small") %>% #count(Tag)
  #mutate(Tag = factor(Tag,levels = c("Singleton", "Polymorphic", "Major", "Shared"))) %>%
  mutate(Tag = factor(Tag,levels = c("Shared", "Major", "Polymorphic", "Singleton"))) %>%
  ggplot(aes(x=name,y=value,fill=type,alpha=Tag)) + 
  scale_alpha_manual(
    values = c(
      "Shared" = 0.4,        # 진함
      "Major" = 0.6,        # 진함
      "Polymorphic" = 0.8,        # 진함
      "Singleton" = 1    # 연함 (원하면 0도 가능)
    )) + 
  geom_bar(stat = 'identity') + 
  labs(x="Sample",y="Number of Small Variant",title = "GRCh38") + 
  scale_fill_manual(values = variant_colors) + 
  theme_step1() + 
  theme(axis.text.x = element_blank(),
        legend.title = element_blank())


graph_variant_count_CHM13 %>% as.data.frame() %>% #head()
  pivot_longer(3:ncol(graph_variant_count_CHM13)) %>% 
  filter(value != "NA") %>% mutate(value = as.numeric(value)) %>%
  mutate(size = str_split_fixed(Scale_SVtype,"_",2)[,1]) %>%
  mutate(type = str_split_fixed(Scale_SVtype,"_",2)[,2]) %>%
  filter(size == "Small") %>% #count(Tag)
  #mutate(Tag = factor(Tag,levels = c("Singleton", "Polymorphic", "Major", "Shared"))) %>%
  mutate(Tag = factor(Tag,levels = c("Shared", "Major", "Polymorphic", "Singleton"))) %>%
  ggplot(aes(x=name,y=value,fill=type,alpha=Tag)) + 
  scale_alpha_manual(
    values = c(
      "Shared" = 0.4,        # 진함
      "Major" = 0.6,        # 진함
      "Polymorphic" = 0.8,        # 진함
      "Singleton" = 1    # 연함 (원하면 0도 가능)
    )) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = variant_colors) + 
  labs(x="Sample",y="Number of Small Variant",title = "CHM13v2") + 
  theme_step1() + 
  theme(axis.text.x = element_blank(),
        legend.title = element_blank())


graph_variant_count_CHM13 %>% as.data.frame() %>% #head()
  pivot_longer(3:ncol(graph_variant_count_CHM13)) %>% 
  filter(value != "NA") %>% mutate(value = as.numeric(value)) %>%
  mutate(size = str_split_fixed(Scale_SVtype,"_",2)[,1]) %>%
  mutate(type = str_split_fixed(Scale_SVtype,"_",2)[,2]) %>%
  filter(size != "Small") %>% #count(Tag)
  #mutate(Tag = factor(Tag,levels = c("Singleton", "Polymorphic", "Major", "Shared"))) %>%
  mutate(type = factor(type,levels = c("DEL", "INS", "SNV"))) %>%
  mutate(Tag = factor(Tag,levels = c("Shared", "Major", "Polymorphic", "Singleton"))) %>%
  ggplot(aes(x=name,y=value,fill=type,alpha=Tag)) + 
  scale_alpha_manual(
    values = c(
      "Shared" = 0.4,        # 진함
      "Major" = 0.6,        # 진함
      "Polymorphic" = 0.8,        # 진함
      "Singleton" = 1    # 연함 (원하면 0도 가능)
    )) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = variant_colors) + 
  labs(x="Sample",y="Number of Structure Variant",title = "CHM13v2") + 
  theme_step1() + 
  theme(axis.text.x = element_blank(),
        legend.title = element_blank())





graph_variant_count_GRCh38 %>% as.data.frame() %>% #head()
  pivot_longer(3:ncol(graph_variant_count_GRCh38)) %>% 
  filter(value != "NA") %>% mutate(value = as.numeric(value)) %>%
  mutate(size = str_split_fixed(Scale_SVtype,"_",2)[,1]) %>%
  mutate(type = str_split_fixed(Scale_SVtype,"_",2)[,2]) %>%
  filter(size != "Small") %>% #count(Tag)
  #mutate(Tag = factor(Tag,levels = c("Singleton", "Polymorphic", "Major", "Shared"))) %>%
  mutate(Tag = factor(Tag,levels = c("Shared", "Major", "Polymorphic", "Singleton"))) %>%
  ggplot(aes(x=name,y=value,fill=type,alpha=Tag)) + 
  scale_alpha_manual(
    values = c(
      "Shared" = 0.4,        # 진함
      "Major" = 0.6,        # 진함
      "Polymorphic" = 0.8,        # 진함
      "Singleton" = 1    # 연함 (원하면 0도 가능)
    )) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = variant_colors) + 
  labs(x="Sample",y="Number of Structure Variant",title = "GRCh38") + 
  theme_step1() + 
  theme(axis.text.x = element_blank(),
        legend.title = element_blank())


####### bar

df <- tribble(
  ~ref,     ~type,   ~count,
  "CHM13",  "INDEL", 10128266,
  "CHM13",  "MNP",     545598,
  #"CHM13",  "OTHER",    29669,
  "CHM13",  "SNP",   17758700,
  "GRCh38", "INDEL",  8945965,
  "GRCh38", "MNP",     398308,
  #"GRCh38", "OTHER",    27839,
  "GRCh38", "SNP",   16550729
) %>%
  mutate(
    #type = factor(type, levels = c("SNP","INDEL","MNP","OTHER")),
    type = factor(type, levels = c("SNP","INDEL","MNP")),
    ref  = factor(ref,  levels = c("CHM13","GRCh38"))
  )


head(df)

ref_cols <- c(
  "CHM13"  = "#3B82F6",  # blue
  "GRCh38" = "#F97316"   # orange
)
head(df)


df %>% 
  ggplot(aes(x = type, y = count, fill = ref)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.72,
    color = "black",
    linewidth = 0.25) +
  scale_fill_manual(values = ref_cols) +
  scale_y_continuous(
    labels = scales::label_number(scale = 1e-6),   # M 단위
    expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL,y = expression("Variant count (" * 10^6 * ")")) +
  theme_step1() +
  theme(
    legend.position = c(0.98, 0.98),        # plot 내부 오른쪽 위
    legend.justification = c(1, 1),
    legend.title = element_blank())


