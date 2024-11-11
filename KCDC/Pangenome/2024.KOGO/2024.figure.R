### KOGO final figure 
library(tidyverse)
if(T){
  extrafont::font_import(pattern = "Arial", prompt = F)
  extrafont::loadfonts()
}


  #setwd('~/Dropbox/SMC_AD_WGS_paper/Data/SMC_cwas_results_20240416/')
  # Set colors
  
  # Function definition
theme_step1 <- function(base_size = 11, base_family = "",
                          base_line_size = base_size / 22,
                          base_rect_size = base_size / 22) {
    theme(title = element_text(family = 'Arial', size = 18, color = 'black'), text = element_text(family = 'Arial', size = 16, color = 'black'),
          axis.title = element_text(family = 'Arial', size = 18, color = 'black'), axis.text = element_text(family = 'Arial', size = 16, color = 'black'), 
          #panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          #panel.background = element_rect(fill = "white", colour = NA), axis.line = element_line(colour = "black", size = rel(1)),
          #legend.background = element_rect(color = 'black'), legend.title = element_text(family = 'Arial', size = 16),
          legend.text = element_text(family = 'Arial', size = 14),
          #legend.direction = "vertical", 
          #legend.box = c("horizontal", "vertical"),
          legend.spacing.x = unit(0.1, 'cm'),
          #plot.margin = unit(c(0.25, 1, 1, 0.5), 'cm'),
          axis.title.y = element_text(margin = margin(r = 10, unit = "pt")))
}
  
resize_heights <- function(g, heights = rep(1, length(idpanels))){
    idpanels <- unique(g$layout[grepl("panel",g$layout$name), "t"])
    g$heights <- unit.c(g$heights)
    g$heights[idpanels] <- unit.c(do.call(unit, list(heights, 'null')))
    g
}


###

library(tidyverse)


df <- read.table("~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/bam.stats.merge.txt",header = T)
head(df)

platform_factor <- c("Illumina","PacBio","ONT")


df %>% filter(!(rname %in% c("chrX","chrY","chrM"))) %>% #group_by(batch,ID) %>%  head()
  mutate(batch = ifelse(batch == "Nanopore","PromethION",batch)) %>% #head()
  mutate(batch = factor(batch, levels=platform_factor)) %>%
  group_by(batch,ID) %>% #head()
  summarise(`mean Depth(X)` = sum(meandepth * endpos)/sum(endpos),
            `Coverage(%)` = sum(covbases)/sum(endpos)*100,
            `mean BaseQ` = sum(meanbaseq * endpos)/sum(endpos),
            `mean MapQ` = sum(meanmapq * endpos)/sum(endpos)) -> a

a %>% pivot_longer(3:6) %>% #head()
  ggplot(aes(x=batch,y=value,fill=batch)) + 
  geom_boxplot() + 
  facet_wrap(~name,scales = "free",nrow = 1) + 
  theme_step1()

a %>% pivot_longer(3:6) %>% #head()
  ggplot(aes(x=batch,y=value,fill=batch)) + 
  geom_boxplot() + 
  theme_step1() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x =  element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom") + 
  facet_wrap(~name,scales = "free",nrow = 1) 


head(a)
a %>% group_by(batch) %>% summarise(sd(`mean Depth(X)`),mean(`mean Depth(X)`)) -> a

### VCF
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/99.kogo/2024.kogo/vcf.pass.stats/rmMXY.stats")
df <- read.table("filter.info.txt",header = T)
ref <- read.table("~/Desktop/KCDC/pangenome/00.datacheck/KBA.Long_Revio_Nanopore_short.IDmatchinagtable.txt",header = T) %>% na.omit()
head(ref)
#colnames(ref)[1:4]<-c("Revio","PromethION","NovaSeq6000","KBAv1")
colnames(ref)[1:4]<-c("PacBio","ONT","Illumina","KBAv1")
ref %>% select(-KBAv2) %>% pivot_longer(cols = 1:3,names_to = "platform",values_to = "ID") ->ref

head(df)
df %>% mutate(ID = str_split_fixed(file,"\\.",2)[,1]) %>% filter(ID == "Pangenome") -> a
df %>% mutate(ID = str_split_fixed(file,"\\.",2)[,1]) %>% filter(ID != "Pangenome") %>% 
  mutate(ID = str_replace(ID,".sorted","")) %>%
  #  mutate(platform = ifelse(str_detect(file,"pbmm2"),"Revio",ifelse(str_detect(file,"bwa"),"NovaSeq6000","PromethION"))) -> df
  mutate(platform = ifelse(str_detect(file,"pbmm2"),"PacBio",ifelse(str_detect(file,"bwa"),"Illumina","ONT"))) -> df

head(df)
platform_factor <- c("Illumina","PacBio","ONT")

library(ggpattern)
library(ggplot2)

pdf(file = "~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/vcf.stats.pdf", width = 15, height = 3.5)
df %>% left_join(ref) %>%  
  select(KBAv1,platform,records,SNPs,indels,multiallelic_sites,ts_tv) %>% #head
  rename(`multi-allelic sites` = multiallelic_sites) %>% 
  rename(`ts/tv ratio` = ts_tv) %>%
  mutate(platform = factor(platform, levels=c(platform_factor))) %>%
  pivot_longer(cols = records:`ts/tv ratio`, names_to = "type", values_to = "value") %>% #head()
  mutate(type = factor(type, levels=c("records","SNPs","indels","multi-allelic sites","ts/tv ratio"))) %>%
  ggplot(aes(x = platform, y = value, fill = platform)) + 
  geom_boxplot() + 
  theme_step1() + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size=18)) + 
  facet_wrap(~type,scales = 'free',nrow = 1)
dev.off()



#### concordace test
setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/99.kogo/2024.kogo/concordance/concordance/snp/")
flist <- list()

flist = grep(list.files("./"),pattern = "by_sample", value=TRUE)
flist
head(flist)
df <- NULL

#revio <- read.table("concordance_revio_nanopore.by_sample.txt",header = T)
#wgs <- read.table("concordance_wgs_nanopore.by_sample.txt",header = T)
#nanopore <- read.table("concordance_wgs_nanopore.by_sample.txt",header = T)

#revio$batch <- "Revio"
#wgs$batch <- "NovaSeq6000"
#nanopore$batch <- "PromethION"


df <- NULL
flist

for (i in flist) {
  data <- read.table(i,header = T)
  data <- cbind(data, TP=(data$REF.REF + data$ALT_1.ALT_1 + data$ALT_2.ALT_2))
  data <- cbind(data, FN=(data$ALT_1.REF + data$ALT_2.REF + data$ALT_2.ALT_1))
  data <- cbind(data, FP=(data$REF.ALT_1 + data$REF.ALT_2 + data$ALT_1.ALT_2))
  
  data$Sensitivity <- data$TP/(data$TP+data$FN)
  data$Precision <- data$TP/(data$TP+data$FP)
  data$Accuracy <- data$TP/(data$TP+data$FP+data$FN)
  
  a <- data %>% #inner_join(maf1) %>% 
    select(sample,Sensitivity,Precision,Accuracy) %>% 
    pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% 
    mutate(ID = str_replace(i,".by_sample.txt",""))
  df <- rbind(df,a)
  #mutate(title = str_remove(str_remove(i,".txt"),"QulityMetrix_"))
  df <- rbind(df,a)
}
head(df)
head(df$variant)
df$variant <- "SNP"
df_snp <- df



setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/99.kogo/2024.kogo/concordance/concordance/indel/")
flist <- list()

flist = grep(list.files("./"),pattern = "by_sample", value=TRUE)
flist
head(flist)
df <- NULL


for (i in flist) {
  data <- read.table(i,header = T)
  data <- cbind(data, TP=(data$REF.REF + data$ALT_1.ALT_1 + data$ALT_2.ALT_2))
  data <- cbind(data, FN=(data$ALT_1.REF + data$ALT_2.REF + data$ALT_2.ALT_1))
  data <- cbind(data, FP=(data$REF.ALT_1 + data$REF.ALT_2 + data$ALT_1.ALT_2))
  
  data$Sensitivity <- data$TP/(data$TP+data$FN)
  data$Precision <- data$TP/(data$TP+data$FP)
  data$Accuracy <- data$TP/(data$TP+data$FP+data$FN)
  
  a <- data %>% #inner_join(maf1) %>% 
    select(sample,Sensitivity,Precision,Accuracy) %>% 
    pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% 
    mutate(ID = str_replace(i,".by_sample.indel.txt",""))
  df <- rbind(df,a)
  #mutate(title = str_remove(str_remove(i,".txt"),"QulityMetrix_"))
  df <- rbind(df,a)
}
head(df)
head(df$variant)
df$variant <- "INDEL"
df_indel <- df


df <- rbind(df_snp,df_indel)
head(df)
tail(df)


pdf(file = "~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/concordance.pdf", width = 4, height = 4)
df %>% #mutate(ID = str_replace(ID,"concordance_","")) %>%
  mutate(ID = ifelse(ID == "concordance_revio_nanopore","PacBio vs ONT",ifelse(ID == "concordance_wgs_revio","Illumina vs PacBio",'Illumina vs ONT'))) %>%
  mutate(ID = factor(ID,levels=c("Illumina vs PacBio","Illumina vs ONT","PacBio vs ONT"))) %>% 
  mutate(type = ifelse(type == "Accuracy",'ACC',ifelse(type == 'Precision','PPV','TPR'))) %>%
  mutate(variant = factor(variant, levels=c("SNP","INDEL"))) %>%
  ggplot(aes(x=type,y=Val,fill=type)) + 
  geom_boxplot() + 
  #theme_step1() + 
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  ) + 
  facet_wrap(variant~ID)
dev.off()

df %>% #mutate(ID = str_replace(ID,"concordance_","")) %>%
  mutate(ID = ifelse(ID == "concordance_revio_nanopore","Revio vs PromethION",ifelse(ID == "concordance_wgs_revio","NovaSeq6000 vs Revio",'NovaSeq6000 vs PromethION'))) %>%
  mutate(ID = factor(ID,levels=c("NovaSeq6000 vs Revio","NovaSeq6000 vs PromethION","Revio vs PromethION"))) %>%
  mutate(variant = factor(variant, levels=c("SNP","INDEL"))) %>% #head()
  group_by(variant,ID,type) %>%  #head()
  summarise(mean = mean(Val)) %>% #head()
  pivot_wider(names_from = type, values_from = mean ) -> a



df %>% #mutate(ID = str_replace(ID,"concordance_","")) %>%
  mutate(ID = ifelse(ID == "concordance_revio_nanopore","PacBio vs ONT",ifelse(ID == "concordance_wgs_revio","Illumina vs PacBio",'Illumina vs ONT'))) %>%
  mutate(ID = factor(ID,levels=c("Illumina vs PacBio","Illumina vs ONT","PacBio vs ONT"))) %>% 
  mutate(type = ifelse(type == "Accuracy",'ACC',ifelse(type == 'Precision','PPV','TPR'))) %>%
  #mutate(variant = factor(variant, levels=c("SNP","INDEL"))) %>%
  filter(variant == "SNP") %>%
  ggplot(aes(x=type,y=Val,fill=type)) + 
  geom_boxplot() + 
  #  labs(title="test") +
  ylim(c(0.980,1)) + 
  theme_step1() + 
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size=18)) +
  facet_wrap(~ID) -> p1

df %>% #mutate(ID = str_replace(ID,"concordance_","")) %>%
  mutate(ID = ifelse(ID == "concordance_revio_nanopore","PacBio vs ONT",ifelse(ID == "concordance_wgs_revio","Illumina vs PacBio",'Illumina vs ONT'))) %>%
  mutate(ID = factor(ID,levels=c("Illumina vs PacBio","Illumina vs ONT","PacBio vs ONT"))) %>% 
  mutate(type = ifelse(type == "Accuracy",'ACC',ifelse(type == 'Precision','PPV','TPR'))) %>%
  #mutate(variant = factor(variant, levels=c("SNP","INDEL"))) %>%
  filter(variant == "INDEL") %>%
  ggplot(aes(x=type,y=Val,fill=type)) + 
  geom_boxplot() + 
  ylim(c(0.980,1)) + 
  theme_step1() + 
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size=18)) + 
  facet_wrap(~ID) -> p2


pdf(file = "~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/concordance.pdf", width = 15, height = 3)
cowplot::plot_grid(p1,p2,labels = c("A","B"),label_size = 20)
dev.off()


####### varinat info

setwd("/Users/ksmpooh/Desktop/KCDC/pangenome/99.kogo/2024.kogo/concordance/stats")

revio_snp <- read_table("Pangenome.Revio_kchip.same_17sample.pbmm2_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlySNP_withtag.vcf.gz.AF_stat.txt.query",col_names = F)
revio_indel <- read_table("Pangenome.Revio_kchip.same_17sample.pbmm2_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlyindels_withtag.vcf.gz.AF_stat.txt.query",col_names = F)

wgs_snp <- read_table("Pangenome.wgs.same_17sample.bwa_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlySNP_withtag.vcf.gz.AF_stat.txt.query",col_names = F)
wgs_indel <- read_table("Pangenome.wgs.same_17sample.bwa_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlyindels_withtag.vcf.gz.AF_stat.txt.query",col_names = F)


nanopore_snp <- read_table("Pangenome.nanopore.same_17sample.minimap2_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlySNP_withtag.vcf.gz.AF_stat.txt.query",col_names = F)
nanopore_indel <- read_table("Pangenome.nanopore.same_17sample.minimap2_hg38.DeepVariant.GLnexus_Jointcalling.sampleID.rmdup.onlyindels_withtag.vcf.gz.AF_stat.txt.query",col_names = F)

head(wgs_indel)
head(revio_snp)
head(revio_indel)
head(wgs_snp)
head(wgs_indel)
head(nanopore_snp)
head(nanopore_indel)

colnames(revio_snp) <- c("ID","QUAL","AC")
colnames(revio_indel) <- c("ID","QUAL","AC")
colnames(wgs_snp) <- c("ID","QUAL","AC")
colnames(wgs_indel) <- c("ID","QUAL","AC")
colnames(nanopore_snp) <- c("ID","QUAL","AC")
colnames(nanopore_indel) <- c("ID","QUAL","AC")

library(ggplot2)
library(ggVennDiagram)
library(VennDiagram)
#platform_factor <- c("NovaSeq6000","Revio","PromethION")
platform_factor <- c("Illumina","PacBio","ONT")
x=list(A = wgs_snp$ID, B = revio_snp$ID, C = nanopore_snp$ID)
x=list(Illumina = wgs_snp$ID, PacBio = revio_snp$ID, ONT = nanopore_snp$ID)


ggVennDiagram(x, label = "count",set_size = TRUE) + 
  scale_fill_gradient(low = "lightblue", high ="white") + 
  theme(legend.position = 'none')


vd <- VennDiagram::venn.diagram(x, filename = NULL,)
grid::grid.draw(vd)
draw.triple.venn()

wgs_snp_only<-wgs_snp %>% filter(!(ID %in% revio_snp$ID)) %>% filter(!(ID %in% nanopore_snp$ID)) %>% dim()
revio_snp_only<-revio_snp %>% filter(!(ID %in% wgs_snp$ID)) %>% filter(!(ID %in% nanopore_snp$ID)) %>% dim()
nanopore_snp_only<-nanopore_snp %>% filter(!(ID %in% wgs_snp$ID)) %>% filter(!(ID %in% revio_snp$ID)) %>% dim()

wgs_revio_snp <- wgs_snp %>% filter(ID %in% revio_snp$ID) %>% filter(!(ID %in% nanopore_snp$ID)) %>% dim()
wgs_nanopore_snp <- wgs_snp %>% filter(ID %in% nanopore_snp$ID) %>% filter(!(ID %in% revio_snp$ID)) %>% dim()
revio_nanopore_snp <- revio_snp %>% filter(ID %in% nanopore_snp$ID) %>% filter(!(ID %in% wgs_snp$ID)) %>% dim()

all_snp <- revio_snp %>% filter(ID %in% nanopore_snp$ID) %>% filter(ID %in% wgs_snp$ID) %>% dim()



x=c()
library(venneuler)
v <-venneuler::venneuler(c(A=wgs_snp_only[1],B=revio_snp_only[1],C=nanopore_snp_only[1],"A&B"=wgs_revio_snp[1],"A&C"=wgs_nanopore_snp[1],"C&B"=revio_nanopore_snp[1],"A&B&C"=all_snp[1]))
plot(v)
#v$labels <- c("")  
#platform_factor <- c("NovaSeq6000","Revio","PromethION")
#x=list(A = wgs_indel$ID, B = revio_snp$ID, C = nanopore_snp$ID)
#x=list(NovaSeq6000 = wgs_indel$ID, Revio = revio_indel$ID, PromethION = nanopore_indel$ID)

x=list(Illumina = wgs_snp$ID, PacBio = revio_snp$ID, ONT = nanopore_snp$ID)

ggVennDiagram(x, label = "count",set_size = 0) + 
  scale_fill_gradient(low = "white", high ="lightblue") + 
  #scale_fill_distiller(palette = "RdBu") + 
  theme(legend.position = 'none') -> p1



x=list(Illumina = wgs_indel$ID, PacBio = revio_indel$ID, ONT = nanopore_indel$ID)

ggVennDiagram(x, label = "count",set_size = 0) + 
  scale_fill_gradient(low = "white", high ="lightblue") + 
  #scale_fill_distiller(palette = "RdBu") + 
  theme(legend.position = 'none') -> p2



cowplot::plot_grid(p1,p2,labels = c("A","B"))


x=list(Illumina = wgs_snp$ID, PacBio = revio_snp$ID, ONT = nanopore_snp$ID)

ggVennDiagram(x, label = "count",set_size = 0,edge_size = 4) + 
  scale_fill_distiller(palette = "RdBu") + 
  theme(legend.position = 'none') -> p1


x=list(Illumina = wgs_indel$ID, PacBio = revio_indel$ID, ONT = nanopore_indel$ID)

ggVennDiagram(x, label = "count",set_size = 0,edge_size = 4) + 
  scale_fill_distiller(palette = "RdBu") + 
    theme(legend.position = 'none') -> p2

library(scales)
show_col(hue_pal()(3))

pdf(file = "~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/venn.pdf", width = 6, height = 4)
cowplot::plot_grid(p1,p2,labels = c("A","B"),align = "hv",vjust = 6)
dev.off()



x=list(Illumina = wgs_snp$ID, PacBio = revio_snp$ID, ONT = nanopore_snp$ID)
venn.diagram(
  x = x,
  filename = '~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/venn.snp.dia.png',
  output=T,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  #fill = myCol,
  col=c("#F8766D", '#00BA38', '#619CFF'),
  fill = c(alpha("#F8766D",1), alpha('#00BA38',1), alpha('#619CFF',1)),
  
  # Numbers
  cex = .4,
  fontfamily = "sans", 
  
  # Set names
  cat.cex = 0.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "Arial",
  rotation = 1
)



x=list(Illumina = wgs_indel$ID, PacBio = revio_indel$ID, ONT = nanopore_indel$ID)
venn.diagram(
  x = x,
    filename = '~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/venn.indel.dia.png',
  output=T,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  #fill = myCol,
  col=c("#F8766D", '#00BA38', '#619CFF'),
  fill = c(alpha("#F8766D",1), alpha('#00BA38',1), alpha('#619CFF',1)),
  
  # Numbers
  cex = .4,
  fontfamily = "sans", 
  
  # Set names
  cat.cex = 0.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "Arial",
  rotation = 1
)



library(png)
snp <- readPNG("~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/venn.snp.dia.png")
indel <- readPNG("~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/venn.indel.dia.png")
snp

#snp <- readImage("~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/venn.snp.dia.png", resize = NULL, rotate = NULL)
#indel <- readImage("~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/venn.indel.dia.png", resize = NULL, rotate = NULL)


pdf(file = "~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/venn.dia.test.pdf", width = 6, height = 4)
cowplot::plot_grid(image_ggplot(snp),image_ggplot(indel),labels = c("A","B"))
dev.off()


ggplot() + annotation_custom(rasterGrob(snp, width=unit(1,"npc"), height=unit(1,"npc"))) -> p1
ggplot() + annotation_custom(rasterGrob(indel, width=unit(1,"npc"), height=unit(1,"npc"))) -> p2



cowplot::plot_grid(p1,p2,labels = c("A","B"),label_size = 20,vjust = 1.1)




######## indel
revio_indel %>% rbind(nanopore_indel) %>% select(ID) %>% unique() %>% filter(!(ID %in% wgs_indel$ID)) %>% #dim()
  mutate(ref = str_split_fixed(ID,"_",4)[,3],alt = str_split_fixed(ID,"_",4)[,4]) %>%
  mutate(ref_length = str_length(ref),alt_length = str_length(alt)) %>%
  mutate(indel_length = ifelse(ref_length > alt_length,ref_length,alt_length)) -> indel_length

head(indel_length)

indel_length %>% filter(ID %in% nanopore_indel$ID)
#install.packages("ggbreak")
library(tidyverse)
library(ggplot2)
library(ggbreak)
#remotes::install_version("ggplot2", version = "3.4.4", repos = "http://cran.us.r-project.org")
indel_length %>% mutate(group = ifelse(indel_length > 50,">50bp","2~50bp")) %>%
  filter(indel_length %in% c(2:50)) %>%
  ggplot(aes(x=indel_length,fill=group)) +
  geom_histogram() + 
  scale_y_break(c(25000, 70000))
  #scale_y_break(c(25000, 70000), scales = "free") + 
  #scale_y_break(c(80000, 140000), scales = "free")
  #scale_y_continuous(
    #breaks = c(0, 100, 200, 300, 400, 100000, 150000), 
    #labels = c(0, 100, 200, 300, 400, 100000, 150000)) + 
  #coord_cartesian(ylim = c(0, 100000)) + 
  #annotate("segment", x = -Inf, xend = Inf, y = 40, yend = 100, linetype = "dashed", color = "red")
  
indel_length %>% count(indel_length >=50)
indel_length %>% mutate(group = ifelse(indel_length > 50,">50bp","2~50bp")) %>%
  filter(indel_length >=50) %>%
  #filter(indel_length >10) %>%
  mutate(indel_length = ifelse(indel_length < 60, '50~59',
                        ifelse(indel_length < 70, '60~69',
                        ifelse(indel_length < 80, '70~79',
                        ifelse(indel_length < 90, '80~89',
                        ifelse(indel_length < 100, '90~99',
  #                      ifelse(indel_length <= 125,'100~125',
                        #ifelse(indel_length <= 150,'126~150',
                        #ifelse(indel_length <= 175,'151~175',
                        #ifelse(indel_length <= 200,'176~200',
                        ifelse(indel_length <= 200,'100~200',
                        ifelse(indel_length <= 300,'201~300',
                        ifelse(indel_length <= 400,'301~400',
                        ifelse(indel_length <= 500,'401~500','501~')))))))))) %>% #head()
  group_by(indel_length) %>% #head()
  count(indel_length) %>% #head()
  mutate(indel_length = factor(indel_length,levels=c('50~59','60~69','70~79','80~89','90~99','100~200','201~300','301~400','401~500','501~'))) %>%
  #summarise(indel_length = count(indel_length))
  ggplot(aes(x=fct_rev(indel_length),y=n,fill=indel_length)) +
  geom_bar(stat='identity') +
  theme_bw() + 
  theme_step1() +
  labs(x="Range of INDEL length (bp)",y="# of INDEL") + 
  theme(legend.position = 'none') + 
  coord_flip() -> c
  

  #'50~59','60~69','70~79','80~89','90~99','100~200',201~300',301~400',401~500','501+'
  
  
  
library(png)
library(ggplot2)
library(grid)
snp <- readPNG("~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/venn.snp.dia.png")
indel <- readPNG("~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/venn.indel.dia.png")
snp


pdf(file = "~/Desktop/KCDC/pangenome/99.kogo/2024.kogo/plot_data/venn.dia.test.pdf", width = 6, height = 4)
cowplot::plot_grid(image_ggplot(snp),image_ggplot(indel),labels = c("A","B"))
dev.off()


ggplot() + annotation_custom(rasterGrob(snp, width=unit(1,"npc"), height=unit(1,"npc"))) -> p1
ggplot() + annotation_custom(rasterGrob(indel, width=unit(1,"npc"), height=unit(1,"npc"))) -> p2



cowplot::plot_grid(p1,p2,c,labels = c("A","B","C"),label_size = 20,vjust = 1.1,nrow = 1)



revio_snp %>% rbind(nanopore_snp) %>% select(ID) %>% unique() %>% filter(!(ID %in% wgs_snp$ID)) %>% dim()
revio_indel %>% rbind(nanopore_indel) %>% select(ID) %>% unique() %>% filter(!(ID %in% wgs_indel$ID)) %>% dim()
