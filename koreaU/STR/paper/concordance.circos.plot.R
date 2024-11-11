# 설치
#install.packages("circlize")
library(tidyverse)
# 패키지 로드
library(circlize)
ref <- read.table("~/Desktop/KU/@research/STR/db/annovar/output_bed_ref")
anno <- read_table("~/Desktop/KU/@research/STR/db/annovar/Raw.anno.merge.processing.onlyneed.txt")
concordance_rate <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/figure/figure2/f2.concordance_bySTR_simpleSTR.txt")

final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR

final_ref %>% select(MOTIFS,ID) %>% rename(STR_ID = ID) %>% 
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> final_ref_pro


head(ref)
ref %>% unique()-> ref
colnames(ref) <- c("chrom","start","end","STR_ID")
head(a)
head(concordance_rate)
concordance_rate %>% left_join(anno) %>% left_join(ref) -> a
head(a)






data<-a[1:100000,] %>% filter(concordance_rate != 1)
data<-a

#head(data,20)

# Circos 플롯 시작
head(data)
head(ref)
#data %>% left_join(ref) -> data
data %>% filter(chrom %in% c(paste0("chr", c(21,22,15,9)))) -> data
head(data)
dim(data)

# cytoband
# 빨간: Cent
# 푸른색: stalk
cyto <- read_table("~/Desktop/KU/@research/STR/db/annovar/hg38_cytoBand.txt",col_names = F)
head(cyto)
colnames(cyto) <- c("chrom","start","end","cytoband","anno")
cyto %>% mutate(cytoband = paste0(str_replace_all(chrom,"chr",""),cytoband)) -> cyto
head(cyto)
cbind(cyto$start, cyto$end)
# Circos 플롯 초기화
# Circos 플롯 초기화
#circos.initializeWithIdeogram(chromosome.index = c(paste0("chr", c(1:22)),'chrX'),species = "hg38")
# chr21, 22 chr 15 chr9
# Load circlize
library(circlize)

# Initialize with cytobands (e.g., for human chromosomes)
#circos.initializeWithIdeogram(species = "hg19")

# Create a heatmap for genomic data, e.g., mutation frequency

# Finalize plot
circos.clear()

data %>% select(chrom,start,end,concordance_rate) %>% rename(value = concordance_rate)

head(data)
head(DMR_hyper)
data <- a
#circos.genomicHeatmap(data %>% select(chrom,start,end,concordance_rate) %>% rename(value = concordance_rate), col = colorRamp2(c(0, 1), c("blue", "red")))

#### goot
concordance_rate %>% left_join(anno) %>% left_join(ref) %>% left_join(final_ref_pro) -> b
#fill = c("grey", "#F8766D", "#619CFF"),

circos.clear()
circos.initializeWithIdeogram(chromosome.index = c(paste0("chr", c(1:22)),"chrX"),species = "hg38")
bed_list = list(b_cyto_bychr %>% filter(concordance_rate == 1) %>% select(chrom,start,end,concordance_rate) %>% arrange(chrom, start) %>% as.data.frame(),
                b_cyto_bychr %>% filter(concordance_rate != 1) %>% select(chrom,start,end,concordance_rate) %>% arrange(chrom, start) %>% as.data.frame())
head(bed_list)
circos.genomicTrack(bed_list, 
                    panel.fun = function(region, concordance_rate, ...) {
                      i = getI(...)
                      circos.genomicLines(region, concordance_rate, col = i, ...)
                    })
circos.genomicDensity(b %>% select(chrom,start,end,concordance_rate) %>% filter(concordance_rate != 1),col = c("#F8766D"), track.height = 0.1)
circos.genomicDensity(b %>% select(chrom,start,end,concordance_rate) %>% filter(concordance_rate == 1),col = c("#619CFF"), track.height = 0.1)

##### Good



concordance_rate %>% left_join(anno) %>% left_join(ref) %>% group_by(chrom,cytoband) %>% filter(concordance_rate != 1) %>%
  summarise(concordance_rate = mean(concordance_rate)) %>% left_join(cyto) -> b

head(bed_list)
circos.clear()
circos.par("track.height" = 0.1, 
           "cell.padding" = c(0, 0, 0, 0), 
           "track.margin" = c(0.01, 0.01),
           #"ideogram_text_cex" = 0.8,  # cytoband 레이블 크기
           "fontsize" = 10)  # 염색체 이름 크기
circos.initializeWithIdeogram(chromosome.index = c(paste0("chr", c(1:22)),"chrX"),species = "hg38")
circos.genomicTrack(bed_list, 
                    panel.fun = function(region, concordance_rate, ...) {
                      i = getI(...)
                      circos.genomicLines(region, concordance_rate, col = i, ...)
                    })
circos.genomicDensity(data %>% select(chrom,start,end,concordance_rate) %>% filter(concordance_rate != 1),col = c("#FF000080"), track.height = 0.1)
circos.genomicDensity(data %>% select(chrom,start,end,concordance_rate) %>% filter(concordance_rate == 1),col = c("blue"), track.height = 0.1)
#bed_list = list(data %>% filter(concordance_rate == 1) %>% select(chrom,start,end,concordance_rate) %>% as.data.frame(),
#                data %>% filter(concordance_rate != 1) %>% select(chrom,start,end,concordance_rate) %>% as.data.frame())



######


#######
concordance_rate %>% left_join(anno) %>% left_join(ref) %>% left_join(final_ref_pro) -> b
head(b)
b %>% filter(concordance_rate != 1) %>% #count(type)
  group_by(RU.length,MOTIFS,GC,main_cyto_type,type) %>%  
  summarise(concordance_rate = mean(concordance_rate)) -> b
  #mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS))
head(b)
b %>% filter(concordance_rate != 1) %>% #count(type)
  group_by(RU.length,MOTIFS,GC,main_cyto_type) %>%  
  summarise(concordance_rate = mean(concordance_rate)) %>% ungroup()-> b_cyto

b_cyto %>% 
  pivot_wider(names_from = main_cyto_type,values_from = concordance_rate) %>%
  mutate_all(~replace(., is.na(.), -1)) -> b_cyto

head(cyto)
head(b)
b %>% select(concordance_rate,chrom,cytoband,main_cyto_type) %>%mutate(concordance_range = ifelse(concordance_rate != 1,0,1)) %>% #left_join(cyto)
  group_by(concordance_range,chrom,cytoband,main_cyto_type) %>%  #head()
  summarise(concordance_rate = mean(concordance_rate)) %>% left_join(cyto) %>% 
  ungroup()-> b_cyto_bychr
head(b_cyto_bychr)

#b_cyto_bychr %>% mutate_all(~replace(., is.na(.), -1)) -> b_cyto_bychr

head(b_cyto_bychr)

#b %>% mutate(concordance_range = ifelse(concordance_rate != 1,0,1)) %>% write.table("~/Desktop/KU/@research/STR/figure/sup.figure/STR.concordance.INFO.withanno.simpleSTR.txt",col.names = T,row.names = F,quote = F,sep = "\t")
b %>% mutate(concordance_range = ifelse(concordance_rate != 1,0,1)) %>% 
  count(chrom,concordance_range) %>% group_by(chrom) %>% mutate(prop = prop.table(n)) %>% #head()
  ggplot(aes(x=chrom,y=n,fill=concordance_range)) + 
  geom_bar(position = 'fill',stat = 'identity')

#•	PAR1: X (1-2,699,520 bp), Y (1-2,699,520 bp)
#•	PAR2: X (154,931,044-156,030,895 bp), Y (56,887,902-57,886,905 bp)
#•	비의사상염색체 구역: X (2,699,521-154,931,043 bp)
b %>% mutate(concordance_range = ifelse(concordance_rate != 1,0,1)) %>% filter(chrom == "chrX") %>% #head()
  mutate(chrXtype = case_when(
    end <= 2699520 ~ "PAR1",
    start >= 154931044 & end <= 156030895 ~ "PAR2",
    start >= 2699521 & end <= 154931043 ~ "nonPAR",
    TRUE ~ "other")) -> b_chrX

circos.clear()
b_chrX



# good 
b_chrX %>% mutate(par = ifelse(chrXtype == "nonPAR",chrXtype,"PAR")) %>% #count(par)
  count(par,concordance_range) %>%
  ggplot(aes(x=par,y=n,fill=as.factor(concordance_range))) +
  geom_bar(stat = "identity",position = "fill")


# 고려
b_chrX %>% mutate(par = ifelse(chrXtype == "nonPAR",chrXtype,"PAR")) %>% #count(par)
  count(par,type,concordance_range) %>% #head()
  ggplot(aes(x=type,y=n,fill=concordance_range)) +
  geom_bar(stat = "identity",position = "fill") + 
  facet_grid(~par,scales = "free")

head(b)

b %>% mutate(concordance_range = ifelse(concordance_rate != 1,0,1)) %>% 
  count(type,concordance_range) %>% 
  ggplot(aes(x=type,y=n,fill=concordance_range)) +
  geom_bar(stat = "identity",position = "fill")

circos.par("canvas.xlim" = c(0, 1), "canvas.ylim" = c(0, 1),
           "start.degree" = 90, "gap.after" = 270)
# Initialize the ideogram for 

circos.initializeWithIdeogram(chromosome.index = c("chrX"),species = "hg38")
par1_start <- 1;par1_end <- 2699520;nonpar_start <- 2699521;nonpar_end <- 154931043;par2_start <- 154931044;par2_end <- 156030895

# Create a track to show the PAR and non-PAR regions as colored bars
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(par1_start, 0, par1_end, 1, col = "blue", border = "blue")
  circos.rect(nonpar_start, 0, nonpar_end, 1, col = "red", border = "red")
  circos.rect(par2_start,0, par2_end, 1, col = "green", border = "green")
  #circos.text((par1_start + par1_end) / 2, 0.5, "PAR1", cex = 1, col = "black", facing = "inside")
  circos.text(0, -1, "PAR1", cex = 1, col = "black", facing = "inside")
  circos.text((nonpar_start + nonpar_end) / 2, 0.5, "non-PAR", cex = 1, col = "black", facing = "inside")
  circos.text((par2_start + par2_end) / 2, -3, "PAR2", cex = 1, col = "black", facing = "clockwise")
}, track.height = 0.05)
circos.genomicRainfall(b_chrX %>% select(chrom,start,end,GC), pch = 16, cex = 0.3, col = c(rgb(1,0,0,0.2), rgb(0,0,1,0.2)),track.height = 0.1)
#circos.genomicRainfall(b_chrX %>% select(chrom,start,end,GC), pch = 16, cex = 0.3, col = c(rgb(1,0,0,0.2), rgb(0,0,1,0.2)),track.height = 0.1)
circos.track(b_chrX %>% select(chrom,start,end,GC), panel.fun = function(region, GC,...) {
  # Extract the start and GC content for the points
  start <- b_chrX$start[CELL_META$sector.numeric.index]
  GC <- b_chrX$GC[CELL_META$sector.numeric.index]
  
  # Plot points at the start position with their respective GC values
  circos.points(start, GC,...)
}, ylim = c(min(b_chrX$GC), max(b_chrX$GC)))

# Add other data tracks or custom visualizations if needed

# Clear the Circos plot after drawing
circos.clear()
  
  
##
circos.clear()
circos.initializeWithIdeogram(chromosome.index = c(paste0("chr", c(1:22)),"chrX"),species = "hg38")
circos.genomicTrack(b_chrX, 
                    panel.fun = function(region, GC, ...) {
                      # numeric.column is automatically passed to `circos.genomicPoints()`
                      circos.genomicPoints(region, GC, ...)
                    })
#circos.genomicHeatmap(b %>% select(chrom,start,end,concordance_rate), col = colorRamp2(c(0, 1), c("blue","red")), 
#side = "inside", border = "white")
head(b)
head(b_cyto)
#bed1 = generateRandomBed(nr = 500)
#bed2 = generateRandomBed(nr = 500)
#bed_list_t = list(bed1, bed2)
#head(bed_list_t)
#circos.genomicRainfall(bed_list, pch = 16, cex = 0.3, col = c(rgb(1,0,0,0.2), rgb(0,0,1,0.2)),track.height = 0.1)
head(b)
head(b_cyto_bychr)
concordance_rate %>% left_join(anno) %>% left_join(ref) %>% left_join(final_ref_pro) -> b

circos.genomicDensity(b %>% select(chrom,start,end,concordance_rate) %>% filter(concordance_rate != 1),col = c("#FF000080"), track.height = 0.1)
circos.genomicDensity(b %>% select(chrom,start,end,concordance_rate) %>% filter(concordance_rate == 1),col = c("blue"), track.height = 0.1)
bed_list = list(b_cyto_bychr %>% filter(concordance_rate == 1) %>% select(chrom,start,end,concordance_rate) %>% arrange(chrom, start) %>% as.data.frame(),
                b_cyto_bychr %>% filter(concordance_rate != 1) %>% select(chrom,start,end,concordance_rate) %>% arrange(chrom, start) %>% as.data.frame())
head(bed_list)
circos.genomicTrack(bed_list, 
                    panel.fun = function(region, concordance_rate, ...) {
                      i = getI(...)
                      circos.genomicLines(region, concordance_rate, col = i, ...)
                    })
circos.genomicDensity(b %>% select(chrom,start,end,concordance_rate) %>% filter(concordance_rate != 1),col = c("#FF000080"), track.height = 0.1)
circos.genomicDensity(b %>% select(chrom,start,end,concordance_rate) %>% filter(concordance_rate == 1),col = c("blue"), track.height = 0.1)


bed = generateRandomBed(nr = 200,nc = 2)
head(bed)


##

head(b_cyto)
circos.clear()

col_meth = colorRamp2(c(0, 0.5, 1), c("white", "blue", "black"))
circos.heatmap(b_cyto %>% select(acen:stalk), split = as.factor(b_cyto$RU.length), col = col_meth, 
               bg.lwd = 1, bg.lty = 1,
               track.height = 0.12,show.sector.labels = TRUE)

library(RColorBrewer)
#unique(b$type)
# Initialize Circos plot with split by RU.length
circos.initialize(factors = factor(b_cyto$RU.length), xlim = c(0, 5))
#circos.initialize()
# Create a list of columns that will be used for the heatmap (acen to stalk)
heatmap_columns <- data.frame(acen = b_cyto$acen, gneg = b_cyto$gneg, gpos = b_cyto$gpos, gvar = b_cyto$gvar, stalk = b_cyto$stalk)
head(heatmap_columns)
# Define color palette for heatmap values
col_fun <- colorRamp2(c(-1,0, 1), c("white","blue", "red"))
b_cyto$RU.length
length(b_cyto$RU.length)
nrow(heatmap_columns)
# Draw the heatmap for acen to stalk values
circos.heatmap(as.matrix(heatmap_columns), split = b_cyto$RU.length, col = col_fun,show.sector.labels = TRUE)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 20) { # the last sector
    cn = colnames(heatmap_columns)
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                1:n - 0.5, cn, 
                cex = 0.5, adj = c(0, 0.5), facing = "inside")
  }
}, bg.border = NA)
circos.clear()
#circos.heatmap()
# Add GC content as a separate track (bar plot)
circos.trackPlotRegion(factors = data$RU.length, ylim = c(0, 1), panel.fun = function(x, y) {
  circos.barplot(data$GC[CELL_META$sector.numeric.index], col = "green", border = "black")
}, track.height = 0.05)

# Finalize and clear the plot
circos.clear()



