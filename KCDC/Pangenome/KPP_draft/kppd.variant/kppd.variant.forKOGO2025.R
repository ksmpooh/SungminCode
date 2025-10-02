library(tidyverse)

hprc <- read_table("main/KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.hprcV1.filTag.filtered.maincontig.tab")
kpp <- read_table("main/KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.kpp132.filTag.filtered.maincontig.tab")
colnames(hprc) <- c("CHROM","POS","TYPE","REF","ALT","AC","AN","NS","AF","MAF","F_MISSING")
colnames(kpp) <- c("CHROM","POS","TYPE","REF","ALT","AC","AN","NS","AF","MAF","F_MISSING")
header <- c("CHROM","POS","node","REF","ALT","FILTER","AC","AN","NS","AF","MAF","F_MISSING")

colnames(hprc) <- header
colnames(kpp) <- header

hprc_snp <- read_table("KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.hprcV1.filTag.filtered.maincontig.multi.snps.filtered.tab")
kpp_snp <- read_table("KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.kpp132.filTag.filtered.maincontig.multi.snps.filtered.tab")

colnames(hprc_snp) <- header
colnames(kpp_snp) <- header



hprc_indel <- read_table("KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.hprcV1.filTag.filtered.maincontig.multi.indels.filtered.tab")
kpp_indel <- read_table("KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.kpp132.filTag.filtered.maincontig.multi.indels.filtered.tab")

colnames(hprc_indel) <- header
colnames(kpp_indel) <- header

#KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.hprcV1.filTag.filtered.maincontig.multi.indels.filtered.tab
#KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.kpp132.filTag.filtered.maincontig.multi.indels.filtered.tab
hprc_snp %>% select(-13)%>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT)) %>% select(ID) %>% unique() -> hprc_snp_ID
kpp_snp %>% select(-13)%>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT)) %>% select(ID) %>% unique() -> kpp_snp_ID

hprc_indel %>% select(-13)%>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT)) %>% select(ID) %>% unique() -> hprc_indel_ID
kpp_indel %>% select(-13) %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT)) %>% select(ID) %>% unique() -> kpp_indel_ID


kpp_snp %>% select(-13)%>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT)) %>% select(ID) %>% unique() -> kpp_snp_ID
F_MISSING
kpp_snp %>% select(-13) %>% filter(F_MISSING < 0.01) %>% select(CHROM,POS) %>% dim()
kpp_indel %>% select(-13) %>% filter(F_MISSING < 0.01) %>% select(CHROM,POS) %>% dim()


kpp_only_snp <- read_table("/BDATA/smkim/pangenome/kppd.variant.check/korean.specific/korean.specifit.snp.F_missingUnder0.01.txt",col_names = header)
kpp_only_indel <- read_table("/BDATA/smkim/pangenome/kppd.variant.check/korean.specific/korean.specifit.indel.F_missingUnder0.01.txt",col_names = header)

kpp_only_snp %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT)) %>% filter(!(ID %in% hprc_snp_ID$ID)) -> kpp_only_snp_withoutHPRC
kpp_only_indel %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT)) %>% filter(!(ID %in% hprc_indel_ID$ID)) -> kpp_only_indel_withoutHPRC


colnames(hprc_indel) <- header
colnames(kpp_indel) <- header

hprc_only_snp <- read_table("/BDATA/smkim/pangenome/kppd.variant.check/korean.specific/KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.hprcV1.filTag.filtered.maincontig.multi.snps.filtered.extracted.tab",col_names = header) %>% filter(CHROM != "#")
hprc_only_indel <- read_table("/BDATA/smkim/pangenome/kppd.variant.check/korean.specific/KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.hprcV1.filTag.filtered.maincontig.multi.indels.filtered.extracted.tab",col_names = header) %>% filter(CHROM != "#")

#hprc_only_snp %>% filter

kpp_only_snp %>% mutate(MAF_bin = cut(MAF, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE, right = FALSE)) %>% count(MAF_bin)
kpp_only_indel %>% mutate(AF_bin = cut(AF, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE, right = FALSE)) %>% count(AF_bin) %>%
  write.table("/BDATA/smkim/pangenome/kppd.variant.check/korean.specific/korean.specific.indel.bin.txt",col.names = T,row.names = F,quote = F,sep = "\t")

kpp_only_snp %>% mutate(AF_bin = cut(AF, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE, right = FALSE)) %>% count(AF_bin) %>%
  write.table("/BDATA/smkim/pangenome/kppd.variant.check/korean.specific/korean.specific.snp.bin.txt",col.names = T,row.names = F,quote = F,sep = "\t")


kpp_only_snp %>% mutate(MAF_bin = cut(MAF, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE, right = FALSE)) %>% count(MAF_bin)
kpp_only_indel_withoutHPRC %>% mutate(AF_bin = cut(AF, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE, right = FALSE)) %>% count(AF_bin) %>%
  write.table("/BDATA/smkim/pangenome/kppd.variant.check/korean.specific/korean.specific.indel.withoutHPRC.bin.txt",col.names = T,row.names = F,quote = F,sep = "\t")

kpp_only_snp_withoutHPRC %>% mutate(AF_bin = cut(AF, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE, right = FALSE)) %>% count(AF_bin) %>%
  write.table("/BDATA/smkim/pangenome/kppd.variant.check/korean.specific/korean.specific.snp.withoutHPRC.bin.txt",col.names = T,row.names = F,quote = F,sep = "\t")


234863/13219042
282039/5281700

hprc_mrg <- read_table("/BDATA/smkim/pangenome/kppd.variant.check/korean.specific/Hprc_only_variants_no_cut.annotated.rmNA.tsv")
kpp_mrg <- read_table("/BDATA/smkim/pangenome/kppd.variant.check/korean.specific/KPP_only_variants_no_cut.annotated.rmNA.tsv")
hprc_mrg %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT)) %>% filter(ID %in% hprc_snp_ID$ID) -> hprc_mrg_snp
kpp_mrg %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT)) %>% filter(ID %in% kpp_snp_ID$ID) -> kpp_mrg_snp
#kpp_mrg %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT)) %>% filter(ID %in% kpp_indel_ID$ID) -> kpp_mrg_indel


hprc_mrg_snp %>% filter(F_MISSING < 0.5) %>% count(Full_gene) %>% dim()
hprc_mrg_snp %>% count(Full_gene)

hprc_mrg_snp %>% count(F)

kpp_mrg_snp %>% filter(F_MISSING < 0.5) %>% count(Full_gene) %>% dim()
kpp_mrg_snp %>% filter(F_MISSING < 0.5) %>% filter(!(ID %in% hprc_snp_ID$ID)) %>% count(Full_gene) %>% dim()
kpp_mrg_snp %>% filter(F_MISSING < 0.5) %>% filter(!(ID %in% hprc_snp_ID$ID)) %>% count(Full_gene) %>% arrange(-n) -> kpp_mrg_top10
kpp_mrg_top10 <- kpp_mrg_top10[1:10,]
kpp_mrg_snp %>% filter(!(ID %in% hprc_snp_ID$ID)) %>% 
  mutate(AF_bin = cut(AF, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE, right = FALSE)) %>% count(Full_gene,AF_bin) %>% #dim()
  filter(Full_gene %in% kpp_mrg_top10$Full_gene) -> kpp_mrg_snp_top10gene_info


hprc_mrg_snp %>% select(ID,AF,Full_gene) %>% rename(AF_HPRC = AF) %>%
  full_join(kpp_mrg_snp %>% select(ID,AF,Full_gene) %>% rename(AF_KPP = AF)) %>%
  filter(Full_gene %in% kpp_mrg_top10$Full_gene) %>%
  group_by(Full_gene) %>% summarise(AF_KPP = mean(AF_KPP),AF_HPRC = mean(AF_HPRC)) -> a
kpp_mrg_snp %>% select(ID,AF,Full_gene) %>% rename(AF_KPP = AF)


hprc_mrg_snp %>% filter(Full_gene %in% kpp_mrg_top10$Full_gene)


kpp_mrg_snp %>% select(ID,AF,Full_gene) %>% rename(AF_KPP = AF) %>%
  left_join(hprc_mrg_snp %>% select(ID,AF,Full_gene) %>% rename(AF_HPRC = AF)) %>% 
  filter(Full_gene %in% kpp_mrg_top10$Full_gene) -> a

9741037 - 9463295
write.table(kpp_mrg_snp_top10gene_info,"/BDATA/smkim/pangenome/kppd.variant.check/korean.specific/kpp.mrg.top10.txt",col.names = T,row.names = F,quote = F,sep = "\t")


kpp_mrg_snp_top10gene_info <- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/korean.variant/kpp.mrg.top10.txt")

head(kpp_mrg_snp_top10gene_info)
dim(kpp_mrg_snp_top10gene_info)

kpp_mrg_snp_top10gene_info %>% 
  ggplot(aes(x=Full_gene,y=n,fill=AF_bin)) + 
  geom_col(position = "dodge")


kpp_mrg_snp %>% filter(F_MISSING < 0.5) %>% filter(!(ID %in% hprc_snp_ID$ID)) %>% count(Full_gene) %>% arrange(-n)
kpp_mrg_snp %>% mutate()
kpp_mrg_snp %>% count(Full_gene) %>% dim()
kpp_mrg_snp %>% filter(!(ID %in% hprc_snp_ID$ID)) %>% count(Full_gene) %>% dim()
mrg <- read_table("/BDATA/smkim/pangenome/kppd.variant.check/korean.specific/mrg.list",col_names = F)
colnames(mrg) <- c("CHROM","POS","node","REF","ALT","FILTER","AC","AN","NS","AF","MAF","F_MISSING","Gene","Gene2")

table(kpp_snp_ID$ID %in% hprc_snp_ID$ID)
table(kpp_indel_ID$ID %in% hprc_indel_ID$ID)
table(hprc_indel_ID$ID %in% kpp_indel_ID$ID)


## CPC
header <- c("CHROM","POS","node","REF","ALT","FILTER","AC","AN","NS","AF","MAF","F_MISSING")
cpc_snp <- read_table("/BDATA/smkim/pangenome/kppd.variant.check/CPC/CPC.Phase1.CHM13v2.contig_change.maincontig.liftoverhsToGRCh38.noGRCh38.filTag.norm.refil.filtered.snps.tab",col_names = header) %>% filter(CHROM != "#")
cpc_indel <- read_table("/BDATA/smkim/pangenome/kppd.variant.check/CPC/CPC.Phase1.CHM13v2.contig_change.maincontig.liftoverhsToGRCh38.noGRCh38.filTag.norm.refil.filtered.indels.tab",col_names = header) %>% filter(CHROM != "#")

cpc_snp %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT))->cpc_snp
cpc_indel %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT))->cpc_indel


### 20250804

hprc_annotate_snp <- read_table("/BDATA/smkim/pangenome/kppd.variant.check/Medical/KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.hprcV1.filTag.filtered.maincontig.multi.snps.filtered.annotated.tsv")  %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT))
kpp_annotate_snp <- read_table("/BDATA/smkim/pangenome/kppd.variant.check/Medical/KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.kpp132.filTag.filtered.maincontig.multi.snps.filtered.annotated.tsv")  %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT))
hprc_annotate_indel <- read_table("/BDATA/smkim/pangenome/kppd.variant.check/Medical/KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.hprcV1.filTag.filtered.maincontig.multi.indels.filtered.annotated.tsv")  %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT))
kpp_annotate_indel <- read_table("/BDATA/smkim/pangenome/kppd.variant.check/Medical/KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.kpp132.filTag.filtered.maincontig.multi.indels.filtered.annotated.tsv")  %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT))


#hprc_annotate_snp < read.table("/BDATA/smkim/pangenome/kppd.variant.check/Medical/KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.hprcV1.filTag.filtered.maincontig.multi.snps.filtered.annotated.tsv",header = T,sep = "\t")  %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT))
#kpp_annotate_snp <- read.table("/BDATA/smkim/pangenome/kppd.variant.check/Medical/KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.kpp132.filTag.filtered.maincontig.multi.snps.filtered.annotated.tsv",header = T)  %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT))
#hprc_annotate_indel <- read.table("/BDATA/smkim/pangenome/kppd.variant.check/Medical/KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.hprcV1.filTag.filtered.maincontig.multi.indels.filtered.annotated.tsv",header = T)  %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT))
#kpp_annotate_indel <- read.table("/BDATA/smkim/pangenome/kppd.variant.check/Medical/KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.kpp132.filTag.filtered.maincontig.multi.indels.filtered.annotated.tsv",header = T)  %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT))

hprc_annotate_snp <- read_table("/BDATA/smkim/pangenome/kppd.variant.check/Medical/noNA/KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.hprcV1.filTag.filtered.maincontig.multi.snps.filtered.annotated.tsv")  %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT))
kpp_annotate_snp <- read_table("/BDATA/smkim/pangenome/kppd.variant.check/Medical/noNA/KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.kpp132.filTag.filtered.maincontig.multi.snps.filtered.annotated.tsv")  %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT))
hprc_annotate_indel <- read_table("/BDATA/smkim/pangenome/kppd.variant.check/Medical/noNA/KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.hprcV1.filTag.filtered.maincontig.multi.indels.filtered.annotated.tsv")  %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT))
kpp_annotate_indel <- read_table("/BDATA/smkim/pangenome/kppd.variant.check/Medical/noNA/KPPD132.GRCh38.noAT.norm.wave.ref.sort.dedupPosOpt.merge.fix.kpp132.filTag.filtered.maincontig.multi.indels.filtered.annotated.tsv")  %>% mutate(ID = paste0(CHROM,":",POS,":",REF,":",ALT))


kpp_annotate_snp %>% select(-13) %>%  na.omit()  %>% #filter(F_MISSING < 0.01) %>%
mutate(onlyKPP = ifelse(ID %in% kpp_mrg_snp$ID, 1,0)) %>% #count(onlyKPP)
mutate(noHPRC = ifelse(ID %in% hprc_annotate_snp$ID,0,1)) %>% count(onlyKPP,noHPRC)

kpp_annotate_indel %>% select(-13) %>%  na.omit()  %>% #filter(F_MISSING < 0.01) %>%
  mutate(onlyKPP = ifelse(ID %in% kpp_mrg_indel$ID, 1,0)) %>% #count(onlyKPP)
  mutate(noHPRC = ifelse(ID %in% hprc_annotate_indel$ID,0,1)) %>% count(onlyKPP,noHPRC)


kpp_mrg_snp  
hprc_snp
#hprc_annotate_snp %>% 

#kpp_annotate_snp %>% 
#snp: 2421547
#indel: 1089707


  
  
  
  




head(hprc_snp)
head(kpp_snp)

head(kpp_only_snp) # F_MISSING < 0.01
kpp_annotate_snp %>% filter(Full_gene != "NA") -> kpp_annotate_snp_noNA
kpp_annotate_indel %>% filter(Full_gene != "NA") -> kpp_annotate_indel_noNA



#hprc_mrg <- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/korean.variant/hprc.mrg.list",col_names = F) %>% na.omit()
#kpp_mrg <- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/korean.variant/kpp.mrg.list",col_names = F) %>% na.omit()
head(hprc_mrg)
head(kpp_mrg)
max(kpp_mrg$X1)

table(hprc_mrg$X2 %in% kpp_mrg$X2)
table(kpp_mrg$X2 %in% hprc_mrg$X2)

top10 <- kpp_mrg[order(-kpp_mrg$X1), ][1:10, ]
top10
colnames(hprc_mrg) <- c("")
colnames(kpp_mrg) <- c("")
kpp_mrg$type <- "KPPD"
hprc_mrg$type <- "HPRC"

kpp_mrg %>% rbind(hprc_mrg) %>% filter(X2 %in% top10$X2) -> a
a_order <- a %>%
  filter(type == "KPPD") %>%
  arrange(desc(X1)) %>%
  pull(X2)

a_order
a$X2 <- factor(a$X2, levels = a_order)

ggplot(a, aes(x = X2, y = X1, fill = type)) +
  geom_col(position = "dodge") +
  labs(#title = "Top 10 Gene",
       x = "Gene",
       y = "# of Variant") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
kpp_mrg %>% rbind(hprc_mrg) %>% filter(X2 %in% top10$X2) %>%
ggplot(aes(x = X2, y = X1, fill = type)) +
  geom_col(position = "dodge") +
  labs(#title = "Top 10 Gene",
    x = "Gene",
    y = "# of Variant") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))


kpp_mrg %>% rbind(hprc_mrg) %>% filter(X2 %in% top10$X2)
kpp_mrg %>% rbind(hprc_mrg) %>% pivot_wider(names_from = type,values_from = X1)

kpp_mrg %>% filter(X2 %in% top10$X2)

kpp_snp_bin <- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/korean.variant/korean.specific.snp.bin.txt")
kpp_indel_bin <- read_table("/Users/ksmpooh/Desktop/KCDC/pangenome/KPPD/korean.variant/korean.specific.indel.bin.txt")


kpp_snp_bin <- factor(kpp_snp_bin$AF_bin, levels = kpp_snp_bin$AF_bin)
ggplot(kpp_snp_bin, aes(x = AF_bin, y = n)) +
  geom_col(fill = "skyblue") +
  labs(title = "MAF ºóµµ ±¸°£º° Count",
       x = "AF Bin",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

head(kpp_snp_bin)
kpp_snp_bin %>% summarise(sum(n))

kpp_indel_bin %>% summarise(sum(n))
