#### variant check
library(tidyverse)
library(stringr)
library(ggvenn)
library(ggVennDiagram)

setwd("~/Desktop/KCDC/HLA_seq/99.Final.vcf/00.VCF/variant/")

t1_snp <- read.table("HLA.Seq.Merged_afterMapping.long_pbmm2.short_trimmed_bwamem2_sort_dedup.hg19_HLAregion.Deepvariant_Variantcalling.GLnexus_Jointcalling.updateID_setID.norm_setID.Variant_SNPS.txt")
t1_indel <- read.table("HLA.Seq.Merged_afterMapping.long_pbmm2.short_trimmed_bwamem2_sort_dedup.hg19_HLAregion.Deepvariant_Variantcalling.GLnexus_Jointcalling.updateID_setID.norm_setID.Variant_INDELS.txt")

t2_snp <- read.table("HLA.LongShort_merged.bwamem2_align.sort.Deepvariant_Variantcalling.GLnexus_Jointcalling.updateID_setID.norm_setID.Variant_SNPS.txt")
t2_indel <- read.table("HLA.Longread.Seq.Q20.Deepvariant_Variantcalling.GLnexus_Jointcalling.sampleIDCheck.updateID.setID.norm_setID.Variant_INDELS.txt")

t0_snp <- read.table("concat_sampleOrder_012.norm_setID.Variant_SNPS.txt")
t0_indel <- read.table("concat_sampleOrder_012.norm_setID.Variant_INDELS.txt")

t0_inter_snp <- read.table("sampleOrder_0002.norm_setID.Variant_SNPS.txt")
t0_inter_indel <- read.table("sampleOrder_0002.norm_setID.Variant_INDELS.txt")



t0 <- rbind(t0_snp,t0_indel)
t0_inter <- rbind(t0_inter_snp,t0_inter_indel)
t1 <- rbind(t1_snp,t1_indel)
t2 <- rbind(t2_snp,t2_indel)


ggvenn::ggvenn(
  list(Theme0 = t0$V5,
       Theme1 = t1$V5,
       Theme2 = t2$V5
  ),
  stroke_size = 0.3, set_name_size = 5,text_size = 5
)

ggvenn::ggvenn(
  list(Theme0 = t0_snp$V5,
       Theme1 = t1_snp$V5,
       Theme2 = t2_snp$V5
  ),
  stroke_size = 0.3, set_name_size = 5,text_size = 5
)


ggvenn::ggvenn(
  list(Theme0 = t0_indel$V5,
       Theme1 = t1_indel$V5,
       Theme2 = t2_indel$V5
  ),
  stroke_size = 0.3, set_name_size = 5,text_size = 5
)



ggvenn::ggvenn(
  list(Theme0_inter = t0_inter$V5,
       Theme1 = t1$V5,
       Theme2 = t2$V5
  ),
  stroke_size = 0.3, set_name_size = 5,text_size = 5
)

ggvenn::ggvenn(
  list(Theme0_inter = t0_inter_snp$V5,
       Theme1 = t1_snp$V5,
       Theme2 = t2_snp$V5
  ),
  stroke_size = 0.3, set_name_size = 5,text_size = 5
)


ggvenn::ggvenn(
  list(Theme0_inter = t0_inter_indel$V5,
       Theme1 = t1_indel$V5,
       Theme2 = t2_indel$V5
  ),
  stroke_size = 0.3, set_name_size = 5,text_size = 5
)

t0_inter_snp %>% filter(!(V5 %in% t1_snp$V5)) %>%
  filter(!(V5 %in% t2_snp$V5))
