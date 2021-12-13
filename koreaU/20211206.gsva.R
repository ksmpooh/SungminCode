### 영상강의
# 20211206
#https://youtu.be/HVEMMq0H33A 

GSVA : gene set variation analysis

# fisher test 기반 gene ontology dㄱ계산
## geneprofiler

library(gProfileR2)
gs <-c("MCM2","MCM3","MCM4","MCM5","MCM6")
res.gs <- gprofiler2::gost(gs)
res.gs$result


res.gs$result %>% arrange(p_value) %>% slice(1:30) %>%
  mutate(pathway = paste(term_name,term_id,sep="\n"),
         pathway = factor(pathway,levels=rev(pathway))) %>%
  ggplot(aes(term_name,-log10(p_value))) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  theme_minimal(base_size = 7)


res.gs <- gprofiler2::gost(gs,sources = "REAC")
res.gs <- gprofiler2::gost(gs,sources = c("REAC","KEGG")

                           