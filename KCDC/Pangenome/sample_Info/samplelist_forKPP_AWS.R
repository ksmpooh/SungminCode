df <- data.frame(
  raw_ID = sprintf("PG24%04d", c(1:55, 57:73, 75:202)),
  stringsAsFactors = FALSE
)

set.seed(123)  # 재현성 필요 없으면 삭제

n <- nrow(df)

df$AWS_ID <- paste0(
  "KPP",
  sprintf("%04d", sample(0:9999, n, replace = FALSE))
)

head(df)

#df %>% writexl::write_xlsx("~/Downloads/Pangenome.opendata.IDmatching.20260114.xlsx")


df$AWS_ID %>% unique() %>% length()

####
