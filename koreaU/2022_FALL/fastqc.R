library(fastqcr)
setwd("~/Desktop/KU/2022_FALL/genetics/fastqc_report/")
qc.file <- list.files("~/Desktop/KU/2022_FALL/genetics/fastqc_report/","*zip")
qc.file <- system.file("fastqc_results", "S1_fastqc.zip",  package = "fastqcr")
qc.file
# [1] "/Library/Frameworks/R.framework/Versions/3.3/Resources/library/fastqcr/fastqc_results/S1_fastqc.zip"
#We start by reading the output using the function qc_read(), which returns a list of tibbles containing the data for specified modules:
  
  # Read all modules
qc <- qc_read(qc.file)
# Elements contained in the qc o

qc <- qc_read(qc.file[1])

qc_plot(qc,"Per sequence GC content")



setwd("~/Desktop/KU/2022_FALL/genetics/fastqc_report/")
qc.files <- list.files("~/Desktop/KU/2022_FALL/genetics/fastqc_report/","*zip")
qc.files
library(stringr)
qc.files
qc.names <- str_replace_all(qc.files,"_fastqc.zip","")
for (i in c(1,2,5,6)) {
  qc.names[i] <- str_replace(qc.names[i],"SRR","CASE_SRR")
}

for (i in c(3,4,7,8)) {
  qc.names[i] <- str_replace(qc.names[i],"SRR","CONTROL_SRR")
}

qc.names <- c("CT1_1","CT1_2","C1_1","C1_2","CT2_1","CT2_2","C2_1","C2_2")

par(mfrow=c(1,3))
qc.files
qc.names
# read all modules in all files
qc <- qc_read_collection(qc.files, sample_names = qc.names)#> Reading: /Users/kassambara/Documents/R/MyPackages/fastqcr/inst/fastqc_results/S1_fastqc.zip#> Reading: /Users/kassambara/Documents/R/MyPackages/fastqcr/inst/fastqc_results/S2_fastqc.zip#> Reading: /Users/kassambara/Documents/R/MyPackages/fastqcr/inst/fastqc_results/S3_fastqc.zip#> Reading: /Users/kassambara/Documents/R/MyPackages/fastqcr/inst/fastqc_results/S4_fastqc.zip#> Reading: /Users/kassambara/Documents/R/MyPackages/fastqcr/inst/fastqc_results/S5_fastqc.zip
# Plot per sequence GC content
qc_plot_collection(qc, "Per sequence GC content")
# Per base sequence quality
qc_plot_collection(qc, "Per base sequence quality")
# Per sequence quality scores
qc_plot_collection(qc, "Per sequence quality scores")
# Per base sequence content
qc_plot_collection(qc, "Per base sequence content")
# Sequence duplication levels
qc_plot_collection(qc, "Sequence duplication levels")

Adapter Content


head(qc)
qc_plot_collection(qc, "Adapter Content")
qc_plot_collection(qc, "Overrepresented sequences")

qc_plot_collection(qc, "Summary")
qc_plot_collection(qc, "Basic Statistics")
"Summary"
qc$summary
qc$basic_statistics
