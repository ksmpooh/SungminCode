# Rscript --vanilla miss-het.plot.R [missing.imiss] [het.het] [PDF.output] [-c call_rate] [-h het_min het_max]
# example:
# Rscript --vanilla miss-het.plot.R MISS.imiss HET.het out.pdf
# Rscript --vanilla miss-het.plot.R MISS.imiss HET.het out.pdf -c 98
# Rscript --vanilla miss-het.plot.R MISS.imiss HET.het out.pdf -h 15 20
# Rscript --vanilla miss-het.plot.R MISS.imiss HET.het out.pdf -c 98 -h 15 20

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: Rscript --vanilla miss-het.plot.R [missing.imiss] [het.het] [PDF.output] [-c call_rate] [-h het_min het_max]")
}

# required positional args
miss_file <- args[1]
het_file  <- args[2]
pdf_out   <- args[3]

# default values
call_rate <- 0.97
het_min <- 15
het_max <- 18

# parse optional flags from 4th arg onward
i <- 4
while (i <= length(args)) {
    if (args[i] == "-c") {
        if ((i + 1) > length(args)) {
            stop("Option -c requires one numeric value. Example: -c 98")
        }

        call_rate_input <- as.numeric(args[i + 1])
        if (is.na(call_rate_input)) {
            stop("call_rate must be numeric. Example: 95 or 0.95")
        }

        if (call_rate_input > 1) {
            call_rate <- call_rate_input / 100
        } else {
            call_rate <- call_rate_input
        }

        i <- i + 2

    } else if (args[i] == "-h") {
        if ((i + 2) > length(args)) {
            stop("Option -h requires two numeric values. Example: -h 15 18")
        }

        het_min_input <- as.numeric(args[i + 1])
        het_max_input <- as.numeric(args[i + 2])

        if (is.na(het_min_input) || is.na(het_max_input)) {
            stop("het_min and het_max must be numeric. Example: -h 15 18")
        }

        het_min <- het_min_input
        het_max <- het_max_input

        i <- i + 3

    } else {
        stop(paste("Unknown option:", args[i]))
    }
}

if (het_min >= het_max) {
    stop("het_min must be smaller than het_max")
}

print("Read : missing.imiss")
miss <- read.table(miss_file, header = TRUE)

print("Read : het.het")
het <- read.table(het_file, header = TRUE)

# missing threshold
miss_cutoff <- 1 - call_rate

# pdf filename without .pdf
pdf_base <- sub("\\.pdf$", "", basename(pdf_out), ignore.case = TRUE)

# labels for title
call_rate_label <- sprintf("%.0f", call_rate * 100)

plot_title <- paste0(
    "Missing vs heterozygosity : ", pdf_base,
    " | CR>=", call_rate_label, "% | HET=", het_min, "-", het_max
)

miss <- cbind(miss, CR = ((1 - miss$F_MISS) * 100))
het  <- cbind(het,  HET = ((het$N.NM. - het$O.HOM.) / het$N.NM.) * 100)

lowSample <- merge(miss, het, by = "FID")

fail_idx <- lowSample$F_MISS > miss_cutoff | lowSample$HET < het_min | lowSample$HET > het_max

pdf(pdf_out, height = 7, width = 10)
plot(lowSample$HET, lowSample$F_MISS,
     xlab = "heterozygosity rate",
     ylab = "missing rate",
     main = plot_title,
     col  = rgb(0, 0, 1, 0.7),
     cex  = 1,
     pch  = 16)

abline(v = het_min, col = rgb(1, 0, 0, 1), lty = 3, lwd = 2)
abline(v = het_max, col = rgb(1, 0, 0, 1), lty = 3, lwd = 2)
abline(h = miss_cutoff, col = rgb(1, 0, 0, 1), lty = 3, lwd = 2)

points(lowSample$HET[fail_idx],
       lowSample$F_MISS[fail_idx],
       col = rgb(1, 0, 0, 0.7),
       cex = 1,
       pch = 16)

dev.off()

rmList <- lowSample[fail_idx, ]

write.table(rmList[, c(1, 2)],
            "rmLQSamples.txt",
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)