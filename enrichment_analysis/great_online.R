#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rGREAT)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript enrichment_analysis/great_online.R <bed_file> <genome> <output_dir>")
}

bed_file <- args[1]
genome <- args[2]
output_dir <- args[3]

if (!file.exists(bed_file)) {
  stop(paste("BED file not found:", bed_file))
}

if (!(genome %in% c("hg38", "mm10"))) {
  stop("Genome must be either hg38 or mm10")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

run_name <- tools::file_path_sans_ext(basename(bed_file))

bed_df <- read.table(
  bed_file,
  sep = "\t",
  header = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

if (ncol(bed_df) < 3) {
  stop("BED file must have at least 3 columns: chrom, start, end")
}

bed_df <- bed_df[, 1:3]
colnames(bed_df) <- c("chrom", "start", "end")

bed_df$start <- as.integer(bed_df$start)
bed_df$end <- as.integer(bed_df$end)

if (any(is.na(bed_df$start)) || any(is.na(bed_df$end))) {
  stop("BED start/end columns contain non-numeric values")
}

invalid_rows <- bed_df$end <= bed_df$start

if (any(invalid_rows)) {
  warning(paste("Removing", sum(invalid_rows), "invalid rows where end <= start"))
  bed_df <- bed_df[!invalid_rows, , drop = FALSE]
}

if (nrow(bed_df) == 0) {
  stop("No valid BED intervals remain after filtering")
}

# BED is 0-based half-open; GRanges is 1-based closed
regions <- GRanges(
  seqnames = bed_df$chrom,
  ranges = IRanges(
    start = bed_df$start + 1,
    end = bed_df$end
  )
)

message("Submitting GREAT job for: ", run_name)
message("Genome: ", genome)
message("Number of regions: ", length(regions))

job <- submitGreatJob(
  regions,
  genome = genome
)

tables <- getEnrichmentTables(job)

if (length(tables) == 0) {
  stop("No enrichment tables returned by GREAT")
}

saveRDS(job, file.path(output_dir, paste0(run_name, ".great_job.rds")))

for (table_name in names(tables)) {
  safe_name <- gsub("[^A-Za-z0-9]+", "_", table_name)
  output_file <- file.path(output_dir, paste0(run_name, ".", safe_name, ".csv"))

  write.csv(
    tables[[table_name]],
    output_file,
    row.names = FALSE
  )
}

metadata_file <- file.path(output_dir, paste0(run_name, ".metadata.txt"))

writeLines(
  c(
    paste("bed_file:", bed_file),
    paste("genome:", genome),
    paste("run_name:", run_name),
    paste("n_regions:", length(regions)),
    paste("tables:", paste(names(tables), collapse = ", "))
  ),
  con = metadata_file
)

message("Finished GREAT analysis for: ", run_name)
message("Results written to: ", output_dir)