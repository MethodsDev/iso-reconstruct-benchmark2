#!/usr/bin/env Rscript

# Load required libraries
library(Isosceles)
library(stringr)
library(dplyr)
library(tibble)
library(optparse)

# Define command line options
option_list <- list(
  make_option(c("-b", "--bam_file"), type = "character", default = NULL, help = "Path to BAM file", metavar = "character"),
  make_option(c("-i", "--gtf_file_id"), type = "character", default = NULL, help = "Path to GTF file for ID task", metavar = "character"),
  make_option(c("-q", "--gtf_file_quant"), type = "character", default = NULL, help = "Path to GTF file for Quant task", metavar = "character"),
  make_option(c("-f", "--genome_fasta_file"), type = "character", default = NULL, help = "Path to genome FASTA file", metavar = "character"),
  make_option(c("-n", "--ncpu"), type = "integer", default = 1, help = "Number of CPUs to use", metavar = "integer"),
  make_option(c("-m", "--mode"), type = "character", default = "pacbio", help = "Mode: 'pacbio', 'ont', or 'Both'", metavar = "character"),
  make_option(c("-t", "--task"), type = "character", default = "ID", help = "Task: 'ID', 'Quant', or 'Both'", metavar = "character")
)

# Parse command line options
opt <- parse_args(OptionParser(option_list = option_list))

# Assign command line options to variables
bam_file <- opt$bam_file
gtf_file_id <- opt$gtf_file_id
gtf_file_quant <- opt$gtf_file_quant
genome_fasta_file <- opt$genome_fasta_file
ncpu <- opt$ncpu
mode <- opt$mode
task <- opt$task

# Validate mode
if (!mode %in% c("pacbio", "ont", "Both")) {
  stop("Invalid mode. Please choose either 'pacbio', 'ont', or 'Both'.")
}

# Validate task
if (!task %in% c("ID", "Quant", "Both")) {
  stop("Invalid task. Please choose either 'ID', 'Quant', or 'Both'.")
}

# Common parameters
bam_files <- c(Sample = bam_file)
min_intron_length <- 30
max_intron_length <- 5e6
known_intron_motifs <- c('GT-AG')
rescue_annotated_introns <- TRUE
bin_size <- 50
min_read_count <- 1
min_relative_expression <- 0
extend_spliced_transcripts <- 100
chunk_size <- 1000000
em.maxiter <- 250
em.conv <- 0.01
use_length_normalization <- TRUE

# Process BAM file
bam_parsed <- bam_to_read_structures(bam_files = bam_files)

# Prepare transcripts function
prepare_transcripts_func <- function(gtf_file) {
  if (mode == "ont") {
    gtf_to_intron_bed(gtf_file, 'intron.bed')
    intron_bed_file <- 'intron.bed'
    known_intron_granges <- rtracklayer::import(intron_bed_file)
  } else {
    known_intron_granges <- NULL
  }

  prepare_transcripts_args <- list(
    gtf_file = gtf_file,
    genome_fasta_file = genome_fasta_file,
    bam_parsed = bam_parsed,
    min_intron_length = min_intron_length,
    max_intron_length = max_intron_length,
    known_intron_motifs = known_intron_motifs,
    rescue_annotated_introns = rescue_annotated_introns,
    known_intron_granges = known_intron_granges,
    bin_size = bin_size
  )

  if (mode == "pacbio") {
    prepare_transcripts_args$min_bam_splice_read_count <- 2
    prepare_transcripts_args$min_bam_splice_fraction <- 0.1
  }

  return(do.call(prepare_transcripts, prepare_transcripts_args))
}

# Generate TCC and export GTF function
generate_tcc_and_export_gtf <- function(transcript_data, run_mode, output_gtf) {
  se_tcc <- bam_to_tcc(
    bam_files = bam_files,
    transcript_data = transcript_data,
    run_mode = run_mode,
    min_read_count = min_read_count,
    min_relative_expression = min_relative_expression,
    extend_spliced_transcripts = extend_spliced_transcripts,
    chunk_size = chunk_size,
    ncpu = ncpu
  )

  se_transcript <- tcc_to_transcript(
    se_tcc = se_tcc,
    em.maxiter = em.maxiter,
    em.conv = em.conv,
    use_length_normalization = use_length_normalization,
    ncpu = ncpu
  )
  export_gtf(se_transcript, output_gtf)

  return(se_transcript)
}

# Conditional execution based on task
if (task %in% c("ID", "Both")) {
  transcript_data_id <- prepare_transcripts_func(gtf_file_id)

  se_transcript_id_loose <- generate_tcc_and_export_gtf(transcript_data_id, 'de_novo_loose', 'Isosceles_de_novo_loose.gtf')
  se_transcript_id_strict <- generate_tcc_and_export_gtf(transcript_data_id, 'de_novo_strict', 'Isosceles_de_novo_strict.gtf')
}

if (task %in% c("Quant", "Both")) {
  transcript_data_quant <- prepare_transcripts_func(gtf_file_quant)

  se_transcript_quant <- generate_tcc_and_export_gtf(transcript_data_quant, 'strict', 'Isosceles_strict.gtf')

  # Process counts
  counts <- assay(se_transcript_quant, "counts")

  if (!is.data.frame(counts)) {
    counts <- as.data.frame(counts)
  }

  transcript_mapping <- list()
  gtf_lines <- readLines('Isosceles_strict.gtf')

  for (line in gtf_lines) {
    if (startsWith(line, "#")) next
    fields <- strsplit(line, "\t")[[1]]
    attributes <- fields[9]
    transcript_id_match <- str_match(attributes, 'transcript_id "([^"]+)"')[,2]
    compatible_tx_match <- str_match(attributes, 'compatible_tx "([^"]+)"')[,2]
    if (!is.na(transcript_id_match) && !is.na(compatible_tx_match)) {
      transcript_mapping[[transcript_id_match]] <- compatible_tx_match
    }
  }

  counts <- counts %>%
    rownames_to_column(var = "transcript_id") %>%
    rename(count = Sample)

  counts$transcript_id <- sapply(counts$transcript_id, function(x) {
    if (x %in% names(transcript_mapping)) {
      return(transcript_mapping[[x]])
    } else {
      return(x)
    }
  })

  write.table(counts, 'Isosceles_quant.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}