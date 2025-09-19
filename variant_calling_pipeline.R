#!/usr/bin/env Rscript

# Variant Calling Pipeline: BAM to VCF
# Description: Pipeline for variant calling from BAM files using multiple callers

# Load required libraries
suppressPackageStartupMessages({
  library(Rsamtools)
  library(VariantAnnotation)
  library(GenomicRanges)
  library(BSgenome)
  library(parallel)
  library(argparser)
  library(logging)
})

# Setup logging
basicConfig()

# Function to validate input files
validate_inputs <- function(bam_file, ref_genome, output_dir) {
  if (!file.exists(bam_file)) {
    stop("BAM file not found: ", bam_file)
  }
  
  if (!file.exists(ref_genome)) {
    stop("Reference genome not found: ", ref_genome)
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    loginfo("Created output directory: %s", output_dir)
  }
  
  # Check BAM index
  bai_file <- paste0(bam_file, ".bai")
  if (!file.exists(bai_file)) {
    loginfo("Creating BAM index...")
    indexBam(bam_file)
  }
  
  return(TRUE)
}

# Function to get basic BAM statistics
get_bam_stats <- function(bam_file) {
  bf <- BamFile(bam_file)
  open(bf)
  
  # Get header info
  header <- scanBamHeader(bf)
  sequences <- header[[1]]$targets
  
  # Basic counts
  total_reads <- countBam(bf)$records
  
  close(bf)
  
  stats <- list(
    total_reads = total_reads,
    sequences = length(sequences),
    seq_names = names(sequences)[1:min(5, length(sequences))]
  )
  
  return(stats)
}

# Function to run FreeBayes variant calling
run_freebayes <- function(bam_file, ref_genome, output_vcf, 
                          min_coverage = 10, min_qual = 20, 
                          ploidy = 2, threads = 1) {
  
  loginfo("Starting FreeBayes variant calling...")
  
  # Construct FreeBayes command
  fb_cmd <- sprintf(
    "freebayes -f %s -C %d -q %d -p %d --vcf %s %s",
    ref_genome, min_coverage, min_qual, ploidy, output_vcf, bam_file
  )
  
  # Run FreeBayes
  exit_code <- system(fb_cmd, wait = TRUE)
  
  if (exit_code != 0) {
    stop("FreeBayes failed with exit code: ", exit_code)
  }
  
  loginfo("FreeBayes completed successfully")
  return(output_vcf)
}

# Function to run BCFtools variant calling
run_bcftools <- function(bam_file, ref_genome, output_vcf,
                         min_coverage = 10, min_qual = 20) {
  
  loginfo("Starting BCFtools variant calling...")
  
  # BCFtools mpileup and call
  temp_bcf <- tempfile(fileext = ".bcf")
  
  mpileup_cmd <- sprintf(
    "bcftools mpileup -f %s -d %d -q %d -O b -o %s %s",
    ref_genome, 1000, min_qual, temp_bcf, bam_file
  )
  
  call_cmd <- sprintf(
    "bcftools call -mv -O v -o %s %s",
    output_vcf, temp_bcf
  )
  
  # Run commands
  exit_code1 <- system(mpileup_cmd, wait = TRUE)
  exit_code2 <- system(call_cmd, wait = TRUE)
  
  if (exit_code1 != 0 || exit_code2 != 0) {
    stop("BCFtools failed")
  }
  
  # Cleanup
  unlink(temp_bcf)
  
  loginfo("BCFtools completed successfully")
  return(output_vcf)
}

# Function to filter VCF
filter_vcf <- function(input_vcf, output_vcf, min_dp = 10, 
                       min_qual = 30, max_missing = 0.1) {
  
  loginfo("Filtering VCF file...")
  
  # Read VCF
  vcf <- readVcf(input_vcf, "hg38")
  
  # Apply filters
  dp_filter <- geno(vcf)$DP >= min_dp
  qual_filter <- qual(vcf) >= min_qual
  
  # Combine filters
  keep_variants <- qual_filter & 
    rowSums(is.na(dp_filter)) / ncol(dp_filter) <= max_missing
  
  # Subset VCF
  filtered_vcf <- vcf[keep_variants, ]
  
  # Write filtered VCF
  writeVcf(filtered_vcf, output_vcf)
  
  loginfo("Filtered %d variants from %d total", 
          sum(keep_variants), length(keep_variants))
  
  return(output_vcf)
}

# Function to annotate variants
annotate_variants <- function(vcf_file, output_file) {
  
  loginfo("Annotating variants...")
  
  # Read VCF
  vcf <- readVcf(vcf_file, "hg38")
  
  # Add basic annotations
  info(vcf)$AC <- geno(vcf)$GT %>% 
    apply(1, function(x) sum(x == "1/1") * 2 + sum(x == "0/1"))
  
  info(vcf)$AF <- info(vcf)$AC / (ncol(vcf) * 2)
  
  # Write annotated VCF
  writeVcf(vcf, output_file)
  
  loginfo("Variant annotation completed")
  return(output_file)
}

# Function to generate summary report
generate_report <- function(vcf_files, output_file) {
  
  loginfo("Generating summary report...")
  
  report <- data.frame(
    caller = character(),
    total_variants = numeric(),
    snvs = numeric(),
    indels = numeric(),
    mean_quality = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(vcf_files)) {
    if (file.exists(vcf_files[i])) {
      vcf <- readVcf(vcf_files[i], "hg38")
      
      # Calculate statistics
      total_vars <- nrow(vcf)
      snvs <- sum(nchar(ref(vcf)) == 1 & nchar(alt(vcf)) == 1)
      indels <- total_vars - snvs
      mean_qual <- mean(qual(vcf), na.rm = TRUE)
      
      report <- rbind(report, data.frame(
        caller = names(vcf_files)[i],
        total_variants = total_vars,
        snvs = snvs,
        indels = indels,
        mean_quality = round(mean_qual, 2)
      ))
    }
  }
  
  # Write report
  write.csv(report, output_file, row.names = FALSE)
  
  loginfo("Summary report saved to: %s", output_file)
  return(report)
}

# Main pipeline function
variant_calling_pipeline <- function(bam_file, ref_genome, output_dir,
                                     sample_name = "sample",
                                     callers = c("freebayes", "bcftools"),
                                     threads = 1) {
  
  loginfo("Starting variant calling pipeline for: %s", sample_name)
  
  # Validate inputs
  validate_inputs(bam_file, ref_genome, output_dir)
  
  # Get BAM statistics
  bam_stats <- get_bam_stats(bam_file)
  loginfo("BAM file contains %d reads across %d sequences", 
          bam_stats$total_reads, bam_stats$sequences)
  
  # Initialize results list
  vcf_files <- list()
  
  # Run FreeBayes if requested
  if ("freebayes" %in% callers) {
    fb_vcf <- file.path(output_dir, paste0(sample_name, "_freebayes.vcf"))
    vcf_files[["freebayes"]] <- run_freebayes(bam_file, ref_genome, fb_vcf)
  }
  
  # Run BCFtools if requested
  if ("bcftools" %in% callers) {
    bc_vcf <- file.path(output_dir, paste0(sample_name, "_bcftools.vcf"))
    vcf_files[["bcftools"]] <- run_bcftools(bam_file, ref_genome, bc_vcf)
  }
  
  # Filter VCFs
  filtered_vcfs <- list()
  for (caller in names(vcf_files)) {
    filtered_vcf <- file.path(output_dir, 
                              paste0(sample_name, "_", caller, "_filtered.vcf"))
    filtered_vcfs[[caller]] <- filter_vcf(vcf_files[[caller]], filtered_vcf)
  }
  
  # Generate summary report
  report_file <- file.path(output_dir, paste0(sample_name, "_summary.csv"))
  summary_report <- generate_report(filtered_vcfs, report_file)
  
  loginfo("Pipeline completed successfully!")
  
  return(list(
    vcf_files = filtered_vcfs,
    summary = summary_report,
    bam_stats = bam_stats
  ))
}

# Command line interface
if (!interactive()) {
  
  # Parse arguments
  p <- arg_parser("Variant calling pipeline from BAM to VCF")
  p <- add_argument(p, "bam", help = "Input BAM file")
  p <- add_argument(p, "ref", help = "Reference genome FASTA file")
  p <- add_argument(p, "output", help = "Output directory")
  p <- add_argument(p, "--sample", default = "sample", help = "Sample name")
  p <- add_argument(p, "--callers", nargs = "+", 
                    default = c("freebayes", "bcftools"),
                    help = "Variant callers to use")
  p <- add_argument(p, "--threads", type = "integer", default = 1,
                    help = "Number of threads")
  
  args <- parse_args(p)
  
  # Run pipeline
  results <- variant_calling_pipeline(
    bam_file = args$bam,
    ref_genome = args$ref,
    output_dir = args$output,
    sample_name = args$sample,
    callers = args$callers,
    threads = args$threads
  )
  
  cat("Variant calling completed successfully!\n")
  cat("Results saved to:", args$output, "\n")
}

# Example usage:
# Rscript variant_calling_pipeline.R sample.bam reference.fasta output_dir --sample my_sample