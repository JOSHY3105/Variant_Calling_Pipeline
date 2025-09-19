# Variant Calling Pipeline

Simple R script to find genetic variants from BAM files.

## What it does

This pipeline takes your aligned sequencing data in BAM format and identifies genetic variants, then saves the results as VCF files. It uses two different variant calling methods, FreeBayes and BCFtools, to ensure we capture all the important variants with good confidence. The script automatically handles quality filtering and generates summary statistics to help you understand your results.

## Requirements

Before running the pipeline, you'll need to install several R packages. You can do this by running `install.packages(c("Rsamtools", "VariantAnnotation", "GenomicRanges", "BSgenome", "parallel", "argparser", "logging"))` in your R console. Additionally, you'll need to have FreeBayes, BCFtools, and samtools installed on your system and available in your PATH.

## Usage

The basic way to run the script is `Rscript variant_calling_pipeline.R your_file.bam reference.fasta results_folder`. If you want more control, you can add options like sample name and number of threads. For example, `Rscript variant_calling_pipeline.R input.bam reference.fasta output --sample my_sample --threads 4` will run the analysis using 4 CPU threads and label your sample as "my_sample".

## Output

When the pipeline finishes, you'll find several files in your output directory. There are raw variant files from both calling methods, filtered versions that contain only high-quality variants, and a summary report with useful statistics about what was found. The filtered files are usually what you want to use for downstream analysis.

## Quality Control

The script automatically applies sensible quality filters to clean up the results. It keeps only variants that are supported by at least 10 reads, have a quality score of 30 or higher, and are present in at least 90% of your data. These settings work well for most projects, but you can modify them in the script if needed.

## Troubleshooting

If you see an error about missing BAM index, don't worry because the script will create one automatically. Memory errors usually mean you need to use fewer threads or work with smaller genomic regions. If the script can't find FreeBayes or BCFtools, make sure they're properly installed and available in your system PATH. Remember that your BAM file should be sorted, which is typically done during the alignment process, and your reference genome should match the one used for the original alignment.

## Additional Notes

Processing time depends on your file size and system resources, so larger datasets will naturally take longer to complete. The pipeline is designed to be robust and will provide detailed logging information to help you track progress and identify any issues that might arise during processing.