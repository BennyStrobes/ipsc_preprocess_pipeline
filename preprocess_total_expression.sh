#!/bin/bash
#SBATCH --time=20:00:00 --mem=12G

preprocess_total_expression_dir="$1"
exon_file="$2"
bam_dir="$3"
visualize_total_expression_dir="$4"
metadata_input_file="$5"
covariate_dir="$6"
fastqc_dir="$7"

if false; then
Rscript preprocess_total_expression.R $preprocess_total_expression_dir $exon_file $bam_dir
fi
Rscript prepare_covariate_files.R $preprocess_total_expression_dir $metadata_input_file $covariate_dir $fastqc_dir


Rscript visualize_processed_total_expression.R $preprocess_total_expression_dir $visualize_total_expression_dir $covariate_dir
