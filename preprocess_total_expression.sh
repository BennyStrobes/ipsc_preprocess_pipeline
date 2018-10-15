#!/bin/bash
#SBATCH --time=20:00:00 --mem=12G --partition=broadwl

preprocess_total_expression_dir="$1"
exon_file="$2"
bam_dir="$3"
visualize_total_expression_dir="$4"
metadata_input_file="$5"
covariate_dir="$6"
fastqc_dir="$7"
mixutre_hmm_cell_line_grouping_dir="$8"
gencode_gene_annotation_file="$9"
gsea_data_dir="${10}"

Rscript preprocess_total_expression.R $preprocess_total_expression_dir $exon_file $bam_dir
date


python preprocess_total_expression_by_cell_lines.py $preprocess_total_expression_dir
date


Rscript prepare_covariate_files.R $preprocess_total_expression_dir $metadata_input_file $covariate_dir $fastqc_dir
date


python pca_gene_set_enrichment.py $covariate_dir $gencode_gene_annotation_file $gsea_data_dir


python prepare_covariate_files_helper.py $covariate_dir
date


Rscript visualize_processed_total_expression.R $preprocess_total_expression_dir $visualize_total_expression_dir $covariate_dir $mixutre_hmm_cell_line_grouping_dir
date

