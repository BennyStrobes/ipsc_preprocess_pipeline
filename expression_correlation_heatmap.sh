#!/bin/bash
#SBATCH --time=20:00:00 --mem=12G --partition=broadwl

expression_file_cross_data_sets="$1"
target_region_file="$2"
visualize_total_expression_dir="$3"



if false; then
python organize_expression_correlation.py $expression_file_cross_data_sets $target_region_file $visualize_total_expression_dir
fi

Rscript visualize_expression_correlation.R $visualize_total_expression_dir