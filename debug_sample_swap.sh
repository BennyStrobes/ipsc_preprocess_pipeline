#!/bin/bash
#SBATCH --time=01:30:00 --mem=20G --partition=broadwl

input_file="$1"
output_dir="$2"
rna_seq_sample="$3"
genotype_dir="$4"

python debug_sample_swap.py $input_file $output_dir $rna_seq_sample $genotype_dir
Rscript visualize_sample_swap.R $output_dir$rna_seq_sample"_ref_all_genotypes.txt" $output_dir$rna_seq_sample"_tot_all_genotypes.txt" $output_dir $rna_seq_sample
