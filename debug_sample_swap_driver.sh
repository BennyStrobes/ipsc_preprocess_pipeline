#!/bin/bash
#SBATCH --time=15:30:00 --mem=20G


processed_allelic_counts_dir="$1"
sample_swap_check_dir="$2"
genotype_dir="$3"
sample_names="$4"



while read standard_id_fastq sequencer_id; do
    standard_id=${standard_id_fastq::${#standard_id_fastq}-16}
    echo $standard_id
    sh debug_sample_swap.sh $processed_allelic_counts_dir"allelic_counts_gene_mapped.txt" $sample_swap_check_dir $standard_id $genotype_dir
done<$sample_names
