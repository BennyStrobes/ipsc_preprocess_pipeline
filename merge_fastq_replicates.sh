#!/bin/bash
#SBATCH --time=01:00:00


fastq_round_1_input_dir="$1"
fastq_round_2_input_dir="$2"
lane_design_round_1_file="$3"
lane_design_round_2_file="$4"
fastq_dir="$5"

python merge_fastq_replicates.py $fastq_round_1_input_dir $fastq_round_2_input_dir $lane_design_round_1_file $lane_design_round_2_file $fastq_dir