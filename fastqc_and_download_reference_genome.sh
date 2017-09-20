#!/bin/bash
#SBATCH --time=12:00:00 --mem=12G
fastq_input_dir="$1"
fastqc_dir="$2"
genome_dir="$3"

date
# Script provided by John Blischak (https://github.com/jdblischak/midway-subread-pipeline)
Rscript download-genome.R $genome_dir


# Script provided by John Blischak (https://github.com/jdblischak/midway-subread-pipeline)
Rscript download-exons.R $genome_dir
date

# Requires fastqc,
fastqc $fastq_input_dir*fastq.gz -o $fastqc_dir
date

# Multiqc combines all of our fastqc files into one organized one!
multiqc --force --outdir $fastqc_dir $fastqc_dir $fastq_input_dir