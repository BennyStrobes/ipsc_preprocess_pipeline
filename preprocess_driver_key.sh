#!/bin/bash
#SBATCH --time=06:00:00 --mem=12G


##Testing
################################################################################################################################################################
################################################################################################################################################################
#INPUT FILES
################################################################################################################################################################
################################################################################################################################################################

#  Directory containing fastq files from sequencing round 1
fastq_round_1_input_dir="/project2/gilad/reem/heart/fastq/"
#  Directory containing fastq files from sequencing round 2
fastq_round_2_input_dir="/project2/gilad/reem/heart2/fastq/"

# File used to map sequencing core names to (cell line, time step) names for the first round of sequencing. Produced by Katie and Reem
lane_design_round_1_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/LaneDesign2forR.csv"  
# File used to map sequencing core names to (cell line, time step) names for the secomd round of sequencing. Produced by Katie and Reem
lane_design_round_2_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/LaneDesignResequencing2forR.csv"  

#  Genotype file created by Bryce Van de Geijn. It's format is vcf-ish.
#  Downloaded from "http://eqtl.uchicago.edu/jointLCL/" under the link "genotypes of 120 YRI individuals" on September 13, 2017
genotype_input="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/genotypesYRI.gen.txt.gz"

#  Heterozygous probabilities (from impute 2) genotype information contained in the following directory
#  Directory contains:
#       1. A file containing an ordered list of samples in the genotype files (YRI_samples.txt) 
#       2. One file for each chromosome of the file format "chr1.hg19.impute2.gz" that contains heterozygous probabilities for all snps in that chromosome
#  Bryce Van de Geijn sent me these files
heterozygous_site_input_dir="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/heterozygous_probability_genotypes/"

#  Metadata/covariates compiled by Katie/Reem
metadata_input_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/pilotmetadata.csv"

#  Gencode hg19 gene annotation file
#  Downloaded from "https://www.gencodegenes.org/releases/19.html" on September 13, 2017
gencode_gene_annotation_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/gencode.v19.annotation.gtf.gz"








################################################################################################################################################################
################################################################################################################################################################
#OUTPUT DIRECTORIES (The following scripts assume these directories all exist)
################################################################################################################################################################
################################################################################################################################################################

#  Directory containing merged (across sequencing rounds) fastq files for each sample
fastq_input_dir="/project2/gilad/bstrober/ipsc_differentiation/preprocess/fastq/"

#  Directory containing fastqc results. One output file for every fastq input file
fastqc_dir="/project2/gilad/bstrober/ipsc_differentiation/preprocess/fastqc_data/"

#  Directory containing reference genome (GRCh37)
genome_dir="/project2/gilad/bstrober/ipsc_differentiation/preprocess/genome/"

#  Directory containing bams. One bam file for each sample.
bam_dir="/project2/gilad/bstrober/ipsc_differentiation/preprocess/bam/"


#  Directory containing processed counts, quantile normalized expression data
preprocess_total_expression_dir="/project2/gilad/bstrober/ipsc_differentiation/preprocess/processed_total_expression/"

#  Directory containing covariate information (covariates, PCs)
covariate_dir="/project2/gilad/bstrober/ipsc_differentiation/preprocess/covariates/"

#  Directory containing plots/figures related to exploratory analysis of the total expression data (preprocess_total_expression_dir)
visualize_total_expression_dir="/project2/gilad/bstrober/ipsc_differentiation/preprocess/visualize_total_expression/"

#  Directory containing various changes to $genotype_input so that it is ammendable to the software pipelines used to process allelic counts (ie WASP and GATK ASEReader)
genotype_dir="/project2/gilad/bstrober/ipsc_differentiation/preprocess/genotype/"

#  Directory to contain various intermediate files developed by the WASP pipeline (many derivatives of the initial bams..)
wasp_intermediate_dir="/project2/gilad/bstrober/ipsc_differentiation/preprocess/wasp_intermediate_files/"

#  Directory containing raw allelic counts
raw_allelic_counts_dir="/project2/gilad/bstrober/ipsc_differentiation/preprocess/raw_allelic_counts/"

#  Directory containing processed allelic counts
processed_allelic_counts_dir="/project2/gilad/bstrober/ipsc_differentiation/preprocess/processed_allelic_counts/"

#  Directory containing plots/figures related to exploratory analysis of the total expression data (preprocess_total_expression_dir)
visualize_allelic_counts_dir="/project2/gilad/bstrober/ipsc_differentiation/preprocess/visualize_allelic_counts/"


sample_swap_check_dir="/project2/gilad/bstrober/ipsc_differentiation/preprocess/debug_sample_swap/"





#############################################################################################################################
# Scripts to preprocess total expresssion data
#############################################################################################################################

# Run the following 4 steps in series.
##########################################
##PART 1
# Concatenate all fastq files for each sample into one "merged" fastq file. Samples have more than one fastq file initially because there are multiple sequencing rounds (to increase read depth)
# Also, create file called "$fastq_dir/fastq_mapping.txt" that contains mapping from merged fastq file to all of the corresponding replicate fastq file names. First column is merged fastq file name, second column is ','-seperated list of original files.
# Both merged fastq files and fastq_mapping.txt will be saved in output directory, $fastq_dir
# Note: This code really isn't the best b/c it very manually parses the fastq files based on the lane_design files. So caution should be taken in applying merge_fastq_replicates to new situations.
# Takes about 20 minutes to run
sh merge_fastq_replicates.sh $fastq_round_1_input_dir $fastq_round_2_input_dir $lane_design_round_1_file $lane_design_round_2_file $fastq_input_dir"2"


##Part 2
# Run "sbatch fastqc_and_download_reference_genome.sh $fastq_input_dir $fastqc_dir $genome_dir" to completion.
# This script runs fastqc on each of the fastq files, as well as downloads the reference genome
# Takes about 8 hours to run
if false; then
sbatch fastqc_and_download_reference_genome.sh $fastq_input_dir $fastqc_dir $genome_dir
fi

##PART 3
# Run "sh submit-subread.sh $fastq_input_dir $bam_dir $genome_dir". This shell submits (sbatch) a job for each sample (fastq file).
# Each job aligns fastq files using subread and creates and bam in $bam_dir
# Each job takes under 2 hours to run
if false; then
sh submit-subread.sh $fastq_input_dir $bam_dir $genome_dir
fi


# PART 4
# Run sbatch preprocess_total_expression.sh $preprocess_total_expression_dir $lane_design_file $exon_file $bam_dir $visualize_total_expression_dir $fastqc_dir
# This script:
#    1. Processes the aligned read counts and creates quantile normalized expression data (preprocess_total_expression.R)
#    2. Prepares covariate files (prepare_covariate_files.R)
#    3. Also does some exploratory visualization analysis of the expression data  (visualize_processed_total_expression.R)
#  Takes about 4 hours to run
exon_file=$genome_dir"exons.saf"
if false; then
sh preprocess_total_expression.sh $preprocess_total_expression_dir $exon_file $bam_dir $visualize_total_expression_dir $metadata_input_file $covariate_dir $fastqc_dir
fi



#############################################################################################################################
# Preprocess allelic counts
#############################################################################################################################

#  wasp_maping_pipeline_part1.sh
#  This includes:
#       1. WASP Mapping Pipeline Step 1: Create text based SNP files. 
#             Completed through 'create_text_based_snp_files.py'
#             WASP description: The text-based input files have three space-delimited columns
#                  (position, ref_allele, alt_allele), and one input file per chromosome.
#                  The filenames must contain the name of the chromosome (e.g. chr2).
#             We use file names 1.snps.txt
# Takes about an hour to run

# file created in merge_fastq_replicates.sh that contains the name of each of our samples. We are going to use this to filter genotypes
sample_names=$fastq_input_dir"fastq_mapping.txt"
if false; then
sbatch wasp_mapping_pipeline_part1.sh $genotype_input $heterozygous_site_input_dir $genotype_dir $genome_dir $sample_names
fi

# Genotype file that is in VCF format.
# Made by wasp_mapping_piepline_part1.sh
vcf_file=$genotype_dir"YRI_genotype.vcf.gz"


# wasp_mapping_pipeline_part2.sh (run in parallel for each sample)
#       2. WASP Mapping Pipeline Step 2: Initial Mapping of fastq files
#             WASP description: Map the fastq files using your favorite mapper/options and filter for
#                  quality using a cutoff of your choice
#             Then sort and index those bams.

#LOOP THROUGH sample_names and get:
if false; then
while read standard_id_fastq sequencer_id; do
    standard_id=${standard_id_fastq::${#standard_id_fastq}-9}
    sbatch wasp_mapping_pipeline_part2.sh $standard_id $genotype_dir $fastq_input_dir $wasp_intermediate_dir $genome_dir $vcf_file $raw_allelic_counts_dir
done<$sample_names
fi

# Now that we have run the WASP mapping pipeline and GATK ASEReadCounter, we now have one allelic count file per sample
# process_and_organize_allelic_counts.sh will:
######### 1. Merge all of the sample specific count files into one table
######### 2. Map the sites to protein coding genes and remove sites that don't lie on a protein-coding gene
######### 3. For various heterozygous probability thresholds, place NA for (sample,site) pairs that have het. prob less than specified threshold
######### 4. Apply various filters for sites based on number of samples that we have mapped read counts to a het. site (etc)
if false; then
sh process_and_organize_allelic_counts.sh $raw_allelic_counts_dir $processed_allelic_counts_dir $genotype_dir $preprocess_total_expression_dir $gencode_gene_annotation_file $visualize_allelic_counts_dir
fi



# debug_sample_swap_driver.sh checks to make sure that every RNA-seq sample (has the correct label) and is paired correctly with its corresponding genotype
# We will test this by looping through each rna-seq sample, and for each sample:
##### Compare fraction of heterozygous sites that show bi-allelic expression for every possible genotype. 
##### There should only be one genotype that results in a high fraction. Make sure this genotype corresponds to the rna-seq sample
##### Takes about 8 hours to run.
if false; then
sbatch debug_sample_swap_driver.sh $processed_allelic_counts_dir $sample_swap_check_dir $genotype_dir $sample_names
fi

