#!/bin/bash
#SBATCH --time=20:30:00 --mem=20G --partition=broadwl

genotype_input="$1"
heterozygous_site_input_dir="$2"
genotype_dir="$3"
genome_dir="$4"
sample_names="$5"
impute2_genotype_dir="$6"
chrom_info_file="$7"
snp2h5_dir="$8"

# Takes about 40 min to run
####################################################################
#  A. WASP Mapping Pipeline Step 1: Create text based SNP-site files. (1 for each chromosome)
#  These files only contain genotype sites (not genotype information)
####################################################################
python create_text_based_snp_files.py $genotype_input $genotype_dir
# Zip up sites file
gzip $genotype_dir*snps.txt

########################################################
# B. Make pseudo-VCF into acceptable vcf file format (Make vcf genotype file and vcf heterozygous site file)
########################################################
# About 10 min
python convert_genotype_to_vcf.py $genotype_input $genotype_dir $heterozygous_site_input_dir $sample_names
# Zip up
bgzip -c $genotype_dir"YRI_genotype.vcf" > $genotype_dir"YRI_genotype.vcf.gz"
tabix -p vcf $genotype_dir"YRI_genotype.vcf.gz"


#########################################################
# C. Get Reference genome in correct format for ASE Mapping
#########################################################
# Zip up each of reference genomes
for chrom_num in {1..22}
do
    gunzip -c $genome_dir"Homo_sapiens.GRCh37.dna_sm.chromosome."$chrom_num".fa.gz" > $genome_dir"Homo_sapiens.GRCh37.dna_sm.chromosome."$chrom_num".fa"
done

gunzip -c $genome_dir"Homo_sapiens.GRCh37.dna_sm.chromosome.MT.fa.gz" > $genome_dir"Homo_sapiens.GRCh37.dna_sm.chromosome.MT.fa"
gunzip -c $genome_dir"Homo_sapiens.GRCh37.dna_sm.chromosome.X.fa.gz" > $genome_dir"Homo_sapiens.GRCh37.dna_sm.chromosome.X.fa"
gunzip -c $genome_dir"Homo_sapiens.GRCh37.dna_sm.chromosome.Y.fa.gz" > $genome_dir"Homo_sapiens.GRCh37.dna_sm.chromosome.Y.fa"

# Concatenate reference genomes into 1
word=""
for C in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT
do
    word=$word$genome_dir"Homo_sapiens.GRCh37.dna_sm.chromosome."$C".fa "
done

# Now sort reference genomes so they are in correct format for ASE mapping
cat $word > $genome_dir"Homo_sapiens.GRCh37.dna_sm.fa"
samtools faidx $genome_dir"Homo_sapiens.GRCh37.dna_sm.fa"
java -jar picard.jar CreateSequenceDictionary R=$genome_dir"Homo_sapiens.GRCh37.dna_sm.fa" O=$genome_dir"Homo_sapiens.GRCh37.dna_sm.dict"


#########################################################
# D. Convert impute2 genotype information to h5 format using WASP's snp2h5 script
#########################################################

all_genotyped_samples_file=$genotype_dir"all_genotyped_samples.txt"
used_samples_file=$genotype_dir"used_samples.txt"

python get_genotype_sample_names.py $genotype_dir"YRI_genotype.vcf" $used_samples_file
python get_genotype_sample_names.py $genotype_dir"YRI_het_prob_genotype_all_samples.vcf" $all_genotyped_samples_file



$snp2h5_dir"snp2h5" --chrom $chrom_info_file \
    --format impute \
    --geno_prob $genotype_dir"geno_probs.h5" \
    --snp_index $genotype_dir"snp_index.h5" \
    --snp_tab $genotype_dir"snp_tab.h5" \
    --haplotype $genotype_dir"haps.h5" \
    --samples $all_genotyped_samples_file \
    $impute2_genotype_dir"chr"*".hg19.impute2.gz" \
    $impute2_genotype_dir"chr"*".hg19.impute2_haps.gz"

#########################################################
# E. Convert fasta information to h5 format using WASP's fasta2h5 script
#########################################################


$snp2h5_dir"fasta2h5" --chrom $chrom_info_file \
    --seq $genotype_dir"seq.h5" \
    $genome_dir"Homo_sapiens.GRCh37.dna_sm.name_edit.chr"*".fa.gz"
