#!/bin/bash
#SBATCH --time=01:30:00 --mem=20G --partition=broadwl

genotype_input="$1"
heterozygous_site_input_dir="$2"
genotype_dir="$3"
genome_dir="$4"
sample_names="$5"

# Takes about 40 min to run
####################################################################
#  1. WASP Mapping Pipeline Step 1: Create text based SNP-site files. (1 for each chromosome)
#  These files only contain genotype sites (not genotype information)
####################################################################
python create_text_based_snp_files.py $genotype_input $genotype_dir
# Zip up sites file
gzip $genotype_dir*snps.txt

########################################################
# 2. Make pseudo-VCF into acceptable vcf file format (Make vcf genotype file and vcf heterozygous site file)
########################################################
# About 10 min
python convert_genotype_to_vcf.py $genotype_input $genotype_dir $heterozygous_site_input_dir $sample_names
# Zip up
bgzip -c $genotype_dir"YRI_genotype.vcf" > $genotype_dir"YRI_genotype.vcf.gz"
tabix -p vcf $genotype_dir"YRI_genotype.vcf.gz"


#########################################################
# 3. Get Reference genome in correct format for ASE Mapping
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

