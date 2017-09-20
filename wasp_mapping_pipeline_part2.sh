#!/bin/bash
#SBATCH --mem=30G --time=03:30:00 --partition=broadwl


standard_id="$1"
genotype_dir="$2"
fastq_dir="$3"
wasp_intermediate_dir="$4"
genome_dir="$5"
vcf_file="$6"
preprocess_allelic_counts_dir="$7"

date

##################################################################################################################
# WASP Mapping Pipeline Step 2: Initial Mapping of fastq files
#             WASP description: Map the fastq files using your favorite mapper/options and filter for
#                  quality using a cutoff of your choice
#             Then sort and index those bams.
# This part takes around 25 minutes

echo "WASP STEP 2: Initial Mapping of fastq file (with subread)"
fq_file=$fastq_dir$standard_id".fastq.gz"
Rscript run-subread.R $fq_file $wasp_intermediate_dir $genome_dir


#  Sort and index the bams
samtools sort -f $wasp_intermediate_dir$standard_id.bam $wasp_intermediate_dir$standard_id.sort.bam
samtools index $wasp_intermediate_dir$standard_id.sort.bam




##################################################################################################################
# WASP Mapping Pipeline Step 3: Use find_intersecting_snps.py to identify reads that may have mapping biases
# find_intersecting_snps.py is a WASP script
# This part takes around 10 minutes
echo "WASP STEP 3: find_intersecting_snps.py"

# Run script
python find_intersecting_snps.py \
            --is_sorted \
            --output_dir $wasp_intermediate_dir \
            --snp_dir $genotype_dir \
            $wasp_intermediate_dir$standard_id.sort.bam



##################################################################################################################
# WASP Mapping Pipeline Step 4: Map the PREFIX.remap.fq.gz using the same mapping arguments used in WASP Mapping Pipeline Step 2.
# This part takes around 10 minutes
echo "WASP STEP 4: Second mapping"

fq_file=$wasp_intermediate_dir$standard_id".sort.remap.fastq.gz"

Rscript run-subread.R $fq_file $wasp_intermediate_dir $genome_dir

#  Sort and index the bams
samtools sort -f $wasp_intermediate_dir$standard_id.sort.remap.bam $wasp_intermediate_dir$standard_id.sort.remap.sort.bam
samtools index $wasp_intermediate_dir$standard_id.sort.remap.sort.bam



#################################################################################################################
# WASP Mapping Pipeline Step 5: Use filter_remapped_reads.py to filter out reads where one or more of the allelic 
#      versions of the reads fail to map back to the same location as the original read
#  This part takes around 3 minutes

echo "WASP STEP 5: filter_remapped_reads.py"
python filter_remapped_reads.py \
        $wasp_intermediate_dir$standard_id.sort.to.remap.bam \
        $wasp_intermediate_dir$standard_id.sort.remap.sort.bam \
        $wasp_intermediate_dir$standard_id.keep.bam



#################################################################################################################
# WASP Mapping Pipeline Step 6: Merge bams that we plan to use. Then sort and index
echo "WASP STEP 6: Merge bams"
#  This part takes around 5 minutes

samtools merge $wasp_intermediate_dir$standard_id.keep.merge.bam \
                $wasp_intermediate_dir$standard_id.keep.bam \
                $wasp_intermediate_dir$standard_id.sort.keep.bam

samtools sort -f $wasp_intermediate_dir$standard_id.keep.merge.bam $wasp_intermediate_dir$standard_id.keep.merge.sort.bam
samtools index $wasp_intermediate_dir$standard_id.keep.merge.sort.bam



#################################################################################################################
# WASP Mapping Pipeline Step 7: Filter duplicate reads.
echo "WASP STEP 7: Filter duplicate reads using rmdup.py"
#  This part takes around 5 minutes
# Use rmdup.py for single end reads and rmdup_pe.py for paired end!
python rmdup.py $wasp_intermediate_dir$standard_id.keep.merge.sort.bam $wasp_intermediate_dir$standard_id.wasp_corrected.bam




#################################################################################################################
# GATK ASEReadCounter (Using WASP Output)

echo "Adding uninformative groupnames to bam file (so GATK accepts the bams)"
#  This is only done so GATK accepts the bams. The group names mean nothing!
java -jar picard.jar AddOrReplaceReadGroups \
      I=$wasp_intermediate_dir$standard_id.wasp_corrected.bam  \
      O=$wasp_intermediate_dir$standard_id.wasp_corrected2.bam \
      RGID=4 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=20

echo "Sorting one more time before running GATK"
samtools sort -f $wasp_intermediate_dir$standard_id.wasp_corrected2.bam $wasp_intermediate_dir$standard_id.wasp_corrected3.bam
samtools index $wasp_intermediate_dir$standard_id.wasp_corrected3.bam


echo "Reordering bam file using picard: reorder sam (ensures correct format for GATK)"
java -jar picard.jar ReorderSam \
    I=$wasp_intermediate_dir$standard_id.wasp_corrected3.bam \
    O=$wasp_intermediate_dir$standard_id.wasp_corrected_gatk_ready.bam \
    R=$genome_dir"Homo_sapiens.GRCh37.dna_sm.fa" \
    CREATE_INDEX=TRUE


echo "Running GATK"
 java -jar GenomeAnalysisTK.jar \
   -R $genome_dir"Homo_sapiens.GRCh37.dna_sm.fa" \
   -T ASEReadCounter \
   -o $preprocess_allelic_counts_dir$standard_id"_wasp_corrected.txt" \
   -I $wasp_intermediate_dir$standard_id.wasp_corrected_gatk_ready.bam \
   -sites $vcf_file \
   --outputFormat "TABLE"



echo "DONE!"
date

