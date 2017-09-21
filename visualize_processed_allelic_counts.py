import numpy as np 
import os
import sys
import pdb
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def get_total_counts(input_data):
    row, col = input_data.shape
    input_mat = input_data[1:,1:]
    total_counts = np.zeros((row-1,col-1))
    for i in range(row-1):
        for j in range(col-1):
            if input_mat[i,j] == 'Nan':
                total_counts[i,j] = float('Nan')
            else:
                total_counts[i,j] = float(input_mat[i,j].split('_')[1])
    return total_counts

def get_ref_counts(input_data):
    row, col = input_data.shape
    input_mat = input_data[1:,1:]
    total_counts = np.zeros((row-1,col-1))
    for i in range(row-1):
        for j in range(col-1):
            if input_mat[i,j] == 'Nan':
                total_counts[i,j] = float('Nan')
            else:
                total_counts[i,j] = float(input_mat[i,j].split('_')[0])
    return total_counts

def number_expressed_het_snps_per_individual_boxplot(input_file, output_file, het_thresh):
    input_data = np.loadtxt(input_file,dtype=str)
    total_counts = get_total_counts(input_data)
    num_sites, num_samples = total_counts.shape
    read_count_threshs = [2,6,10,15,30]
    plotting_array = []
    fig = plt.figure()
    for read_count_thresh in read_count_threshs:
        nested_array = []
        for sample_i in range(num_samples):
            num_sites = len(np.where(total_counts[:,sample_i] > read_count_thresh)[0])
            nested_array.append(num_sites)
        plotting_array.append(nested_array)
    plt.boxplot(plotting_array,notch=True, patch_artist=True)
    plt.xticks(range(1, 1+len(read_count_threshs)), np.asarray(read_count_threshs).astype(str))
    plt.xlabel('Reads / het-SNP')
    plt.ylabel('het-SNPs per individual')
    plt.title('Het thresh: ' + str(het_thresh))
    fig.savefig(output_file)

def get_sites(input_file):
    sites = {}
    f = open(input_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ss = data[0].split('_')
        site_id = ss[0] + '_' + ss[1] + '_' + ss[2] + '_' + ss[3] + '_' + ss[4]
        sites[site_id] = 1
    return sites

def het_prob_histogram(het_prob_file, output_file, site_file):
    sites = get_sites(site_file)
    f = open(het_prob_file)
    probs = []
    fig = plt.figure()

    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#'):
            continue
        site_id = data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3] + '_' + data[4]
        if site_id not in sites:
            continue
        probz = np.asarray(data[9:]).astype(float)
        for prob in probz:
            probs.append(prob)
    plt.hist(probs, bins=80)
    plt.xlabel('Heterozygous probability (from impute2)')
    plt.ylabel('Number of sites')
    fig.savefig(output_file)


def percent_biallelic_het_snps_scatter(input_file, output_file):
    input_data = np.loadtxt(input_file,dtype=str)
    total_counts = get_total_counts(input_data)
    ref_counts = get_ref_counts(input_data)
    percent_biallelic = []
    num_sites, num_samples = total_counts.shape
    fig = plt.figure()

    for i in range(num_samples):
        total = 0
        biallelic = 0
        for j in range(num_sites):
            if np.isnan(total_counts[j,i]):
                continue
            if total_counts[j,i] < 10:
                continue
            total = total + 1
            if np.isnan(ref_counts[j,i]):
                print('ERORORRO')
                pdb.set_trace()
            if ref_counts[j,i] != total_counts[j,i]:
                biallelic = biallelic+1
        if 1.0*biallelic/total < .7:
            print(input_data[0,i+1] + '\t' + str(biallelic) + '\t' + str(total))

        percent_biallelic.append(1.0*biallelic/total)
    plt.scatter(range(num_samples), percent_biallelic)
    plt.xlabel('Sample Num')
    plt.ylabel('Percent Biallelic')
    fig.savefig(output_file)


input_dir = sys.argv[1]  # Input directory
output_dir = sys.argv[2]  # Directory to save output images
genotype_dir = sys.argv[3]

het_prob_file = genotype_dir + 'YRI_het_prob_genotype.vcf'

het_prob_histogram(het_prob_file, output_dir + 'het_prob_histogram.png', input_dir + 'allelic_counts_gene_mapped_het_prob_999.txt')
# percent_biallelic_het_snps_scatter( input_dir + 'allelic_counts_gene_mapped_het_prob_9.txt', output_dir + 'percent_biallelic_het_snps_scatter.png')







het_threshs = [.55, .6, .65, .7, .75, .8, .85, .9, .95,.99,.999]  # Various heterozygous thresholds we picked
#het_threshs= [.55]
# Perform analysis for each het_thresh
for het_thresh in het_threshs:
    # Input file for specified heterozygous threshold
    input_file = input_dir + 'allelic_counts_gene_mapped_het_prob_' + str(het_thresh).split('.')[1] + '.txt'

    #number_expressed_het_snps_per_individual_boxplot(input_file, output_dir + 'number_expressed_het_snps_per_individual_boxplot_' +  str(het_thresh).split('.')[1] + '.png', het_thresh)
