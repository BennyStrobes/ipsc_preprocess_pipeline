import numpy as np
import os
import sys
import pdb
import gzip
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


#  Create list of sites (variants) we are interested in 
def get_list_of_sites(input_file):
    f = open(input_file)
    head_count = 0  #  For header
    dicti = {}  # Initialize dictionary to keep track of sites
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # Skip header
            head_count = head_count + 1
            continue
        # Remove gene name from site id identifier
        full_site_id = data[0]
        site_id_info = full_site_id.split('_')
        site_id = site_id_info[0] + '_' + site_id_info[1] + '_' + site_id_info[2] + '_' + site_id_info[3] + '_' + site_id_info[4]
        dicti[site_id] = 1.0
    return dicti


#  Create dictionary that has keys that are sites, and values that are (nested) dictionaries.
#  The nested dictionaries have keys that are cell_line_ids if that cell_line is heterozygous at the site
#  therefore dictionary can be accessed by is_site_heterozygous[siteID][cellLine].
def learn_heterozygous_sites(sites, het_thresh, het_prob_genotype_file):
    dicti = {}  # Initialize dictionary
    f = open(het_prob_genotype_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#CHROM'):  # If sample id line
            ordered_cell_lines = np.asarray(data[9:])
        if line.startswith('#'):  # Skip other headers
            continue
        site_id = data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3] + '_' + data[4]  # create site_id from line
        if site_id not in sites:  # Ignore sites that we do not have in our data
            continue
        het_probs = np.asarray(data[9:]).astype(float)  # Convert heterozygous probs into array
        # Error checking
        if site_id in dicti:
            print('EROROROR')
            pdb.set_trace()
        dicti[site_id] = {}
        #  Loop through each cellLine
        for i, het_prob in enumerate(het_probs):
            i_cell_line = ordered_cell_lines[i]  # Name of ith cell line
            if het_prob >= het_thresh:  # Is ith cell line heterozygous at this site
                dicti[site_id][i_cell_line] = 1.0
    return dicti, ordered_cell_lines



input_file = sys.argv[1]
output_dir = sys.argv[2]
rna_seq_sample = sys.argv[3]
genotype_dir = sys.argv[4]








f = open(input_file)
head_count = 0
sites = {}
for line in f:
    line = line.rstrip()
    data = line.split()
    if head_count == 0:
        head_count = head_count + 1
        samples = data[1:]
        lines = []
        for sample in samples:
            lines.append(sample.split('_')[0])
            sites[sample.split('_')[0]] = {}
        continue
    counts = data[1:]
    full_site_id = data[0]
    site_id_info = full_site_id.split('_')
    site_id = site_id_info[0] + '_' + site_id_info[1] + '_' + site_id_info[2] + '_' + site_id_info[3] + '_' + site_id_info[4]
    for i,count in enumerate(counts):
        i_cell_line = lines[i]
        if count == 'Nan':
            continue
        ref = int(count.split('_')[0])
        total = int(count.split('_')[1])
        if total > 4 and ref != 0 and ref != total:
            if site_id not in sites[i_cell_line]:
                sites[i_cell_line][site_id] = 1
lines = np.asarray(lines)
counts = []
unique_lines = np.unique(lines)
for line in unique_lines:
    counts.append(len(sites[line]))


fig = plt.figure()
plt.scatter(range(len(counts)), counts)
plt.xticks(range(len(counts)), unique_lines)
plt.xticks(rotation=90)
plt.xlabel('Cell Line')
plt.gcf().subplots_adjust(bottom=0.15)
plt.ylabel('Number of unique bi-allelic sites')
fig.savefig(output_dir + 'num_unique_sites.png')



f.close()











full_genotype_file = genotype_dir + 'YRI_het_prob_genotype_all_samples.vcf'


sites = get_list_of_sites(input_file)

is_site_heterozygous, ordered_cell_lines = learn_heterozygous_sites(sites, .95, full_genotype_file)



output_file_ref = output_dir + rna_seq_sample + '_ref_all_genotypes.txt'
output_file_tot = output_dir + rna_seq_sample + '_tot_all_genotypes.txt'


f = open(input_file)
t_ref = open(output_file_ref, 'w')
t_tot = open(output_file_tot, 'w')


head_count = 0
for line in f:
    line = line.rstrip()
    data = line.split()
    if head_count == 0:
        head_count = head_count + 1
        t_ref.write(data[0] + '\t' + '\t'.join(ordered_cell_lines) + '\n')
        t_tot.write(data[0] + '\t' + '\t'.join(ordered_cell_lines) + '\n')

        samples = data[1:]
        for i,val in enumerate(samples):
            if val == rna_seq_sample:
                column_num = i
        continue
    t_ref.write(data[0])
    t_ref.write(data[0])

    counts = data[1:]
    sample_count = counts[column_num]
    full_site_id = data[0]
    site_id_info = full_site_id.split('_')
    site_id = site_id_info[0] + '_' + site_id_info[1] + '_' + site_id_info[2] + '_' + site_id_info[3] + '_' + site_id_info[4]
    # Loop through genotypes
    for i,cell_line in enumerate(ordered_cell_lines):
        if cell_line in is_site_heterozygous[site_id]:
            if sample_count == 'Nan':
                t_ref.write('\t' + '0')
                t_tot.write('\t' + '0')
            else:
                t_ref.write('\t' + sample_count.split('_')[0])
                t_tot.write('\t' + sample_count.split('_')[1])
        else:
            t_ref.write('\tNan')
            t_tot.write('\tNan')
    t_ref.write('\n')
    t_tot.write('\n')
t_ref.close()
t_tot.close()
f.close()

