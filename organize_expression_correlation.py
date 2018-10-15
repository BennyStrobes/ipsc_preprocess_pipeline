import numpy as np 
import os
import sys
import pdb
import gzip


def extract_ipsc_ca_genes(file_name, genes):
	f = gzip.open(file_name)
	head_count = 0
	used = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble_id = data[-1].split('_')[2]
		if ensamble_id in used:
			continue
		used[ensamble_id] = 1
		if ensamble_id not in genes:
			genes[ensamble_id] = 1
		else:
			genes[ensamble_id] = genes[ensamble_id] + 1
	return genes

def extract_ipsc_genes(file_name, genes):
	f = open(file_name)
	used = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble_id = data[0].split('.')[0]
		if ensamble_id in used:
			continue
		used[ensamble_id] = 1
		if ensamble_id not in genes:
			genes[ensamble_id] = 1
		else:
			genes[ensamble_id] = genes[ensamble_id] + 1
	return genes

def extract_gtex_genes(file_name, genes):
	f = open(file_name)
	count = 0
	used = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if count < 3:
			count = count + 1
			continue
		ensamble_id = data[0].split('.')[0]
		if ensamble_id in used:
			continue
		used[ensamble_id] = 1
		if ensamble_id not in genes:
			genes[ensamble_id] = 1
		else:
			genes[ensamble_id] = genes[ensamble_id] + 1
	return genes

def extract_ipsc_time_genes(file_name, genes, converter):
	f = gzip.open(file_name)
	head_count = 0
	used = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		identifier = data[0] + '_' + data[7] + '_' + data[8]
		ensamble_id = converter[identifier]
		if ensamble_id in used:
			continue
		used[ensamble_id] = 1
		if ensamble_id not in genes:
			genes[ensamble_id] = 1
		else:
			genes[ensamble_id] = genes[ensamble_id] + 1
	return genes


def get_list_of_genes_found_in_all_data_sets(expression_file_cross_data_sets, target_region_to_ensamble):
	genes = {}
	f = open(expression_file_cross_data_sets)
	ipsc_ca_used = 0
	ipsc_used = 0
	gtex_used = 0
	ipsc_time_used =0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if data[2] == 'ipsc_ca' and ipsc_ca_used == 0:
			genes = extract_ipsc_ca_genes(data[1], genes)
			ipsc_ca_used = 1
		if data[2] == 'ipsc' and ipsc_used == 0:
			genes = extract_ipsc_genes(data[1], genes)
			ipsc_used = ipsc_used + 1
		if data[2] == 'GTEx' and gtex_used == 0:
			genes = extract_gtex_genes(data[1], genes)
			gtex_used = gtex_used + 1
		if data[2] == 'ipsc_time' and ipsc_time_used == 0:
			genes = extract_ipsc_time_genes(data[1], genes, target_region_to_ensamble)
			ipsc_time_used = ipsc_time_used + 1
	f.close()
	shared_genes = {}
	for geney in genes.keys():
		if genes[geney] == 4:
			shared_genes[geney] = 1
		if genes[geney] > 4:
			print('assumprtioner eoror')
			pdb.set_trace()
	return shared_genes 

def get_ipsc_ca_gene_counts(file_name, genes):
	f = gzip.open(file_name)
	head_count = 0
	used = {}
	count_vec = np.zeros(len(genes))
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble_id = data[-1].split('_')[2]
		if ensamble_id not in genes:
			continue
		used[ensamble_id] = 1
		read_counts = float(data[-3])
		count_vec[genes[ensamble_id]] = read_counts 
	f.close()
	if len(used) != len(genes):
		print('assumption erroror!!')
	return count_vec

def get_ipsc_gene_counts(file_name, genes, sample_name):
	short_sample_name = sample_name.split('_')[0]
	f = open(file_name)
	head_count = 0
	used = {}
	count_vec = np.zeros(len(genes))
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			for i, ele in enumerate(data):
				if ele == short_sample_name:
					index = i
			continue
		ensamble_id = data[0].split('.')[0]
		if ensamble_id not in genes:
			continue
		used[ensamble_id] = 1
		read_counts = float(data[index])
		count_vec[genes[ensamble_id]] = read_counts 
	f.close()
	if len(used) != len(genes):
		print('assumption erroror!!')
	return count_vec

def get_gtex_gene_counts(file_name, genes, sample_name):
	short_sample_name = sample_name.split('_')[0]
	count = 0
	used = {}
	count_vec = np.zeros(len(genes))
	f = open(file_name)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if count < 2:
			count = count + 1
			continue
		if count == 2:
			for i, ele in enumerate(data):
				if ele == short_sample_name:
					index = i
			count = count + 1
			continue
		ensamble_id = data[0].split('.')[0]
		if ensamble_id not in genes:
			continue
		used[ensamble_id] = 1
		read_counts = float(data[index])
		count_vec[genes[ensamble_id]] = read_counts 
	f.close()
	if len(used) != len(genes):
		print('assumption erroror!!')
	return count_vec

def get_ipsc_time_gene_counts(file_name, genes, converter):
	f = gzip.open(file_name)
	head_count = 0
	used = {}
	count_vec = np.zeros(len(genes))
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		identifier = data[0] + '_' + data[7] + '_' + data[8]
		ensamble_id = converter[identifier]
		if ensamble_id not in genes:
			continue
		used[ensamble_id] = 1
		read_counts = float(data[-2])
		count_vec[genes[ensamble_id]] = read_counts 
	f.close()
	if len(used) != len(genes):
		print('assumption erroror!!')
	return count_vec


def fill_in_expression_values(genes, expression_file_cross_data_sets, target_region_to_ensamble):
	gene_counts = {}
	f = open(expression_file_cross_data_sets)
	loop_counter = 0
	for line in f:
		loop_counter = loop_counter + 1
		print(loop_counter)
		line = line.rstrip()
		data = line.split()
		if data[2] == 'ipsc_ca':
			sample_name = data[0]
			counts = get_ipsc_ca_gene_counts(data[1],genes)
			gene_counts[sample_name] = counts
		if data[2] == 'ipsc':
			sample_name = data[0]
			counts = get_ipsc_gene_counts(data[1], genes, sample_name)
			gene_counts[sample_name] = counts
		if data[2] == 'GTEx':
			sample_name = data[0]
			counts = get_gtex_gene_counts(data[1], genes, sample_name)
			gene_counts[sample_name] = counts
		if data[2] == 'ipsc_time':
			sample_name = data[0]
			counts = get_ipsc_time_gene_counts(data[1], genes, target_region_to_ensamble)
			gene_counts[sample_name] = counts
	return gene_counts



def get_mapping_from_target_region_to_ensamble(target_region_file):
	f = open(target_region_file)
	mapper = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		ensamble_id = data[6].split('_')[2]
		identifier = data[0] + '_' + data[7] + '_' + data[8]
		mapper[identifier] = ensamble_id
	return mapper


expression_file_cross_data_sets = sys.argv[1]
target_region_file = sys.argv[2]
visualize_total_expression_dir = sys.argv[3]


target_region_to_ensamble = get_mapping_from_target_region_to_ensamble(target_region_file)

# Get list of genes found in all data sets
genes = get_list_of_genes_found_in_all_data_sets(expression_file_cross_data_sets, target_region_to_ensamble)

ordered_genes = sorted(genes.keys())
gene_dict = {}
for i,geney in enumerate(ordered_genes):
	gene_dict[geney] = i

gene_counts = fill_in_expression_values(gene_dict, expression_file_cross_data_sets, target_region_to_ensamble)

ordered_samples = sorted(gene_counts.keys())

t = open(visualize_total_expression_dir + 'raw_counts_cmp_data_sets.txt', 'w')
t.write('ensamble_id' + '\t' + '\t'.join(np.asarray(ordered_samples)) + '\n')
for index, gene_name in enumerate(ordered_genes):
	t.write(gene_name)
	for sample in ordered_samples:
		county = gene_counts[sample][index]
		t.write('\t' + str(county))
	t.write('\n')
t.close()
