import numpy as np 
import os
import sys
import pdb
import gzip


def extract_ensamble_ids(input_file, cov_num, num_genes):
	f = open(input_file)
	head_count = 0
	all_genes = {}
	gene_list = []
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0].split('_')[0]
		gene_loading = float(data[cov_num])
		all_genes[gene_id] = 1
		gene_list.append((gene_id, abs(gene_loading)))
	f.close()

	gene_list_sorted = sorted(gene_list, key=lambda tup: tup[1], reverse=True)
	counter = 0
	hit_genes = {}
	for tupler in gene_list_sorted:
		if len(hit_genes) == num_genes:
			continue
		geney = tupler[0]
		hit_genes[geney] = 1
	return hit_genes, all_genes

def extract_ensamble_ids_at_late_time(input_file, cov_num, num_genes):
	f = open(input_file)
	head_count = 0
	all_genes = {}
	gene_list = []
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0].split('_')[0]
		time_step = int(data[0].split('_')[1])
		gene_loading = float(data[cov_num])
		all_genes[gene_id] = 1
		if time_step >= 7:
			gene_list.append((gene_id, abs(gene_loading)))
	f.close()

	gene_list_sorted = sorted(gene_list, key=lambda tup: tup[1], reverse=True)
	counter = 0
	hit_genes = {}
	for tupler in gene_list_sorted:
		if len(hit_genes) == num_genes:
			continue
		geney = tupler[0]
		hit_genes[geney] = 1
	return hit_genes, all_genes

def convert_from_ensamble_to_gene_symbol(ensamble_hits, gencode_file):
    f = gzip.open(gencode_file)
    gene_symbol_hits = []
    for line in f:
        line = line.decode('utf-8').rstrip()
        data = line.split()
        if line.startswith('#'):
            continue
        line_ensamble_id = data[9].split('"')[1].split('.')[0]
        line_gene_symbol = data[17].split('"')[1]
        if line_ensamble_id in ensamble_hits:
            gene_symbol_hits.append(line_gene_symbol)
    return np.unique(gene_symbol_hits)


def print_array(file_name, array):
    t = open(file_name,'w')
    for ele in array:
        t.write(ele + '\n')
    t.close()

def sort_gsea(save_file, new_save_file):
    f = open(save_file)
    t = open(new_save_file,'w')
    pairs = []
    for i,line in enumerate(f):
        line = line.rstrip()
        data = line.split()
        if i < 4:
            continue
        pvalue = float(data[6])
        pairs.append((pvalue, line))
    sorted_pairs = sorted(pairs, key=lambda x: x[0])
    for pair in sorted_pairs:
        liner = pair[1]
        t.write(liner + '\n')
    t.close()


def run_gene_set_enrichment_on_cell_line_pc(num_genes, input_file, cov_num, output_file_root, gencode_file, gsea_data_dir):
	ensamble_hits, ensamble_background = extract_ensamble_ids(input_file, cov_num, num_genes)

	gene_symbol_hits = convert_from_ensamble_to_gene_symbol(ensamble_hits, gencode_file)

	gene_symbol_background = convert_from_ensamble_to_gene_symbol(ensamble_background, gencode_file)

	hits_file = output_file_root + '_hit_genes.txt'

	background_file = output_file_root + '_background_genes.txt'


	print_array(hits_file, gene_symbol_hits)
	print_array(background_file, gene_symbol_background)

	genesets = ['h.all.v5.1.symbols.gmt.txt', 'c2.cp.biocarta.v5.1.symbols.gmt.txt', 'c2.cp.kegg.v5.1.symbols.gmt.txt', 'c2.cp.reactome.v5.1.symbols.gmt.txt']

	names = ['hallmark', 'biocarta', 'kegg', 'reactome']
	for i, val in enumerate(genesets):
		name = names[i]
		geneset_file = gsea_data_dir + val
		save_file = output_file_root + '_' + name + '_gsea_output.txt'
		os.system('gsea ' + hits_file + ' ' + background_file + ' ' + geneset_file + ' ' + save_file)

		new_save_file = output_file_root + '_' + name + '_gsea_sorted_output.txt'
		sort_gsea(save_file, new_save_file)
		# Remove un-sorted file
		os.system('rm ' + save_file)

	# Remove some unnecessary files
	os.system('rm ' + hits_file)
	os.system('rm ' + background_file)

def run_gene_set_enrichment_on_late_time_cell_line_pc(num_genes, input_file, cov_num, output_file_root, gencode_file, gsea_data_dir):
	ensamble_hits, ensamble_background = extract_ensamble_ids_at_late_time(input_file, cov_num, num_genes)

	gene_symbol_hits = convert_from_ensamble_to_gene_symbol(ensamble_hits, gencode_file)

	gene_symbol_background = convert_from_ensamble_to_gene_symbol(ensamble_background, gencode_file)

	hits_file = output_file_root + '_hit_genes.txt'

	background_file = output_file_root + '_background_genes.txt'


	print_array(hits_file, gene_symbol_hits)
	print_array(background_file, gene_symbol_background)

	genesets = ['h.all.v5.1.symbols.gmt.txt', 'c2.cp.biocarta.v5.1.symbols.gmt.txt', 'c2.cp.kegg.v5.1.symbols.gmt.txt', 'c2.cp.reactome.v5.1.symbols.gmt.txt']

	names = ['hallmark', 'biocarta', 'kegg', 'reactome']
	for i, val in enumerate(genesets):
		name = names[i]
		geneset_file = gsea_data_dir + val
		save_file = output_file_root + '_' + name + '_gsea_output.txt'
		os.system('gsea ' + hits_file + ' ' + background_file + ' ' + geneset_file + ' ' + save_file)

		new_save_file = output_file_root + '_' + name + '_gsea_sorted_output.txt'
		sort_gsea(save_file, new_save_file)
		# Remove un-sorted file
		os.system('rm ' + save_file)

	# Remove some unnecessary files
	os.system('rm ' + hits_file)
	os.system('rm ' + background_file)





covariate_dir = sys.argv[1]
gencode_file = sys.argv[2]
gsea_data_dir = sys.argv[3]
num_genes = 200
input_file = covariate_dir + 'cell_line_ignore_missing_gene_loadings_9.txt'


for cov_num in range(9):
	print(cov_num)
	cov_num_temp = cov_num + 1
	output_file_root = covariate_dir + 'gene_set_enrichment_on_late_time_cell_line_pc_' + str(cov_num_temp)
	run_gene_set_enrichment_on_late_time_cell_line_pc(num_genes, input_file, cov_num_temp, output_file_root, gencode_file, gsea_data_dir)


for cov_num in range(9):
	print(cov_num)
	cov_num_temp = cov_num + 1
	output_file_root = covariate_dir + 'gene_set_enrichment_on_cell_line_pc_' + str(cov_num_temp)
	run_gene_set_enrichment_on_cell_line_pc(num_genes, input_file, cov_num_temp, output_file_root, gencode_file, gsea_data_dir)
