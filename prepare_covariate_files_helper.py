import numpy as np 
import os
import sys
import pdb





def get_sample_names(file_name):
    f = open(file_name)
    arr = []
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        arr.append(data[0])
    return np.asarray(arr)



def get_cell_line_pcs(cell_line_pc_file, pc_num):
    f = open(cell_line_pc_file)
    head_count = 0
    dicti = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        dicti[data[0]] = float(data[pc_num])
    return dicti



covariate_dir = sys.argv[1]





cell_line_pc_file = covariate_dir + 'cell_line_ignore_missing_principal_components_9.txt'

global_pc_file = covariate_dir + 'principal_components_10.txt'

output_file = covariate_dir + 'cell_line_ignore_missing_principal_components_9_time_time_per_sample.txt'



sample_names = get_sample_names(global_pc_file)

cell_line_pcs = {}
for pc_num in range(1,10):
    cell_line_pcs[pc_num] = get_cell_line_pcs(cell_line_pc_file, pc_num)

t = open(output_file,'w')
t.write('Sample_id')
for pc_num in range(1,10):
    t.write('\t' + 'PC' + str(pc_num))
t.write('\n')

for sample_name in sample_names:
    t.write(sample_name)
    time_step = float(sample_name.split('_')[1])
    cell_line = sample_name.split('_')[0]
    for pc_num in range(1,10):
        t.write('\t' + str(time_step*cell_line_pcs[pc_num][cell_line]))
    t.write('\n')
t.close()

