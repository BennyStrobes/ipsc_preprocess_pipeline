import numpy as np
import os
import sys
import pdb

# Create dictionary-list of all cell lines
def extract_list_of_cell_lines(lane_design_round_1_file, lane_design_round_2_file):
    cell_lines = {}  # Initialize list
    #  Loop through first lane design file
    f = open(lane_design_round_1_file)
    for line in f:
        line = line.rstrip()
        data = line.split(',')
        # Loop through every element in the row
        for ele in data:
            if len(ele) > 2:  # If element is a string of length greater than 2, then it is a cell line.
                cell_lines[ele] = 1  # Add to list of cell lines
    f.close()
    #  Loop through second lane design file
    f = open(lane_design_round_2_file)
    for line in f:
        line = line.rstrip()
        data = line.split(',')
        # Loop through every element in the row
        for ele in data:
            if len(ele) > 2:  # If element is a string of length greater than 2, then it is a cell line.
                cell_lines[ele] = 1  # Add to list of cell lines

    # Quick sanity check. Check to make sure we got 11 unique cell lines:
    if len(cell_lines) != 11:
        print('ERROR: Incorrect number of cell lines extracted')
        pdb.set_trace()
    return cell_lines

# Extract array of tuples from fastq_round_1_input_dir where each tuple represents (lane,sample_num) corresponding to this sample_id
def lane_sample_num_extraction_round_1(lane_design_round_1_file, cell_line, time_step):
    f = open(lane_design_round_1_file)
    sample_num = 1  # Index lines of lane_design file
    output_list = []  # Output list
    for line in f:
        line = line.rstrip() 
        data = line.split(',')
        for index in range(8): # Loop through 8 lanes in file
            cell_line_i = data[index*2]  # cell_line in ith lane
            time_step_i = data[index*2 + 1]  # time step of ith lane
            lane_num = index + 1
            if cell_line == cell_line_i and time_step_i == time_step:  # if this lane,sample_num contains our cellLine_timeStep
                output_list.append((lane_num,sample_num))
        sample_num = sample_num + 1    
    return output_list

# Create dictionary mapping from every sample id (cellLine_timeStep) to an array where each element of the array is the the complete path of one of the fastq files corresponding to that sample 
def create_sample_id_to_fastqs_mapping(mapping, sample_id, cell_line, time_step, fastq_round_1_input_dir, fastq_round_2_input_dir, lane_design_round_1_file, lane_design_round_2_file):
    # Extract array of tuples from fastq_round_1_input_dir where each tuple represents (lane,sample_num) corresponding to this sample_id
    lane_sample_num_list1 = lane_sample_num_extraction_round_1(lane_design_round_1_file, cell_line, str(time_step))
    match_list_1 = False
    # now loop through all fastq files in fastq_round_1_input_dir, and add path to mapping dictionary if fastq comes from this sample
    for fastq_file in os.listdir(fastq_round_1_input_dir):
        if fastq_file.endswith('fastq.gz') == False:  # skip files that don't end with fastq.gz
            continue
        if fastq_file.startswith('Undetermined'):  # Skip the "undetermined fastq files"
            continue
        lane = fastq_file[9:10]  # Extract lane of current fasq_file
        sample_num = fastq_file[11:13].split('_')[0]  # Extract sample_num of current fastq file
        fastq_tuple = (int(lane),int(sample_num))  # Put lane and sample_num into a tuple data structure
        if fastq_tuple in lane_sample_num_list1:  # Check to see if this fastq file corresponds to our sample id
            match_list_1 = True  # we got a match
            mapping[sample_id].append(fastq_round_1_input_dir + fastq_file)  # Add fastq file to dictionary
    if match_list_1 == False:  # Error checking
        print('Missing hit!!!')
        pdb.set_trace()

    # Now do the same for fastq_round_2_input_dir
    lane_sample_num_list2 = lane_sample_num_extraction_round_1(lane_design_round_2_file, cell_line, str(time_step))
    match_list_2 = False
    # now loop through all fastq files in fastq_round_1_input_dir, and add path to mapping dictionary if fastq comes from this sample
    for fastq_file in os.listdir(fastq_round_2_input_dir):
        if fastq_file.endswith('fastq.gz') == False:  # skip files that don't end with fastq.gz
            continue
        if fastq_file.startswith('Undetermined'):  # Skip the "undetermined fastq files"
            continue
        if fastq_file.startswith('YG-172S'):  # Traditional sequener labeling
            lane = fastq_file[9:10]  # Extract lane of current fasq_file
            sample_num = fastq_file[11:13].split('_')[0]
            fastq_tuple = (int(lane),int(sample_num))
            if fastq_tuple in lane_sample_num_list2:  # we got a match
                match_list_2 = True
                mapping[sample_id].append(fastq_round_2_input_dir + fastq_file)
        elif fastq_file.startswith('YG-RE-172S'):  # Alternative sequencer labeling
            # Sometimes there is a substring 'day', and other times there is a substring 'Day' (Unclear why), but we deal with it here
            if 'day' in fastq_file:
                cell_line_i = fastq_file.split('pl8-')[1].split('day')[0]
                time_step_i = fastq_file.split('day')[1].split('_')[0]
            elif 'Day' in fastq_file:
                cell_line_i = fastq_file.split('pl8-')[1].split('Day')[0]
                time_step_i = fastq_file.split('Day')[1].split('_')[0]
            else:  # Just a sanity check
                print('missing word day')
                pdb.set_trace()
            if cell_line_i + '_' + time_step_i == sample_id:  # we got a match
                match_list_2 = True
                mapping[sample_id].append(fastq_round_2_input_dir + fastq_file)
    return mapping



# Command line args:
fastq_round_1_input_dir = sys.argv[1]  # Directory containing fastq files from first sequencing round
fastq_round_2_input_dir = sys.argv[2]  # Directory containing fastq files from second sequencing round
lane_design_round_1_file = sys.argv[3]  # File containing information for mapping from sequencer name to cellLine_timeStep (sample) notation for the first sequencing round
lane_design_round_2_file = sys.argv[4]  # File containing information for mapping from sequencer name to cellLine_timeStep (sample) notation for the first sequencing round
fastq_dir = sys.argv[5]  # Output file

# NOTE: This whole script is kinda confusing/messy because of the mislabeled samples, contaminated samples, and change in sequencer name notation in the second round.
#  Because of this, this script is definitely not generalizable and will need to be carefully evaluated (to make sure its still correct) if any future changes occur.


# First create dictionary-list of all cell lines
cell_lines = extract_list_of_cell_lines(lane_design_round_1_file, lane_design_round_2_file)
#  Cell line 19128 is contaminated. Throw it out:
cell_lines.pop('19128', None)

# Now create dictionary mapping from every sample id (cellLine_timeStep) to an array where each element of the array is the the complete path of one of the fastq files corresponding to that sample 
# So loop through all sample ids
mapping = {} # Initialize mapping
for cell_line in cell_lines.keys():
    for time_step in range(16):
        if cell_line == '18870' and time_step == 4:  # Don't have this time point
            continue
        if cell_line == '19108' and time_step == 4:  # Don't have this time point
            continue
        if cell_line == '18855' and time_step == 2:  # Don't have this time point
            continue
        if cell_line == '18520' and time_step == 2:  # Don't have this time point
            continue
        sample_id = cell_line + '_' + str(time_step)
        mapping[sample_id] = []  # Initialize sampleID in mapping
        mapping = create_sample_id_to_fastqs_mapping(mapping, sample_id, cell_line, time_step, fastq_round_1_input_dir, fastq_round_2_input_dir, lane_design_round_1_file, lane_design_round_2_file)


# Print mapping to $fasq_dir"fastq_mapping.txt"
t = open(fastq_dir + 'fastq_mapping.txt','w')
for key in sorted(mapping.keys()):
    cell_line = key.split('_')[0]
    time_step = key.split('_')[1]
    if cell_line == '19238':
        word = '18499' + '_' + time_step
    else:
        word = key
    file_name = word + '_merged.fastq.gz'
    t.write(file_name + '\t' + ','.join(mapping[key]) +  '\n')

    # Now cat fastq files into one merged file
    stringer = 'cat '
    for stringy in mapping[key]:
        stringer = stringer + stringy + ' '
    stringer = stringer + '> ' + fastq_dir + file_name
    os.system(stringer)


t.close()