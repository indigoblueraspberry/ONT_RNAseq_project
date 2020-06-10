# Following the results (*_uniq_junc_reads.txt) of get_uniq_junc_reads_of_novel_mul_ex_trans.py, 
# this script summarizes and normalizes all novel multi-exon transcripts with unique junctions in all samples.


import pandas as pd
import re
import os


# get the sequencing depth for all samples by processing STAR Log.final.out files

all_depth = {}
files_dir = 'D:\\MCGDYY\\ont_project\\star_stats\\'
for file in os.listdir(files_dir):
	if re.search('Log.final.out', file) is not None:
		sample = file[:6]
		file_full_dir = files_dir + file
		df = open(file_full_dir, 'r')
		for i, line in enumerate(df):
			if i == 5:
				line = line.strip('\n')
				depth = line.split('\t')[1]
				all_depth[sample] = int(depth)


# initialize a table with columns of all novel transcripts with unique junctions

uniq_junc_novel_trans = pd.read_csv('D:\\MCGDYY\\ont_project\\lists\\list_of_novel_trans_with_uniq_junc.csv')
uniq_junc_novel_trans = uniq_junc_novel_trans['transcript']
cols = []
for i in uniq_junc_novel_trans:
	cols.append(i + '_N')
	cols.append(i + '_T')
	cols.append(i + '_value')
	cols.append(i + '_status')
sum_table = pd.DataFrame(columns = cols)


# take the outputs from get_uniq_junc_reads_of_novel_mul_ex_trans.py,
# for a novel transcript with multiple unique junctions, average all the unique junction reads

res_dir = 'D:\\MCGDYY\\ont_project\\novel_trans_with_uniq_junc\\'
for sample in os.listdir(res_dir):
	sample_path = res_dir + sample
	sample_id = sample[:4] # patient ID
	sample_type = sample[4:6] # -1 or -2
	sample_df = pd.read_csv(sample_path, header = 0, sep = '\t')

	trans_counter = {}
	for i in range(len(sample_df)):
		transcript = sample_df.loc[i, 'transcript']
		uniq_reads = sample_df.loc[i, 'uniq_reads']
			
		if transcript not in trans_counter:	# first time to be discovered
			trans_counter[transcript] = 1

			if sample_type == '-1':
				sum_table.loc[sample_id, transcript + '_T'] = uniq_reads
			else:
				sum_table.loc[sample_id, transcript + '_N'] = uniq_reads

		else:	# the TST has more than one unique junctions
			trans_counter[transcript] += 1

			if sample_type == '-1':
				sum_table.loc[sample_id, transcript + '_T'] = \
				(sum_table.loc[sample_id, transcript + '_T'] + uniq_reads) / trans_counter[transcript]
			else:
				sum_table.loc[sample_id, transcript + '_N'] = \
				(sum_table.loc[sample_id, transcript + '_N'] + uniq_reads) / trans_counter[transcript]

sum_table = sum_table.fillna(value = 0)
# sum_table.to_csv('D:\\MCGDYY\\ont_project\\quantification\\sum_novel_mul_trans_uniq_junc.csv')


# normalize read counts in sum_table for total sequencing depth

for sample_id in sum_table.index:
	for col in sum_table.columns:

		if re.search('_N', col):
			depth = all_depth[sample_id + '-2']
			sum_table.loc[sample_id, col] = sum_table.loc[sample_id, col] * 1000000000 / depth
		if re.search('_T', col):
			depth = all_depth[sample_id + '-1']
			sum_table.loc[sample_id, col] = sum_table.loc[sample_id, col] * 1000000000 / depth


sum_table.to_csv('D:\\MCGDYY\\ont_project\\quantification\\sum_norm_novel_mul_trans_uniq_junc.csv')