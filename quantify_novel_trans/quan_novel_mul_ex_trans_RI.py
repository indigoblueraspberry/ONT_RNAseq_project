import pandas as pd
import re
import os

novel_trans_RI = pd.read_csv('D:/MCGDYY/ont_project/IR_counts.txt', sep = '\t', index_col = 0, header = 0)

cols = []
for i in novel_trans_RI.index:
	cols.append(i + '_N')
	cols.append(i + '_T')
	cols.append(i + '_value')
	cols.append(i + '_status')
sum_table = pd.DataFrame(columns = cols)

for i in novel_trans_RI.columns:
	patient = re.findall("/(.*?)-", i)[0][-4:]
	kind = re.findall("-(.*?)_", i)[0]
	
	if patient != '1438' and patient != '1483':	# these two patients have missing data
		for k in novel_trans_RI.index:
			counts = novel_trans_RI.loc[k, i]
			if kind == '1':
				sum_table.loc[patient, k + '_T'] = counts
			else:
				sum_table.loc[patient, k + '_N'] = counts

# normalization

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

for sample_id in sum_table.index:
		for col in sum_table.columns:
			if re.search('_N', col):
				depth = all_depth[sample_id + '-2']
				sum_table.loc[sample_id, col] = sum_table.loc[sample_id, col] * 1000000000 / depth
			if re.search('_T', col):
				depth = all_depth[sample_id + '-1']
				sum_table.loc[sample_id, col] = sum_table.loc[sample_id, col] * 1000000000 / depth

sum_table = sum_table.fillna(0)
sum_table.to_csv('D:/MCGDYY/ont_project/quantification/sum_norm_novel_mul_trans_RI.csv')