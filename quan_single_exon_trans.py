# Used to quantify single-exon transcripts based on the results from rsem

import pandas as pd
import os
import re

# create a table of novel single-exon transcripts
cols = []
novel_sig_ex_list = []
sig_ex_file = open('D:\\MCGDYY\\ont_project\\gtf_files\\sig_ex.gtf', 'r')
for line in sig_ex_file:
	info = line.split('\t')[-1]
	transcript = re.findall(r"transcript_id \"(.*?)\";", info)[0]
	if re.search('ENST', transcript) is None:
		novel_sig_ex_list.append(transcript)
		cols.append(transcript + '_N')
		cols.append(transcript + '_T')
		cols.append(transcript + '_value')
		cols.append(transcript + '_status')
sum_table = pd.DataFrame(columns = cols)

# create a function to get quantification information from the rsem files
def get_quan(file_path, sample_type):	# sample_type = N or T
	df = pd.read_csv(file_path, sep = '\t', header = 0, index_col = 0)
	for i in novel_sig_ex_list:
		quan = df.loc[i, 'TPM']
		sum_table.loc[sample_id, i + '_' + sample_type] = quan

# iterate through rsem quantification results
rsem_dir = 'D:\\MCGDYY\\ont_project\\rsem\\'
for file in os.listdir(rsem_dir):
	file_path = rsem_dir + file
	sample_id = file[:4]

	if re.search('-2.isoforms.results', file):
		get_quan(file_path, 'N')
	if re.search('-1.isoforms.results', file):
		get_quan(file_path, 'T')

	print(file)

sum_table.to_csv('D:\\MCGDYY\\ont_project\\quantification\\sum_novel_single_exon_trans.csv')