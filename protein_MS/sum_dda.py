import pandas as pd
import os
import re

# annotate the novel peptides identified by MS

uniport = open('D:\\MCGDYY\\ont_project\\MS\\uniport.fa', 'r').read()
novel_ref_all = open('D:\\MCGDYY\\ont_project\\MS\\novel_peptides.fa', 'r').read()

sum_table = pd.DataFrame()
dda_dir = 'D:\\MCGDYY\\ont_project\\MS\\new_DDA\\res_pep\\'
for file_name in os.listdir(dda_dir):
	sample = file_name.split('.')[0]
	print(sample)
	full_path = dda_dir + file_name
	df = pd.read_excel(full_path)
	for i in df.index:
		trans_name = df.loc[i, r'Master Protein Accessions']
		if trans_name == 'sp' or pd.isnull(trans_name):	# look for novel transcripts
			seq = df.loc[i, 'Sequence']

			search_uni = re.findall(seq, uniport)
			if len(search_uni) == 0:
				search_nov_ref = re.findall(seq, novel_ref_all)
				if len(search_nov_ref) == 1:

					with open('D:\\MCGDYY\\ont_project\\MS\\novel_peptides.fa', 'r') as novel_ref:
						while True:
							title = novel_ref.readline()
							if not title:
								break
							ref_seq = novel_ref.readline()
							if re.search(seq, ref_seq) is not None:
								real_name = title.strip('\n')[1:]
								sum_table.loc[real_name, sample] = 1
								df.loc[i, r'Master Protein Accessions'] = real_name
								break
	# df.to_csv('D:\\MCGDYY\\ont_project\\MS\\new_DDA\\res_pep_anno\\' + sample + '.csv', index = None)

# sum_table.to_csv('D:\\MCGDYY\\ont_project\\MS\\new_DDA\\sum_dda.csv', index = 1)


######################## refine the sum_table

# build a table of strand information
strand_table = pd.DataFrame(columns = ['strand'])
nov_list = pd.read_csv('D:\\MCGDYY\\ont_project\\gtf_files\\novel_mul_ex_trans_exons.gtf', comment = '#', header = None, sep = '\t')
for i in nov_list.index:
	info = nov_list.iloc[i, 8]
	trans = re.findall(r'transcript_id \"(.*?)\"', info)[0]
	if trans not in strand_table.index:
		strand = nov_list.iloc[i, 6]
		strand_table.loc[trans, 'strand'] = strand

# remove the FPs
for i in sum_table.index:
	trans = i.split('_')[0]
	strand = i.split('_')[2]
	if strand == 'f':
		strand = '+'
	else:
		strand = '-'

	try:
		real_strand = strand_table.loc[trans, 'strand']
		if strand != real_strand:
			sum_table = sum_table.drop(index = i)
	except KeyError:
		sum_table = sum_table.drop(index = i)
sum_table = sum_table.fillna(0)
# sum_table.to_csv('D:\\MCGDYY\\ont_project\\MS\\new_DDA\\sum_dda.csv', index = 1)

# merge fragments into one peptide
refined = pd.DataFrame()
for i in sum_table.index:
	trans = i.split('_')[0]
	for k in sum_table.columns:
		value = sum_table.loc[i, k]
		if value == 1:
			refined.loc[trans, k] = 1
refined = refined.fillna(0)
# refined.to_csv('D:\\MCGDYY\\ont_project\\MS\\new_DDA\\refined_dda.csv', index = 1)