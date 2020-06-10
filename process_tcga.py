import pandas as pd
import re
import os

#####################################################################################
######### PART 1: might wanna use multiple threads for this part ####################
#####################################################################################

# access the information of novel transcripts with unique junctions
novel_trans_uniq_junc = pd.read_csv('D:/MCGDYY/ont_project/novel_trans_uniq_junc.csv')
novel_trans_uniq_junc['count'] = 0
for i in novel_trans_uniq_junc.index:
	start = novel_trans_uniq_junc.loc[i, 'start']
	end = novel_trans_uniq_junc.loc[i, 'end']
	novel_trans_uniq_junc.loc[i, 'start1'] = start - 10
	novel_trans_uniq_junc.loc[i, 'start2'] = start + 10
	novel_trans_uniq_junc.loc[i, 'end1'] = end - 10
	novel_trans_uniq_junc.loc[i, 'end2'] = end + 10

# create a summary table for novel transcripts with unqiue junctions
list_novel_trans_uniq_junc = set(novel_trans_uniq_junc['transcript'])
cols = []
for i in list_novel_trans_uniq_junc:
	cols.append(i + '_N')
	cols.append(i + '_T')
	cols.append(i + '_value')
	cols.append(i + '_status')
sum_table = pd.DataFrame(columns = cols)
main_df = pd.DataFrame(columns = cols)

# iterate through TCGA junction files and quantify the novel transcripts
main_path = 'D:\\MCGDYY\\ont_project\\TCGA\\test\\'
for junc_file in os.listdir(main_path):
	sample_name = junc_file[:-5]
	file_path = main_path + junc_file

	with open(file_path, 'r') as lc:
		for line in lc:
			line = line.strip('\n')
			chromosome = line.split('\t')[0]
			pos1 = int(line.split('\t')[1])
			pos2 = int(line.split('\t')[2])
			count = int(line.split('\t')[4])
			strand = line.split('\t')[5]

			# +- 10 bp is deemed as the same junction
			for i in novel_trans_uniq_junc.index:
				if chromosome == novel_trans_uniq_junc.loc[i, 'chromosome'] and \
				strand == novel_trans_uniq_junc.loc[i, 'strand'] and \
				novel_trans_uniq_junc.loc[i, 'start1'] <= pos1 <= novel_trans_uniq_junc.loc[i, 'start2'] and \
				novel_trans_uniq_junc.loc[i, 'end1'] <= pos2 <= novel_trans_uniq_junc.loc[i, 'end2']:
					novel_trans_uniq_junc.loc[i, 'count'] += count

		# novel_trans_uniq_junc.to_csv('D:/MCGDYY/ont_project/TCGA/per_intron_counts/sample_quan.csv', index=False)

		# average counts of all junctions in each transcript
		for k in list_novel_trans_uniq_junc:
			all_counts = novel_trans_uniq_junc[novel_trans_uniq_junc['transcript'] == k]['count']
			mean_counts = all_counts.mean()
			# append counts to the sum_table
			sum_table.loc[sample_name, k + '_T'] = mean_counts

		# sum_table.to_csv('D:/MCGDYY/ont_project/TCGA/sum_counts/sample_sum.csv')

#####################################
######### PART 2 ####################
#####################################

# concat all sample_sum.csv files

main_dir = 'D:\\MCGDYY\\ont_project\\TCGA\\sum_counts\\'
for i in os.listdir(main_dir):
	full_path = main_dir + i
	df = pd.read_csv(full_path, index_col = 0)
	main_df = pd.concat((main_df, df), axis = 0, sort = False)

main_df.to_csv('D:/MCGDYY/ont_project/TCGA/TCGA_junc_counts.csv')


# normalize the counts to sequencing depth

# access the sequencing depth information
depth_info = {}
dep_file = open('D:/MCGDYY/ont_project/TCGA/depth.txt')
for i in dep_file:
	name = i.strip('\n').split('\t')[0]
	depth = int(i.strip('\n').split('\t')[1])
	depth_info[name] = depth

# normalization
for trans in list_novel_trans_uniq_junc:
	col_name = trans + '_T'
	for sample in main_df.index:
		main_df.loc[sample, col_name] = main_df.loc[sample, col_name] * 1000000000 / depth_info[sample]

	median = main_df[col_name].median()
	for sample in main_df.index:
		if main_df.loc[sample, col_name] > median:
			main_df.loc[sample, trans + '_status'] = 'H'
		else:
			main_df.loc[sample, trans + '_status'] = 'L'

#################################################################################################################################
########################################### PART 3: novel transcripts with retained introns #####################################
#################################################################################################################################

novel_trans_RI = pd.read_csv('D:/MCGDYY/ont_project/TCGA/TCGA_RI_counts.txt', sep = '\t', comment = '#', index_col = 0)

cols = []
for i in novel_trans_RI.index:
	cols.append(i + '_N')
	cols.append(i + '_T')
	cols.append(i + '_value')
	cols.append(i + '_status')
sum_table = pd.DataFrame(columns = cols)

all_depth = pd.read_csv('D:/MCGDYY/ont_project/TCGA/depth.txt', sep = '\t', index_col = 0, names = ['depth'])
for i in novel_trans_RI.columns[5:]:
	patient = re.findall("files/(.*?).bam", i)[0]
	depth = all_depth.loc[patient, 'depth']
	for k in novel_trans_RI.index:
		sum_table.loc[patient, k + '_T'] = novel_trans_RI.loc[k, i] * 1000000000 / depth

# define H or L
for i in sum_table.columns:
	if '_T' in i:
		trans = i.split('_')[0]
		if trans not in list_novel_trans_uniq_junc:	# if has both uniq junc and retained introns, choose uniq junc
			median = sum_table[i].median()

			for k in sum_table.index:
				counts = sum_table.loc[k, i]
				if counts >= median:
					sum_table.loc[k, trans + '_status'] = 'H'
				else:
					sum_table.loc[k, trans + '_status'] = 'L'

sum_table = sum_table.fillna(0)

########################################################################################
########################################### PART 4 #####################################
########################################################################################

final = pd.concat((main_df, sum_table), axis = 1, sort = False)
new_index = [i[:-4] for i in final.index]
final.index = new_index
clinical = pd.read_csv('D:/MCGDYY/ont_project/TCGA/TCGA_clinical_info.csv', index_col = 0)
final = pd.concat((clinical, final), axis = 1, sort = False, join = 'inner')

final.to_csv('D:\\MCGDYY\\ont_project\\TCGA\\TCGA_novel_trans_uniq_junc_RI.csv')