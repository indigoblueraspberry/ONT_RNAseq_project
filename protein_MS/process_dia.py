import pandas as pd
import re

novel_mul_ex_trans = []
df = open('D:\\MCGDYY\\ont_project\\lists\\list_of_all_novel_mul_ex_trans.txt', 'r')
for i in df:
	i = i.strip('\n')
	novel_mul_ex_trans.append(i)

# batch1
df = pd.read_csv('D:\\MCGDYY\\ont_project\\MS\\DIA_res\\batch1_processed.csv')
sum_df = pd.DataFrame()
for i in df.columns:
	if 'PG' in i:
		sample = i.split('_')[1][:4]
		kind = i.split('_')[1][4]

		for k in df.index:
			trans = df.loc[k, 'Protein']
			if trans in novel_mul_ex_trans:
				exp = df.loc[k, i]
				if '?' not in exp:
					sum_df.loc[sample, trans + '_pro_' + kind] = 1

# batch2
df = pd.read_csv('D:\\MCGDYY\\ont_project\\MS\\DIA_res\\batch2_processed.csv')
for i in df.columns:
	if 'Overall' in i:
		sample = i.split('_')[0][:4]
		kind = i.split('_')[0][4]

		for k in df.index:
			trans = df.loc[k, 'Protein']
			if trans in novel_mul_ex_trans:
				exp = df.loc[k, i]
				if '?' not in exp:
					sum_df.loc[sample, trans + '_pro_' + kind] = 1

sum_df = sum_df.fillna(0)
sum_df = sum_df.sort_index(axis = 0)
sum_df = sum_df.sort_index(axis = 1)

# remove unbalanced groups
for col in sum_df.columns:
	group1 = 0
	group2 = 0
	for sample in sum_df.index:
		status = sum_df.loc[sample, col]
		if status == 0:
			group1 += 1
		else:
			group2 += 1

	if min(group1, group2) / max(group1, group2) < 0.5:
		sum_df = sum_df.drop(col, axis = 1)

# clinical info
surv_info = pd.read_csv('D:\\MCGDYY\\ont_project\\clinical_info.csv', index_col = 0)
surv_info.index = surv_info.index.astype(str)
sum_df = pd.concat((surv_info, sum_df), axis = 1, join = 'inner')

sum_df.to_csv('D:\\MCGDYY\\ont_project\\MS\\DIA_res\\DIA_sum.csv', index = True)

