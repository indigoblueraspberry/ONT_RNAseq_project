import pandas as pd
import math
import statistics


raw_df = pd.read_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\NovelQuant_res.csv', index_col = 0)
surv_info = pd.read_csv('D:\\MCGDYY\\ont_project\\prognosis\\clinical_info.csv', index_col = 0)
surv_info.index = surv_info.index.astype('str')
final = surv_info

cols = []
for i in raw_df.index:
	cols.append(i + '_N')
	cols.append(i + '_T')
	cols.append(i + '_value')
	cols.append(i + '_status')
sum_table = pd.DataFrame(columns = cols)

for k in raw_df.columns:
	patient = k.split('-')[0]
	kind = k.split('-')[1]
	if patient != '1438' and patient != '1483':	# these two patients have missing data
		for i in raw_df.index:
			exp = raw_df.loc[i, k]
			if kind == '1':
				sum_table.loc[patient, i + '_T'] = exp
			else:
				sum_table.loc[patient, i + '_N'] = exp

for k in sum_table.columns:
	if '_status' in k:
		trans = k[:-7]
		group1 = 0
		group2 = 0
		median = sum_table[trans + '_T'].median()
		for i in sum_table.index:
			N_exp = sum_table.loc[i, trans + '_N']
			T_exp = sum_table.loc[i, trans + '_T']
			sum_table.loc[i, trans + '_value'] = math.log2( ( T_exp + 1 ) / ( N_exp + 1 ))
			if sum_table.loc[i, trans + '_T'] > median:
				sum_table.loc[i, trans + '_status'] = 1
				group1 += 1
			else:
				sum_table.loc[i, trans + '_status'] = 0
				group2 += 1

		if min(group1, group2) / max(group1, group2) >= 0.5:
			final = pd.concat((final, sum_table[trans + '_N'], sum_table[trans + '_T'], sum_table[trans + '_value'], sum_table[trans + '_status']), 
			axis = 1, join = 'inner')

final.to_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\t_exp_median.csv')