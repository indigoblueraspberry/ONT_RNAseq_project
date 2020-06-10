import pandas as pd
import statistics
import math

surv_info = pd.read_csv('D:\\MCGDYY\\ont_project\\prognosis\\clinical_info.csv', index_col = 0)
final = surv_info
	
####################################################### PART 1: novel transcripts that include intronic regions

RI = pd.read_csv('D:\\MCGDYY\\ont_project\\quantification\\sum_norm_novel_mul_trans_RI.csv', index_col = 0)
RI_list = []

for col in RI.columns:
	if '_status' in col:
		trans = col[:-7]
		RI_list.append(trans)
		group1 = 0
		group2 = 0
		median = RI[trans + '_T'].median()

		for sample in RI.index:
			N_exp = RI.loc[sample, trans + '_N']
			T_exp = RI.loc[sample, trans + '_T']
			# calculate log2FC
			RI.loc[sample, trans + '_value'] = math.log2((T_exp + 1) / (N_exp + 1))

			if T_exp > median:
				RI.loc[sample, trans + '_status'] = 1
				group1 += 1
			else:
				RI.loc[sample, trans + '_status'] = 0
				group2 += 1

		if min(group1, group2) / max(group1, group2) >= 0.5:
			final = pd.concat((final, RI[trans + '_N'], RI[trans + '_T'], RI[trans + '_value'], RI[trans + '_status']), 
				axis = 1, join = 'inner')


########################################## PART 2: novel transcripts with unique junctions

junc = pd.read_csv('D:\\MCGDYY\\ont_project\\quantification\\sum_norm_novel_mul_trans_uniq_junc.csv', index_col = 0)

for col in junc.columns:
	if '_status' in col:
		trans = col[:-7]
		if trans not in RI_list:	# if the transcript is already identified by intronic regions, then use RI quantification
			group1 = 0
			group2 = 0
			median = junc[trans + '_T'].median()

			for sample in junc.index:
				N_exp = junc.loc[sample, trans + '_N']
				T_exp = junc.loc[sample, trans + '_T']
				junc.loc[sample, trans + '_value'] = math.log2((T_exp + 1) / (N_exp + 1))

				if T_exp > median:
					junc.loc[sample, trans + '_status'] = 1
					group1 += 1
				else:
					junc.loc[sample, trans + '_status'] = 0
					group2 += 1

			# remove unbalanced groups
			if min(group1, group2) / max(group1, group2) >= 0.5:
				final = pd.concat((final, junc[trans + '_N'], junc[trans + '_T'], junc[trans + '_value'], junc[trans + '_status']), 
				axis = 1, join = 'inner')


final.to_csv('D:/MCGDYY/ont_project/quantification/t_exp_median.csv')

