import pandas as pd
import re
import statistics

########################################## PART 1: novel transcripts with unique junctions

mul = pd.read_csv('D:\\MCGDYY\\ont_project\\quantification\\sum_norm_novel_mul_trans_uniq_junc.csv', index_col = 0)

sum_table = pd.DataFrame(columns = ['frequency'])

for col in mul.columns:
	if re.search('_status', col) is not None:
		trans = col[:-7]
		# initialize the frequency with 0
		sum_table.loc[trans, 'frequency'] = 0
		median = mul[trans + '_T'].median()

		for sample in mul.index:
			N_exp = mul.loc[sample, trans + '_N']
			T_exp = mul.loc[sample, trans + '_T']

			if T_exp > median:
				mul.loc[sample, trans + '_status'] = 'H'
			else:
				mul.loc[sample, trans + '_status'] = 'L'

			if N_exp == 0 and T_exp > 0:
				sum_table.loc[trans, 'frequency'] += 1
			if N_exp > 0 and T_exp/N_exp > 10:
				sum_table.loc[trans, 'frequency'] += 1
				
####################################################### PART 2: novel transcripts that include intronic regions

# filter
sample_size = mul.shape[0]
mul_tsts_filtered = sum_table[sum_table['frequency'] > sample_size * 0.3]

# merge tst info with clinical info
surv_info = pd.read_csv('D:\\MCGDYY\\ont_project\\clinical_info.csv', index_col = 0)
merged = surv_info
for tst in mul_tsts_filtered.index:
	merged = pd.concat((merged, mul[tst + '_N'], mul[tst + '_T'], mul[tst + '_value'], mul[tst + '_status']), axis = 1, join = 'inner')

merged = pd.concat((risk_df, merged), axis = 1, join = 'inner')

merged.to_csv('D:/MCGDYY/ont_project/tsts_10x_30_q3_risk_median.csv')
