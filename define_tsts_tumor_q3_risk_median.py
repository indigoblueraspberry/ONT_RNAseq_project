# 1: calculate median tumor expression, and category each patient into H: > Q3 and L: <= Q3
# 2: calculate risk score based on the two significant TSTs (from R results), and category patients (H: > median, L: < median)

import pandas as pd
import re
import statistics


surv_info = pd.read_csv('D:\\MCGDYY\\ont_project\\clinical_info.csv', index_col = 0)
mul = pd.read_csv('D:\\MCGDYY\\ont_project\\quantification\\sum_norm_novel_mul_trans_uniq_junc.csv', index_col = 0)


sum_table = pd.DataFrame(columns = ['frequency'])

for col in mul.columns:
	if re.search('_status', col) is not None:
		trans = col[:-7]
		# initialize the frequency with 0
		sum_table.loc[trans, 'frequency'] = 0
		q3 = mul[trans + '_T'].quantile(0.75)

		for sample in mul.index:
			N_exp = mul.loc[sample, trans + '_N']
			T_exp = mul.loc[sample, trans + '_T']

			if T_exp > q3:
				mul.loc[sample, trans + '_status'] = 'H'
			else:
				mul.loc[sample, trans + '_status'] = 'L'

			if N_exp == 0 and T_exp > 0:
				sum_table.loc[trans, 'frequency'] += 1
			if N_exp > 0 and T_exp/N_exp > 10:
				sum_table.loc[trans, 'frequency'] += 1
				
# filter
sample_size = mul.shape[0]
mul_tsts_filtered = sum_table[sum_table['frequency'] > sample_size * 0.3]

# calculate risk score based on 2 TSTs
tst1 = '93d07f14-5252-481c-9988-60b2a9653a5f-1_T'
tst2 = '547bd10a-0658-4b4c-85bd-76c16dd03a31_T'
coef1 = 0.006806
coef2 = 0.016079

risk_df = pd.DataFrame(columns = ['risk_score', 'risk'])
for sample in mul.index:
	T_exp1 = mul.loc[sample, tst1]
	T_exp2 = mul.loc[sample, tst2]
	each_score = T_exp1 * coef1 + T_exp2 * coef2
	risk_df.loc[sample, 'risk_score'] = each_score

median_risk = risk_df['risk_score'].median()
for sample in risk_df.index:
	if risk_df.loc[sample, 'risk_score'] > median_risk:
		risk_df.loc[sample, 'risk'] = 'H'
	else:
		risk_df.loc[sample, 'risk'] = 'L'

# merge tst info with clinical info
merged = surv_info
for tst in mul_tsts_filtered.index:
	merged = pd.concat((merged, mul[tst + '_N'], mul[tst + '_T'], mul[tst + '_value'], mul[tst + '_status']), axis = 1, join = 'inner')

merged = pd.concat((risk_df, merged), axis = 1, join = 'inner')

merged.to_csv('D:/MCGDYY/ont_project/tsts_10x_30_q3_risk_median.csv')
