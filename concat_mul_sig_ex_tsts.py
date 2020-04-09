# takes results from summarize_novel_tsts_by_uniq_junc.py (sum_norm_novel_tsts_uniq_junc.csv) and 
# quan_single_exon_trans.py (sum_novel_single_exon_tsts.csv) to concat multiple and single -exon TSTs, and define status of each TST

import pandas as pd
import re
import statistics

include_sig = False	# decide here if want to include single-exon transcripts in analysis
surv_info = pd.read_csv('D:\\MCGDYY\\ont_project\\clinical_info.csv', index_col = 0)
mul = pd.read_csv('D:\\MCGDYY\\ont_project\\quantification\\sum_norm_novel_mul_trans_uniq_junc.csv', index_col = 0)
sig = pd.read_csv('D:\\MCGDYY\\ont_project\\quantification\\sum_novel_single_exon_trans.csv', index_col = 0)


def cal_freq(df):	# df is either mul or sig
	sum_table = pd.DataFrame(columns = ['frequency'])

	for col in df.columns:
		if re.search('_status', col) is not None:
			trans = col[:-7]
			median = df[trans + '_T'].median()
			mean = df[trans + '_T'].mean()
			q1 = df[trans + '_T'].quantile(0.25)
			q3 = df[trans + '_T'].quantile(0.75)
			# initialize the frequency with 0
			sum_table.loc[trans, 'frequency'] = 0

			for sample in df.index:
				N_exp = df.loc[sample, trans + '_N']
				T_exp = df.loc[sample, trans + '_T']
				df.loc[sample, trans + '_status'] = 0

				if N_exp == 0 and T_exp > 0:
					sum_table.loc[trans, 'frequency'] += 1
					# df.loc[sample, trans + '_status'] = 1
				if N_exp > 0 and T_exp/N_exp > 10:
					sum_table.loc[trans, 'frequency'] += 1
					# df.loc[sample, trans + '_status'] = 1

				# define the status
				if T_exp > q3:
					df.loc[sample, trans + '_status'] = 'H'
				elif q1 < T_exp <= q3:
					df.loc[sample, trans + '_status'] = 'M'
				else:
					df.loc[sample, trans + '_status'] = 'L'

	return sum_table


mul_tsts = cal_freq(mul)
sig_tsts = cal_freq(sig)

# filter
sample_size = mul.shape[0]
mul_tsts_filtered = mul_tsts[mul_tsts['frequency'] > sample_size * 0.3]
sig_tsts_filtered = sig_tsts[sig_tsts['frequency'] > sample_size * 0.3]

# merge tst info with clinical info
merged = surv_info
for tst in mul_tsts_filtered.index:
	merged = pd.concat((merged, mul[tst + '_N'], mul[tst + '_T'], mul[tst + '_status']), axis = 1, join = 'inner')

if include_sig:
	for tst in sig_tsts_filtered.index:
		merged = pd.concat((merged, sig[tst + '_N'], sig[tst + '_T'], sig[tst + '_status']), axis = 1, join = 'inner')	

merged.to_csv('D:/MCGDYY/ont_project/tsts_10x_30_quantiles.csv')
