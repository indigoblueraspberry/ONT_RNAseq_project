# differentially expressed analysis on novel transcripts

import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

exp_df = pd.read_csv('D:\\MCGDYY\\ont_project\\quantification\\t_exp_median.csv')
sum_df = pd.DataFrame(columns = ['log2FC', 'wilcoxon_p', 'adj_p'])

for i in exp_df.columns:
	if '_value' in i:
		trans_name = i.split('_')[0]
		log2FC = exp_df[i].mean()
		sum_df.loc[trans_name, 'log2FC'] = log2FC
		wc_p = stats.wilcoxon(exp_df[i])[1]
		sum_df.loc[trans_name, 'wilcoxon_p'] = wc_p
# multiple comparison correction
fdr_bh = multipletests(sum_df['wilcoxon_p'], method='fdr_bh')[1]
sum_df['adj_p'] = fdr_bh
DE_trans = sum_df[(abs(sum_df['log2FC']) > 0.58) & (sum_df['adj_p'] < 0.05)]
DE_list = DE_trans.index

DE_df = exp_df
for i in exp_df.columns[15:]:
	trans_name = i.split('_')[0]
	if trans_name not in DE_list:
		DE_df = DE_df.drop(columns = [i])

DE_df.to_csv('D:\\MCGDYY\\ont_project\\quantification\\t_exp_median_DEonly.csv', index = False)