import pandas as pd
import re
from scipy import stats
from statsmodels.stats.multitest import multipletests


clinical_trans = pd.read_csv('D:\\MCGDYY\\ont_project\\prognosis\\union_os_rfs.csv')
exp_df = pd.read_csv('D:\\MCGDYY\\ont_project\\quantification\\t_exp_median.csv')

for i in clinical_trans.index:
	# rename transcripts
	trans = clinical_trans.loc[i, 'transcript']
	pat = re.compile('trans_(.*?)_status')
	trans = re.match(pat, trans)[1]
	trans = trans.replace('_', '-')
	trans = trans + '_value'
	clinical_trans.loc[i, 'transcript'] = trans

	FC_col = exp_df[trans]
	wc_p = stats.wilcoxon(FC_col)[1]
	clinical_trans.loc[i, 'Wilcoxon_p'] = wc_p
	mean_log2FC = FC_col.mean()
	clinical_trans.loc[i, 'mean_log2FC'] = mean_log2FC

# adjust p-value for multiple comparison using FDRs
fdr_bh_p = multipletests(clinical_trans['Wilcoxon_p'], method='fdr_bh')[1]
clinical_trans['BH_FDR_p'] = fdr_bh_p
# define up or down 
for i in clinical_trans.index:
	mean_log2FC = clinical_trans.loc[i, 'mean_log2FC']
	fdr_bh_p = clinical_trans.loc[i, 'BH_FDR_p']
	if fdr_bh_p < 0.05:
		if mean_log2FC > 0.58:
			clinical_trans.loc[i, 'status'] = 'up'
		elif mean_log2FC < -0.58:
			clinical_trans.loc[i, 'status'] = 'down'
		else:
			clinical_trans.loc[i, 'status'] = 'no_diff'
	else:
		clinical_trans.loc[i, 'status'] = 'no_diff'

# create tables for boxplots
def boxplot_table(para):	# para is either os or rfs
	sub_df = clinical_trans[clinical_trans['HR_' + para] != '*']
	FC_df = pd.DataFrame()
	for i in sub_df.index:
		order = sub_df.loc[i, 'order_' + para]
		trans = sub_df.loc[i, 'transcript']
		status = sub_df.loc[i, 'status']
		log2FC_values = exp_df[trans]
		for k in log2FC_values:
			cur_len = len(FC_df)
			FC_df.loc[cur_len, 'order'] = order
			FC_df.loc[cur_len, 'transcript'] = trans
			FC_df.loc[cur_len, 'log2FC'] = k
			FC_df.loc[cur_len, 'status'] = status
	FC_df.to_csv('D:\\MCGDYY\\ont_project\\prognosis\\FC_boxplot_' + para + '.csv', index = False)

boxplot_table('os')
boxplot_table('rfs')