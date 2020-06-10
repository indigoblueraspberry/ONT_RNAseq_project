import pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
import re

# load quantification of novel transcripts with uniq junc and RI whose unbalanced groups have been removed 
sum_table = pd.read_csv('D:/MCGDYY/ont_project/quantification/t_exp_median.csv', index_col = 0)

####################################### Wilcoxon signed-rank test on log2FC

final = pd.DataFrame()
for i in sum_table.columns:
	if '_value' in i:
		FC_col =sum_table[i]
		trans = i.split('_')[0]
		avg = FC_col.mean()
		final.loc[trans, 'log2FC_mean'] = avg
		wc_p = wilcoxon(FC_col)[1]	# default two-sided
		final.loc[trans, 'Wilcoxon_p'] = wc_p
final.dropna(axis = 0, inplace = True)	# remove none-expression
fdr_bh_p = multipletests(final['Wilcoxon_p'], method='fdr_bh')[1]	# adjust p-value for multiple comparison using FDRs
final['fdr_bh_p'] = fdr_bh_p
DE_novel = final[(abs(final['log2FC_mean']) > 0.58) & (final['fdr_bh_p'] < 0.05)]

######################################## DE novel lncRNA candidates

# associated with prognoses and intersection by lncRNA prediction of 3 software
novel_lncRNA = pd.read_csv('D:\\MCGDYY\\ont_project\\prog_assoc_trans\\novel_lncRNA\\intersec_novel_lncRNA.csv')
novel_lncRNA = list(novel_lncRNA['Transcript'])
DE_novel_lncRNA = pd.DataFrame()
for i in DE_novel.index:
	if i in novel_lncRNA:
		DE_novel_lncRNA = pd.concat((DE_novel_lncRNA, DE_novel.loc[i].to_frame().T), axis = 0)
# DE_novel_lncRNA.to_csv('D:\\MCGDYY\\ont_project\\prog_assoc_trans\\novel_lncRNA\\DE_novel_lncRNA.csv')

######################################## DE novel mRNA candidates

# find the mRNA supported by MS and associated with prognoses
union_os_rfs = pd.read_csv('D:\\MCGDYY\\ont_project\\lists\\union_os_rfs.csv', index_col = 0)
ms = pd.read_csv('D:/MCGDYY/ont_project/MS/new_DDA/refined_dda.csv', index_col = 0)
novel_mRNA = pd.merge(ms, union_os_rfs, how = 'inner', left_on = ms.index, right_on = 'union_os_rfs')
novel_mRNA = list(novel_mRNA['union_os_rfs'])
DE_novel_mRNA = pd.DataFrame()
for i in DE_novel.index:
	if i in novel_mRNA:
		DE_novel_mRNA = pd.concat((DE_novel_mRNA, DE_novel.loc[i].to_frame().T), axis = 0)
# DE_novel_mRNA.to_csv('D:\\MCGDYY\\ont_project\\prog_assoc_trans\\novel_mRNA\\DE_novel_mRNA.csv')




# # extract sequence
# lnc_seq = ''
# with open('D:\\MCGDYY\\ont_project\\lncRNA\\novel_uniq_RI_real_seq.fa') as seq_file:
# 	while True:		
# 		line = seq_file.readline()
# 		if not line:
# 			break
# 		if line.startswith('>'):
# 			trans = re.findall('>(.*?)_', line)[0]
# 			if trans in DE_novel_lnc.index:
# 				lnc_seq += line
# 				line = seq_file.readline()
# 				lnc_seq += line

# # with open('D://MCGDYY//ont_project//lncRNA/DE_novel_lncRNA.fa', 'w') as w:
# #     w.write(lnc_seq)


