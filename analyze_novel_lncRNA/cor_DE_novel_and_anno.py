import pandas as pd
from scipy import stats

novel_counts = pd.read_csv('D:\\MCGDYY\\ont_project\\quantification\\t_exp_median.csv', index_col = 0)
all_counts_df = pd.read_csv('D:\\MCGDYY\\ont_project\\quantification\\all_counts.csv', index_col = 0)
novel_t2g = pd.read_csv('D:\\MCGDYY\\ont_project\\lists\\all_novel_list.txt', index_col = 0, names = ['gene'], sep = '\t')


def correlate(novel_type, anno_type):	# either 'lncRNA' or 'mRNA'

	DE_novel = pd.read_csv('D:\\MCGDYY\\ont_project\\prog_assoc_trans\\novel_' + novel_type + '\\DE_novel_' + novel_type + '.csv', index_col = 0)
	DE_anno = pd.read_csv('D:\\MCGDYY\\ont_project\\prog_assoc_trans\\novel_' + novel_type + '\\DE_anno_' + anno_type + '.csv', index_col = 0)
	anno_t2g = pd.read_csv('D:\\MCGDYY\\ont_project\\lists\\anno_' + anno_type + '_trans.txt', index_col = 1, names = ['gene', 'xx', 'xxx'], sep = '\t')

	DE_novel_list = list(DE_novel.index)
	clean_novel = pd.DataFrame()
	for i in DE_novel_list:
		clean_novel = pd.concat((clean_novel, novel_counts[i + '_T']), axis = 1)
	clean_novel.columns = [i.split('_')[0] for i in clean_novel.columns]

	DE_anno_list = list(DE_anno.index)
	clean_anno = pd.DataFrame()
	all_counts = all_counts_df.T
	for i in DE_anno_list:
		clean_anno = pd.concat((clean_anno, all_counts[i]), axis = 1)
	clean_anno = clean_anno[clean_anno.index.str.contains('-1')]
	clean_anno.index = [i.split('-')[0] for i in clean_anno.index]
	clean_anno = clean_anno.drop(index = '1848')	# 1848 has no clinical info

	# calculate Spearman's correlation
	sum_table = pd.DataFrame()
	sig_table = pd.DataFrame()
	for i in clean_novel.columns:
		g1 = novel_t2g.loc[i, 'gene']
		for k in clean_anno.columns:
			g2 = anno_t2g.loc[k, 'gene']
			rho, pval = stats.spearmanr(clean_novel[i], clean_anno[k])
			sum_table.loc[k + '_rho', i] = rho
			sum_table.loc[k + '_pval', i] = pval

			if g1 != g2:
				if abs(rho) > 0.8 and pval < 0.05:
					sig_table.loc[k, i] = str(rho) + ', ' + str(pval)
	return sum_table, sig_table


sum1, sig1 = correlate(novel_type = 'mRNA', anno_type = 'lncRNA')
sum2, sig2 = correlate(novel_type = 'lncRNA', anno_type = 'mRNA')