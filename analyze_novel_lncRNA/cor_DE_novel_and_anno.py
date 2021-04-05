import pandas as pd
from scipy import stats

novel_counts = pd.read_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\t_exp_median.csv', index_col = 0)
all_counts_df = pd.read_csv('D:\\MCGDYY\\ont_project\\quantification\\all_counts.csv', index_col = 0)
novel_t2g = pd.read_csv('D:\\MCGDYY\\ont_project\\lists\\all_novel_list.txt', index_col = 0, names = ['gene'], sep = '\t')
DE_novel = pd.read_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\lncRNA\\intersec_novel_lncRNA.csv', index_col = 0)
DE_anno = pd.read_csv('D:\\MCGDYY\\ont_project\\lncRNA\\DE_anno_mRNA.csv', index_col = 0)
anno_t2g = pd.read_csv('D:\\MCGDYY\\ont_project\\lists\\anno_mRNA_trans.txt', index_col = 1, names = ['gene', 'xx', 'xxx'], sep = '\t')

DE_novel_list = list(i.split('_')[0] for i in DE_novel.index)
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

		# remove comparison inside the same gene
		if g1 != g2:
			if abs(rho) > 0.8 and pval < 0.05:
				sig_table.loc[k, i] = str(rho) + ', ' + str(pval)
sum_table.to_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\lncRNA\\all_cor_novel_anno.csv')				
sig_table.to_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\lncRNA\\raw_sig_cor_novel_anno.csv')

# annotate results
novel_names = pd.read_csv('D:\\MCGDYY\\ont_project\\lists\\all_novel_list.txt', index_col = 0, header = None, sep = '\t')
trans_2_gene = pd.read_csv('D:\\MCGDYY\\ont_project\\lists\\all_trans_to_gene_list.txt', index_col = 0, header = None, sep = '\t')
final_table = pd.DataFrame(columns = ['novel', 'parent_gene_ensembl', 'parent_gene', 'anno_ensembl', 'anno', 'rho', 'p-value'])

for i in sig_table.columns:
	for k in sig_table.index:
		if pd.isnull(sig_table.loc[k, i]) is False:
			cur_len = len(final_table)
			final_table.loc[cur_len, 'novel'] = i
			parent_gene_ensembl = novel_names.loc[i, 1]
			final_table.loc[cur_len, 'parent_gene_ensembl'] = parent_gene_ensembl
			try:
				parent_gene = trans_2_gene.loc[parent_gene_ensembl, 2].values[0]
			except AttributeError:
				parent_gene = trans_2_gene.loc[parent_gene_ensembl, 2]
			final_table.loc[cur_len, 'parent_gene'] = parent_gene
			final_table.loc[cur_len, 'anno_ensembl'] = k
			final_table.loc[cur_len, 'anno'] = trans_2_gene[trans_2_gene[1] == k].iloc[0, 2]
			final_table.loc[cur_len, 'rho'] = sig_table.loc[k, i].split(',')[0]
			final_table.loc[cur_len, 'p-value'] = sig_table.loc[k, i].split(',')[1]
final_table.to_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\lncRNA\\compiled_sig_cor_novel_anno.csv', index = False)