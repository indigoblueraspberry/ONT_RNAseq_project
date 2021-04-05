import pandas as pd

coexp_df = pd.read_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\lncRNA\\compiled_sig_cor_novel_anno_0.6.csv')
targetscan_df = pd.read_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\lncRNA\\GRN\\processed_targetscan.txt', header = None, sep = '\t')
mirtarbase_df = pd.read_excel('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\lncRNA\\GRN\\hsa_MTI.xlsx', header = 0)
mirtarbase_df = mirtarbase_df[mirtarbase_df['Support Type'] == 'Functional MTI']
# coexp_df = coexp_df[abs(coexp_df['rho']) > 0.7] 
sum_df = ''

for i in targetscan_df.index:
	new_name = targetscan_df.loc[i, 0].split('_')[0]
	targetscan_df.loc[i, 3] = new_name

for i in coexp_df.index:
	novel = coexp_df.loc[i, 'novel']
	anno = coexp_df.loc[i, 'anno']
	anno = anno.split('-')[0]
	tmp_1 = targetscan_df[targetscan_df[3] == novel]
	for mirna in tmp_1[1]:
		tmp_2 = mirtarbase_df[(mirtarbase_df['miRNA'] == mirna) & (mirtarbase_df['Target Gene'] == anno)]
		if len(tmp_2) != 0:
			print(coexp_df.loc[i, 'rho'])
			sum_df += novel + '\t' + mirna + '\n'
			sum_df += mirna + '\t' + anno + '\n'

with open('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\lncRNA\\GRN\\lncRNA_miRNA_mRNA_pairs.txt', 'w') as writer:
	writer.write(sum_df)


