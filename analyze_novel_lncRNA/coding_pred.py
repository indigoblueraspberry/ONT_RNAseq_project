# predict coding potentials of novel transcripts with CNCI, CPC2 and PLEK
# also integrate protein MS data to confirm

import pandas as pd

CNCI_df = pd.read_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\lncRNA\\CNCI_out.txt', sep = '\t')
CPC2_df = pd.read_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\lncRNA\\CPC2_out.txt', sep = '\t')
PLEK_df = pd.read_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\lncRNA\\PLEK_out.txt', sep = '\t', header = None)
sum_table = pd.DataFrame(columns = ['CNCI', 'CPC2', 'PLEK'])

for i in CNCI_df.index:
	name = CNCI_df.loc[i, r'Transcript ID']
	res = CNCI_df.loc[i, 'index']
	sum_table.loc[name, 'CNCI'] = res

for i in CPC2_df.index:
	name = CPC2_df.loc[i, '#ID']
	res = CPC2_df.loc[i, 'label']
	sum_table.loc[name, 'CPC2'] = res

for i in PLEK_df.index:
	name = PLEK_df.iloc[i, 2][1:]
	res = PLEK_df.iloc[i, 0]
	if res == 'Coding':
		sum_table.loc[name, 'PLEK'] = 'coding'
	else:
		sum_table.loc[name, 'PLEK'] = 'noncoding'

sum_table = sum_table.fillna('noncoding')

for i in sum_table.index:
	trans = i.split('_')[0]
	sum_table.loc[i, 'Transcript'] = trans
	gene = i.split('_')[1]
	sum_table.loc[i, 'Gene'] = gene

def sum_data(feature):	# feature is either coding or noncoding

	# only CPC2
	count1 = len(sum_table[(sum_table['CPC2'] == feature) & (sum_table['CNCI'] != feature) & (sum_table['PLEK'] != feature)])
	print(feature + ' in CPC2 only: ' + str(count1))
	# only CNCI
	count2 = len(sum_table[(sum_table['CNCI'] == feature) & (sum_table['CPC2'] != feature) & (sum_table['PLEK'] != feature)])
	print(feature + ' in CNCI only: ' + str(count2))
	# only PLEK
	count3 = len(sum_table[(sum_table['PLEK'] == feature) & (sum_table['CPC2'] != feature) & (sum_table['CNCI'] != feature)])
	print(feature + ' in PLEK only: ' + str(count3))

	# CNCI & CPC2
	count4 = len(sum_table[(sum_table['CNCI'] == feature) & (sum_table['CPC2'] == feature)])
	print(feature + ' intersection between CNCI and CPC2: ' + str(count4))
	# CNCI & PLEK
	count5 = len(sum_table[(sum_table['CNCI'] == feature) & (sum_table['PLEK'] == feature)])
	print(feature + ' intersection between CNCI and PLEK: ' + str(count5))
	# CPC2 & PLEK
	count6 = len(sum_table[(sum_table['CPC2'] == feature) & (sum_table['PLEK'] == feature)])
	print(feature + ' intersection between CPC2 and PLEK: ' + str(count6))
	
	# CNCI & CPC2 & PLEK
	final = sum_table[(sum_table['CPC2'] == feature) & (sum_table['PLEK'] == feature) & (sum_table['CNCI'] == feature)]
	print(feature + ' intersection among three: ' + str(len(final)))

	# CPC2 all
	count7 = len(sum_table[sum_table['CPC2'] == feature])
	print(feature + ' CPC2: ' + str(count7))

	# CPC2 all
	count8 = len(sum_table[sum_table['CNCI'] == feature])
	print(feature + ' CNCI: ' + str(count8))

	# CPC2 all
	count9 = len(sum_table[sum_table['PLEK'] == feature])
	print(feature + ' PLEK: ' + str(count9))

	# final.to_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\lncRNA\\intersec_novel_lncRNA.csv')

sum_data('noncoding')