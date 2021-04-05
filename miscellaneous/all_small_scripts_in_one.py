# an all-in-one script...plz run each part separately as the parts might share same variables

"""
# PART 1

##########################################################################################
###### produce a plot for distributions of annotated and novel transcript length #########
##########################################################################################

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

anno = pd.DataFrame(columns = ['name', 'length'])
novel = pd.DataFrame(columns = ['name', 'length'])

with open('D:\\MCGDYY\\ont_project\\flair_out\\flair.collapse.isoforms.stringent.fa', 'r') as all_seq:
	while True:
		line = all_seq.readline().strip('\n')
		if not line:
			break

		if line.startswith('>'):
			trans_name = line.split('_')[0][1:]
			line = all_seq.readline().strip('\n')
			trans_len = len(line)
			tmp = pd.DataFrame({'name': [trans_name], 'length': [trans_len]})

			if trans_name.startswith('ENST'):
				anno = anno.append(tmp)
			else:
				novel = novel.append(tmp)
anno['length'] = anno['length'].astype('int')
novel['length'] = novel['length'].astype('int')
t_test_p = stats.ttest_ind(anno['length'], novel['length'])[1]	# p-value is wayyyy less than 0.001

plot = sns.kdeplot(anno['length'], shade = True, color = "b")
plot = sns.kdeplot(novel['length'], shade = True, color = "r")
plt.legend(labels = ['annotated transcripts', 'novel transcripts'], fontsize = 20)
plt.xlabel('Transcript length (bp)', fontsize = 30)
plt.ylabel('Density', fontsize = 30)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plot.spines['right'].set_visible(False)
plot.spines['top'].set_visible(False)
plt.text(0.5, 0.5, 'P < 0.001',  transform=plot.transAxes, fontsize = 20)
plt.show(plot)

"""

"""
# PART2

########################################################################################
######################################## summarize AS events ###########################
########################################################################################

def process_res(in_file_path, as_type):
	anno = []
	novel = []
	with open(in_file_path, 'r') as in_file:
		line = in_file.readline()
		while True:
			if not line:
				break

			line = in_file.readline().strip('\n')
			trans_names = line.split('\t')[-1]
			trans_names = trans_names.split(',')
			for trans in trans_names:
				if trans.startswith('ENST'):
					if trans not in anno:
						anno.append(trans)
				else:
					if trans not in novel:
						novel.append(trans)
	print(as_type + ' in annotated transcripts: ' + str(len(anno)))
	print(as_type + ' in novel transcripts: ' + str(len(novel)))

process_res('D:\\MCGDYY\\ont_project\\flair_out\\flair.diffsplice.alt3.events.quant.tsv', 'alt3')
process_res('D:\\MCGDYY\\ont_project\\flair_out\\flair.diffsplice.alt5.events.quant.tsv', 'alt5')
process_res('D:\\MCGDYY\\ont_project\\flair_out\\flair.diffsplice.es.events.quant.tsv', 'es')
process_res('D:\\MCGDYY\\ont_project\\flair_out\\flair.diffsplice.ir.events.quant.tsv', 'ir')

"""

"""

# PART 3

############################################################################################
########## explore sequencing depth vs. identified annotated and novel transcripts #########
############################################################################################

import pandas as pd
import matplotlib.pyplot as plt

count_table = pd.read_csv('D:\\MCGDYY\\ont_project\\flair_out\\counts.stringent.tsv', sep = '\t', index_col = 0)
new_cols = []
for i in count_table.columns:
	k = i.split('_')[0] + '_' + i.split('_')[1]
	new_cols.append(k)
count_table.columns = new_cols
sum_table = pd.DataFrame(columns = ['matched_name', 'depth', 'anno_trans', 'novel_trans'], index = new_cols)
sum_table = sum_table.fillna(0)

# renaming and sequencing depth
match_names = {'N12_N': 'A_N',
				'N12_T': 'A_T',
				'N12_M': 'A_M',
				'N23_N': 'B_N',
				'N23_T': 'B_T',
				'N23_M': 'B_M',
				'N25_N': 'C_N',
				'N25_T': 'C_T',
				'N25_M': 'C_M'}
match_dep = {'N12_N': 637960,
				'N12_T': 1066841,
				'N12_M': 363883,
				'N23_N': 117215,
				'N23_T': 166795,
				'N23_M': 706471,
				'N25_N': 939831,
				'N25_T': 1264000,
				'N25_M': 805493}

for trans in count_table.index:
	for sample in count_table.columns:
		if sum_table.loc[sample, 'matched_name'] == 0:
			sum_table.loc[sample, 'matched_name'] = match_names[sample]
			sum_table.loc[sample, 'depth'] = match_dep[sample]
		count = count_table.loc[trans, sample]
		if count > 0:
			if 'ENST' in trans:
				sum_table.loc[sample, 'anno_trans'] += 1
			else:
				sum_table.loc[sample, 'novel_trans'] += 1

ax = plt.subplot(111)
ax.plot(sum_table['depth']/1000, sum_table['anno_trans']/1000, 'o', color = 'blue')
ax.plot(sum_table['depth']/1000, sum_table['novel_trans']/1000, 'o', color = 'red')
plt.legend(labels = ['annotated transcripts', 'novel transcripts'])
plt.xlabel('Pass read number (k)')
plt.ylabel('Identified transcript number (k)')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.show(ax)

"""

"""
# PART 4

############################################################################################
########## compare exon numbers of identified annotated and novel transcripts #########################
############################################################################################

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import re
from scipy import stats

sum_table = pd.DataFrame(columns = ['exon_number'])
with open('D:\\MCGDYY\\ont_project\\flair_out\\flair.collapse.isoforms.stringent.gtf', 'r') as df:
	for line in df:
		info = line.split('\t')[-1]
		trans_name = re.findall(r'transcript_id \"(.*?)\"', info)[0]
		exon_num = int(re.findall(r'exon_number \"(.*?)\";', info)[0]) + 1	# it's 0-based
		if trans_name not in sum_table.index:
			sum_table.loc[trans_name, 'exon_number'] = exon_num
		else:
			if exon_num > sum_table.loc[trans_name, 'exon_number']:
				sum_table.loc[trans_name, 'exon_number'] = exon_num
sum_table['name'] = sum_table.index
anno = sum_table[sum_table['name'].str.contains('ENST')]
novel = sum_table[-sum_table['name'].str.contains('ENST')]
ttest_p = stats.ttest_ind(anno['exon_number'], novel['exon_number'])[1]	# p-value is way less than 0.001

plot = sns.kdeplot(anno['exon_number'], shade = True, color = "b")
plot = sns.kdeplot(novel['exon_number'], shade = True, color = "r")
plt.legend(labels = ['annotated transcripts', 'novel transcripts'], fontsize = 20)
plt.xlabel('Exon number', fontsize = 30)
plt.ylabel('Density', fontsize = 30)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plot.spines['right'].set_visible(False)
plot.spines['top'].set_visible(False)
plt.text(0.5, 0.5, 'P < 0.001',  transform=plot.transAxes, fontsize = 20)
plt.show(plot)

"""

"""
# PART 5
# Summary of NovelQuant results

import pandas as pd
import re

RI_df = pd.read_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\NovelQuant_res\\retained_introns.gtf', sep = '\t', header = None)
RI_list = []
for i in RI_df[8]:
	trans = re.findall(r"transcript_id \"(.*?)\";", i)[0]
	if trans not in RI_list:
		RI_list.append(trans)

print(len(RI_list))
print(len(RI_df))
print(len(RI_df[RI_df[5] == 'IR']))
print(len(RI_df[RI_df[5] == 'CE']))
print(len(RI_df[RI_df[5] == 'E3E']))
print(len(RI_df[RI_df[5] == 'E5E']))

UJ_df = pd.read_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\NovelQuant_res\\uniq_eej.gtf', sep = '\t', header = 0)
UJ_list = []
for i in UJ_df['transcript']:
	if i not in UJ_list:
		UJ_list.append(i)

print(len(UJ_list))
print(len(UJ_df))

intersec = []
for i in RI_list:
	if i in UJ_list:
		intersec.append(i)
print(len(intersec))

all_trans = []
for i in RI_list:
	if i not in all_trans:
		all_trans.append(i)
for i in UJ_list:
	if i not in all_trans:
		all_trans.append(i)
print(len(all_trans))
"""

"""
# PART 6
# Summary of SUPPA results

import pandas as pd

sum_table = pd.DataFrame(columns = ['annotated', 'novel'])

# RI

df = pd.read_csv('D:\\MCGDYY\\ont_project\\SUPPA\\novel_AS_RI_strict.ioe', sep = '\t')
anno_list = []
novel_list = []
for i in df.index:
	info = df.loc[i, 'alternative_transcripts']
	for k in info.split(','):
		if 'ENST' in k:
			if k not in anno_list:
				anno_list.append(k)
		else:
			if k not in novel_list:
				novel_list.append(k)
sum_table.loc['RI', 'annotated'] = len(anno_list)
sum_table.loc['RI', 'novel'] = len(novel_list)

# SE

df = pd.read_csv('D:\\MCGDYY\\ont_project\\SUPPA\\novel_AS_SE_strict.ioe', sep = '\t')
anno_list = []
novel_list = []
for i in df.index:
	sub_info = df.loc[i, 'alternative_transcripts'].split(',')
	all_info = df.loc[i, 'total_transcripts'].split(',')
	for k in all_info:
		if k not in sub_info:
			if 'ENST' in k:
				anno_list.append(k)
			else:
				novel_list.append(k)
sum_table.loc['SE', 'annotated'] = len(anno_list)
sum_table.loc['SE', 'novel'] = len(novel_list)

# A3

df = pd.read_csv('D:\\MCGDYY\\ont_project\\SUPPA\\novel_AS_A3_strict.ioe', sep = '\t')
anno_list = []
novel_list = []
for i in df.index:
	all_info = df.loc[i, 'alternative_transcripts'].split(',')
	for k in all_info:
		if 'ENST' in k:
			anno_list.append(k)
		else:
			novel_list.append(k)
sum_table.loc['A3', 'annotated'] = len(anno_list)
sum_table.loc['A3', 'novel'] = len(novel_list)

# A5

df = pd.read_csv('D:\\MCGDYY\\ont_project\\SUPPA\\novel_AS_A5_strict.ioe', sep = '\t')
anno_list = []
novel_list = []
for i in df.index:
	all_info = df.loc[i, 'alternative_transcripts'].split(',')
	for k in all_info:
		if 'ENST' in k:
			anno_list.append(k)
		else:
			novel_list.append(k)
sum_table.loc['A5', 'annotated'] = len(anno_list)
sum_table.loc['A5', 'novel'] = len(novel_list)

# AF

df = pd.read_csv('D:\\MCGDYY\\ont_project\\SUPPA\\novel_AS_AF_strict.ioe', sep = '\t')
anno_list = []
novel_list = []
for i in df.index:
	all_info = df.loc[i, 'alternative_transcripts'].split(',')
	for k in all_info:
		if 'ENST' in k:
			anno_list.append(k)
		else:
			novel_list.append(k)
sum_table.loc['AF', 'annotated'] = len(anno_list)
sum_table.loc['AF', 'novel'] = len(novel_list)

# AL

df = pd.read_csv('D:\\MCGDYY\\ont_project\\SUPPA\\novel_AS_AL_strict.ioe', sep = '\t')
anno_list = []
novel_list = []
for i in df.index:
	all_info = df.loc[i, 'alternative_transcripts'].split(',')
	for k in all_info:
		if 'ENST' in k:
			anno_list.append(k)
		else:
			novel_list.append(k)
sum_table.loc['AL', 'annotated'] = len(anno_list)
sum_table.loc['AL', 'novel'] = len(novel_list)

# MX

df = pd.read_csv('D:\\MCGDYY\\ont_project\\SUPPA\\novel_AS_MX_strict.ioe', sep = '\t')
anno_list = []
novel_list = []
for i in df.index:
	all_info = df.loc[i, 'alternative_transcripts'].split(',')
	for k in all_info:
		if 'ENST' in k:
			anno_list.append(k)
		else:
			novel_list.append(k)
sum_table.loc['MX', 'annotated'] = len(anno_list)
sum_table.loc['MX', 'novel'] = len(novel_list)

sum_table.to_csv('D:\\MCGDYY\\ont_project\\SUPPA\\sum.csv')
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

novel_df = pd.read_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\NovelQuant_res.csv', index_col = 0)
master_df = pd.DataFrame(columns = ['CDO1', 'CYP', 'CDO1-novel', 'CYP-novel'])
for sample in os.listdir('D:\\MCGDYY\\ont_project\\rsem\\'):
    if 'isoforms' in sample:
        full_path = 'D:\\MCGDYY\\ont_project\\rsem\\' + sample
        sample = sample.split('.')[0]
        each_df = pd.read_csv(full_path, sep = '\t', index_col = 0)
        CDO1 = each_df.loc['ENST00000250535.5', 'TPM']
        master_df.loc[sample, 'CDO1'] = CDO1
        CYP = each_df.loc['ENST00000301141.10', 'TPM']
        master_df.loc[sample, 'CYP'] = CYP
        
        master_df.loc[sample, 'CDO1-novel'] = novel_df.loc['b5082e5d-a8e4-4986-841d-9f416d0d57fe', sample]
        master_df.loc[sample, 'CYP-novel'] = novel_df.loc['9c1c526a-89d8-4fed-afb7-8a1eb93189a3', sample]

for i in master_df.index:
    if '-1' in i:
        master_df.loc[i, 'type'] = 'T'
    else:
        master_df.loc[i, 'type'] = 'N'

# groups = master_df.groupby('type')
# for name, group in groups:
#     if name == 'N':
#         col = 'b'
#     else:
#         col = 'r'
#     plot = plt.scatter(np.log2(group['CDO1'].astype(np.float64)), np.log2(group['CDO1-novel'].astype(np.float64)), 
#             label = name, c = col)
# plt.legend(labels = ['peritumor', 'tumor'], fontsize = 15)
# plt.xlabel(r'log2(CDO1 expression in TPM)', fontsize = 15)
# plt.ylabel(r'log2(CDO1-novel expression)', fontsize = 15)   
# plt.show(plot)  

groups = master_df.groupby('type')
for name, group in groups:
    if name == 'N':
        col = 'b'
    else:
        col = 'r'
    plot = plt.scatter(np.log2(group['CYP'].astype(np.float64)), np.log2(group['CYP-novel'].astype(np.float64)), 
            label = name, c = col)
plt.legend(labels = ['peritumor', 'tumor'], fontsize = 15)
plt.xlabel(r'log2(CYP2A6 expression in TPM)', fontsize = 15)
plt.ylabel(r'log2(CYP2A6-novel expression)', fontsize = 15)
plt.show(plot)