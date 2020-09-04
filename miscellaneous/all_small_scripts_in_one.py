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
plt.legend(labels = ['annotated transcripts', 'novel transcripts'])
plt.xlabel('Transcript length (bp)')
plt.ylabel('Density')
plt.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0))
plot.spines['right'].set_visible(False)
plot.spines['top'].set_visible(False)
plt.text(0.65, 0.7, 'p-value < 0.001',  transform=plot.transAxes)
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
plt.legend(labels = ['annotated transcripts', 'novel transcripts'])
plt.xlabel('Exon number')
plt.ylabel('Density')
plot.spines['right'].set_visible(False)
plot.spines['top'].set_visible(False)
plt.text(0.65, 0.7, 'p-value < 0.001',  transform=plot.transAxes)
plt.show(plot)

"""