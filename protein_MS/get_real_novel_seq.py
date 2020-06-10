# get the strand-specific sequence of novel transcripts

import pandas as pd
import re
from Bio.Seq import Seq

# get the list of novel transcripts with unique junctions and retained junctions
novel_uniq_RI_list = []
novel_uniq_RI_file = open('D:/MCGDYY/ont_project/lists/list_of_novel_trans_uniq_junc_RI.txt', 'r')
for i in novel_uniq_RI_file:
	i = i.strip('\n')
	novel_uniq_RI_list.append(i)

# get strandness info
gtf = pd.read_csv('D:\\MCGDYY\\ont_project\\flair_out\\flair.collapse.isoforms.stringent.gtf', header = None, sep = '\t')
strand_info = pd.DataFrame(columns = ['strand'])
trans_list = []
for i in gtf.index:
	info = gtf.iloc[i, 8]
	trans = re.findall(r"transcript_id \"(.*?)\";", info)[0]
	if 'ENST' not in trans and trans not in trans_list:
		trans_list.append(trans)
		gene = re.findall(r"gene_id \"(.*?)\";", info)[0]
		new_name = trans + '_' + gene
		strand = gtf.iloc[i, 6]
		strand_info.loc[new_name, 'strand'] = strand

# reverse complement fa
rev_fa = ''
filtered = ''
with open('D:\\MCGDYY\\ont_project\\flair_out\\flair.collapse.isoforms.stringent.fa', 'r') as fa:
	while True:
		title = fa.readline().strip('\n')
		if not title:
			break

		if title.startswith('>') and 'ENST' not in title:
			transcript = title.split('_')[0][1:]
			strand = strand_info.loc[title[1:], 'strand']
			rev_fa += title + '\n'		
			seq = fa.readline().strip('\n')

			if strand == '+':
				rev_fa += seq + '\n'
			else:
				seq = str(Seq(seq).reverse_complement())
				rev_fa += seq + '\n'

			if transcript in novel_uniq_RI_list and len(seq) > 200:
				filtered += title + '\n'
				filtered += seq + '\n'

with open('D:/MCGDYY/ont_project/flair_out/all_novel_real_seq.fa', 'w') as w:
    w.write(rev_fa)

with open('D:/MCGDYY/ont_project/flair_out/novel_uniq_RI_real_seq.fa', 'w') as w:
    w.write(filtered)