import pandas as pd
import re

df = pd.read_csv('D:\\MCGDYY\\ont_project\\gtf_files\\gencode.v32.GRCh38.novel.trans.gtf', sep = '\t', comment = '#',
	names = ['chr', 'source', 'type', 'start', 'end', 'none1', 'strand', 'none2', 'info'])

anno_genes = df[df['type'] == 'gene'].reset_index(drop = True)

novel_trans = df[df['type'] == 'exon']
novel_trans = novel_trans[-novel_trans['info'].str.contains('ENST')].reset_index(drop = True)
# append two columns of gene and transcript names
for i in range(len(novel_trans)):
	info = novel_trans.loc[i, 'info']
	novel_trans.loc[i, 'gene'] = re.findall(r"gene_id \"(.*?)\";", info)[0]
	novel_trans.loc[i, 'transcript'] = re.findall(r"transcript_id \"(.*?)\";", info)[0]
novel_trans = novel_trans.sort_values(by = ['chr', 'transcript', 'start']).reset_index(drop = True)

# check on the novel intergenic transcripts
intergenic = novel_trans[-novel_trans['gene'].str.contains('ENSG')]
intergenic_list = intergenic[['gene', 'transcript']].drop_duplicates(subset = 'transcript', keep = 'first')

intergenic_list['intergenic'] = 'yes'
for i in intergenic_list.index:
	trans_name = intergenic_list.loc[i, 'transcript']
	that_trans = intergenic[intergenic['transcript'] == trans_name].reset_index(drop = True)
	chr1 = str(that_trans.loc[0, 'chr'])
	start1 = that_trans.loc[0, 'start']
	end1 = that_trans.loc[len(that_trans) - 1, 'end']

	for k in anno_genes.index:
		chr2 = str(anno_genes.loc[k, 'chr'])
		start2 = anno_genes.loc[k, 'start']
		end2 = anno_genes.loc[k, 'end']

		if chr1 == chr2:
			if end1 <= start2 or start1 >= end2:
				pass
			else:
				info = anno_genes.loc[k, 'info']
				gene_name = re.findall(r"gene_id \"(.*?)\";", info)[0]
				intergenic_list.loc[i, 'intergenic'] = gene_name

intergenic_list.to_csv('D:/MCGDYY/ont_project/haha.txt', sep = '\t', header=True, index= False)