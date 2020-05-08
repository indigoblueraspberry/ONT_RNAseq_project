import pandas as pd
import csv

intron_gtf = pd.read_csv('D:\\MCGDYY\\ont_project\\gtf_files\\introns.gtf', sep = '\t', chunksize = 10000, header= None,\
						names = ['chromosome', 'source', 'type', 'start', 'end', 'none1', 'strand', 'none2', 'info'])
intron_gtf = pd.concat(intron_gtf, ignore_index = True)
# retained only annotated introns
anno_intron = intron_gtf[intron_gtf['source'] != 'FLAIR']

# get all exons of novel transcripts
flair_ex_gtf = pd.read_csv('D:\\MCGDYY\\ont_project\\flair_out\\flair.collapse.isoforms.stringent.gtf', header = None, sep = '\t',
						names = ['chromosome', 'source', 'type', 'start', 'end', 'none1', 'strand', 'none2', 'info'])
ex_novel = flair_ex_gtf[-flair_ex_gtf['info'].str.contains('ENST')]
# sort the df
ex_novel.sort_values(by = 'chromosome', inplace = True)
ex_novel.reset_index(drop = True, inplace = True)

# iteration
retained_introns = pd.DataFrame(columns = ['chromosome', 'source', 'type', 'start', 'end', 'none1', 'strand', 'none2', 'info'])
for i in anno_intron.index:
	in_chr = anno_intron.loc[i, 'chromosome']
	in_start = anno_intron.loc[i, 'start']
	in_end = anno_intron.loc[i, 'end']

	tmp_df = ex_novel[ex_novel['chromosome'] == str(in_chr)]
	for k in tmp_df.index:
		ex_start = tmp_df.loc[k, 'start']
		ex_end = tmp_df.loc[k, 'end']

		if in_start >= ex_start and in_end <= ex_end:
			line = anno_intron.loc[i]
			line['info'] =  tmp_df.loc[k, 'info']
			retained_introns = retained_introns.append(line)

retained_introns.to_csv('D:\\MCGDYY\\ont_project\\retained_introns.gtf', 
	header = None, index = False, sep = '\t', quoting = csv.QUOTE_NONE)