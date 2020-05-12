import pandas as pd

# get all introns of annotated transcripts
anno_trans_introns = pd.read_csv('D:\\MCGDYY\\ont_project\\gtf_files\\anno_trans_introns.gtf', 
	sep = '\t', chunksize = 10000, header= None, comment = '#',
	names = ['chromosome', 'source', 'type', 'start', 'end', 'none1', 'strand', 'none2', 'info'])
anno_trans_introns = pd.concat(anno_trans_introns, ignore_index = True)

# get all exons of novel transcripts
novel_trans_exons = pd.read_csv('D:\\MCGDYY\\ont_project\\gtf_files\\novel_trans_exons.gtf', 
	header = None, sep = '\t', comment = '#',
	names = ['chromosome', 'source', 'type', 'start', 'end', 'none1', 'strand', 'none2', 'info'])

# check if intronic regions are included in the novel transcripts
retained_introns = ''
for i in anno_trans_introns.index:
	in_chr = anno_trans_introns.loc[i, 'chromosome']
	in_start = anno_trans_introns.loc[i, 'start']
	in_end = anno_trans_introns.loc[i, 'end']
	in_strand = anno_trans_introns.loc[i, 'strand']

	# select specific chromosome and nearby sites to reduce iteration and speed up the script
	tmp_df = novel_trans_exons[novel_trans_exons['chromosome'] == str(in_chr)]
	tmp_df = tmp_df[tmp_df['start'] - in_start < 10000]
	for k in tmp_df.index:
		ex_strand = tmp_df.loc[k, 'strand']
		# make sure they are both on the same strand
		if in_strand == ex_strand:
			ex_start = tmp_df.loc[k, 'start']
			ex_end = tmp_df.loc[k, 'end']

			# full intron retaintion
			if in_start >= ex_start and in_end <= ex_end:
				retained_introns += str(in_chr) + '\t' + '*' + '\t' + 'intron' + '\t' + str(in_start) + '\t' + str(in_end) + '\t' \
				+ 'IR' + '\t' + in_strand + '\t' + '*' + '\t' + tmp_df.loc[k, 'info'] + '\n'

			# alternative 5' splice site
			if in_start <= ex_end and in_end >= ex_end:
				retained_introns += str(in_chr) + '\t' + '*' + '\t' + 'intron' + '\t' + str(in_start) + '\t' + str(ex_end) + '\t' \
				+ 'ALT5' + '\t' + in_strand + '\t' + '*' + '\t' + tmp_df.loc[k, 'info'] + '\n'

			# alternative 3' splice site
			if in_start <= ex_start and in_end >= ex_start:
				retained_introns += str(in_chr) + '\t' + '*' + '\t' + 'intron' + '\t' + str(ex_start) + '\t' + str(in_end) + '\t' \
				+ 'ALT3' + '\t' + in_strand + '\t' + '*' + '\t' + tmp_df.loc[k, 'info'] + '\n'

with open('D:\\MCGDYY\\ont_project\\gtf_files\\retained_introns.gtf', 'w') as w:
	w.write(retained_introns)