# Used to collect junction reads for novel multiple-exon TSTs that have unique junctions

import pandas as pd
import re
import csv
import os

############################################ Step 1: identify unique junctions

# intron.gtf is too big for pd.read_csv
intron_gtf = pd.read_csv('D:\\MCGDYY\\ont_project\\gtf_files\\introns.gtf', sep = '\t', chunksize = 10000, header= None,\
						names = ['chromosome', 'source', 'type', 'start', 'end', 'none1', 'strand', 'none2', 'info'])
intron_gtf = pd.concat(intron_gtf, ignore_index = True)

# append two columns of gene and transcript names
gene = []
trans = []
for entry in intron_gtf['info']:
    each_gene = re.findall(r"gene_id \"(.*?)\";", entry)[0]
    each_trans = re.findall(r"transcript_id \"(.*?)\";", entry)[0]
    gene.append(each_gene)
    trans.append(each_trans)

intron_gtf['gene'] = gene
intron_gtf['transcript'] = trans

# find unique junctions that can identify all specific transcripts, dedup based on chr, start and end.
uniq_junc = intron_gtf.drop_duplicates(subset = ['chromosome', 'start', 'end'], keep = False)
# find uniq junctions for novel transcripts
uniq_junc_novel_trans = uniq_junc[-uniq_junc['transcript'].str.contains('ENST')]
uniq_junc_novel_trans = uniq_junc_novel_trans.reset_index(drop = True)

# create a list of novel transcripts with unique junctions
dedup_trans = uniq_junc_novel_trans[['gene', 'transcript']].drop_duplicates(subset = 'transcript', keep = 'first')
dedup_trans.to_csv('D:\\MCGDYY\\ont_project\\list_of_novel_trans_with_uniq_junc.csv', index = False)


##################################################### Step 2: process STAR SJ.out.tab to extract junction reads

star_junc_main_dir = 'D:\\MCGDYY\\ont_project\\star_junctions\\'
for file in os.listdir(star_junc_main_dir):
	sample_id = file[:6]
	star_junc_path = star_junc_main_dir + file
	star_junc = pd.read_csv(star_junc_path, sep = '\t', header = None,\
							names = ['chromosome', 'start', 'end', 'strand', 'canonical', 'annotated', 'uniq_reads', 'multi_reads', 'max_overhang'])
	merged = pd.merge(uniq_junc_novel_trans, star_junc, how= 'inner', on = ['chromosome', 'start', 'end'], sort = True)

	# quoting = csv.QUOTE_NONE to avoid "" in the text
	merged.to_csv('D:\\MCGDYY\\ont_project\\novel_trans_with_uniq_junc\\' + sample_id + '_uniq_junc_reads.txt', sep='\t', header=True, index=False, quoting = csv.QUOTE_NONE)


