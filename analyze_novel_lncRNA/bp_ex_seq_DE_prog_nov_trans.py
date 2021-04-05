import pandas as pd

# create tables for box plots of log2FC of DE and clinically relevant novel transcripts 

exp_df = pd.read_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\t_exp_median_DEonly.csv')
DE_prog_novel_trans_list = []

def create_bp_input(prog):	# prog = os or rfs

	sig_prog_df = pd.read_csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\prognosis\\sig_res_' + prog + '.csv')
	box_plot_df = 'order	transcript	log2FC	status' + '\n'

	for i in sig_prog_df.index:
		trans_name = sig_prog_df.loc[i, 'transcript']
		trans_name = trans_name.replace('_', '-')
		trans_name = trans_name.replace('trans-', '')
		trans_name = trans_name.replace('-status', '')

		avg_log2FC = exp_df[trans_name + '_value'].mean()
		if avg_log2FC > 0.58:
			status = 'up'
		if avg_log2FC < -0.58:
			status = 'down'
		order = sig_prog_df.loc[i, 'order']

		for k in exp_df[trans_name + '_value']:
			box_plot_df += str(order) + '\t' + trans_name + '\t' + str(k) + '\t' + status + '\n'
		if trans_name not in DE_prog_novel_trans_list:
			DE_prog_novel_trans_list.append(trans_name)

	with open('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\prognosis\\FC_boxplot_' + prog + '.txt', 'w') as w:
		w.write(box_plot_df)

create_bp_input('os')
create_bp_input('rfs')


# extract sequence of DE and clinically relevant novel transcripts

extracted_seq = ''
for i in DE_prog_novel_trans_list:
	with open('D:\\MCGDYY\\ont_project\\lncRNA\\lncRNA_pred\\novel_uniq_RI_real_seq_200bp.fa', 'r') as master_seq:
		line = master_seq.readline()
		while line:
			if line.startswith('>'):
				if i in line:
					extracted_seq += line
					line = master_seq.readline()
					extracted_seq += line
					break
			line = master_seq.readline()

with open('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\coding_pred\\DE_clinical_novel_trans_real_seq.fa', 'w') as w:
	w.write(extracted_seq)


# find the parent gene names for the DE and clinically relevant novel transcripts

all_novel_list = pd.read_csv('D:\\MCGDYY\\ont_project\\lists\\all_novel_list.txt', index_col = 0, sep = '\t', header = None)
DE_prog_parent_genes = ''
for i in DE_prog_novel_trans_list:
	parent_gene = all_novel_list.loc[i, 1]
	DE_prog_parent_genes += i + '\t' + parent_gene + '\n'
with open('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\prognosis\\DE_prog_parent_genes.txt', 'w') as w:
	w.write(DE_prog_parent_genes)
