import pandas as pd
import re
from scipy import stats
from statsmodels.stats.multitest import multipletests
from Bio.Seq import Seq

def processor(para):	# para = either 'os' or 'rfs'

	sig_res = pd.read_csv('D:\\MCGDYY\\ont_project\\prognosis\\sig_res_' + para + '.csv')	# sig_res_os.csv or sig_res_rfs.csv
	sum_table = pd.read_csv('D:\\MCGDYY\\ont_project\\quantification\\t_exp_median.csv')

	for i in sig_res.index:
		ori_name = sig_res.loc[i, 'transcript']
		name = ori_name.replace('_', '-')
		pat = re.compile('trans-(.*?)-status')
		name = re.match(pat, name)[1]
		name = name + '_value'
		sig_res.loc[i, 'transcript'] = name

		FC_col = sum_table[name]
		wc_p_val = stats.wilcoxon(FC_col)[1]
		sig_res.loc[i, 'Wilcoxon_p'] = wc_p_val
		avg_log2FC = sum_table[name].mean()
		sig_res.loc[i, 'mean_log2FC'] = avg_log2FC

	sig_res.index = sig_res['transcript']
	fdr_bh_p = multipletests(sig_res['Wilcoxon_p'], method='fdr_bh')[1]	# adjust p-value for multiple comparison using FDRs
	sig_res['BH_FDR_p'] = fdr_bh_p

	final_df = pd.DataFrame()
	for i in sig_res.index:
		tmp_df = pd.DataFrame(columns = ['order', 'transcript', 'log2FC', 'status'])
		fdr_bh_p = sig_res.loc[i, 'BH_FDR_p']
		avg_log2FC = sig_res.loc[i, 'mean_log2FC']
		# according to BH FDR p-value and log2FC to determine the status
		if fdr_bh_p < 0.05:
			if avg_log2FC > 0.58:
				sig_res.loc[i, 'status'] = 'up'
			elif avg_log2FC < -0.58:
				sig_res.loc[i, 'status'] = 'down'
			else:
				sig_res.loc[i, 'status'] = 'not_sig'
		else:
			sig_res.loc[i, 'status'] = 'not_sig'

		tmp_df['log2FC'] = sum_table[i]
		tmp_df['order'] = sig_res.loc[i, 'order']
		tmp_df['transcript'] = i
		tmp_df['status'] = sig_res.loc[i, 'status']
		final_df = pd.concat((final_df, tmp_df), axis = 0)

	# sig_res.to_csv('D:\\MCGDYY\\ont_project\\prognosis\\sig_res_' + para + '_detailed.csv', index = None)
	# final_df.to_csv('D:\\MCGDYY\\ont_project\\prognosis\\FC_boxplot_' + para + '.csv', index = None)

	return sig_res


# extract sequence of prognosis associated novel transcripts

call_os = processor('os')
os_list = list(call_os['transcript'])
call_rfs = processor('rfs')
rfs_list = list(call_rfs['transcript'])
union = set(os_list + rfs_list)
union = [i.split('_')[0] for i in union]

nt_fa = ''
aa_fa = ''
with open('D:\\MCGDYY\\ont_project\\prog_assoc_trans\\novel_uniq_RI_real_seq_200bp.fa') as seq_file:
	while True:		
		line = seq_file.readline().strip('\n')
		if not line:
			break
		if line.startswith('>'):
			trans = re.findall('>(.*?)_', line)[0]
			if trans in union:
				nt_name = line
				nt_fa += nt_name + '\n'
				nt_seq = seq_file.readline().strip('\n')
				nt_fa += nt_seq + '\n'

				for i in range(3):
					aa_name = nt_name + '_' + str(i + 1)
					aa_fa += aa_name + '\n'
					aa_seq = Seq(nt_seq[i:]).translate()
					aa_fa += str(aa_seq) + '\n'

with open('D:\\MCGDYY\\ont_project\\prog_assoc_trans\\prog_assoc_novel_nt.fa', 'w') as w:
    w.write(nt_fa)

with open('D:\\MCGDYY\\ont_project\\prog_assoc_trans\\prog_assoc_novel_aa.fa', 'w') as w:
    w.write(aa_fa)