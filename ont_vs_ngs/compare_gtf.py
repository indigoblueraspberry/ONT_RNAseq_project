# find overlapping and unique transcripts from two gtf files generated from gffcompare
# remember to sort gtf files based on chromosome and coordinates

import pandas as pd
import  re


# # convert flair gtf output to a standard gtf format

# flair_gtf = pd.read_csv('D:\\MCGDYY\\ont_project\\flair_out\\flair.collapse.isoforms.stringent.gtf', sep = '\t', header = None)
# flair_sort = flair_gtf.sort_values(by = [0, 3])
# flair_list = []
# for i in flair_sort.index:
# 	info = flair_sort.loc[i, 8]
# 	trans = re.findall(r'transcript_id "(.*?)";', info)[0]
# 	if trans not in flair_list:
# 		flair_list.append(trans)

# flair_compiled = ''
# for i in flair_list:
# 	sub = flair_gtf[flair_gtf[8].str.contains(i)].sort_values(by = 3)
# 	flair_compiled += str(sub.iloc[0, 0]) + '\tFLAIR\ttranscript\t' + str(sub.iloc[0, 3]) + '\t' + str(sub.iloc[len(sub) - 1, 4]) + '\t' + i + '\n'
# 	for k in sub.index:
# 		line = '\t'.join(str(j) for j in sub.loc[k].to_list())
# 		flair_compiled += line + '\n'

# with open('D:\\MCGDYY\\ont_project\\ont_vs_ngs\\assem\\flair_out_compiled.gtf', 'w') as w:
# 	w.write(flair_compiled)

# convert gtf to a list of dicts

def convert_gtf(gtf_path):

	trans_list = []
	trans_dict = {}
	ex_count = 0

	line = True
	with open(gtf_path, 'r') as input_gtf:
		while line:
			line = input_gtf.readline()
			if line.startswith('#') is False and len(line) != 0:
				entry_type = line.split('\t')[2]

				if entry_type == 'transcript':
					trans_list.append(trans_dict)
					chrom = line.split('\t')[0]
					trans_start = int(line.split('\t')[3])
					trans_end = int(line.split('\t')[4])
					trans_type = re.findall(r'class_code "(.*?)";', line)[0]
					if trans_type == '=':
						trans_type = 'annotated'
						trans_id = re.findall(r'cmp_ref "(.*?)";', line)[0]
					elif trans_type == 'u':
						trans_type = 'intergenic'
						trans_id = 'intergenic'
					else:
						trans_type = 'novel_genic'
						trans_id = re.findall(r'cmp_ref "(.*?)";', line)[0]

					trans_dict = {'trans_type': trans_type, 'trans_id': trans_id, 'chrom': chrom, 'trans_start': trans_start, 'trans_end': trans_end}
					ex_count = 0

				if entry_type == 'exon':
					ex_count += 1
					ex_start = int(line.split('\t')[3])
					ex_end = int(line.split('\t')[4])
					trans_dict['ex_start' + str(ex_count)] = ex_start
					trans_dict['ex_end' + str(ex_count)] = ex_end

	# the first dict in the list is empty
	trans_list.pop(0)
	return trans_list

def compare(in_gtf_a, in_gtf_b, bp_diff = 5):

	# create a copy of list a for modification later, otherwise modification during a for loop will cause a problem
	in_gtf_a_cp = in_gtf_a.copy()
	# retain only transcripts from a
	intersec_one = []
	# retain transcript info from both a and b
	intersec_both = []
	total = len(in_gtf_a)
	counter = 0
	counter_list = []

	for entry_a in in_gtf_a:
		len_a = len(entry_a)
		chrom_a = entry_a['chrom']
		trans_start_a = entry_a['trans_start']
		trans_end_a = entry_a['trans_end']
		tried = False
		
		for entry_b in in_gtf_b:
			# try to broadly match two entries
			if chrom_a == entry_b['chrom']:
				tried = True			
				if len_a == len(entry_b):
					# try to fine match two entries
					# two INTRON chains less than diff bp are the same transcript
					# two free ends of terminal exons less than 100 bp are considered the same, this also works for single-exon transcript
					go = True
					for i in range(5, len_a):
						diff = abs(list(entry_a.values())[i] - list(entry_b.values())[i])
						print(diff)
						if i == 5 and diff <= 100:
							go = True
						elif 5 < i < len_a - 1 and diff <= bp_diff:	# base pair difference default is 5 bp
							go = True
						elif i == len_a - 1 and diff <= 100:
							go = True
						else:
							go = False
							break
					# found the match!
					if go:
						intersec_one.append(entry_a)
						intersec_both.append(str(entry_a) + str(entry_b))
						in_gtf_a_cp.remove(entry_a)
						in_gtf_b.remove(entry_b)	# the for loop having list b will break anyway, so no need to use a copy of list b
						break
			
			# if you have tried within the same chromosome, and then the new chromosomes are different, 
			# means you have passed the current chromosome, so no match			
			else:
				if tried:
					break

		# visualize with an odometer
		counter += 1
		percent = round(100*counter/total)		
		if percent not in counter_list:
			 counter_list.append(percent)
			 print(percent)

	return intersec_one, intersec_both, in_gtf_a_cp, in_gtf_b

def write_it_out(in_dict, out_path):

	holder = ''
	for i in in_dict:
		line = '\t'.join(str(k) for k in i.values())
		holder += line + '\n'

	with open(out_path, 'w') as w:
		w.write(holder)


in_gtf_st = convert_gtf('D:\\MCGDYY\\ont_project\\test\\sclp.sub.gtf')
print('list st done')
in_gtf_sclp = convert_gtf('D:\\MCGDYY\\ont_project\\test\\sclp.sub.gtf')
print('list sclp done')
in_gtf_flair = convert_gtf('D:\\MCGDYY\\ont_project\\test\\sclp.sub.gtf')
print('list flair done')

st_sclp, st_sclp_detail, uniq_st, uniq_sclp = compare(in_gtf_st, in_gtf_sclp)
print('st and sclp done')
ngs_flair, ngs_flair_detail, uniq_ngs, uniq_flair = compare(st_sclp, in_gtf_flair)
print('ngs and flair done')

write_it_out(st_sclp, 'D:\\MCGDYY\\ont_project\\test\\st_sclp.txt')
write_it_out(uniq_st, 'D:\\MCGDYY\\ont_project\\test\\uniq_st.txt')
write_it_out(uniq_sclp, 'D:\\MCGDYY\\ont_project\\test\\uniq_sclp.txt')
lines = ''
for line in st_sclp_detail:
	lines += line + '\n'
with open('D:\\MCGDYY\\ont_project\\test\\st_sclp_detail.txt', 'w') as w:
	w.write(lines)

write_it_out(ngs_flair, 'D:\\MCGDYY\\ont_project\\test\\ngs_flair.txt')
write_it_out(uniq_ngs, 'D:\\MCGDYY\\ont_project\\test\\uniq_ngs.txt')
write_it_out(uniq_flair, 'D:\\MCGDYY\\ont_project\\test\\uniq_flair.txt')
lines = ''
for line in ngs_flair_detail:
	lines += line + '\n'
with open('D:\\MCGDYY\\ont_project\\test\\ngs_flair_detail.txt', 'w') as w:
	w.write(lines)