# find overlapping and unique transcripts from two gtf files
# remember to sort gtf files based on chromosome and coordinates

import  re


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

					if 'cmp_ref' in line:	# for scallop whose transcript id is in cmp_ref
						trans_id = re.findall(r'cmp_ref "(.*?)";', line)[0]
					else:
						trans_id = re.findall(r'transcript_id "(.*?)";', line)[0]
					trans_dict = {'trans_id': trans_id, 'chrom': chrom, 'trans_start': trans_start, 'trans_end': trans_end}
					ex_count = 0

				if entry_type == 'exon':
					ex_count += 1
					ex_start = int(line.split('\t')[3])
					ex_end = int(line.split('\t')[4])
					trans_dict['ex_start' + str(ex_count)] = ex_start
					trans_dict['ex_end' + str(ex_count)] = ex_end

	# the first dict is empty
	trans_list.pop(0)
	return trans_list


in_gtf_a = convert_gtf('D:\\MCGDYY\\ont_project\\test\\sclp.sub.gtf')
print('list a done')
in_gtf_b = convert_gtf('D:\\MCGDYY\\ont_project\\test\\st.sub.gtf')
print('list b done')
intersec = []

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
				# try to fine match two entries, two entries less than 10 bp are the same transcript
				for i in range(3, len_a):
					diff = abs(list(entry_a.values())[i] - list(entry_b.values())[i])
					if diff <= 10:
						go = True
					else:
						go = False
						break
				# found the match!
				if go:
					intersec.append(str(entry_a) + str(entry_b))
					in_gtf_a.remove(entry_a)
					in_gtf_b.remove(entry_b)
					break
		
		# if you have tried within the same chromosome, and then the new chromosomes are different, 
		# means you have passed the current chromosome, so no match			
		else:
			if tried:
				break

unique_a = in_gtf_a
unique_b = in_gtf_b