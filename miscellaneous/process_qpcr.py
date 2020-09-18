import pandas as pd

sample_list = ['1122', '1123', '1124', '1130', '1131', '1133',
				'1135', '1136', '1137', '1138', '1144', '1146',
				'1147', '1149', '1151', '1152', '1153', '1159']
targets = ['F', 'BA']
ref = '18s'

raw_df = pd.read_excel('D:\\MCGDYY\\ont_project\\QPCR\\sum_table.xlsx', names = ['sample', 'target', 'Ct'])
# remove rows with undetermined in Ct
raw_df = raw_df[raw_df['Ct'] != 'Undetermined']

# build the dictionary for reference
ref_df = raw_df[raw_df['target'] == '18s']
ref_dict = {}
for i in ref_df.index:
	sample = ref_df.loc[i, 'sample']
	Ct = ref_df.loc[i, 'Ct']
	ref_dict[sample] = Ct

# iterate the raw data frame
sum_df = pd.DataFrame(columns = ['target', 'sample', 'ctl_dCt', 'trt_dCt', 'ctl_2^-dCt', 'trt_2^-dCt'])
count = len(sum_df)
for target in targets:
	for sample in sample_list:
		sum_df.loc[count, 'target'] = target
		sum_df.loc[count, 'sample'] = sample
		for kind in ['N', 'T']:
			temp_df = raw_df[(raw_df['sample'] == sample + kind) & (raw_df['target'] == target)]
			mean_Ct = temp_df['Ct'].mean()
			delta_Ct = mean_Ct - ref_dict[sample + kind]
			if kind == 'N':
				sum_df.loc[count, 'ctl_dCt'] = delta_Ct
				sum_df.loc[count, 'ctl_2^-dCt'] = 2**-delta_Ct
			if kind == 'T':
				sum_df.loc[count, 'trt_dCt'] = delta_Ct
				sum_df.loc[count, 'trt_2^-dCt'] = 2**-delta_Ct
		count += 1

sum_df.to_csv('D:\\MCGDYY\\ont_project\\QPCR\\final.csv')

