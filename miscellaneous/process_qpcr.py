import pandas as pd

sample_list = ['1111', '1112', '1113', '1114', '1118', '1119',
				'1122', '1123']
targets = ['E']
ref = 'GAPDH'

raw_df = pd.read_excel('D:\\MCGDYY\\ont_project\\QPCR\\primerE.xlsx', names = ['sample', 'target', 'Ct'])
# remove rows with undetermined in Ct
raw_df = raw_df[raw_df['Ct'] != 'Undetermined']

# iterate the raw data frame
sum_df = pd.DataFrame(columns = ['target', 'sample', 'ctl_ref_Ct', 'ctl_target_Ct', 'ctl_2^-dCt', 'trt_ref_Ct', 'trt_target_Ct', 'trt_2^-dCt'])
count = len(sum_df)
for target in targets:
	for sample in sample_list:
		sum_df.loc[count, 'target'] = target
		sum_df.loc[count, 'sample'] = sample
		for kind in ['N', 'T']:
			target_df = raw_df[(raw_df['sample'] == sample + kind) & (raw_df['target'] == target)]
			ref_df = raw_df[(raw_df['sample'] == sample + kind) & (raw_df['target'] == ref)]
			target_Ct = target_df['Ct'].mean()
			ref_Ct = ref_df['Ct'].mean()
			delta_Ct = target_Ct - ref_Ct
			if kind == 'N':
				sum_df.loc[count, 'ctl_ref_Ct'] = ref_Ct
				sum_df.loc[count, 'ctl_target_Ct'] = target_Ct
				sum_df.loc[count, 'ctl_2^-dCt'] = 2**-delta_Ct
			if kind == 'T':
				sum_df.loc[count, 'trt_ref_Ct'] = ref_Ct
				sum_df.loc[count, 'trt_target_Ct'] = target_Ct
				sum_df.loc[count, 'trt_2^-dCt'] = 2**-delta_Ct
		count += 1

sum_df.to_csv('D:\\MCGDYY\\ont_project\\QPCR\\final.csv')

