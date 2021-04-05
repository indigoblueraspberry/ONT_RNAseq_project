import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats


raw_df = pd.read_excel('D:\\MCGDYY\\ont_project\\QPCR\\raw_all_E_BA_18S.xlsx', names = ['sample', 'target', 'Ct'])
# remove rows with undetermined in Ct
raw_df = raw_df[raw_df['Ct'] != 'Undetermined']

targets = ['E']
ref = 'GAPDH'
sample_list = []
for sample in raw_df['sample']:
	sample = sample.split('-')[0]
	sample_list.append(sample)
sample_list = set(sample_list)


# iterate the raw data frame
sum_df = pd.DataFrame(columns = ['target', 'sample', 'ctl_ref_Ct', 'ctl_target_Ct', 'ctl_2^-dCt', 'trt_ref_Ct', 'trt_target_Ct', 'trt_2^-dCt'])
count = len(sum_df)
for target in targets:
	for sample in sample_list:
		sum_df.loc[count, 'target'] = target
		sum_df.loc[count, 'sample'] = sample
		for kind in ['-1', '-2']:
			target_df = raw_df[(raw_df['sample'] == sample + kind) & (raw_df['target'] == target)]
			ref_df = raw_df[(raw_df['sample'] == sample + kind) & (raw_df['target'] == ref)]
			target_Ct = target_df['Ct'].mean()
			ref_Ct = ref_df['Ct'].mean()
			delta_Ct = target_Ct - ref_Ct
			if kind == '-2':
				sum_df.loc[count, 'ctl_ref_Ct'] = ref_Ct
				sum_df.loc[count, 'ctl_target_Ct'] = target_Ct
				sum_df.loc[count, 'ctl_2^-dCt'] = 2**-delta_Ct
			if kind == '-1':
				sum_df.loc[count, 'trt_ref_Ct'] = ref_Ct
				sum_df.loc[count, 'trt_target_Ct'] = target_Ct
				sum_df.loc[count, 'trt_2^-dCt'] = 2**-delta_Ct
		count += 1

sum_df = sum_df.dropna()
# sum_df.to_csv(D:\\MCGDYY\\ont_project\\QPCR\\perimerE_junc1\\tumor_and_normal\\filtered\\compiled.csv', index = False)


stats_res = stats.wilcoxon(sum_df['ctl_2^-dCt'], sum_df['trt_2^-dCt'])
# plot_res = plt.boxplot([sum_df['ctl_2^-dCt'], sum_df['trt_2^-dCt']])

print(stats_res)
