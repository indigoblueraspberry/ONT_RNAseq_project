import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def process(in_file):

	anno_len = []
	anno_ex = []
	novel_len = []
	novel_ex = []

	with open(in_file, 'r') as df:
		line = True
		while line:
			line = df.readline()
			items = line.split('\t')
			loc = items[0]
			try:
				trans_len = int(items[4]) - int(items[3])
				trans_ex = (len(items) - 5) / 2
			except IndexError:
				continue

			if loc == 'annotated':
				anno_len.append(trans_len)
				anno_ex.append(trans_ex)
			else:
				novel_len.append(trans_len)
				novel_ex.append(trans_ex)

	return [anno_len, anno_ex, novel_len, novel_ex]

uniq_ngs = process('D:\\MCGDYY\\ont_project\\ont_vs_ngs\\assem\\st_sclp.txt')
uniq_ont = process('D:\\MCGDYY\\ont_project\\ont_vs_ngs\\assem\\st_sclp.txt')
ngs_ont = process('D:\\MCGDYY\\ont_project\\ont_vs_ngs\\assem\\st_sclp.txt')

"""
1	uniq_ont_len	2	uniq_ont_ex

3	intersec_len	4	intersec_ex

5	uniq_ngs_len	6	uniq_ngs_ex

"""
fig, axes = plt.subplots(3, 2, sharex = 'col', sharey = 'row')

# 1
sns.kdeplot(ax = axes[0, 0], data = np.log10(uniq_ont[0]), shade = True, color = "b")
sns.kdeplot(ax = axes[0, 0], data = np.log10(uniq_ont[2]), shade = True, color = "r")

# 2
sns.kdeplot(ax = axes[0, 1], data = np.log10(uniq_ont[1]), shade = True, color = "b")
sns.kdeplot(ax = axes[0, 1], data = np.log10(uniq_ont[3]), shade = True, color = "r")

# 3
sns.kdeplot(ax = axes[1, 0], data = np.log10(ngs_ont[0]), shade = True, color = "b")
sns.kdeplot(ax = axes[1, 0], data = np.log10(ngs_ont[2]), shade = True, color = "r")

# 4
sns.kdeplot(ax = axes[1, 1], data = np.log10(ngs_ont[1]), shade = True, color = "b")
sns.kdeplot(ax = axes[1, 1], data = np.log10(ngs_ont[3]), shade = True, color = "r")

# 5
sns.kdeplot(ax = axes[2, 0], data = np.log10(uniq_ngs[0]), shade = True, color = "b")
sns.kdeplot(ax = axes[2, 0], data = np.log10(uniq_ngs[2]), shade = True, color = "r")

# 6
sns.kdeplot(ax = axes[2, 1], data = np.log10(uniq_ngs[1]), shade = True, color = "b")
sns.kdeplot(ax = axes[2, 1], data = np.log10(uniq_ngs[3]), shade = True, color = "r")

axes[1, 0].set_ylabel('Density')
axes[2, 0].set_xlabel('Log10(transcript length in bp)')
axes[2, 1].set_xlabel('Exon number')

plt.show()