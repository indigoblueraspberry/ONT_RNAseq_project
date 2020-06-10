# This script is used to extract intron coordinates from a gtf file.
# It additionally outputs a file to collect all single-exon transcripts.

import pandas as pd
import re

gtf_path = 'D:\\MCGDYY\\ont_second_run\\gencode.v32.GRCh38.novel.trans.gtf'
out_intron_gtf_path = 'D:\\MCGDYY\\ont_second_run\\introns.gtf'
out_sig_ex_trans_path = 'D:\\MCGDYY\\ont_second_run\\sig_ex.gtf'

# extract lines with exons

gtf = pd.read_csv(gtf_path, sep = '\t', header = None, comment = '#')
exons = gtf[gtf[2] == 'exon']

transcripts = []
for i in exons[8]:
    trans_name = re.findall(r"transcript_id \"(.*?)\";", i)[0]
    transcripts.append(trans_name)
exons[9] = transcripts
# sort the exon table by chromosome, transcript name and then the start position
exons = exons.sort_values(by = [0, 9, 3]).reset_index(drop = True)


# iterate through exons to distinguish between single and multiple -exon transcripts
sig = open(out_sig_ex_trans_path, 'a')
mul = open(out_intron_gtf_path, 'a')
multi = False
meter = 1

try:
    for i in range(len(exons)):
        chromosome = exons.loc[i, 0]
        source = exons.loc[i, 1]
        start_exon = exons.loc[i, 3]
        end_exon = exons.loc[i, 4]
        strand = exons.loc[i, 6]
        info = exons.loc[i, 8]
        this_trans = exons.loc[i, 9]
        next_trans = exons.loc[i + 1, 9]

        # for single-exon transcripts
        if this_trans != next_trans:
            if multi == False:
                sig_ex_trans = str(chromosome) + '\t' + source + '\t' + 'exon' + '\t' + str(start_exon) + '\t' + str(end_exon) + \
                            '\t' + '.' + '\t' + strand + '\t' + '.' + '\t' + info + '\n'
                sig.write(sig_ex_trans)

        # for multiple-exon transcripts
        if this_trans == next_trans:
            multi = True
            start_intron = end_exon + 1
            end_intron = exons.loc[i + 1, 3] - 1
            multi_ex_trans =  str(chromosome) + '\t' + source + '\t' + 'intron' + '\t' + str(start_intron) + '\t' + str(end_intron) + \
                                '\t' + '.' + '\t' + strand + '\t' + '.' + '\t' + info + '\n'
            mul.write(multi_ex_trans)
        else:   # for the last exon in a multiple-exon transcript
            multi = False

        print(str(meter) + '/' + str(len(exons)) + ' done.')
        meter += 1

# for the last line in the gtf file
except KeyError:
    if multi == False:   # if the last line is single_exon transcript
        sig_ex_trans = str(chromosome) + '\t' + source + '\t' + 'exon' + '\t' + str(start_exon) + '\t' + str(end_exon) + \
                            '\t' + '.' + '\t' + strand + '\t' + '.' + '\t' + info + '\n'
        sig.write(sig_ex_trans)
    else:   # if multiple
        print(str(meter) + '/' + str(len(exons)) + ' done.')

sig.close()
mul.close()