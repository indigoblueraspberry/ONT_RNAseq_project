# 6-frame-shift translation of all novel transcripts and further process into fragments stopping at the stop codon

from Bio.Seq import Seq

def fragmenting (aa_seq, pep_name):
    fragment = ''
    collector = ''
    frag_count = 1
    
    for aa in aa_seq:
        if aa != '*':
            fragment += aa
        else:
            # only retain fragments > 6 aa
            if len(fragment) >= 6:
                # only retain fragments containing K or R
                if 'K' in fragment or 'R' in fragment:
                    collector += pep_name + '_' + str(frag_count) + '\n'
                    collector += fragment + '\n'
                    frag_count += 1
                    
            fragment = ''

        
    return collector

novel = ''
all_fragments = ''

with open('D:\\MCGDYY\\ont_project\\flair_out\\flair.collapse.isoforms.stringent.fa', 'r') as all_trans:
    lines = all_trans.read().splitlines()
    
    for i in range(len(lines)):
        # look for novel transcripts
        if i % 2 == 0 and 'ENST' not in lines[i]:
            trans_name = lines[i]
            
            for k in range(3):
                # forward translation
                rna_seq_f = Seq(lines[i+1][k:])
                peptide_seq_f = rna_seq_f.translate()
                peptide_name_f =  trans_name + '_f_' + str(k+1)             
                
                novel += peptide_name_f + '\n'
                novel += str(peptide_seq_f) + '\n'
                # fragment each transcript 
                all_fragments += fragmenting(peptide_seq_f, peptide_name_f)            
                
                # reverse translation
                rna_seq_r = Seq(lines[i+1]).reverse_complement()[k:]
                peptide_seq_r = rna_seq_r.translate()
                peptide_name_r = trans_name + '_r_' + str(k+1)
                
                novel += peptide_name_r + '\n'
                novel += str(peptide_seq_r) + '\n'
                # fragment each transcript 
                all_fragments += fragmenting(peptide_seq_r, peptide_name_r)        
            

with open('D:\\MCGDYY\\ont_project\\MS\\novel_peptides.fa', 'w') as w:
    w.write(novel)


with open('D:\\MCGDYY\\ont_project\\MS\\novel_peptide_fragments.fa', 'w') as w:
    w.write(all_fragments)  