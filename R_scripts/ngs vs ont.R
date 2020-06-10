library(stringr)
library(ggplot2)
library(clusterProfiler) # for GO analysis
library(org.Hs.eg.db) # human gene annotation db
library(pathview) # view individual pathway


############################################### gene

# load ont gene count by featureCount
ONT_gene <- read.csv(file = 'D:\\MCGDYY\\ont_project\\ont_vs_ngs\\ont_gene_counts.txt', 
                     sep = '\t', header = TRUE, skip = 1, stringsAsFactors = FALSE)
ONT_gene <- ONT_gene[,-(2:6)]
colnames(ONT_gene)[2:length(ONT_gene)] <- str_split_fixed(colnames(ONT_gene)[2:length(ONT_gene)], "\\.", n = Inf )[,8]

# load NGS gene quantification by RSEM
NGS_gene_N <- read.csv(file = 'D:\\MCGDYY\\ont_project\\ont_vs_ngs\\N12_N.genes.results', 
                       sep = '\t', header = TRUE, stringsAsFactors = FALSE )
NGS_gene_T <- read.csv(file = 'D:\\MCGDYY\\ont_project\\ont_vs_ngs\\N12_T.genes.results', 
                       sep = '\t', header = TRUE, stringsAsFactors = FALSE )
NGS_gene_M <- read.csv(file = 'D:\\MCGDYY\\ont_project\\ont_vs_ngs\\N12_M.genes.results', 
                       sep = '\t', header = TRUE, stringsAsFactors = FALSE )
NGS_gene <- data.frame(cbind(NGS_gene_N$gene_id, NGS_gene_N$FPKM, NGS_gene_T$FPKM, NGS_gene_M$FPKM),
                       stringsAsFactors = FALSE)
colnames(NGS_gene) <- c('Geneid', 'NGS_N', 'NGS_T', 'NGS_M')
NGS_gene$NGS_N <- as.numeric(NGS_gene$NGS_N)
NGS_gene$NGS_T <- as.numeric(NGS_gene$NGS_T)
NGS_gene$NGS_M <- as.numeric(NGS_gene$NGS_M)

# merge ONT and NGS data
master_gene <- merge(NGS_gene, ONT_gene, by = 'Geneid')

ggplot(mapping = aes(log(master_gene$NGS_N +1 ), log(master_gene$N12_N_prim + 1))) +
  geom_point() +
  xlab('log(TPM+1)') +
  ylab('log(counts+1)') +
  #coord_cartesian(xlim = c(0,11), ylim = c(0,11)) +
  theme(axis.text = element_text(color = "black", 
                                 size = 15)) +
  theme_classic()

cor.test(log(master_gene$NGS_N +1 ), log(master_gene$N12_N_prim + 1), method = 'pearson')
cor.test(master_gene$NGS_N, master_gene$N12_N_prim, method = 'spearman')





# transcript

# this ont input is produced by FLAIR
ONT_trans <- read.csv(file = 'D:/MCGDYY/ont_stuff/new_stuff/counts_matrix.tsv', 
                      sep = '\t', header = TRUE, stringsAsFactors = FALSE)
ONT_trans <- mutate(ONT_trans, 
                    trans = str_split_fixed(ONT_trans$id, '_', 2)[,1],
                    gene = str_split_fixed(ONT_trans$id, '_', 2)[,2])

# retain only trans having > 3 reads in at least one sample
#ONT_trans <- subset(ONT_trans, ONT_N >= 3 | ONT_T >= 3 | ONT_M >= 3)

NGS_trans_N <- read.csv(file = 'D:/MCGDYY/ont_stuff/new_stuff/N12_N.isoforms.results', 
                  sep = '\t', header = TRUE, stringsAsFactors = FALSE )
NGS_trans_T <- read.csv(file = 'D:/MCGDYY/ont_stuff/new_stuff/N12_T.isoforms.results', 
                  sep = '\t', header = TRUE, stringsAsFactors = FALSE )
NGS_trans_M <- read.csv(file = 'D:/MCGDYY/ont_stuff/new_stuff/N12_M.isoforms.results', 
                  sep = '\t', header = TRUE, stringsAsFactors = FALSE )
NGS_trans <- data.frame(cbind(NGS_trans_N$transcript_id, NGS_trans_N$TPM, NGS_trans_T$TPM, NGS_trans_M$TPM),
                       stringsAsFactors = FALSE)
colnames(NGS_trans) <- c('trans', 'NGS_N', 'NGS_T', 'NGS_M')
NGS_trans$NGS_N <- as.numeric(NGS_trans$NGS_N)
NGS_trans$NGS_T <- as.numeric(NGS_trans$NGS_T)
NGS_trans$NGS_M <- as.numeric(NGS_trans$NGS_M)

master_trans <- merge(NGS_trans, ONT_trans, by = 'trans')

ggplot(mapping = aes(log(master_trans$NGS_M +1 ), log(master_trans$N12_M_batch1 + 1))) +
  geom_point() +
  xlab('log(TPM+1)') +
  ylab('log(counts+1)') +
  #coord_cartesian(xlim = c(0,11), ylim = c(0,11)) +
  theme(axis.text = element_text(color = "black", 
                                 size = 10)) +
  theme_classic()
cor.test(log(master_trans$NGS_M +1 ), log(master_trans$N12_M_batch1 + 1), method = 'pearson')

# GO and KEGG

TST <- subset(ONT_trans, ONT_N == 0 & ONT_T > 0)
TST_of_anno_gene <- subset(TST, grepl('^ENSG', TST$gene))
# remove duplicated genes
#candidate_genes <- TST_of_anno_gene[!duplicated(TST_of_anno_gene$gene),]

candidate_genes <- TST_of_anno_gene[order(-TST_of_anno_gene$ONT_T),]

first100 <- candidate_genes$gene[1:100]
ensembl <- gsub("\\..*", "", candidate_genes$gene)
entrez <- bitr(ensembl, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')

go <- enrichGO(gene = ensembl, OrgDb = org.Hs.eg.db, keyType = 'ENSEMBL', ont = 'ALL')
barplot(go, showCategory = 10)

kegg <- enrichKEGG(entrez$ENTREZID, organism = 'hsa', keyType = 'kegg')
dotplot(kegg, showCategory = 10)

dim(go[go$ONTOLOGY == 'BP'])
dim(go[go$ONTOLOGY == 'MF'])
dim(go[go$ONTOLOGY == 'CC'])

MST <- subset(ONT_trans, ONT_N == 0 & ONT_T == 0 & ONT_M > 0)
MST_of_anno_gene <- subset(MST, grepl('^ENST', MST$trans))
#candidate_genes <- MST_of_anno_gene[!duplicated(MST_of_anno_gene$gene),]
candidate_genes <- MST_of_anno_gene[order(-MST_of_anno_gene$ONT_M),]

ensembl <- gsub("\\..*", "", MST_of_anno_gene$trans)
go <- enrichGO(gene = ensembl, OrgDb = org.Hs.eg.db, keyType = 'ENSEMBL', ont = 'ALL')
barplot(go, showCategory = 10)







