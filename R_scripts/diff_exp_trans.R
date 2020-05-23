library(DESeq2)
library(ggplot2)
library(tximport)
library(clusterProfiler)
library(org.Hs.eg.db)

###################### Differential expressed genes

df <- read.csv('D:\\MCGDYY\\ont_project\\scripts_to_results\\tumor_exp_median\\tsts_10x_30_median.csv', 
               check.names = FALSE, row.names = 1)
colnames(df)[15:length(colnames(df))] <- paste('TST_', colnames(df)[15:length(colnames(df))] , sep = '')
colnames(df)[15:length(colnames(df))] <- gsub('-', '_', colnames(df)[15:length(colnames(df))] )
sample_names <- row.names(df)

# construct all file paths
files <- file.path('D:\\MCGDYY\\ont_project\\rsem', paste0(sample_names, '-1.genes.results'), fsep = '\\')
names(files) <- sample_names
txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
# let 0 gene length be 1 (probably happening in novel genes?),
# otherwise this will raise an error later
txi$length[txi$length == 0] <- 1
coldata <- subset.data.frame(df, 'TST_0f7926c1_b340_4a18_9003_767a964dc9ca_status', subset = TRUE)

# analyzed by DEseq2
dds <- DESeqDataSetFromTximport(txi, coldata, design = ~ TST_0f7926c1_b340_4a18_9003_767a964dc9ca_status)
dds <- DESeq(dds)
res <- results(dds, contrast = c('TST_0f7926c1_b340_4a18_9003_767a964dc9ca_status', 'H', 'L'))

# compile and filter results
all_df <- merge(as.matrix(res), txi, by = 'row.names')
rownames(all_df) <- all_df$Row.names
all_df <- all_df[,-1]
all_df <- all_df[order(all_df$padj),]
# remove novel genes
anno_genes <- subset(all_df, grepl('^ENSG', row.names(all_df)))

threshold <- as.factor(ifelse(anno_genes$padj < 0.05 & abs(anno_genes$log2FoldChange) > 1,
                              ifelse(anno_genes$log2FoldChange > 1 ,'Up','Down'),'Not'))
ggplot(anno_genes, aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  geom_point()

####################################### GO terms

diff_anno_genes <- subset(anno_genes, pvalue < 0.05 & abs(log2FoldChange) > log2(2))
gene_names <- row.names(diff_anno_genes)
gene_names <- gsub("\\..*", "", gene_names)

go_res <- enrichGO(gene = gene_names, OrgDb = org.Hs.eg.db, keyType = 'ENSEMBL', ont = 'MF')
go_plot <- barplot(go_res, showCategory = 10)

############################# kegg

converted_table <- bitr(gene_names, fromType = 'ENSEMBL', toType = c("SYMBOL", "ENTREZID"), 
                        OrgDb = org.Hs.eg.db)
entrez <- converted_table$ENTREZID

kegg_res <- enrichKEGG(gene = entrez, organism = 'hsa')
kegg_plot <- dotplot(kegg_res, showCategory = 10)

















