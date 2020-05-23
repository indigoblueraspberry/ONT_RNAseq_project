library(clusterProfiler)
library(org.Hs.eg.db)


df <- read.csv('D:\\MCGDYY\\ont_stuff\\new_stuff\\TSTs_stringent.txt', sep = '\t')

gene_names <- df[,12]
gene_names <- gsub("\\..*", "", gene_names)

# look for annotated genes
ensembl <- subset(gene_names, grepl('^ENSG', gene_names))

# Biological Processes has the most enriched genes
go_res <- enrichGO(gene = ensembl, OrgDb = org.Hs.eg.db, keyType = 'ENSEMBL', ont = 'BP')
go_plot <- barplot(go_res, showCategory = 10)

# convert ensemble names to entrez ID
converted_table <- bitr(ensembl, fromType = 'ENSEMBL', toType = c("SYMBOL", "ENTREZID"), 
                        OrgDb = org.Hs.eg.db)
entrez <- converted_table$ENTREZID

kegg_res <- enrichKEGG(gene = entrez, organism = 'hsa')
kegg_plot <- dotplot(kegg_res, showCategory = 10)

