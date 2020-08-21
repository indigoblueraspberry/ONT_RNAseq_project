library(DESeq2)
library(tximport)
library(ggplot2)

# iterate through RNAseq quantification

files <- list.files('D:/MCGDYY/ont_project/rsem', pattern = 'isoform', full.names = TRUE)
names(files) <- list.files('D:/MCGDYY/ont_project/rsem', pattern = 'isoform')
txi.rsem <- tximport(files, type = 'rsem', txIn = TRUE, txOut = TRUE)
col_data <- data.frame(condition = factor(rep(c('T', 'N'), 60)))
rownames(col_data) <- colnames(txi.rsem$counts)

dds <- DESeqDataSetFromTximport(txi.rsem, col_data, ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c('condition', 'T', 'N'))

all_df <- as.data.frame(res)
all_df <- na.omit(all_df)
all_df <- subset(all_df, grepl('ENST', row.names(all_df)))
all_df <- all_df[order(all_df$padj),]

# find quantification for DE mRNA (RSEM doesn't work well on quantifying miRNA due to short length)

mRNA_file <- read.csv('D:\\MCGDYY\\ont_project\\lists\\anno_mRNA_trans.txt', sep = '\t', header = FALSE)
all_mRNA <- merge(all_df, mRNA_file, by.x = 0, by.y = 'V2', all = FALSE)
DE_mRNA <- subset(all_mRNA, padj < 0.05 & abs(log2FoldChange) > log2(1.5))
#write.csv(DE_mRNA, 'D:\\MCGDYY\\ont_project\\prog_assoc_trans\\novel_lncRNA\\DE_anno_mRNA.csv', row.names = FALSE)

# find quantification for DE lncRNA

lncRNA_file <- read.csv('D:\\MCGDYY\\ont_project\\lists\\anno_lncRNA_trans.txt', sep = '\t', header = FALSE)
all_lncRNA <- merge(all_df, lncRNA_file, by.x = 0, by.y = 'V2', all = FALSE)
DE_lncRNA <- subset(all_lncRNA, padj < 0.05 & abs(log2FoldChange) > log2(1.5))
#write.csv(DE_lncRNA, 'D:\\MCGDYY\\ont_project\\prog_assoc_trans\\novel_mRNA\\DE_anno_lncRNA.csv', row.names = FALSE)

# plot a volcano plot for DE mRNA

threshold <- as.factor(ifelse(all_mRNA$padj < 0.05 & abs(all_mRNA$log2FoldChange) >= 0.58 ,
                              ifelse(all_mRNA$log2FoldChange >= 0.58 ,'up','down'),'no_diff'))
ggplot(all_mRNA, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point() +
  geom_vline(xintercept = c(-0.58, 0.58),
             linetype = 'dotted') +
  geom_hline(yintercept = 1.3,
             linetype = 'dotted') +
  theme_bw() +
  theme(legend.title=element_blank(),
        axis.text = element_text(color = 'black'),
        panel.grid = element_blank()) +
  xlim(-30, 30) +
  xlab('log2FC') +
  ylab('-log10(adj. p-value)') +
  scale_color_manual(values = c("blue", "green", "red"))

  




