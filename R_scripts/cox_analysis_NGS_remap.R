library(survival)
library(survminer)
library(stringr)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)

#################################################################################
############## volcano plot for all novel transcripts ###########################
#################################################################################

DE_table <- read.csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\sum_DE_all_novel.csv')

# volcano plot
threshold <- as.factor(ifelse(DE_table$adj_p < 0.05 & abs(DE_table$log2FC) >= 0.58,
                              ifelse(DE_table$log2FC >= 0.58 ,'up','down'),'no_diff'))
vol <- ggplot(DE_table, aes(x = log2FC, y = -log10(adj_p), color = threshold)) +
  geom_point() +
  geom_vline(xintercept = c(-0.58, 0.58),
             linetype = 'dotted') +
  geom_hline(yintercept = 1.3,
             linetype = 'dotted') +
  theme_bw() +
  theme(legend.title=element_blank(),
        legend.position = c(0.9,0.35),
        legend.box.background = element_rect(colour = 'black'),
        legend.background = element_blank(),
        legend.text = element_text(size = 10),
        axis.text = element_text(color = 'black', size = 15),
        axis.title = element_text(size = 20),
        panel.grid = element_blank()) +
  xlim(-5, 5) +
  xlab('log2FC') +
  ylab('-log10(adj. p-value)') +
  scale_color_manual(values = c("blue", "green", "red"))


###################################################################################
####################################### analyzing survival ########################
###################################################################################

# categorized based on median expression in tumor

df <- read.csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\t_exp_median_DEonly.csv', 
               check.names = FALSE, row.names = 1)
# R can't have special characters and numbers at the beginning of the col names.
colnames(df)[15:length(colnames(df))] <- paste('trans_', colnames(df)[15:length(colnames(df))] , sep = '')
colnames(df)[15:length(colnames(df))] <- gsub('-', '_', colnames(df)[15:length(colnames(df))] )

# univariable cox

variables <- colnames(df)[endsWith(colnames(df), '_status')]
univ_formulas <- sapply(variables,
                        function(x){as.formula(paste('Surv(OS_time, OS) ~ ', x, sep = ''))})
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = df)})
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value <- signif(x$waldtest["pvalue"], digits = 2)
                         wald.test <- signif(x$waldtest["test"], digits = 2)
                         beta <- x$coef[1]
                         HR <- x$coef[2]
                         HR.confint.lower <- x$conf.int[,"lower .95"]
                         HR.confint.upper <- x$conf.int[,"upper .95"]
                         res <- c(beta, HR, HR.confint.lower, HR.confint.upper, wald.test, p.value)
                         names(res)<-c("beta", "HR", "95%_CI_lower", "95%_CI_upper", "wald.test", "p.value")
                         return(res)
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))

# multivariable cox

sig_res <- subset(res, res[,6]<0.05)
sig_var <- str_c(row.names(sig_res), collapse = '+')
multi_formula <- as.formula(paste('Surv(OS_time, OS) ~ ', sig_var, sep = ''))
multi_models <- coxph(multi_formula, data = df)
summary(multi_models)

# post processing

sig_res <- as.data.frame(sig_res)
sig_res <- mutate(sig_res, transcript = row.names(sig_res))
sig_res <- mutate(sig_res, log2HR = log2(HR))
sig_res <- sig_res[order(sig_res$log2HR, decreasing = TRUE),]
sig_res_os <- mutate(sig_res, order = 1:nrow(sig_res))
#write.csv(sig_res_os, 'D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\prognosis\\sig_res_os.csv', row.names = FALSE)


###################################################################################
####################################### analyzing recurrence ######################
###################################################################################

# univariable cox

variables <- colnames(df)[endsWith(colnames(df), '_status')]
univ_formulas <- sapply(variables,
                        function(x){as.formula(paste('Surv(RFS_time, RFS) ~ ', x, sep = ''))})
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = df)})
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value <- signif(x$waldtest["pvalue"], digits = 2)
                         wald.test <- signif(x$waldtest["test"], digits = 2)
                         beta <- x$coef[1]
                         HR <- x$coef[2]
                         HR.confint.lower <- x$conf.int[,"lower .95"]
                         HR.confint.upper <- x$conf.int[,"upper .95"]
                         res <- c(beta, HR, HR.confint.lower, HR.confint.upper, wald.test, p.value)
                         names(res)<-c("beta", "HR", "95%_CI_lower", "95%_CI_upper", "wald.test", "p.value")
                         return(res)
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))

# multivariable cox

sig_res <- subset(res, res[,6]<0.05)
sig_var <- str_c(row.names(sig_res), collapse = '+')
multi_formula <- as.formula(paste('Surv(RFS_time, RFS) ~ ', sig_var, sep = ''))
multi_models <- coxph(multi_formula, data = df)
summary(multi_models)

# post processing

sig_res <- as.data.frame(sig_res)
sig_res <- mutate(sig_res, transcript = row.names(sig_res))
sig_res <- mutate(sig_res, log2HR = log2(HR))
sig_res <- sig_res[order(sig_res$log2HR, decreasing = TRUE),]
sig_res_rfs <- mutate(sig_res, order = 1:nrow(sig_res))
#write.csv(sig_res_rfs, 'D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\prognosis\\sig_res_rfs.csv', row.names = FALSE)

# output the union of os and rfs
union_os_rfs <- merge(sig_res_os, sig_res_rfs, by = 'transcript', all = TRUE, suffixes = c('_os', '_rfs'))
union_os_rfs[is.na(union_os_rfs)] <- '*'
#write.csv(union_os_rfs, 'D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\prognosis\\union_os_rfs.csv', row.names = FALSE)

##################################################
######################### plotting os ############
##################################################

# use KM to plot H group vs. L group

# CDO1
fit <- survfit(Surv(RFS_time, RFS) ~ trans_b5082e5d_a8e4_4986_841d_9f416d0d57fe_status, data = df)
ggsurvplot(fit, pval = TRUE, xlab = 'Days', ylab = 'Recurrence Free Survival (%)',
           legend.labs = c('CDO1-novel low', 'CDO1-novel high'),
           legend.title = '',  palette  = c('blue', 'red'), legend = c(0.7, 0.9),
           font.x = 30, font.y = 30, font.legend = 25, font.tickslab = 25)

# CYP2A6
fit <- survfit(Surv(RFS_time, RFS) ~ trans_9c1c526a_89d8_4fed_afb7_8a1eb93189a3_status, data = df)
ggsurvplot(fit, pval = TRUE, xlab = 'Days', ylab = 'Recurrence Free Survival (%)',
           legend.labs = c('CYP2A6-novel low', 'CYP2A6-novel high'),
           legend.title = '',  palette  = c('blue', 'red'), legend = c(0.7, 0.9),
           font.x = 30, font.y = 30, font.legend = 25, font.tickslab = 25)

fit <- survfit(Surv(OS_time, OS) ~ trans_9c1c526a_89d8_4fed_afb7_8a1eb93189a3_status, data = df)
ggsurvplot(fit, pval = TRUE, xlab = 'Days', ylab = 'Overall Survival (%)',
           legend.labs = c('CYP2A6-novel low', 'CYP2A6-novel high'),
           legend.title = '',  palette  = c('blue', 'red'), legend = c(0.7, 0.2),
           font.x = 30, font.y = 30, font.legend = 25, font.tickslab = 25)

# plot uni cox log2HR distribution

# sig_res_os <- read.csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\prognosis\\sig_res_os.csv')
# novel_gene <- read.csv('D:\\MCGDYY\\ont_project\\lists\\all_novel_list.txt', sep = '\t', header = FALSE)
# gene_symbol <- read.csv('D:\\MCGDYY\\ont_project\\lists\\all_trans_to_gene_list.txt', sep = '\t', header = FALSE)
# gene_symbol <- unique(gene_symbol[c('V1', 'V3')])
# sig_res_os <- left_join(sig_res_os, novel_gene, by = c('transcript' = 'V1'))
# sig_res_os <- left_join(sig_res_os, gene_symbol, by = c('V2' = 'V1'))

p1 <- ggplot(data = sig_res_os, aes(x = order, y = log2HR)) +
  geom_point(size = 1) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlim(0, dim(sig_res_os)[1] + 1)

# # plot p-value of uni cox
# 
# p2 <- ggplot(sig_res_os, aes(x = order, y = 1, fill = p.value)) +
#       geom_tile() +
#       scale_fill_gradientn(colours = c("red", "white")) +
#       theme(panel.grid = element_blank(),
#             panel.background = element_blank(),
#             axis.text = element_blank(),
#             axis.ticks = element_blank(),
#             panel.border = element_rect(fill = NA),
#             legend.direction = 'horizontal',
#             legend.position = 'none',
#             legend.title = element_blank()) +
#       xlab('') +
#       ylab('p-value') +
#       xlim(0, dim(sig_res_os)[1] + 1)

# jump out of R and use python to process union_os_rfs.csv and get FC_info
fc_info <- read.csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\prognosis\\FC_boxplot_os.txt', header = TRUE, sep = '\t')

# plot log2FC and lable up and down regulations

p3 <- ggplot(data = fc_info, aes(x = order, y = log2FC, group = order, fill = status)) +
  geom_boxplot(size = 0.1, outlier.size = 0.1) +
  scale_fill_manual(breaks=c("up", "down"), values=c("blue", "red")) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none') +
  xlim(0, dim(sig_res_os)[1] + 1)

ggarrange(p1, p3, ncol = 1, nrow = 2, heights = c(1, 1), align = 'v')

##################################################
######################### plotting rfs ###########
##################################################

# use KM to plot H group vs. L group

fit <- survfit(Surv(RFS_time, RFS) ~ trans_e420b9f9_dcca_496c_9035_d0678183caeb_status, data = df)
ggsurvplot(fit, pval = TRUE)

# plot uni cox log2HR distribution

p4 <- ggplot(data = sig_res_rfs, aes(x = order, y = log2HR)) +
  geom_point(size = 1) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlim(0, dim(sig_res_rfs)[1] + 1)

# plot p-value of uni cox

# p5 <- ggplot(sig_res_rfs, aes(x = order, y = 1, fill = p.value)) +
#   geom_tile() +
#   scale_fill_gradientn(colours = c("red", "white")) +
#   theme(panel.grid = element_blank(),
#         panel.background = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title.x = element_blank(),
#         panel.border = element_rect(fill = NA),
#         legend.direction = 'horizontal',
#         legend.position = 'none',
#         legend.title = element_blank()) +
#   ylab('p-value') +
#   xlim(0, dim(sig_res_rfs)[1] + 1)

# jump out of R and use python to process union_os_rfs.csv and get FC_info
fc_info <- read.csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\prognosis\\FC_boxplot_rfs.txt', header = TRUE, sep = '\t')

# plot log2FC and lable up and down regulations

p6 <- ggplot(data = fc_info, aes(x = order, y = log2FC, group = order, fill = status)) +
  geom_boxplot(size = 0.1, outlier.size = 0.1) +
  scale_fill_manual(breaks=c("up", "down"), values=c("blue", "red")) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none') +
  xlim(0, dim(sig_res_rfs)[1] + 1)

ggarrange(p4, p6, ncol = 1, nrow = 2, heights = c(1,1), align = 'v')

ggarrange(vol, 
          ggarrange(p4, p6, p1, p3, ncol = 1, nrow = 4, heights = c(1,1,1,1), align = 'v'), ncol = 2, widths = c(1, 2))

####################################################################################
## enrichment analysis of novel prognostic mRNAs               #####################
####################################################################################

parent_genes <- read.csv('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\lncRNA\\intersec_novel_mRNA.csv',header = TRUE)

# enrichment analysis

ensembl <- gsub("\\..*", "", parent_genes$Gene)
entrez <- bitr(ensembl, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')

go <- enrichGO(gene = ensembl, OrgDb = org.Hs.eg.db, keyType = 'ENSEMBL', ont = 'BP')
barplot(go, showCategory = 5)

kegg <- enrichKEGG(entrez$ENTREZID, organism = 'hsa', keyType = 'kegg')
dotplot(kegg, showCategory = 5)