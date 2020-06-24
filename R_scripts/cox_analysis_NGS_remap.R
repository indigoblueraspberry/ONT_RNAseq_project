library(survival)
library(survminer)
library(stringr)
library(ggplot2)
library(ggpubr)

###################################################################################
####################################### analyzing survival ########################
###################################################################################

# categorized based on median expression in tumor

df <- read.csv('D:\\MCGDYY\\ont_project\\quantification\\t_exp_median.csv', 
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
#write.csv(sig_res_os, 'D:\\MCGDYY\\ont_project\\prognosis\\sig_res_os.csv', row.names = FALSE)

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
#write.csv(sig_res_rfs, 'D:\\MCGDYY\\ont_project\\prognosis\\sig_res_rfs.csv', row.names = FALSE)

# output the union of os and rfs
union_os_rfs <- merge(sig_res_os, sig_res_rfs, by = 'transcript', all = TRUE, suffixes = c('_os', '_rfs'))
union_os_rfs[is.na(union_os_rfs)] <- '*'
#write.csv(union_os_rfs, 'D:\\MCGDYY\\ont_project\\prognosis\\union_os_rfs.csv', row.names = FALSE)

##################################################
######################### plotting os ############
##################################################

# use KM to plot H group vs. L group

fit <- survfit(Surv(OS_time, OS) ~ trans_8a16c5d1_9028_4073_a601_1a40f8d5d134_status, data = df)
ggsurvplot(fit, pval = TRUE)

# plot uni cox log2HR distribution

p1 <- ggplot(data = sig_res_os, aes(x = order, y = log2HR)) +
      geom_point(shape = 19) +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())

# plot p-value of uni cox

p2 <- ggplot(sig_res_os, aes(x = order, y = 1, fill = p.value)) +
      geom_tile() +
      scale_fill_gradientn(colours = c("red", "white")) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_rect(fill = NA),
            legend.direction = 'horizontal',
            legend.position = 'bottom',
            legend.title = element_blank()) +
      xlab('') +
      ylab('p-value') +
      coord_fixed(ratio = 10)

# jump out of R and use python to process union_os_rfs.csv and get FC_info
fc_info <- read.csv('D:\\MCGDYY\\ont_project\\prognosis\\FC_boxplot_os.csv', header = TRUE)

# plot log2FC and lable up and down regulations

p3 <- ggplot(data = fc_info, aes(x = order, y = log2FC, group = order, fill = status)) +
  geom_boxplot() +
  scale_fill_manual(breaks=c("up","no_diff","down"), values=c("blue", "green", 'red')) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'right')

ggarrange(p1, p2, p3, ncol = 1, nrow = 3)

##################################################
######################### plotting rfs ###########
##################################################

# use KM to plot H group vs. L group

fit <- survfit(Surv(RFS_time, RFS) ~ trans_8a16c5d1_9028_4073_a601_1a40f8d5d134_status, data = df)
ggsurvplot(fit, pval = TRUE)

# plot uni cox log2HR distribution

p4 <- ggplot(data = sig_res_rfs, aes(x = order, y = log2HR)) +
  geom_point(shape = 19) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# plot p-value of uni cox

p5 <- ggplot(sig_res_rfs, aes(x = order, y = 1, fill = p.value)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("red", "white")) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.direction = 'horizontal',
        legend.position = 'bottom',
        legend.title = element_blank()) +
  xlab('') +
  ylab('p-value') +
  coord_fixed(ratio = 10)

# jump out of R and use python to process union_os_rfs.csv and get FC_info
fc_info <- read.csv('D:\\MCGDYY\\ont_project\\prognosis\\FC_boxplot_rfs.csv', header = TRUE)

# plot log2FC and lable up and down regulations

p6 <- ggplot(data = fc_info, aes(x = order, y = log2FC, group = order, fill = status)) +
  geom_boxplot() +
  scale_fill_manual(breaks=c("up","no_diff","down"), values=c("blue", "green", 'red')) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'right')

ggarrange(p4, p5, p6, ncol = 1, nrow = 3)