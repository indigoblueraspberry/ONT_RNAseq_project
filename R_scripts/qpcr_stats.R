library(readxl)
library(stringr)
library(FSA)
library(ggplot2)
library(Rmisc)

df <- read_xlsx('D:/MCGDYY/ont_stuff/MS/qpcr_quan.xlsx', sheet = 'clean', col_names = TRUE)

for (i in unique(df$Target_Name)){
  
  each_gene <- subset(df, df$Target_Name == i)
  each_gene$Sample_Name <- str_split_fixed(each_gene$Sample_Name, pattern = '', n = 4)[,4]
  each_gene$Sample_Name <- as.factor(each_gene$Sample_Name)
  
  #kruskal-wallis
  KW <- kruskal.test(each_gene$dCT_18s_quan ~ each_gene$Sample_Name)
  DT <- dunnTest(each_gene$dCT_18s_quan ~ each_gene$Sample_Name, method = 'bh')
  
  print('############################################################')
  print(i)
  print('############################################################')
  print(KW)
  print(DT)

sum <- summarySE(data = each_gene, measurevar = 'dCT_18s_quan', groupvars = 'Sample_Name')

ggplot(sum, aes(x = Sample_Name,
                y = dCT_18s_quan)) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymin = dCT_18s_quan - se, ymax = dCT_18s_quan + se), 
                width = 0.5)+
  scale_x_discrete(limits = c('N', 'T', 'M')) +
  theme_classic() +
  xlab('') +
  ylab('Quantification') +
  theme(axis.text = element_text(color = "black", size = 10))

ggsave(paste(i, '.png', sep = ''), path = 'D:/MCGDYY/ont_stuff/MS/', device = 'png' )

}

