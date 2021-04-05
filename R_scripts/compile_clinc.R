df <- readxl::read_excel('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\final.xlsx')

df$HBV[df$HBV == 1] <- 'yes'
df$HBV[df$HBV == 0] <- 'no'
df$microsatellite[df$microsatellite == '有'] <- 'yes'
df$microsatellite[df$microsatellite == '无'] <- 'no'
df$PVTT[df$PVTT == '有'] <- 'yes'
df$PVTT[df$PVTT == '无'] <- 'no'
df$envelop[df$envelop == '完整'] <- 'complete'
df$envelop[df$envelop == '无'] <- 'no'
df$envelop[df$envelop == '不完整'] <- 'incomplete'
df$gender <- ifelse(df$gender == '男', 'male', 'female')
df$age <- ifelse(df$age >= 65, 'old', 'young')
df$smoke[df$smoke == 1] <- 'yes'
df$smoke[df$smoke == 0] <- 'no'
df$alcohol[df$alcohol == 1] <- 'yes'
df$alcohol[df$alcohol == 0] <- 'no'
df$size <- ifelse(df$size >= 5, 'large', 'small')
df$vascular_invasion[df$vascular_invasion == 1] <- 'yes'
df$vascular_invasion[df$vascular_invasion == 0] <- 'no'
df$differentiation[df$differentiation == 1] <- 'low'
df$differentiation[df$differentiation == 2] <- 'high'
df$cirrhosis[df$cirrhosis == '结节型'] <- 'yes'
df$cirrhosis[df$cirrhosis == '无'] <- 'no'

char_list <- c('age', 'gender', 'size', 'HBV', 'PVTT', 
               'vascular_invasion', 'differentiation', 'cirrhosis', 'BCLC')

result <- sapply(char_list, function(x){
  data_table <- table(df[[x]], df$`b5082e5d-a8e4-4986-841d-9f416d0d57fe_status`)
  p <- fisher.test(data_table)$p.value
  res <- c(x, p)
  print(x)
  print(data_table)
  print(signif(p, digits = 3))
  print('##############################################')
  return(res)
  })



# b5082e5d-a8e4-4986-841d-9f416d0d57fe is CDO1
fisher.test(table(df$HBV, df$`b5082e5d-a8e4-4986-841d-9f416d0d57fe_status`))
# 9c1c526a-89d8-4fed-afb7-8a1eb93189a3 is CYP2A6
