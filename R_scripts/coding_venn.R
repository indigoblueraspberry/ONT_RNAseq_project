library(VennDiagram)
library(readxl)

nc <- read_excel('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\lncRNA\\data_4_venn.xlsx',
                   sheet = 'noncoding')
draw.triple.venn(
  area1 = nc$CPC2,
  area2 = nc$CNCI,
  area3 = nc$PLEK,
  n12 = nc$CNCI_CPC2,
  n23 = nc$CNCI_PLEK,
  n13 = nc$CPC2_PLEK,
  n123 = nc$CPC2_CNCI_PLEK,
  category = c('CPC2', 'CNCI', 'PLEK'),
  fill = c("light blue", "#99ff99", "pink"),
  cex = 3,
  cat.cex = 3
)

c <- read_excel('D:\\MCGDYY\\ont_project\\NovelQuant_pipeline\\lncRNA\\data_4_venn.xlsx',
                 sheet = 'coding')
draw.triple.venn(
  area1 = c$CPC2,
  area2 = c$CNCI,
  area3 = c$PLEK,
  n12 = c$CNCI_CPC2,
  n23 = c$CNCI_PLEK,
  n13 = c$CPC2_PLEK,
  n123 = c$CPC2_CNCI_PLEK,
  category = c('CPC2', 'CNCI', 'PLEK'),
  fill = c("light blue", "#99ff99", "pink"),
  cex = 3,
  cat.cex = 3
)