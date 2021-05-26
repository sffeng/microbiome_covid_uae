df_level6 <- read.csv('data/2020-11-09 qiime andreas/level-6.csv',row.names = "index")
df_level5 <- read.csv('data/2020-11-09 qiime andreas/level-5.csv',row.names = "index")
df_level4 <- read.csv('data/2020-11-09 qiime andreas/level-4.csv',row.names = "index")
df_level3 <- read.csv('data/2020-11-09 qiime andreas/level-3.csv',row.names = "index")
df_level2 <- read.csv('data/2020-11-09 qiime andreas/level-2.csv',row.names = "index")

df_kraken6 <- read.delim('data/2020-11-15 kraken/kraken_level_6.csv',sep='\t',row.names="X")
df_kraken5 <- read.delim('data/2020-11-15 kraken/kraken_level_5.csv',sep='\t',row.names="X")
df_kraken4 <- read.delim('data/2020-11-15 kraken/kraken_level_4.csv',sep='\t',row.names="X")
df_kraken3 <- read.delim('data/2020-11-15 kraken/kraken_level_3.csv',sep='\t',row.names="X")
df_kraken2 <- read.delim('data/2020-11-15 kraken/kraken_level_2.csv',sep='\t',row.names="X")

library(readxl)
df_level6_mira <- read.csv('data/2021-03-14 Mira/level-6.csv')


df_picrust3a <- read.csv('data/2021-04-13 Andreas PiCrust/PiCrust_Categorize_by_Function_KEGG_L3a.csv')
df_picrust2 <- read.csv('data/2021-04-13 Andreas PiCrust/PiCrust_categorize_by_Function_KEGG_L2.csv')
