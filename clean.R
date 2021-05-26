# level 6
df6_counts <- df_level6[,1:358] #indices found by hand
df6_other <- df_level6[,-(1:358)]

df <- t(df6_counts)
df6_counts <- df

# level 5
df5_counts <- df_level5[,1:138] #indices found by hand
df5_other <- df_level5[,-(1:138)]

df <- t(df5_counts)
df5_counts <- df

# level 4
df4_counts <- t(df_level4[,1:82]) #indices found by hand
df4_other <- df_level4[,-(1:82)]

# level 3
df3_counts <- t(df_level3[,1:33]) #indices found by hand
df3_other <- df_level3[,-(1:33)]

# level 2
df2_counts <- t(df_level2[,1:20]) #indices found by hand
df2_other <- df_level2[,-(1:20)]

# Remove S003J-0057 from the kraken data, not present in qiime data (found by hand)
###df6_kraken_counts <- cbind(df_kraken6[-126,],df6_other)
df6_kraken_counts <- t(as.matrix(df_kraken6[-126,]))
df6_kraken_other <- df6_other
rownames(df6_kraken_other) <- colnames(df6_kraken_counts)

df5_kraken_counts <- t(as.matrix(df_kraken5[-126,]))
df5_kraken_other <- df5_other
rownames(df5_kraken_other) <- colnames(df5_kraken_counts)

df4_kraken_counts <- t(as.matrix(df_kraken4[-126,]))
df4_kraken_other <- df4_other
rownames(df4_kraken_other) <- colnames(df4_kraken_counts)

df3_kraken_counts <- t(as.matrix(df_kraken3[-126,]))
df3_kraken_other <- df3_other
rownames(df3_kraken_other) <- colnames(df3_kraken_counts)

df2_kraken_counts <- t(as.matrix(df_kraken2[-126,]))
df2_kraken_other <- df2_other
rownames(df2_kraken_other) <- colnames(df2_kraken_counts)

## Cleaning data from Mira, 2021-03-14

index <- df_level6_mira$index
df <- df_level6_mira[,-1]
df6mira_counts <- df[,1:358] #indices found by hand
df6mira_other <- df[,-(1:358)]

df <- t(df6mira_counts)
df6mira_counts <- df

##### Cleaning data from Andreas, 2021-04-13 ####
library(stringr)

# level 2
df <- as.matrix(df_picrust2[-43,-1])
dimnames(df)[[1]] <- df_picrust2$X.OTU.ID[-43]
dimnames(df)[[2]] <- str_replace(dimnames(df)[[2]],'J\\.','J-')
#dfother <- merge(t(df),df2_other,by='row.names',all.x=T) # bring it in from above
dfother <- data.frame(t(df))
dfother$Covid <- ifelse(str_detect(rownames(dfother),'^BQ'),
                        'positive','negative')
dfother$gender <-  NA
dfother <- dfother[c('gender','Covid')]
# rownames(dfother) <- dfother$Row.names
# dfother <- subset(dfother,select=-c(Row.names))
pc2_counts <- as.matrix(df)
pc2_other <- dfother

# level 3a
df <- as.matrix(df_picrust3a[,-c(1,135:137)])
dimnames(df)[[1]] <- df_picrust3a$X.OTU.ID
dimnames(df)[[2]] <- str_replace(dimnames(df)[[2]],'J\\.','J-')
dfother <- data.frame(t(df))
dfother$Covid <- ifelse(str_detect(rownames(dfother),'^BQ'),
                        'positive','negative')
dfother$gender <-  NA
dfother <- dfother[c('gender','Covid')]
pc3a_counts <- as.matrix(df)
pc3a_other <- dfother