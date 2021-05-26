# Qiime2 differential analysis adjusted for agegroup and Gender

cts <- df6mira_counts
oth <- df6mira_other
# cts <- cts[,!is.na(oth$Age)]
# oth <- oth[!is.na(oth$Age),]
res6 <- analyze_Covid(cts,oth,pthresh=.01,fname="MIM_level6_DESeq_p_ordered.csv")


#res5 <- analyze_adjusted_Covid(df1,df2,pthresh=.01,fname="level6_adjDESeq_p_ordered.csv")


# Qiime2 differential analysis on Covid
res6 <- analyze_Covid(df6_counts,df6_other,pthresh=.01,fname="level6_DESeq_p_ordered.csv")
res5 <- analyze_Covid(df5_counts,df5_other,pthresh=.01,fname="level5_DESeq_p_ordered.csv")
res4 <- analyze_Covid(df4_counts,df4_other,pthresh=.01,fname="level4_DESeq_p_ordered.csv")
res3 <- analyze_Covid(df3_counts,df3_other,pthresh=.01,fname="level3_DESeq_p_ordered.csv")
res2 <- analyze_Covid(df2_counts,df2_other,pthresh=.01,fname="level2_DESeq_p_ordered.csv")
write.csv(res6[res6$baseMean>1000,],"level6_DESeq_p_ordered_basemean1000.csv")
write.csv(res5[res5$baseMean>1000,],"level5_DESeq_p_ordered_basemean1000.csv")
write.csv(res4[res4$baseMean>1000,],"level4_DESeq_p_ordered_basemean1000.csv")
write.csv(res3[res3$baseMean>1000,],"level3_DESeq_p_ordered_basemean1000.csv")
write.csv(res2[res2$baseMean>1000,],"level2_DESeq_p_ordered_basemean1000.csv")
saveplotMA(res6,"level6_DESeq_plot.pdf")
saveplotMA(res5,"level5_DESeq_plot.pdf")
saveplotMA(res4,"level4_DESeq_plot.pdf")
saveplotMA(res3,"level3_DESeq_plot.pdf")
saveplotMA(res2,"level2_DESeq_plot.pdf")

# KRAKEN differential analysis covid
res6k <- analyze_Covid(df6_kraken_counts,df6_kraken_other,pthresh=.01,fname="KR-level6_DESeq_p_ordered.csv")
write.csv(res6k[res6k$baseMean>1000,],"KR-level6_DESeq_p_ordered_basemean1000.csv")
saveplotMA(res6k,"KR-level6_DESeq_plot.pdf")
res5k <- analyze_Covid(df5_kraken_counts,df5_kraken_other,pthresh=.01,fname="KR-level5_DESeq_p_ordered.csv")
write.csv(res5k[res5k$baseMean>1000,],"KR-level5_DESeq_p_ordered_basemean1000.csv")
saveplotMA(res5k,"KR-level5_DESeq_plot.pdf")
res4k <- analyze_Covid(df4_kraken_counts,df4_kraken_other,pthresh=.01,fname="KR-level4_DESeq_p_ordered.csv")
write.csv(res4k[res4k$baseMean>1000,],"KR-level4_DESeq_p_ordered_basemean1000.csv")
saveplotMA(res4k,"KR-level4_DESeq_plot.pdf")
res3k <- analyze_Covid(df3_kraken_counts,df3_kraken_other,pthresh=.01,fname="KR-level3_DESeq_p_ordered.csv")
write.csv(res3k[res3k$baseMean>1000,],"KR-level3_DESeq_p_ordered_basemean1000.csv")
saveplotMA(res3k,"KR-level3_DESeq_plot.pdf")
res2k <- analyze_Covid(df2_kraken_counts,df2_kraken_other,pthresh=.01,fname="KR-level2_DESeq_p_ordered.csv")
write.csv(res2k[res2k$baseMean>1000,],"KR-level2_DESeq_p_ordered_basemean1000.csv")
saveplotMA(res2k,"KR-level2_DESeq_plot.pdf")


# Kraken Viral Load analysis including control group -- 
res <- analyze_viral_load(df6_kraken_counts,df6_kraken_other,pthresh=.01,fname="KR_level6_DESeq_p_ordered_VL.csv")
write.csv(res[res$baseMean>1000,],"KR_level6_DESeq_p_ordered_basemean1000_VL.csv")
saveplotMA(res,"KR-level6_DESeq_plot_VL.pdf")

res <- analyze_viral_load(df5_kraken_counts,df5_kraken_other,pthresh=.01,fname="KR_level5_DESeq_p_ordered_VL.csv")
write.csv(res[res$baseMean>1000,],"KR_level5_DESeq_p_ordered_basemean1000_VL.csv")
saveplotMA(res,"KR-level5_DESeq_plot_VL.pdf")

res <- analyze_viral_load(df4_kraken_counts,df4_kraken_other,pthresh=.01,fname="KR_level4_DESeq_p_ordered_VL.csv")
write.csv(res[res$baseMean>1000,],"KR_level4_DESeq_p_ordered_basemean1000_VL.csv")
saveplotMA(res,"KR-level4_DESeq_plot_VL.pdf")

res <- analyze_viral_load(df3_kraken_counts,df3_kraken_other,pthresh=.01,fname="KR_level3_DESeq_p_ordered_VL.csv")
write.csv(res[res$baseMean>1000,],"KR_level3_DESeq_p_ordered_basemean1000_VL.csv")
saveplotMA(res,"KR-level3_DESeq_plot_VL.pdf")

res <- analyze_viral_load(df2_kraken_counts,df2_kraken_other,pthresh=.01,fname="KR_level2_DESeq_p_ordered_VL.csv")
write.csv(res[res$baseMean>1000,],"KR_level2_DESeq_p_ordered_basemean1000_VL.csv")
saveplotMA(res,"KR-level2_DESeq_plot_VL.pdf")

## KRAKEN Not including control group
res <- analyze_viral_load(df6_kraken_counts,df6_kraken_other,pthresh=.01,fname="KR_level6_DESeq_p_ordered_VL.csv", useControl = FALSE,viralType = "High")
write.csv(res[res$baseMean>1000,],"KR_level6_DESeq_p_ordered_basemean1000_VL.csv")
saveplotMA(res,"KR-level6_DESeq_plot_VL.pdf")

res <- analyze_viral_load(df5_kraken_counts,df5_kraken_other,pthresh=.01,fname="KR_level5_DESeq_p_ordered_VL.csv", useControl = FALSE,viralType = "High")
write.csv(res[res$baseMean>1000,],"KR_level5_DESeq_p_ordered_basemean1000_VL.csv")
saveplotMA(res,"KR-level5_DESeq_plot_VL.pdf")

res <- analyze_viral_load(df4_kraken_counts,df4_kraken_other,pthresh=.01,fname="KR_level4_DESeq_p_ordered_VL.csv", useControl = FALSE)
write.csv(res[res$baseMean>1000,],"KR_level4_DESeq_p_ordered_basemean1000_VL.csv")
saveplotMA(res,"KR-level4_DESeq_plot_VL.pdf")

res <- analyze_viral_load(df3_kraken_counts,df3_kraken_other,pthresh=.01,fname="KR_level3_DESeq_p_ordered_VL.csv", useControl = FALSE)
write.csv(res[res$baseMean>1000,],"KR_level3_DESeq_p_ordered_basemean1000_VL.csv")
saveplotMA(res,"KR-level3_DESeq_plot_VL.pdf")

res <- analyze_viral_load(df2_kraken_counts,df2_kraken_other,pthresh=.01,fname="KR_level2_DESeq_p_ordered_VL.csv", useControl = FALSE)
write.csv(res[res$baseMean>1000,],"KR_level2_DESeq_p_ordered_basemean1000_VL.csv")
saveplotMA(res,"KR-level2_DESeq_plot_VL.pdf")



# Viral Load analysis including control group -- Qiime2
## Not including control group
res <- analyze_viral_load(df6_counts,df6_other,pthresh=.01,fname="q2_level6_DESeq_p_ordered_VL.csv", useControl = FALSE,viralType = "High")
write.csv(res[res$baseMean>1000,],"q2_level6_DESeq_p_ordered_basemean1000_VL.csv")
saveplotMA(res,"q2-level6_DESeq_plot_VL.pdf")

res <- analyze_viral_load(df5_counts,df5_other,pthresh=.01,fname="q2_level5_DESeq_p_ordered_VL.csv", useControl = FALSE,viralType = "High")
write.csv(res[res$baseMean>1000,],"q2_level5_DESeq_p_ordered_basemean1000_VL.csv")
saveplotMA(res,"q2-level5_DESeq_plot_VL.pdf")

res <- analyze_viral_load(df4_counts,df4_other,pthresh=.01,fname="q2_level4_DESeq_p_ordered_VL.csv", useControl = FALSE)
write.csv(res[res$baseMean>1000,],"q2_level4_DESeq_p_ordered_basemean1000_VL.csv")
saveplotMA(res,"q2-level4_DESeq_plot_VL.pdf")

res <- analyze_viral_load(df3_counts,df3_other,pthresh=.01,fname="q2_level3_DESeq_p_ordered_VL.csv", useControl = FALSE)
write.csv(res[res$baseMean>1000,],"q2_level3_DESeq_p_ordered_basemean1000_VL.csv")
saveplotMA(res,"q2-level3_DESeq_plot_VL.pdf")

res <- analyze_viral_load(df2_counts,df2_other,pthresh=.01,fname="q2_level2_DESeq_p_ordered_VL.csv", useControl = FALSE)
write.csv(res[res$baseMean>1000,],"q2_level2_DESeq_p_ordered_basemean1000_VL.csv")
saveplotMA(res,"q2-level2_DESeq_plot_VL.pdf")


# KRAKEN differential analysis Fibre
res6k <- analyze_Fibre(df6_kraken_counts,df6_kraken_other,pthresh=.01,fname="KR-level6_DESeq_p_ordered_Fibre.csv")
write.csv(res6k[res6k$baseMean>1000,],"KR-level6_DESeq_p_ordered_basemean1000_Fibre.csv")
saveplotMA(res6k,"KR-level6_DESeq_plot_Fibre.pdf")


# yesCovid <- df$Covid=="positive"
# noCovid <- !yesCovid
# 
# pthresh <- .01
# 
# pv <- numeric(0)
# sig_i <- numeric(0)
# 
# for (i in 1:ncol(df)) {
#   if (is.numeric(df[[i]])) {
#     w.test <- wilcox.test(df[[i]][yesCovid],df[[i]][noCovid])  
#     if (w.test$p.value < pthresh){
#       pv <- c(pv,w.test$p.value)
#       sig_i <- c(sig_i,i)
#     }
#   }
#   
# }

# Picrust
res <- analyze_Covid(pc2_counts,pc2_other)
write.csv(res,"Picrust_2_DESeq.csv")
write.csv(res[res$baseMean>1000,],"Picrust_2_DESeq_p_ordered_basemean1000.csv")
saveplotMA(res,"Picrust_2_DESeq_plot.pdf")

res <- analyze_Covid(pc3a_counts,pc3a_other)
write.csv(res,"Picrust_3a_DESeq.csv")
write.csv(res[res$baseMean>1000,],"Picrust_3a_DESeq_p_ordered_basemean1000.csv")
saveplotMA(res,"Picrust_3a_DESeq_plot.pdf")
