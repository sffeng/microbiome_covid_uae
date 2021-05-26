library(DESeq2)

saveplotMA <- function(res, fname="plotMA.pdf"){
  pdf(fname)
  plotMA(res)
  dev.off()
}

analyze_adjusted_Covid <- function(df_counts,
                          df_other,
                          pthresh=.05,
                          fname="deseq_p_ordered.csv"){
  cts <- as.matrix(df_counts)
  df_other$AgeGroups <- cut(df_other$Age,c(0, 25, 35, 51, Inf))
  
  coldata <- df_other[,c("Covid","gender","AgeGroups")]
  
  if (!all(rownames(coldata)==colnames(cts))){
    print("ERROR IN ROWNAMES AND COLNAMSEgoo")
  }
  dds <- DESeqDataSetFromMatrix(countData=cts,
                                colData = coldata,
                                design = ~ Covid + gender + AgeGroups)
  
  dds <- DESeq(dds,test="LRT",reduced=~ gender + AgeGroups)
  res <- results(dds,alpha=pthresh)
  resOrdered <- res[order(res$pvalue),]
  (nSig <- sum(res$padj < 0.01, na.rm=TRUE))
  write.csv(resOrdered,fname)
  resOrdered
}


analyze_Covid <- function(df_counts,
                          df_other,
                          pthresh=.05,
                          fname="deseq_p_ordered.csv"){
  cts <- as.matrix(df_counts)
  coldata <- df_other[,c("Covid","gender")]
  
  if (!all(rownames(coldata)==colnames(cts))){
    print("ERROR IN ROWNAMES AND COLNAMES")
  }
  dds <- DESeqDataSetFromMatrix(countData=cts,
                                colData = coldata,
                                design = ~ Covid)
  
  dds <- DESeq(dds)
  res <- results(dds,alpha=pthresh)
  resOrdered <- res[order(res$pvalue),]
  (nSig <- sum(res$padj < 0.01, na.rm=TRUE))
  write.csv(resOrdered,fname)
  resOrdered
}



analyze_Fibre <- function(df_counts,
                          df_other,
                          pthresh=.05,
                          fname="deseq_p_ordered.csv"){
  cts <- as.matrix(df_counts)
  coldata <- df_other[,c("Fibre","gender")]
  
  if (!all(rownames(coldata)==colnames(cts))){
    print("ERROR IN ROWNAMES AND COLNAMES")
  }
  dds <- DESeqDataSetFromMatrix(countData=cts,
                                colData = coldata,
                                design = ~ Fibre)
  
  dds <- DESeq(dds)
  res <- results(dds,alpha=pthresh)
  resOrdered <- res[order(res$pvalue),]
  (nSig <- sum(res$padj < 0.01, na.rm=TRUE))
  write.csv(resOrdered,fname)
  resOrdered
}



analyze_viral_load <- function(df_counts,
                          df_other,
                          pthresh=.05,
                          fname="deseq_viral_load_p_ordered.csv",
                          useControl=TRUE,
                          viralType="Normal"){
  
  if (!useControl){
    print("Removing Control Group")
    df_counts <- df_counts[,1:78]
    df_other <- df_other[1:78,]
  }
  
  if (viralType=="High"){
    I <- df_other$viral.loaad_stool %in% c(0,3)
    df_counts <- df_counts[,I]
    df_other <- df_other[I,]
  }
  
  cts <- as.matrix(df_counts)
  df_other$hasViral <- as.factor(df_other$viral.loaad_stool>0)
  

  
  coldata <- df_other[,c("hasViral","gender")]
  
  if (!all(rownames(coldata)==colnames(cts))){
    print("ERROR IN ROWNAMES AND COLNAMES")
  }
  dds <- DESeqDataSetFromMatrix(countData=cts,
                                colData = coldata,
                                design = ~ hasViral)
  
  dds <- DESeq(dds)
  res <- results(dds,alpha=pthresh)
  resOrdered <- res[order(res$pvalue),]
  (nSig <- sum(res$padj < 0.01, na.rm=TRUE))
  write.csv(resOrdered,fname)
  resOrdered
}



# custom function to transpose while preserving names
transpose_df <- function(df) {
  t_df <- data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  t_df <- t_df %>%
    tibble::rownames_to_column(.data = .) %>%
    tibble::as_tibble(.)
  return(t_df)
}
