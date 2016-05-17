#' Train a random forest model with caret
#'
#' @param gene_set A vector containing gene names to subset
#' @param yvals Feature in the experimental design data frame to use as data labels
#' @param exp.design Experimental design data frame
#' @param counts.mat A matrix of gene count data
#'
#' @return caret.rf A an object containing a caret Randomforest
#' @export
#'
train_caret_rf <- function(gene_set, yvals, exp.design, counts.mat){
  # Find rows to keep
  rows_to_keep <- rownames(counts.mat) %in% gene_set
  # Subset the count matrix to the genes of interest
  counts.subset.mat <- counts.mat[rows_to_keep,]
  # Transpose the matrix
  counts.subset.transpose.mat <- t(counts.subset.mat)
  # Use caret to train a RF classifier
  caret.rf <- caret::train(x = counts.subset.transpose.mat,
                           y = yvals,
                           method="rf",
                           trControl=trainControl(method="cv",number=5),
                           prox=TRUE,
                           allowParallel=TRUE,
                           ntree=5000)
  return(caret.rf)
}
