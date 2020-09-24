#' Calculate reference set mean expression
#'
#' In this case, we have a matrix of expression values, and a character vector of samples in the reference set.
#' What we want to do is calculate the rowMeans from just the samples in the reference set
#'
#' @param expr.mat A matrix of expression values
#' @param reference_set The reference set
#'
#' @return A transformed matrix
#' @export
#'
calculate_reference_set_mean_expression <- function(expr.mat, reference_set){
  return(apply(expr.mat[,reference_set], 1, mean))
}


#' Calculate reference set coefficient of variation
#'
#' In this case, we have a matrix of expression values, and a character vector of samples in the reference set
#' What we want to do is calculate the coefficients of variation from just the samples in the reference set
#'
#' @param expr.mat A matrix of expression values
#' @param reference_set The reference set
#'
#' @return A transformed matrix
#' @export
#'
calculate_reference_set_coefficient_of_variation <- function(expr.mat, reference_set){
  row_mean <- apply(expr.mat[,reference_set], 1, mean)
  row_sd <- apply(expr.mat[,reference_set], 1, stats::sd)
  return(row_sd / row_mean * 100)
}

#' Transform expression matrix to long data frame
#'
#' This function takes a matrix (columns are samples, rows are genes),
#' and converts it into a 'long' data frame suitable for plotting
#' It also handles sorting the data by chromosomal coordinate.
#'
#' @param expression.mat A matrix of expression values.
#' @param genes.bed A bed-like file containing gene coordinates.
#'
#' @return expression.long.df
#' @export
expression_matrix_to_long_df <- function(expression.mat, genes.bed){
  # Convert back to data frame and set rownames -
  # This assumes the rownames of the matrix are already set as the gene names
  tmp.df <- as.data.frame(expression.mat)
  tmp.df$gene <- rownames(expression.mat)
  # Join the expression data with the sorted BED gene list -
  #  this has the effect of sorting by gene coordinate
  expression.df <- dplyr::inner_join(genes.bed, tmp.df, by = 'gene')
  expression.df$sequence <- seq(1, nrow(expression.df))
  expression.long.df <-
    expression.df %>%
    tidyr::gather(.data$library, .data$cnv_est, -.data$start,
                  -.data$stop, -.data$chrom, -.data$gene, -.data$sequence)
  return(expression.long.df)
}

#' Transform long data frame to sorted matrix
#'
#' @param long.df A tidy-formatted data frame of CNA data.
#'
#' @return ret.mat The re-sorted matrix.
#' @export
#'
long_df_to_ranked_matrix <- function(long.df){
  # Convert back to matrix, with rownames as the gene, and columns as the library name
  # Spread to go back to wide, and remove unnecessary columns
  tmp.wide.df <-
    long.df %>%
    dplyr::select(-.data$chrom, -.data$start, -.data$stop) %>%
    tidyr::spread(key = .data$library, value = .data$cnv_est) %>%
    dplyr::arrange(sequence)
  # Extract the rownames, and convert back to matrix
  ranked_genes <- tmp.wide.df$gene
  # Remove unnecessary columns
  ret.mat <- as.matrix(tmp.wide.df[,c(-1,-2)])
  # Set rownames
  rownames(ret.mat) <- ranked_genes
  return(ret.mat)
}

#' Transform data by taking the log2(count+1)
#'
#' @param x A matrix containing un-transformed expression (in some count format)
#'
#' @return A matrix with transformed values
#' @export
#'
#' @examples
#' \dontrun{log2_plus1_transform(matrix(c(1, 2, 4, 8)))}
log2_plus1_transform <- function(x){
  log2_plus_one.mat <- log2(x + 1)
  return(log2_plus_one.mat)
}

#' Rescale all genes by reference expression levels
#'
#' @param expr.mat a matrix of expression values
#' @param reference_means a named vector (With values for all genes, which is subset to the same genes retained in the matrix)
#'
#' @return expr.mat
#' @export
#'
rescale_genes_by_reference <- function(expr.mat, reference_means){
  return(expr.mat - reference_means[rownames(expr.mat)])
}

#' Run the CNA workflow
#'
#' @param expr.mat A matrix of gene-expression values
#' @param genes.bed A bed-like file of gene coordinates
#' @param window_size Window size to use for rolling mean
#' @param mean_tpm_threshold TPM threshold to retain genes
#' @param cv_threshold CV threshold to retain genes
#' @param reference_labels Reference labels
#' @param scale_samples_by_zero Scale samples by zero?
#'
#' @return A re-scaled matrix
#' @export
#'
run_cna_workflow <- function(expr.mat, genes.bed, window_size = NA,
                             mean_tpm_threshold = 0, cv_threshold = 0,
                             reference_labels,
                             scale_samples_by_zero = TRUE) {

  # Transform data by taking log2(TPM+1)
  expr.mat <- log2_plus1_transform(expr.mat)

  # Sort by chromosomal coordinate
  #  (This is done by converting back and forth to a long dataframe)
  expr.mat <-
    expr.mat %>%
    expression_matrix_to_long_df(genes.bed) %>%
    long_df_to_ranked_matrix()

  # Generate a reference set for mean gene expression and CV
  #  using a supplied reference set (generally, of copy-neutral samples)
  reference_means <- calculate_reference_set_mean_expression(expr.mat, reference_labels)
  reference_cvs <- calculate_reference_set_coefficient_of_variation(expr.mat, reference_labels)

  # Threshold by mean expression and coefficient of variation
  #  - retain genes only with mean expression > mean_tpm_threshold
  #  - retain genes only with CV < cv_threshold
  expr.mat <- subset_genes_by_thresholds(expr.mat, mean_tpm_threshold, cv_threshold,
                                         reference_means, reference_cvs)

  # Centre gene expression estimates by subtracting the mean
  #  value for each gene from the reference set
  expr.mat <- rescale_genes_by_reference(expr.mat, reference_means)

  # Calculate a rolling mean of expression values
  if (!is.na(window_size)) {
    expr.mat <- scale_by_rolling_mean(expr.mat, genes.bed, window_size)
  }

  # Scale samples around 0
  if (scale_samples_by_zero) {
    expr.mat <- scale(expr.mat)
  }

  return(expr.mat)
}

#' Scale expression estimates by rolling mean
#'
#' @param expr.mat a matrix of expression values
#' @param genes.bed A bed-like file containing gene coordinates
#' @param window_size Window size to use for rolling mean
#'
#' @return A transformed matrix
#' @export
#'
scale_by_rolling_mean <- function(expr.mat, genes.bed, window_size){
  # Convert to a long data frame - still ordered by gene
  tmp.long.df <-
    expr.mat %>%
    expression_matrix_to_long_df(genes.bed)
  # Split by chromosome and sample and apply rollmean
  tmp.rollmean_by_chrom.long.df <-
  tmp.long.df %>%
    dplyr::group_by(.data$chrom, .data$library) %>%
    dplyr::mutate(cnv_est = zoo::rollmean(x = .data$cnv_est, k = window_size,
                                   align = "center", fill = c(0,0,0)))
  # Ungroup and convert back to matrix
  ret.mat <-
    tmp.rollmean_by_chrom.long.df %>%
    dplyr::ungroup() %>%
    long_df_to_ranked_matrix()
  return(ret.mat)
}

#' Subset expression matrix by thresholds
#'
#' @param expr.mat A matrix of expression values
#' @param mean_tpm_threshold The mean TPM threshold. Genes with values above this threshold will be retained
#' @param cv_threshold The coefficient of variation threshold. Genes with values below this threshold will be retained
#' @param reference_means The reference means, generated by the calculate_reference_set_mean_expression function
#' @param reference_cvs The reference CVs, generated by the calculate_reference_set_coefficient_of_variation function
#'
#' @return expr.mat
#' @export
#'
subset_genes_by_thresholds <- function(expr.mat, mean_tpm_threshold,
                                       cv_threshold, reference_means,
                                       reference_cvs){
  # Figure out which genes pass both thresholds
  pass_threshold <- reference_means >= mean_tpm_threshold & reference_cvs <= cv_threshold
  # Subset the matrix for rows (genes) that pass thresholds
  return(expr.mat[pass_threshold,])
}
