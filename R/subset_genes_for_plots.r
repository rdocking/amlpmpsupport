#' Subset the TPM matrix to the genes of interest
#'
#' @param tpm.mat Matrix containing expression values
#' @param exp.design Experimental design data frame
#' @param gene_list Vector of gene names to extract
#'
#' @return hits_to_plot_long A data frame containing hits for plotting
#' @export
#'
subset_genes_for_plots <- function(tpm.mat, exp.design, gene_list){
  hits_to_plot <- tpm.mat[gene_list,]
  hits_to_plot <- as.data.frame(t(hits_to_plot))
  hits_to_plot$rna_seq_lib <- rownames(hits_to_plot)
  # Convert to long for plotting
  hits_to_plot_long <- tidyr::gather(hits_to_plot,
                              key = gene, value = TPM,
                              -rna_seq_lib)
  # Add the experimental design
  hits_to_plot_long <- dplyr::left_join(hits_to_plot_long,
                                 exp.design, by='rna_seq_lib')
  return(hits_to_plot_long)
}
