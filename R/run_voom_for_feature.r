#' Run voom for a given feature comparison
#'
#' @param yvals Feature in the experimental design data frame to use as data labels
#' @param exp.design Experimental design data frame
#' @param counts.mat Matrix containing count data
#'
#' @return df A data frame containing the differential expression scores for all genes
#' @export
#'
run_voom_for_feature <- function(yvals, exp.design, counts.mat){
  # Set up the design matrix
  design.mat <- model.matrix(~yvals, exp.design)
  # Calculate the normalizing factor for each library
  norm.factor <- edgeR::calcNormFactors(counts.mat)
  # Run voom
  voom.dat <- limma::voom(counts.mat,
                          design.mat,
                          plot = TRUE,
                          lib.size = colSums(counts.mat) * norm.factor)
  # Clear plots
  grDevices::dev.off()
  # Run lmFit and eBayes
  voom.fit <- limma::lmFit(voom.dat, design.mat)
  voom.fit <- limma::eBayes(voom.fit)
  # Extract the list of all the genes
  voom.hits.all <- limma::topTable(voom.fit,
                                   coef = NULL,
                                   number = Inf,
                                   adjust.method = "BH",
                                   p.value = Inf)
  # Add the gene name back as a column for later
  voom.hits.all$gene <- rownames(voom.hits.all)
  return(voom.hits.all)
}

# Add Entrez gene IDs to a data frame of voom hits
#
# @param voom_hits A data frame containing voom hits, with HUGO gene ids
#
# @return df A data frame containing the Entrez gene IDs in a column called 'entrez'
# @export
#
#add_entrez_ids_from_hugo_genes <- function(voom_hits){
#  voom_hits$entrez = as.numeric(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
#                                                      keys=voom_hits$gene,
#                                                      column="ENTREZID",
#                                                      keytype="SYMBOL",
#                                                      multiVals="first"))
#  return(voom_hits)
#}


#' Munge differential expression hits in voom format to IPA-friendly format
#'
#' @param voom_hits A data frame containing voom hits
#'
#' @return df A data frame containing the voom hits in IPA-friendly format
#' @export
#'
munge_voom_to_ipa <- function(voom_hits){
  ipa_hits <- transmute(voom_hits,
                        ID = gene,
                        FOLD = limma_logFC_to_signed_foldchange(logFC),
                        P_VALUE = P.Value,
                        Q_VALUE = adj.P.Val)
  return(ipa_hits)
}

#' Convert limma's logFC value to signed fold-change
#'
#' Note that this needs to convert negative logFC numbers to absolute values to get correct fold-change estimates.
#' limma reports `logFC` as the 'estimate of the log2-fold-change'. For fold-changes in the \emph{negative} direction
#' we need to use the absolute value of the fold-change before applying the sign
#'
#' @param logFC A logFC value reported by limma
#'
#' @return foldchange A signed fold-change value
#' @export
#'
limma_logFC_to_signed_foldchange <- function(logFC){
  return (ifelse(abs(logFC >= 0), 2^logFC, (2^(abs(logFC)) * -1)))
}
