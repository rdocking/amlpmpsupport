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
                          lib.size = colSums(sailfish.genes.counts.mat) * norm.factor)
  # Clear plots
  dev.off()
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
  # Add Entrez IDs for all genes
  voom.hits.all$entrez = as.numeric(AnnotationDbi::mapIds(org.Hs.eg.db,
                                                          keys=voom.hits.all$gene,
                                                          column="ENTREZID",
                                                          keytype="SYMBOL",
                                                          multiVals="first"))
  return(voom.hits.all)
}
