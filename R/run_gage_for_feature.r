#' Run GAGE gene-set enrichment analysis for a given feature
#'
#' @param voom.hits A data frame from the voom analysis, _unthresholded_
#' @param gsets A list of gene sets, constructed by the geneIds function above
#' @param adj_p_threshold Threshold cutoff to use for voom adjusted p-value
#' @param same.dir same.dir argument for the `gage` function
#'
#' @return pathways.merged A data frame containing up- and down-regulated pathways
#' @export
#'
run_gage_for_feature <- function(voom.hits, gsets,
                                 adj_p_threshold = 0.05, same.dir = TRUE){
  # Subset to adjusted p-value threshold <= 0.01
  threshold.hits <- dplyr::filter(voom.hits, adj.P.Val <= adj_p_threshold)
  # Set up a named vector of fold changes
  foldchanges = threshold.hits$logFC
  names(foldchanges) = threshold.hits$entrez
  # Run gage to find the enriched pathways:
  pathways = gage::gage(foldchanges,
                  gsets=gsets,
                  same.dir=same.dir)
  # Inspect up- and down-regulated pathways as a data frame
  pathways.up <- as.data.frame(pathways$greater)
  pathways.up$pathway <- rownames(pathways.up)
  pathways.down <- as.data.frame(pathways$lesser)
  pathways.down$pathway <- rownames(pathways.down)
  # Concatenate and return the data frames
  pathways.merged <- dplyr::bind_rows(pathways.up, pathways.down)
  return(pathways.merged)
}
