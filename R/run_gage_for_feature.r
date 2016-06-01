#' Run GAGE gene-set enrichment analysis for a given feature
#'
#' @param voom.hits A data frame from the voom analysis, _unthresholded_
#' @param gsets A list of gene sets, constructed by the geneIds function above
#' @param same.dir same.dir argument for the `gage` function
#'
#' @return pathways.merged A data frame containing up- and down-regulated pathways
#' @export
#'
run_gage_for_feature <- function(voom.hits, gsets, same.dir = TRUE){
  # Set up a named vector of fold changes
  foldchanges = voom.hits$logFC
  names(foldchanges) = voom.hits$entrez
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
