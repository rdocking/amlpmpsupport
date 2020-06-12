#' Construct a data frame of RNA-SeQC results
#'
#' Given a data frame containing library names and RNA-SeQC paths, Construct a data frame of RNA-SeQC results
#' NOTE that this function makes strong assumptions about the shape of the incoming data frame
#'
#' @param exp.design Experimental design data frame
#'
#' @return df A data frame containing the RNA-SeQC results
#' @export
#'
construct_rna_qc_dataframe <- function(exp.design){
  # Set up an empty data frame
  rna_seqc.df <- data.frame()
  for (lib in exp.design$rna_seq_lib) {
    # Select the row that we want from the experimental design
    data_row <- filter(exp.design, rna_seq_lib == lib)
    # Read the QC data into a temporary df
    tmp.df <- rna_seqc_parser(data_row$rna_seqc_path)
    tmp.df$library_name <- lib
    # Bind the count data to the growing df
    rna_seqc.df <- bind_rows(rna_seqc.df, tmp.df)
    # Arrange by library name
    rna_seqc.df <- arrange(rna_seqc.df, Sample)
  }
  # Return the data frame
  return(rna_seqc.df)
}

#' Construct a data frame of Sailfish gene-level results
#'
#' Given a data frame containing library names and Sailfish paths, Construct a data frame of Sailfish results
#' NOTE that this function makes strong assumptions about the shape of the incoming data frame,
#' and currently only fetches \emph{gene-level} results
#'
#' @param exp.design Experimental design data frame
#'
#' @return df A data frame containing the Sailfish results
#' @export
#'
construct_sailfish_gene_dataframe <- function(exp.design){
  sailfish.genes.df <- data.frame()
  for (lib in exp.design$rna_seq_lib) {
    # Select the row that we want
    data_row <- filter(exp.design, rna_seq_lib == lib)
    # Read the events data into a temporary df
    tmp.df <- sailfish_gene_parser_post_0.7.0(data_row$sailfish_path)
    tmp.df$library_name <- lib
    # Bind the count data to the growing df
    sailfish.genes.df <- bind_rows(sailfish.genes.df, tmp.df)
  }
  return(sailfish.genes.df)
}

#' From a data frame of Sailfish gene-level results, construct an expression matrix
#'
#' Given a data frame containing Sailfish gene-level results, from the construct_sailfish_gene_dataframe
#' function, convert that dataframe to a matrix suitable for downstream analysis
#'
#' @param sailfish.df Data frame containing sailfish results
#' @param metric Expression metric to use (either 'NumReads' or 'TPM')
#'
#' @return df A data frame containing the Sailfish results
#' @export
#'
convert_sailfish_df_to_matrix <- function(sailfish.df, metric = c("NumReads", "TPM")){
  # Select out the relevant columns
  # What I want is just the gene name, the desired metric, and the library names
  # Note that we need to use the standard evaluation versions here because of the arguments
  metric_subset.df <- dplyr::select_(sailfish.df, "Name", metric, "library_name")
  # Spread the data into a wide format...
  metric_subset.wide.df <- tidyr::spread_(metric_subset.df,
                                          key="library_name",
                                          value=metric)
  # Set the rownames as the transcript ID
  rownames(metric_subset.wide.df) <- metric_subset.wide.df$Name
  # Then remove the transcript column and convert to a matrix
  metric_subset.wide.df <- dplyr::select(metric_subset.wide.df, -Name)
  metric_subset.mat <- as.matrix(metric_subset.wide.df)
  return(metric_subset.mat)
}

#' Convert a 'tidy' expression DF to a matrix
#'
#' Assumes the only other columns other than 'key' and 'val' are expression measurements
#'
#' @param df A tidy data frame
#' @param key Column name for the column containing gene names
#' @param val Column name for the column containing sample IDs
#'
#' @return x A matrix containing the expression results
#' @export
#'
spread_expression_df_to_matrix <- function(df, key, value){
  # Quote the input variables
  key <- enquo(key)
  value <- enquo(value)
  # Spread to wide and convert to DF
  tmp.df <- spread(df, key = !!key, value = !!value) %>% as.data.frame()
  # Set the first column as the rownames and return as a matrix
  rownames(tmp.df) <- tmp.df[,1]
  return(as.matrix(tmp.df[,-1]))
}

#' Convert an expression matrix to a tidy DF
#'
#' @param x Expression matrix
#' @param gene_col Column name for the column containing gene names
#' @param id_col Column name for the column containing sample IDs
#' @param expr_col Column name for the column containing expression values
#'
#' @return df A tidy data frame containing the expression results
#' @export
#'
gather_expression_matrix_to_tidy_df <- function(x, gene_col, id_col, expr_col){
  x %>%
    as.data.frame() %>%
    rownames_to_column(gene_col) %>%
    gather(key = !!id_col, value = !!expr_col, -gene_col)
}


#' Pull EnrichR enriched terms for a single pathway source
#'
#' @param enrichr_lst Result list returned from EnrichR
#' @param pathway_source Name of pathway source to query
#' @param cohort_name Cohort name to apply
#' @param padj_thresh Adjusted p-value threshold
#'
#' @return
#' @export
#'
#' @examples
pull_enrichr_enriched_terms <- function(enrichr_lst, pathway_source, cohort_name, padj_thresh){
    enrichr_lst[pathway_source][[1]] %>%
      dplyr::select(Term, Overlap, padj = Adjusted.P.value) %>%
      dplyr::mutate(cohort = cohort_name) %>%
      dplyr::filter(padj <= padj_thresh)
  }

#' Identify enriched terms across cohorts
#'
#' @param long_df Long dataframe containing enrichment results across cohorts
#'
#' @return
#' @export
#'
#' @examples
pull_enrichr_terms_to_plot <- function(long_df){

  ranked_terms.df <-
    long_df %>%
    dplyr::group_by(Term) %>%
    dplyr::summarise(count = n(),
                     mean = mean(padj)) %>%
    dplyr::arrange(-count, mean)

  terms_to_keep <-
    ranked_terms.df %>%
    dplyr::filter(count >= 2) %>%
    dplyr::pull(Term)

  return(list("ranked_terms" = ranked_terms.df,
              "terms_to_keep" = terms_to_keep))
}
