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
    tmp.df <- ifxparsers::rna_seqc_parser(data_row$rna_seqc_path)
    tmp.df$library_name <- lib
    # Bind the count data to the growing df
    rna_seqc.df <- bind_rows(rna_seqc.df, tmp.df)
    # Arrange by library name
    rna_seqc.df <- arrange(rna_seqc.df, Sample)
  }
  # Return the data frame
  return(rna_seqc.df)
}
