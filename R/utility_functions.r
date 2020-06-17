#' Add nested library information to a patient-level data frame
#'
#' @param patients.df Data frame containing patient-level data.
#' @param libraries.df Data frame containing library-level data.
#' @param key Key to join the two data-frames by (default is 'tfl_id').
#'
#' @return patients.w_libs.df A data frame patients with nested library information
#' @export
#'
nest_libraries_by_patient <- function(patients.df, libraries.df, key = 'tfl_id') {

  # Map functions - find all the library names and concatenate them
  select_fun <- function(df) dplyr::select(df, .data$library_name)
  cat_fun <- function(df) stringr::str_c(df, sep = ", ")

  # This code aggregates all the library data by tfl_id, then simplifies down
  #  to a list of libraries associated with each tfl_id and joins it back
  #  to the patient-level data frame
  libraries.df %>%
    dplyr::group_by(.data$tfl_id) %>%
    tidyr::nest_legacy() %>%
    dplyr::mutate(lib_list = purrr::map(.data$data, select_fun),
                  libraries = purrr::map_chr(.data$lib_list, cat_fun)) %>%
    dplyr::right_join(patients.df, by = key)
}

#' Set a function for pretty-printing kables
#'
#' @param df Data frame containing data to format as a table
#' @param ... Additional options to be passed to kable
#'
#' @return A kable table
#' @export
#'
#' @examples pretty_kable(data.frame(x = c(1000, 999), y = c(0.345:0.351)))
pretty_kable <- function(df, ...){
  # This has options I like - comma separator for thousands, small number of sig digits
  knitr::kable(df, digits = 2,
               format.args = list(big.mark = ','),
               ...)
}

#' Set a function for pretty-printing numbers in the text
#'
#' @param x Number to print
#' @param ... Additional options to be passed to prettyNum
#'
#' @return A formatted number
#' @export
#'
#' @examples pretty_num(5.55654)
pretty_num <- function(x, ...){
  # This has options I like - small number of sig digits
  prettyNum(x, digits = 3, ...)
}

#' # Report basic information about a matrix
#'
#' @param mat A matrix
#'
#' @export
#'
#' @examples matrix_glimpse(matrix(runif(1000), nrow = 200, ncol = 500))
matrix_glimpse <- function(mat) {
  print("Dimensions: ")
  print(dim(mat))
  print("Range: ")
  print(range(mat))
  print("Head: ")
  print(mat[1:3,1:3])
}

#' Quick function for mapping to select out list items
#'
#' This is mostly useful in the context of purrr::map pipelines.
#'
#' @param x A vector.
#' @param position The position to select out.
#'
#' @return The vector at the given index.
#' @export
#'
#' @examples select_by_position(c(1,2,3), 2)
select_by_position <- function(x, position) x[position]

#' Convert a small number to a text representation
#'
#' @param x An integer to convert to text.
#' @param capitalize Boolean - whether or not to capitalize text
#'
#' @return A character string.
#' @export
#'
#' @examples
#' num_to_text(5, capitalize = TRUE)
#' num_to_text(23, capitalize = FALSE)
num_to_text <- function(x, capitalize = FALSE){

  # Input checks
  assertthat::is.number(x)
  assertthat::assert_that(x <= 30)

  num_table <- tibble::tribble(
    ~x, ~small,
    1, 'one',
    2, 'two',
    3, 'three',
    4, 'four',
    5, 'five',
    6, 'six',
    7, 'seven',
    8, 'eight',
    9, 'nine',
    10, 'ten',
    11, 'eleven',
    12, 'twelve',
    13, 'thirteen',
    14, 'fourteen',
    15, 'fifteen',
    16, 'sixteen',
    17, 'seventeen',
    18, 'eighteen',
    19, 'nineteen',
    20, 'twenty',
    21, 'twenty-one',
    22, 'twenty-two',
    23, 'twenty-three',
    24, 'twenty-four',
    25, 'twenty-five',
    26, 'twenty-six',
    27, 'twenty-seven',
    28, 'twenty-eight',
    29, 'twenty-nine',
    30, 'thirty') %>%
    dplyr::mutate(big = stringr::str_to_sentence(.data$small))

  if (capitalize) {
    ret <- num_table[x,]$big
  } else {
    ret <- num_table[x,]$small
  }
  return(ret)
}

#' Pull up-regulated genes
#'
#' @param df Data frame containing DE results
#' @param adjpval_threshold Adjusted p-value threshold
#' @param log2fc_threshold Log2-fold-change threshold
#'
#' @return A vector of up-regulated genes
#' @export
upregulated_genes <- function(df, adjpval_threshold, log2fc_threshold){
  dplyr::filter(df,
         padj <= adjpval_threshold,
         log2FoldChange >= log2fc_threshold) %>%
    dplyr::pull(gene_name)
}

#' Pull down-regulated genes
#'
#' @param df Data frame containing DE results.
#' @param adjpval_threshold Adjusted p-value threshold.
#' @param log2fc_threshold Log2-fold-change threshold.
#'
#' @return A vector of down-regulated genes.
#' @export
downregulated_genes <- function(df, adjpval_threshold, log2fc_threshold){
  dplyr::filter(df,
         .data$padj <= adjpval_threshold,
         .data$log2FoldChange <= -1 * log2fc_threshold) %>%
    dplyr::pull(.data$gene_name)
}

#' Group / Tally / Table
#'
#' @param x A data frame.
#' @param ... Unquoted variables to group by.
#'
#' @return A knitr::kable().
#' @export
#'
#' @examples group_tally_table(mtcars, cyl, carb)
group_tally_table <- function(x, ...) {
  x %>%
    dplyr::group_by(...) %>%
    dplyr::tally(wt = NULL) %>%
    knitr::kable()
}

#' Save ggplot to PDF and PNG
#'
#' @param filename_root Root filename (without extension).
#' @param plot ggplot object.
#' @param width plot width in inches.
#' @param height plot height in inches.
#' @param ... parameters passed on to ggsave.
#'
#' @export
ggsave_pdf_and_png <- function(filename_root, plot, width = 4, height = 2.472, ...){

  plot_file_pdf <- fs::path_ext_set(filename_root, 'pdf')
  ggplot2::ggsave(plot_file_pdf, plot = plot,
         device = "pdf", width = width, height = height,
         ...)

  plot_file_png <- fs::path_ext_set(filename_root, 'png')
  ggplot2::ggsave(plot_file_png, plot = plot,
         device = "png", width = width, height = height,
         ...)
}

#' A function to sort TFL IDs into ascending order
#'
#' @param df Data frame containing TFL IDs in 'patient_external_id' column.
#'
#' @return A re-sorted data frame.
#' @export
tfl_split_sort <- function(df){

  # Split identifier
  tmp.df <- df %>%
    tidyr::separate(patient_external_id, into = c('prefix', 'suffix'),
                    sep = '-', remove = FALSE, fill = "right",
                    convert = TRUE)

  tmp.df$prefix <- as.integer(tmp.df$prefix)
  tmp.df$suffix <- as.integer(tmp.df$suffix)

  dplyr::arrange(tmp.df, prefix, suffix) %>%
    dplyr::select(-prefix, -suffix)
}
