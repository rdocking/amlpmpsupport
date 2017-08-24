#' Add nested library information to a patient-level data frame
#'
#' @param patients.df Data frame containing patient-level data
#' @param libraries.df Data frame containing library-level data
#' @param key Key to join the two data-frames by (default is 'tfl_id')
#'
#' @return patients.w_libs.df A data frame patients with nested library information
#' @export
#'
nest_libraries_by_patient <- function(patients.df, libraries.df) {

  # Map functions - find all the library names and concatenate them
  select_fun <- function(df) dplyr::select(df, library_name)
  cat_fun <- function(df) stringr::str_c(df, sep = ", ")

  # This code aggregates all the library data by tfl_id, then simplifies down
  #  to a list of libraries associated with each tfl_id and joins it back
  #  to the patient-level data frame
  libraries.df %>%
    dplyr::group_by(tfl_id) %>%
    tidyr::nest() %>%
    dplyr::mutate(lib_list = purrr::map(data, select_fun),
                  libraries = purrr::map_chr(lib_list, cat_fun)) %>%
    # dplyr::select(-data, -lib_list) %>%
    dplyr::right_join(patients.df, by = 'tfl_id')
}
