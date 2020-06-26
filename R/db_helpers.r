#' Connect to AML Project ORM Database
#'
#' @param dbname Database name
#' @param host Database host
#' @param user Database username
#'
#' @return db_conn A database connection object
#' @export
connect_to_aml_project_orm <- function(dbname = 'aml_project_db',
                                       host = 'karsanlab-db01',
                                       user) {
  # Connect to the database
  db_conn <- dplyr::src_postgres(dbname = dbname,
                               host = host,
                               user = user)
  # Return the DB connection object
  return(db_conn)
}


#' Retrieve DB Libraries
#'
#' @param db_conn Database connection object
#'
#' @return db_libraries.df A data-frame containing library information from the database
#' @export
retrieve_db_libraries <- function(db_conn) {
  # Connect to the relevant tables
  patient.tbl <- dplyr::tbl(db_conn, "patient")
  library.tbl <- dplyr::tbl(db_conn, "library")
  specimen.tbl <- dplyr::tbl(db_conn, "specimen")
  specimen_subset.tbl <- dplyr::tbl(db_conn, "specimen_subset")
  platform_version.tbl <- dplyr::tbl(db_conn, "platform_version")
  platform.tbl <- dplyr::tbl(db_conn, "platform")
  # Start from patient, and query downwards
  db_libraries.query <- patient.tbl %>%
    dplyr::select(.data$id, .data$external_id, .data$is_patient) %>%
    dplyr::rename(patient_id = .data$id) %>%
    dplyr::rename(patient_external_id = .data$external_id) %>%
    # Join to specimen
    dplyr::left_join(specimen.tbl, by = "patient_id") %>%
    dplyr::rename(specimen_id = .data$id) %>%
    dplyr::rename(specimen_external_id = .data$external_id) %>%
    dplyr::select(-.data$patient_id, -.data$meta) %>%
    # Join to specimen_subset
    dplyr::left_join(specimen_subset.tbl, by = "specimen_id") %>%
    dplyr::rename(specimen_subset_id = .data$id) %>%
    dplyr::rename(specimen_subset_external_id = .data$external_id) %>%
    # Join to library
    dplyr::left_join(library.tbl, by = "specimen_subset_id") %>%
    dplyr::rename(library_id = .data$id) %>%
    dplyr::rename(library_name = .data$name) %>%
    # Join to platform and platform_version
    dplyr::left_join(platform_version.tbl,
                     by = c("platform_version_id" = "id")) %>%
    dplyr::rename(platform_version_name = .data$name) %>%
    dplyr::left_join(platform.tbl,
                     by = c("platform_id" = "id")) %>%
    dplyr::rename(platform_name = .data$name) %>%
    dplyr::select(-.data$fields, -.data$platform_id, -.data$platform_version_id,
                  -.data$specimen_subset_id, -.data$specimen_id,
                  -.data$library_id, -.data$sequencing_effort_id,
                  -.data$library_qc_info, -.data$type) %>%
    # Filter out rows with no library name
    dplyr::filter(!is.na(.data$library_name))
  db_libraries.df <- dplyr::collect(db_libraries.query)
  return(db_libraries.df)
}

#' Retrieve DB Paths
#'
#' @param db_conn Database connection object
#'
#' @return db_paths.df A data-frame containing path information from the database
#' @export
retrieve_db_paths <- function(db_conn) {
  # Connect to the relevant tables
  patient.tbl <- dplyr::tbl(db_conn, "patient")
  library.tbl <- dplyr::tbl(db_conn, "library")
  specimen.tbl <- dplyr::tbl(db_conn, "specimen")
  specimen_subset.tbl <- dplyr::tbl(db_conn, "specimen_subset")
  platform_version.tbl <- dplyr::tbl(db_conn, "platform_version")
  platform.tbl <- dplyr::tbl(db_conn, "platform")
  analysis.tbl <- dplyr::tbl(db_conn, "analysis")
  path.tbl <- dplyr::tbl(db_conn, "path")
  filetype.tbl <- dplyr::tbl(db_conn, "filetype")
  method.tbl <- dplyr::tbl(db_conn, "method")
  pipeline_version.tbl <- dplyr::tbl(db_conn, "pipeline_version")
  # Start from patient, and query downwards
  db_paths.query <-
    patient.tbl %>%
    dplyr::select(.data$id, .data$external_id, .data$is_patient) %>%
    dplyr::rename(patient_id = .data$id) %>%
    dplyr::rename(patient_external_id = .data$external_id) %>%
    # Join to specimen
    dplyr::left_join(specimen.tbl, by = "patient_id") %>%
    dplyr::rename(specimen_id = .data$id) %>%
    dplyr::rename(specimen_external_id = .data$external_id) %>%
    dplyr::select(-.data$patient_id, -.data$meta) %>%
    # Join to specimen_subset
    dplyr::left_join(specimen_subset.tbl, by = "specimen_id") %>%
    dplyr::rename(specimen_subset_id = .data$id) %>%
    dplyr::rename(specimen_subset_external_id = .data$external_id) %>%
    # Join to library
    dplyr::left_join(library.tbl, by = "specimen_subset_id") %>%
    dplyr::rename(library_id = .data$id) %>%
    dplyr::rename(library_name = .data$name) %>%
    # Join to platform and platform_version
    dplyr::left_join(platform_version.tbl,
                     by = c("platform_version_id" = "id")) %>%
    dplyr::rename(platform_version_name = .data$name) %>%
    dplyr::left_join(platform.tbl,
                     by = c("platform_id" = "id")) %>%
    dplyr::rename(platform_name = .data$name) %>%
    # Join to analysis
    dplyr::left_join(analysis.tbl,
                     by = c("library_id")) %>%
    dplyr::rename(analysis_id = .data$id) %>%
    # Join to path
    dplyr::left_join(path.tbl,
                     by = c("analysis_id")) %>%
    dplyr::rename(path_id = .data$id) %>%
    # Join to filetype, method, and pipeline version
    dplyr::left_join(filetype.tbl,
                     by = c("filetype_id" = "id")) %>%
    dplyr::rename(filetype_name = .data$name) %>%
    dplyr::left_join(method.tbl,
                     by = c("method_id" = "id")) %>%
    dplyr::rename(method_name = .data$name) %>%
    dplyr::left_join(pipeline_version.tbl,
                     by = c("pipeline_version_id" = "id")) %>%
    dplyr::rename(pipeline_version_name = .data$name) %>%
    # Remove rows that have no path information
    dplyr::filter(!is.na(.data$path)) %>%
    # Select out UUID fields and unused fields
    dplyr::select(-.data$fields, -.data$platform_id, -.data$platform_version_id,
                  -.data$specimen_subset_id, -.data$specimen_id,
                  -.data$library_id, -.data$sequencing_effort_id,
                  -.data$library_qc_info, -.data$type, -.data$description,
                  -.data$path_id, -.data$filetype_id, -.data$method_id,
                  -.data$pipeline_version_id, -.data$analysis_id)
  db_paths.df <- dplyr::collect(db_paths.query)
  return(db_paths.df)
}


#' From a subtype indication, return a simplified disease type
#'
#' @param subtype Curated disease subtype
#'
#' @return disease A simplified disease name
#' @export
subtype_to_disease <- function(subtype) {
  disease <- NA
  if (subtype == 'AML-MDS') {
    disease <- subtype
  } else if (startsWith(subtype, "AML")) {
    disease <- 'AML'
  } else {
    disease <- subtype
  }
  return(disease)
}


#' Retrieve DB Curated Results
#'
#' @param db_conn Database connection object
#'
#' @return db_curated_results.df A data-frame containing curated result information from the database
#' @export
retrieve_db_curated_results <- function(db_conn) {
  # Connect to the relevant tables
  patient.tbl <- dplyr::tbl(db_conn, "patient")
  specimen.tbl <- dplyr::tbl(db_conn, "specimen")
  curated_result.tbl <- dplyr::tbl(db_conn, "curated_result")
  # Starting from the curated results table in the database...
  curated_result.tbl %>%
    # Join to specimen
    dplyr::inner_join(specimen.tbl, by = c("specimen_id" = "id")) %>%
    dplyr::rename(specimen_external_id = .data$external_id) %>%
    # Join to patient
    dplyr::inner_join(patient.tbl, by = c("patient_id" = "id")) %>%
    dplyr::rename(patient_external_id = .data$external_id) %>%
    # Select matching fields
    dplyr::select(.data$patient_external_id, .data$specimen_external_id,
                  .data$subtype, .data$flt3_status,
                  .data$npm1_status, .data$kit_status, .data$cebpa_status,
                  .data$cytogenetic_risk, .data$karyotype, .data$fusions,
                  .data$tier_one_mutation_status, .data$kmt2a_status,
                  .data$dnmt3a_status) %>%
    # Sort by specimen ID
    dplyr::arrange(.data$specimen_external_id) ->
    curated_results.query
    # Collect results into a data frame
    curated_results.df <- dplyr::collect(curated_results.query)
    # Add a simplified disease type
    curated_results.df$disease <- purrr::map_chr(curated_results.df$subtype,
                                                 subtype_to_disease)
    return(curated_results.df)
}

#' Retrieve DB Specimens
#'
#' @param db_conn Database connection object
#'
#' @return db_specimens.df A data-frame containing specimen information from the database
#' @export
retrieve_db_specimens <- function(db_conn) {
  # Connect to the relevant tables
  patient.tbl <- dplyr::tbl(db_conn, "patient")
  specimen.tbl <- dplyr::tbl(db_conn, "specimen")
  specimen_subset.tbl <- dplyr::tbl(db_conn, "specimen_subset")
  # Start from patient, and query downwards
  db_specimens.query <-
    patient.tbl %>%
    dplyr::select(.data$id, .data$external_id, .data$is_patient) %>%
    dplyr::rename(patient_id = .data$id) %>%
    dplyr::rename(patient_external_id = .data$external_id) %>%
    # Join to specimen
    dplyr::left_join(specimen.tbl, by = "patient_id") %>%
    dplyr::rename(specimen_id = .data$id) %>%
    dplyr::rename(specimen_external_id = .data$external_id) %>%
    dplyr::select(-.data$patient_id) %>%
    # Join to specimen_subset
    dplyr::left_join(specimen_subset.tbl, by = "specimen_id") %>%
    dplyr::rename(specimen_subset_id = .data$id) %>%
    dplyr::rename(specimen_subset_external_id = .data$external_id) %>%
    dplyr::select(-.data$specimen_subset_id, -.data$specimen_id)

  db_specimens.df <- dplyr::collect(db_specimens.query)
  return(db_specimens.df)
}

#' Retrieve DB Comments
#'
#' @param db_conn Database connection object
#'
#' @return db_comments.df A data-frame containing curated result comments from the database
#' @export
retrieve_db_comments <- function(db_conn) {
  # Connect to the relevant tables
  patient.tbl <- dplyr::tbl(db_conn, "patient")
  specimen.tbl <- dplyr::tbl(db_conn, "specimen")
  curated_result.tbl <- dplyr::tbl(db_conn, "curated_result")
  curated_result_comment.tbl <- dplyr::tbl(db_conn, "curated_result_comment")
  # Starting from the curated results table in the database...
  curated_result_comment.tbl %>%
    dplyr::inner_join(curated_result.tbl, by = c("curated_result_id" = "id")) %>%
    # Join to specimen
    dplyr::inner_join(specimen.tbl, by = c("specimen_id" = "id")) %>%
    dplyr::rename(specimen_external_id = .data$external_id) %>%
    # Join to patient
    dplyr::inner_join(patient.tbl, by = c("patient_id" = "id")) %>%
    dplyr::rename(patient_external_id = .data$external_id) %>%
    # Select matching fields
    dplyr::select(.data$patient_external_id, .data$specimen_external_id,
                  .data$comment, .data$added_by, .data$added_date) ->
    curated_comments.query
  # Collect results into a data frame
  curated_comments.df <- dplyr::collect(curated_comments.query)
  return(curated_comments.df)
}
