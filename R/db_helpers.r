#' Connect to AML Project ORM Database
#'
#' @param dbname Database name
#' @param host Database host
#' @param user Database username
#'
#' @return db_conn A database connection object
#' @export
#'
#' @examples
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
#'
#' @examples
retrieve_db_libraries <- function(db_conn) {
  # Connect to the relevant tables
  patient.tbl <- tbl(db_conn, "patient")
  library.tbl <- tbl(db_conn, "library")
  specimen.tbl <- tbl(db_conn, "specimen")
  specimen_subset.tbl <- tbl(db_conn, "specimen_subset")
  platform_version.tbl <- tbl(db_conn, "platform_version")
  platform.tbl <- tbl(db_conn, "platform")
  # Start from patient, and query downwards
  patient.tbl %>%
    dplyr::select(id, external_id, is_patient) %>%
    dplyr::rename(patient_id = id,
           patient_external_id = external_id) %>%
    # Join to specimen
    dplyr::left_join(specimen.tbl, by = "patient_id") %>%
    dplyr::rename(specimen_id = id,
           specimen_external_id = external_id) %>%
    dplyr::select(-patient_id, -meta) %>%
    # Join to specimen_subset
    dplyr::left_join(specimen_subset.tbl, by = "specimen_id") %>%
    dplyr::rename(specimen_subset_id = id,
           specimen_subset_external_id = external_id) %>%
    # Join to library
    left_join(library.tbl, by = "specimen_subset_id") %>%
    dplyr::rename(library_id = id,
           library_name = name) %>%
    # Join to platform and platform_version
    dplyr::left_join(platform_version.tbl,
              by = c("platform_version_id" = "id")) %>%
    dplyr::rename(platform_version_name = name) %>%
    dplyr::left_join(platform.tbl,
              by = c("platform_id" = "id")) %>%
    dplyr::rename(platform_name = name) %>%
    dplyr::select(-fields, -platform_id, -platform_version_id, -specimen_subset_id,
                  -specimen_id, -library_id, -sequencing_effort_id,
                  -library_qc_info, -type) ->
    db_libraries.query
  db_libraries.df <- collect(db_libraries.query)
  return(db_libraries.df)
}

#' Retrieve DB Paths
#'
#' @param db_conn Database connection object
#'
#' @return db_paths.df A data-frame containing path information from the database
#' @export
#'
#' @examples
retrieve_db_paths <- function(db_conn) {
  # Connect to the relevant tables
  patient.tbl <- tbl(db_conn, "patient")
  library.tbl <- tbl(db_conn, "library")
  specimen.tbl <- tbl(db_conn, "specimen")
  specimen_subset.tbl <- tbl(db_conn, "specimen_subset")
  platform_version.tbl <- tbl(db_conn, "platform_version")
  platform.tbl <- tbl(db_conn, "platform")
  analysis.tbl <- tbl(db_conn, "analysis")
  path.tbl <- tbl(db_conn, "path")
  filetype.tbl <- tbl(db_conn, "filetype")
  method.tbl <- tbl(db_conn, "method")
  pipeline_version.tbl <- tbl(db_conn, "pipeline_version")
  # Start from patient, and query downwards
  patient.tbl %>%
    dplyr::select(id, external_id, is_patient) %>%
    dplyr::rename(patient_id = id,
                  patient_external_id = external_id) %>%
    # Join to specimen
    dplyr::left_join(specimen.tbl, by = "patient_id") %>%
    dplyr::rename(specimen_id = id,
                  specimen_external_id = external_id) %>%
    dplyr::select(-patient_id, -meta) %>%
    # Join to specimen_subset
    dplyr::left_join(specimen_subset.tbl, by = "specimen_id") %>%
    dplyr::rename(specimen_subset_id = id,
                  specimen_subset_external_id = external_id) %>%
    # Join to library
    left_join(library.tbl, by = "specimen_subset_id") %>%
    dplyr::rename(library_id = id,
                  library_name = name) %>%
    # Join to platform and platform_version
    dplyr::left_join(platform_version.tbl,
                     by = c("platform_version_id" = "id")) %>%
    dplyr::rename(platform_version_name = name) %>%
    dplyr::left_join(platform.tbl,
                     by = c("platform_id" = "id")) %>%
    dplyr::rename(platform_name = name) %>%
    # Join to analysis
    dplyr::left_join(analysis.tbl,
                     by = c("library_id")) %>%
    dplyr::rename(analysis_id = id) %>%
    # Join to path
    dplyr::left_join(path.tbl,
                     by = c("analysis_id")) %>%
    dplyr::rename(path_id = id) %>%
    # Join to filetype, method, and pipeline version
    dplyr::left_join(filetype.tbl,
                     by = c("filetype_id" = "id")) %>%
    dplyr::rename(filetype_name = name) %>%
    dplyr::left_join(method.tbl,
                     by = c("method_id" = "id")) %>%
    dplyr::rename(method_name = name) %>%
    dplyr::left_join(pipeline_version.tbl,
                     by = c("pipeline_version_id" = "id")) %>%
    dplyr::rename(pipeline_version_name = name) %>%
    # Remove rows that have no path information
    filter(!is.na(path)) %>%
    # Select out UUID fields and unused fields
    dplyr::select(-fields, -platform_id, -platform_version_id, -specimen_subset_id,
                  -specimen_id, -library_id, -sequencing_effort_id,
                  -library_qc_info, -type, -description,
                  -path_id, -filetype_id, -method_id, -pipeline_version_id,
                  -analysis_id) ->
    db_paths.query
  db_paths.df <- collect(db_paths.query)
  return(db_paths.df)
}


#' Retrieve DB Curated Results
#'
#' @param db_conn Database connection object
#'
#' @return db_curated_results.df A data-frame containing curated result information from the database
#' @export
#'
#' @examples
retrieve_db_curated_results <- function(db_conn) {
  # Connect to the relevant tables
  patient.tbl <- tbl(db_conn, "patient")
  specimen.tbl <- tbl(db_conn, "specimen")
  curated_result.tbl <- tbl(db_conn, "curated_result")
  # Starting from the curated results table in the database...
  curated_result.tbl %>%
    # Join to specimen
    dplyr::inner_join(specimen.tbl, by = c("specimen_id" = "id")) %>%
    dplyr::rename(specimen_external_id = external_id) %>%
    # Join to patient
    dplyr::inner_join(patient.tbl, by = c("patient_id" = "id")) %>%
    # Select matching fields
    dplyr::select(specimen_external_id, subtype, flt3_status,
                  npm1_status, kit_status, cebpa_status,
                  cytogenetic_risk, karyotype, fusions,
                  tier_one_mutation_status, kmt2a_status,
                  dnmt3a_status) %>%
    # Sort by specimen ID
    arrange(specimen_external_id) ->
    curated_results.query
    # Collect results into a data frame
    curated_results.df <- collect(curated_results.query)
    return(curated_results.df)
}
