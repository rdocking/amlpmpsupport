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
    dplyr::rename(patient_id = id) %>%
    dplyr::rename(patient_external_id = external_id) %>%
    # Join to specimen
    dplyr::left_join(specimen.tbl, by = "patient_id") %>%
    dplyr::rename(specimen_id = id) %>%
    dplyr::rename(specimen_external_id = external_id) %>%
    dplyr::select(-patient_id, -meta) %>%
    # Join to specimen_subset
    dplyr::left_join(specimen_subset.tbl, by = "specimen_id") %>%
    dplyr::rename(specimen_subset_id = id) %>%
    dplyr::rename(specimen_subset_external_id = external_id) %>%
    # Join to library
    left_join(library.tbl, by = "specimen_subset_id") %>%
    dplyr::rename(library_id = id) %>%
    dplyr::rename(library_name = name) %>%
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
    dplyr::rename(patient_id = id) %>%
    dplyr::rename(patient_external_id = external_id) %>%
    # Join to specimen
    dplyr::left_join(specimen.tbl, by = "patient_id") %>%
    dplyr::rename(specimen_id = id) %>%
    dplyr::rename(specimen_external_id = external_id) %>%
    dplyr::select(-patient_id, -meta) %>%
    # Join to specimen_subset
    dplyr::left_join(specimen_subset.tbl, by = "specimen_id") %>%
    dplyr::rename(specimen_subset_id = id) %>%
    dplyr::rename(specimen_subset_external_id = external_id) %>%
    # Join to library
    left_join(library.tbl, by = "specimen_subset_id") %>%
    dplyr::rename(library_id = id) %>%
    dplyr::rename(library_name = name) %>%
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


# From a subtype indication, return a simplified disease type
#' Convert subtype indications to disease status
#'
#' @param subtype Curated disease subtype
#'
#' @return disease A simplified disease name
#' @export
#'
#' @examples
subtype_to_disease <- function(subtype) {
  disease <- NA
  if(subtype == 'AML-MDS'){
    disease <- subtype
  } else if(startsWith(subtype, "AML")){
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
    dplyr::rename(patient_external_id = external_id) %>%
    # Select matching fields
    dplyr::select(patient_external_id, specimen_external_id,
                  subtype, flt3_status,
                  npm1_status, kit_status, cebpa_status,
                  cytogenetic_risk, karyotype, fusions,
                  tier_one_mutation_status, kmt2a_status,
                  dnmt3a_status) %>%
    # Sort by specimen ID
    arrange(specimen_external_id) ->
    curated_results.query
    # Collect results into a data frame
    curated_results.df <- collect(curated_results.query)
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
#'
#' @examples
retrieve_db_specimens <- function(db_conn) {
  # Connect to the relevant tables
  patient.tbl <- tbl(db_conn, "patient")
  specimen.tbl <- tbl(db_conn, "specimen")
  specimen_subset.tbl <- tbl(db_conn, "specimen_subset")
  # Start from patient, and query downwards
  patient.tbl %>%
    dplyr::select(id, external_id, is_patient) %>%
    dplyr::rename(patient_id = id) %>%
    dplyr::rename(patient_external_id = external_id) %>%
    # Join to specimen
    dplyr::left_join(specimen.tbl, by = "patient_id") %>%
    dplyr::rename(specimen_id = id) %>%
    dplyr::rename(specimen_external_id = external_id) %>%
    dplyr::select(-patient_id) %>%
    # Join to specimen_subset
    dplyr::left_join(specimen_subset.tbl, by = "specimen_id") %>%
    dplyr::rename(specimen_subset_id = id) %>%
    dplyr::rename(specimen_subset_external_id = external_id) %>%
    dplyr::select(-specimen_subset_id, -specimen_id) ->
    db_specimens.query
  db_specimens.df <- collect(db_specimens.query)
  return(db_specimens.df)
}

#' Retrieve DB Comments
#'
#' @param db_conn Database connection object
#'
#' @return db_comments.df A data-frame containing curated result comments from the database
#' @export
#'
#' @examples
retrieve_db_comments <- function(db_conn) {
  # Connect to the relevant tables
  patient.tbl <- tbl(db_conn, "patient")
  specimen.tbl <- tbl(db_conn, "specimen")
  curated_result.tbl <- tbl(db_conn, "curated_result")
  curated_result_comment.tbl <- tbl(db_conn, "curated_result_comment")
  # Starting from the curated results table in the database...
  curated_result_comment.tbl %>%
    dplyr::inner_join(curated_result.tbl, by = c("curated_result_id" = "id")) %>%
    # dplyr::rename(specimen_external_id = external_id) %>%
    # Join to specimen
    dplyr::inner_join(specimen.tbl, by = c("specimen_id" = "id")) %>%
    dplyr::rename(specimen_external_id = external_id) %>%
    # Join to patient
    dplyr::inner_join(patient.tbl, by = c("patient_id" = "id")) %>%
    dplyr::rename(patient_external_id = external_id) %>%
    # Select matching fields
    dplyr::select(patient_external_id, specimen_external_id,
                  comment, added_by, added_date) ->
    curated_comments.query
  # Collect results into a data frame
  curated_comments.df <- collect(curated_comments.query)
  return(curated_comments.df)
}
