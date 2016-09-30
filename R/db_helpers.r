#' Connect to AML Project ORM Database
#'
#' @param dbname
#' @param host
#' @param user
#'
#' @return db_conn
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
#' @param variables
#'
#' @return
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
    rename(patient_id = id,
           patient_external_id = external_id) %>%
    # Join to specimen
    left_join(specimen.tbl, by = "patient_id") %>%
    rename(specimen_id = id,
           specimen_external_id = external_id) %>%
    dplyr::select(-patient_id, -meta) %>%
    # Join to specimen_subset
    left_join(specimen_subset.tbl, by = "specimen_id") %>%
    rename(specimen_subset_id = id,
           specimen_subset_external_id = external_id) %>%
    # Join to library
    left_join(library.tbl, by = "specimen_subset_id") %>%
    rename(library_id = id,
           library_name = name) %>%
    # Join to platform and platform_version
    left_join(platform_version.tbl,
              by = c("platform_version_id" = "id")) %>%
    rename(platform_version_name = name) %>%
    left_join(platform.tbl,
              by = c("platform_id" = "id")) %>%
    rename(platform_name = name) %>%
    dplyr::select(-fields, -platform_id, -platform_version_id, -specimen_subset_id,
                  -specimen_id, -library_id, -sequencing_effort_id,
                  -library_qc_info, -type) ->
    db_libraries.query
  db_libraries.df <- collect(db_libraries.query)
  return(db_libraries.df)
}

