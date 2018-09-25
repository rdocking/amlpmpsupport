# This script contains functions for applying stratification models based on molecular markers


#' Apply ELN2017 stratification based on cytogenetics
#'
#' @param cyto_status Cytogenetic status
#' @param npm1 NPM1 status
#' @param flt3_itd FLT3-ITD status
#' @param flt3_itd_support FLT3-ITD support (high or low)
#' @param cebpa CEBPA status
#' @param tp53 TP53 status
#' @param runx1 RUNX1 status
#' @param asxl1 ASXL1 status
#'
#' @return A stratification status and rationale
#' @export
#'
#' @examples
apply_eln2017_cyto <-
  function(cyto_status, npm1, flt3_itd, flt3_itd_support,
           cebpa, tp53, runx1, asxl1){

  # Treat NA values as 'missing'
  if (is.na(cyto_status)) cyto_status = 'missing'
  if (is.na(npm1)) npm1 = 'missing'
  if (is.na(cebpa)) cebpa = 'missing'
  if (is.na(tp53)) tp53 = 'missing'
  if (is.na(asxl1)) asxl1 = 'missing'
  if (is.na(runx1)) runx1 = 'missing'
  if (is.na(flt3_itd)) flt3_itd = 'missing'

  # Set FLT3 status as 'high' or 'low' based on presence and support
  # TODO is to fix this to be something more rigorous
  if (flt3_itd == 'positive' & flt3_itd_support == 'FLT3-ITD_high') {
    flt3_status = 'high'
  } else {
    flt3_status = 'low'
  }

  # Start as intermediate
  status <- 'intermediate'
  rationale <- 'No rules applied'

  # Check cytogenetics, then mutations
  if (cyto_status %in% c('t_15_17', 'cbf_fusion')) {
    status <- 'favourable'
    rationale <- 'Favourable cytogenetics'
  # Poor-risk cytogenetics
  } else if (cyto_status == 'adverse_risk') {
    status <- 'adverse'
    rationale <- 'Adverse-risk cytogenetics'
  # Intermediate-risk cyto and normal karyotype - check mutations
  } else if (cyto_status %in% c('normal_karyotype', 'intermediate_risk', 'missing')) {

    status <- 'intermediate'
    rationale <- 'Intermediate-risk cytogenetics'

    # CEBPA biallelic
    if (cebpa == 'CEBPAbi') {
      status <- 'favourable'
      rationale <- 'Intermediate SVs, biallelic CEBPA'
      # Then FLT3 and NPM1
    } else if (npm1 == 'NPM1fs' & flt3_status == 'low') {
      status <- 'favourable'
      rationale <- 'Intermediate SVs, NPM1+FLT3-'
      # Then FLT3 positive - note the use of a hard-coded cutoff here
    } else if (npm1 == 'NPM1fs' & flt3_status == 'high') {
      status <- 'intermediate'
      rationale <- 'Intermediate SVs, NPM1+FLT3+'
    } else if (flt3_status == 'high') {
      status <- 'adverse'
      rationale <- 'Intermediate SVs, NPM1-FLT3+'
      # Then TP53, RUNX1, ASXL1
    } else if (!tp53 == 'missing') {
      status <- 'adverse'
      rationale <- 'Intermediate SVs, TP53 mutation'
    } else if (!runx1 == 'missing') {
      status <- 'adverse'
      rationale <- 'Intermediate SVs, RUNX1 mutation'
    } else if (!asxl1 == 'missing') {
      status <- 'adverse'
      rationale <- 'Intermediate SVs, ASXL1 mutation'
    }
  }
  # Return the status and rationale
  return(data_frame(status = status, rationale = rationale))
}


#' Apply ELN2015 stratification based on cytogenetics
#'
#' @param cyto_status Cytogenetic status
#' @param npm1 NPM1 status
#' @param flt3_itd FLT3-ITD status
#' @param cebpa CEBPA status
#'
#' @return A stratification status and rationale
#' @export
#'
#' @examples
apply_eln2015_cyto <- function(cyto_status, npm1, flt3_itd, cebpa){

    # Treat NA values as 'missing'
    if (is.na(cyto_status)) cyto_status = 'missing'
    if (is.na(npm1)) npm1 = 'missing'
    if (is.na(cebpa)) cebpa = 'missing'
    if (is.na(flt3_itd)) flt3_itd = 'missing'

    # Convert all values except 'positive' for FLT3-ITD to 'negative'
    if (!flt3_itd == 'positive') flt3_itd = 'negative'

    # Start as intermediate, no rules applied
    status <- 'intermediate-II'
    rationale <- 'Intermediate SVs, NPM1-FLT3-'

    # Check cytogenetics, then mutations
    if (cyto_status %in% c('t_15_17', 'cbf_fusion')) {
      status <- 'favourable'
      rationale <- 'Favourable cytogenetics'
    # Poor-risk cytogenetics
    } else if (cyto_status == 'adverse_risk') {
      status <- 'adverse'
      rationale <- 'Adverse-risk cytogenetics'
    # Intermediate-risk cyto and normal karyotype - check mutations
    } else if (cyto_status %in% c('normal_karyotype', 'intermediate_risk', 'missing')) {

      status <- 'intermediate-II'
      rationale <- 'Intermediate-risk cytogenetics'

      # True intermediates - KMT2A fusions and other alterations
      if (cyto_status %in% 'intermediate_risk') {
        status <- 'intermediate-II'
        rationale <- 'Intermediate-risk cytogenetics'
      # CEBPA biallelic
      } else if (cebpa == 'CEBPAbi') {
        status <- 'favourable'
        rationale <- 'Intermediate SVs, biallelic CEBPA'
      # Then FLT3 and NPM1
      } else if (npm1 == 'NPM1fs' & flt3_itd == 'negative') {
        status <- 'favourable'
        rationale <- 'Intermediate SVs, NPM1+FLT3-'
      # Then FLT3 positive
      } else if (npm1 == 'NPM1fs' & flt3_itd == 'positive') {
        status <- 'intermediate-I'
        rationale <- 'Intermediate SVs, NPM1+FLT3+'
      } else if (npm1 != 'NPM1fs' & flt3_itd == 'positive') {
        status <- 'intermediate-I'
        rationale <- 'Intermediate SVs, NPM1-FLT3+'
      } else if (npm1 != 'NPM1fs' & flt3_itd == 'negative') {
        status <- 'intermediate-I'
        rationale <- 'Intermediate SVs, NPM1-FLT3-'
      }
    }
    # Return the status and rationale
    return(data_frame(status = status, rationale = rationale))
  }
