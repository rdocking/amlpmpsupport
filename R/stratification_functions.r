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
  } else if (cyto_status %in% c('adverse_risk', 'poor_risk')) {
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
  return(tibble::tibble(status = status, rationale = rationale))
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
    } else if (cyto_status %in% c('adverse_risk', 'poor_risk')) {
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
    return(tibble::tibble(status = status, rationale = rationale))
  }


#' Apply ELN2017-RNA stratification
#'
#' @param sv_status Structural variant status
#' @param npm1 NPM1 status
#' @param flt3_itd FLT3-ITD status
#' @param flt3_itd_support FLT3-ITD support (high or low)
#' @param cebpa CEBPA status
#' @param tp53 TP53 status
#' @param runx1 RUNX1 status
#' @param asxl1 ASXL1 status
#' @param expression_outlier Expression outlier status
#' @param debug Print debug output
#'
#' @return A dataframe containing the stratification category and rationale
#' @export
apply_eln2017_rna <- function(sv_status, npm1, flt3_itd, flt3_itd_support, cebpa,
                              tp53, runx1, asxl1, expression_outlier, debug = FALSE){

  if (debug) {
    print(glue::glue("In: SV status: {sv_status} NPM1: {npm1} FLT3: {flt3_itd} {flt3_itd_support}"))
  }

  # Treat NA values as 'missing'
  if (is.na(npm1)) npm1 = 'missing'
  if (is.na(cebpa)) cebpa = 'missing'
  if (is.na(tp53)) tp53 = 'missing'
  if (is.na(runx1)) runx1 = 'missing'
  if (is.na(asxl1)) asxl1 = 'missing'
  if (is.na(flt3_itd)) flt3_itd = 'missing'
  if (is.na(expression_outlier)) expression_outlier = 'missing'

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
  # Check favourable-risk first
  if (sv_status %in% c('t_15_17', 't_8_21', 'inv_16')) {
    status <- 'favourable'
    rationale <- 'Favourable structural variant'

    # MLL translocations and KMT2A PTD
  } else if (sv_status == 't_9_11') {
    status <- 'intermediate'
    rationale <- 'MLL translocation'

  } else if (sv_status == 'KMT2A-PTD positive') {
    status <- 'adverse'
    rationale <- 'KMT2A-PTD positive'

    # Adverse-risk SVs
  } else if (sv_status == 'adverse_risk_fusion') {
    status <- 'adverse'
    rationale <- 'Adverse-risk structural variant'

    # Expression outlier status
  } else if (expression_outlier == 'MECOM high') {
    status <- 'adverse'
    rationale <- 'Outlier MECOM expression'

    # Intermediate-risk SVs - check mutations
  } else if (sv_status == 'intermediate_risk') {

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

  # Debug statements
  if (debug) {
    print(glue::glue("Out: {status} Rationale: {rationale}"))
  }

  # Return the status and rationale
  return(tibble::tibble(status = status, rationale = rationale))
}


#' Update an existing stratification by applying a gene expression score
#'
#' @param original_cat The original stratification category
#' @param expression_score The gene expression signature score
#' @param fav_theshold Threshold for adjusting stratification to 'favourable' category
#' @param adverse_threshold Threshold for adjusting stratification to 'adverse' category
#'
#' @return A revised stratification category and justification
#' @export
apply_expression_signature_reclassification <- function(original_cat, expression_score, fav_threshold, adverse_threshold){

  # Note that the signatures are scaled so that higher-values = worse outcomes
  status <- original_cat
  rationale <- "Second tertile expression"
  if (expression_score <= fav_threshold) {
    status <- 'favourable'
    rationale <- "First tertile expression"
  } else if (expression_score >= adverse_threshold) {
    status <- 'adverse'
    rationale <- "Third tertile expression"
  }
  return(tibble::tibble(status = status, rationale = rationale))

}
