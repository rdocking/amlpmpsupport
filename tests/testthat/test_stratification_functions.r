context("Stratification function checks")
library(amlpmpsupport)
library(tidyverse)

# Examples adapted from http://r-pkgs.had.co.nz/tests.html
#test_that("Test examples", {
#  expect_equal(1, 1)
#  a <- list(1:10, letters)
#  expect_output(str(a), "List of 2")
#  expect_match("Testing is fun!", "Testing")
#})

# Test stratifications for different models
test_that("Testing ELN2017-Cyto", {

  # Set up default inputs
  cyto_status <- NA
  npm1 <- NA
  flt3_itd <- NA
  flt3_itd_support <- NA
  cebpa <- NA
  tp53 <- NA
  runx1 <- NA
  asxl1 <- NA

  # All missing
  df <- apply_eln2017_cyto(cyto_status, npm1, flt3_itd, flt3_itd_support,
                           cebpa, tp53, runx1, asxl1)
  expect_equal(df$status, 'intermediate')
  expect_equal(df$rationale, 'Intermediate-risk cytogenetics')

  # Good-risk cyto
  cyto_status <- 't_15_17'
  df <- apply_eln2017_cyto(cyto_status, npm1, flt3_itd, flt3_itd_support,
                         cebpa, tp53, runx1, asxl1)
  expect_equal(df$status, 'favourable')
  expect_equal(df$rationale, 'Favourable cytogenetics')

  cyto_status <- 'cbf_fusion'
  df <- apply_eln2017_cyto(cyto_status, npm1, flt3_itd, flt3_itd_support,
                           cebpa, tp53, runx1, asxl1)
  expect_equal(df$status, 'favourable')
  expect_equal(df$rationale, 'Favourable cytogenetics')

  # Adverse-risk cyto
  cyto_status <- 'adverse_risk'
  df <- apply_eln2017_cyto(cyto_status, npm1, flt3_itd, flt3_itd_support,
                           cebpa, tp53, runx1, asxl1)
  expect_equal(df$status, 'adverse')
  expect_equal(df$rationale, 'Adverse-risk cytogenetics')

  # Intermediate with no mutations
  cyto_status <- 'intermediate_risk'
  df <- apply_eln2017_cyto(cyto_status, npm1, flt3_itd, flt3_itd_support,
                           cebpa, tp53, runx1, asxl1)
  expect_equal(df$status, 'intermediate')
  expect_equal(df$rationale, 'Intermediate-risk cytogenetics')


#   # CEBPA biallelic
#   if (cebpa == 'CEBPAbi') {
#     status <- 'favourable'
#     rationale <- 'Intermediate SVs, biallelic CEBPA'
#     # Then FLT3 and NPM1
#   } else if (npm1 == 'NPM1fs' & flt3_status == 'low') {
#     status <- 'favourable'
#     rationale <- 'Intermediate SVs, NPM1+FLT3-'
#     # Then FLT3 positive - note the use of a hard-coded cutoff here
#   } else if (npm1 == 'NPM1fs' & flt3_status == 'high') {
#     status <- 'intermediate'
#     rationale <- 'Intermediate SVs, NPM1+FLT3+'
#   } else if (flt3_status == 'high') {
#     status <- 'adverse'
#     rationale <- 'Intermediate SVs, NPM1-FLT3+'
#     # Then TP53, RUNX1, ASXL1
#   } else if (!tp53 == 'missing') {
#     status <- 'adverse'
#     rationale <- 'Intermediate SVs, TP53 mutation'
#   } else if (!runx1 == 'missing') {
#     status <- 'adverse'
#     rationale <- 'Intermediate SVs, RUNX1 mutation'
#   } else if (!asxl1 == 'missing') {
#     status <- 'adverse'
#     rationale <- 'Intermediate SVs, ASXL1 mutation'
#   }
#

})
