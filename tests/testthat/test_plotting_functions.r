context("Plotting function checks")
library(amlpmpsupport)
library(tidyverse)

# Examples adapted from http://r-pkgs.had.co.nz/tests.html
#test_that("Test examples", {
#  expect_equal(1, 1)
#  a <- list(1:10, letters)
#  expect_output(str(a), "List of 2")
#  expect_match("Testing is fun!", "Testing")
#})

# Test plotting utility functions
test_that("Plotting function tests", {

  # Test that the `eng_log_breaks` function returns appropriate breaks
  expect_equal(eng_log_breaks(1), c(1, 2, 5))
  expect_equal(eng_log_breaks(100:1000), c(100, 200, 500, 1000, 2000, 5000))

})
