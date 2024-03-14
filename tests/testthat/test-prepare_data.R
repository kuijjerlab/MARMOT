context("Data manipulation")
library(JDRnet)

test_that("list of correct length is made", {
  test_files <- c("test_file1.Rda", "test_file2.tsv", "test_file3.txt")
  expect_equal(class(.create_omics_list(test_files)), "list")
  expect_equal(length(.create_omics_list(test_files)), length(test_files))
})
