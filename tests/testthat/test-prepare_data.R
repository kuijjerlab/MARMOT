context("Data manipulation")
library(JDRnet)

test_that("incorrect file format throws error", {
  test_files <- c("test_file.fastq")
  expect_error(.load_data(test_files),
    "Invalid file format. Please make sure you provide a text or .RData file.")
})

test_that("list of correct length is made", {
  test_files <- c("test_file1.Rda", "test_file2.tsv", "test_file3.txt")
  expect_equal(class(.create_omics_list(test_files, sep = )), "list")
  expect_equal(length(.create_omics_list(test_files)), length(test_files))
})

test_that("non-overlapping samples throws an error", {
  test_files <- c("test_file2.tsv", "test_file4.txt")
  expect_error(.create_omics_list(test_files),
    "No common samples between the omics. Please ensure omics share at least some samples.")
})