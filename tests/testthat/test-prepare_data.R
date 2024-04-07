library(JDRnet)

test_that("incorrect file format throws error", {
  test_files <- c("test_file.fastq")
  expect_error(suppressWarnings(.load_data(test_files)),
    "Invalid file format: test_file.fastq. Please make sure you provide a supported file format.")
})

test_that("list of correct length is made", {
  test_files <- c("test_file1.Rda", "test_file2.tsv", "test_file3.txt")
  expect_equal(class(suppressWarnings(.create_omics_list(test_files, overlap_samples = FALSE))), "list")
  expect_equal(length(suppressWarnings(.create_omics_list(test_files, overlap_samples = FALSE))), 3)
})

test_that("no overlapping samples throws an error", {
  test_files <- c("test_file2.tsv", "test_file4.txt")
  expect_error(suppressWarnings(.create_omics_list(test_files)),
    "No common samples between the omics. Please ensure omics share at least some samples.")
})

test_that("too few samples gives a warning", {
  test_files <- c("test_file2.tsv", "test_file5.txt")
  expect_warning(.create_omics_list(test_files),
  "Fewer than 15 samples common to all omics. This may cause unpredictable results.")
})

test_that("non-numeric data throws an error", {
  test_files <- c("test_file6.txt")
  expect_error(suppressWarnings(.create_omics_list(test_files)),
    "Data is not numeric. Please make sure all your omics are numeric.")
})