
# Run all tests
testthat::test_dir(file.path("tests", "testthat"))

# Run individual tests
testthat::test_file(file.path("tests", "testthat", "test-1-module.R"))
