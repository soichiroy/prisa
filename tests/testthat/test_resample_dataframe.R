test_that(".ResampleDataFrame performs simple random sampling correctly", {
  set.seed(1234)
  data <- data.frame(id = 1:10, value = rnorm(10))

  # Simple random sampling
  resampled_data <- .ResampleDataFrame(data, cluster_var = NULL)

  # Check that the resampled data has the same number of rows
  expect_equal(nrow(resampled_data), nrow(data))
})

test_that(".ResampleDataFrame performs cluster sampling correctly", {
  set.seed(1111)
  data <- data.frame(cluster = rep(letters[1:3], each = 3), value = rnorm(9))

  # Cluster sampling
  resampled_data <- .ResampleDataFrame(data, cluster_var = "cluster")

  expect_equal(
    length(unique(data$cluster)),
    length(unique(resampled_data$id_groups_tmp))
  )
})
