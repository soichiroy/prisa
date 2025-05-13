test_that(".CheckInput handles valid and invalid inputs correctly", {
  data <- data.frame(cluster = c("A", "A", "B", "B", "C", "C"), value = 1:6)
  options <- list(cluster_var_name = "cluster")

  # Valid input
  expect_silent(.CheckInput(data, options))

  # Invalid input: cluster_var_name not in data
  options_invalid <- list(cluster_var_name = "nonexistent")
  expect_error(
    .CheckInput(data, options_invalid),
    "Cluster variable nonexistent must be in data."
  )
})

test_that(".CheckClusterSize warns when clusters have insufficient observations", {
  data <- data.frame(cluster = c("A", "A", "B", "C"), value = 1:4)
  var_cluster <- "cluster"

  # Valid input: all clusters have sufficient observations
  data_valid <- data.frame(
    cluster = c("A", "A", "B", "B", "C", "C"),
    value = 1:6
  )
  expect_silent(.CheckClusterSize(data_valid, var_cluster))

  # Invalid input: some clusters have insufficient observations
  expect_warning(
    .CheckClusterSize(data, var_cluster),
    "Cluster variable cluster must have at least 2 observations per cluster in the labeled set."
  )
})

