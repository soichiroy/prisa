# Valid input: all clusters have sufficient observations
data_valid <- data.frame(
  cluster = rep(LETTERS[1:6], each = 2),
  value = 1:12
)
test_that(".CheckInput handles valid and invalid inputs correctly", {
  options <- list(cluster_var_name = "cluster")

  # Valid input
  expect_silent(.CheckInput(data_valid, options))

  # Invalid input: cluster_var_name not in data
  options_invalid <- list(cluster_var_name = "nonexistent")
  expect_error(
    .CheckInput(data_valid, options_invalid),
    "Cluster variable nonexistent must be in data."
  )
})

test_that(".CheckClusterSize warns when clusters have insufficient observations", {
  data <- data.frame(cluster = c("A", "A", "B", "C"), value = 1:4)
  var_cluster <- "cluster"

  expect_silent(.CheckClusterSize(data_valid, var_cluster))

  # Invalid input: some clusters have insufficient observations
  warnings <- capture_warnings(.CheckClusterSize(data, var_cluster))
  expect_match(
    warnings,
    "Cluster variable cluster must have at least 2 observations per cluster in the labeled set.",
    all = FALSE
  )
  expect_match(
    warnings,
    "Cluster variable cluster must have at least 5 clusters in the labeled set.",
    all = FALSE
  )
})

