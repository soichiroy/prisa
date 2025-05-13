.CheckInput <- function(data, options) {
  if (!is.null(options$cluster_var_name)) {
    .CheckClusterSize(data, options$cluster_var_name)
  }
  invisible(data)
}

#' Check if clusters have enough observations, otherwise warn
#' 
#' @noRd 
.CheckClusterSize <- function(data, var_cluster) {
  k_min_obs_per_cluster <- 2
  # Check if the cluster variable is present in the data
  if (!var_cluster %in% names(data)) {
    stop(paste("Cluster variable", var_cluster, "must be in data."))
  }

  nobs_per_cluster <- table(data[[var_cluster]])
  # Check if observations within each cluster are greater than k_min_obs_per_cluster
  if (any(nobs_per_cluster < k_min_obs_per_cluster)) {
    warning(paste(
      "Cluster variable",
      var_cluster,
      "must have at least",
      k_min_obs_per_cluster,
      "observations per cluster in the labeled set."
    ))
  }
  invisible(data)
}
