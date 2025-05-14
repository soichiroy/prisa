#' Run bootstrap
#'
#' @importFrom foreach foreach %do% %dopar%
#' @import dplyr
#' @noRd
.RunBootstrap <- function(
  main_model,
  proxy_model,
  data_list,
  point_estimate,
  options,
  args_main_model,
  args_proxy_model
) {
  # Extract option values
  n_boot <- options$n_boot
  use_full <- options$use_full
  use_parallel <- options$use_parallel
  seed_value <- options$seed_value
  n_cores <- options$n_cores

  # Register parallel backend
  .SetupParallelBackend(
    use_parallel = use_parallel,
    n_cores = n_cores,
    seed_value = seed_value
  )

  var_cluster <- options$cluster_var_name
  if (isFALSE(options$debug_mode) && !is.null(var_cluster)) {
    var_cluster <- NULL
    message(
      "Clustering is currently disabled. Use debug_mode = TRUE to enable."
    )
  }

  # Bootstrap for the labeled data to estimate the variance-covariance matrix
  # of the main estimators based on the labeled data.
  boot_estimate_labeled <- foreach(
    i = seq_len(n_boot),
    .combine = rbind,
    .packages = c("dplyr"),
    .inorder = FALSE,
    .export = c(
      ".GetPointEstimates",
      ".ResampleDataFrame"
    )
  ) %dopar%
    {
      data_labeled_resampled <- list(
        dat_labeled = .ResampleDataFrame(data_list$dat_labeled, var_cluster)
      )
      .GetPointEstimates(
        main_model,
        proxy_model,
        data_labeled_resampled, 
        args_main_model,
        args_proxy_model,
        use_label_only = TRUE
      )
    }
  vcov_labeled <- cov(boot_estimate_labeled)

  # Exit the function if boot_full is FALSE (skip the bootstrap for the entire
  # data).
  if (isFALSE(use_full)) {
    # Return vcov estimate
    return(list(
      vcov_labeled = vcov_labeled,
      vcov_full = NULL,
      vcov_main_diff = NULL
    ))
  }

  # Estimate the variance covariance matrix of the biased estimator based on
  # the unlabeled data or the full data
  data_main <- data_list$dat_full

  # Bootstrap for the full (or unlabeled) data
  boot_estimate_full <- foreach(
    i = seq_len(n_boot),
    .combine = rbind,
    .inorder = FALSE,
    .packages = c("dplyr"),
    .export = c(".ResampleDataFrame")
  ) %dopar%
    {
      data_main_resampled <- .ResampleDataFrame(data_main, var_cluster)
      do.call(proxy_model, c(list(data_main_resampled), args_proxy_model))
    }

  # VCOV(tau_proxy_ell - tau_proxy_full)
  boot_estimate_diff <- boot_estimate_labeled[, -1] - boot_estimate_full
  vcov_full <- cov(boot_estimate_diff)

  # Cov(tau_main_ell, tau_proxy_ell - tau_proxy_full)
  vcov_main_diff <-
    as.vector(cov(boot_estimate_labeled[, 1], boot_estimate_diff))

  # Return bootstrapped variance-covariance estimates
  list(
    vcov_labeled = vcov_labeled,
    vcov_full = vcov_full,
    vcov_main_diff = vcov_main_diff
  )
}

#' Resample data frame that allows for cluster sampling
#'
#' @param df A data frame to be resampled.
#' @param cluster_var A string representing the name of the cluster variable.
#'  If NULL, the function will perform simple random sampling.``
#' @return A resampled data frame.
#' @noRd
#' @importFrom dplyr slice_sample across all_of
#' @importFrom tidyr nest unnest
.ResampleDataFrame <- function(df, cluster_var) {
  if (is.null(cluster_var)) {
    return(slice_sample(df, prop = 1, replace = TRUE))
  }

  # Resample the data frame based on the cluster variable
  # This operation can be slow for large data frames
  resampled_df <- df %>%
    group_by(across(all_of(cluster_var))) %>%
    nest() %>%
    ungroup() %>%
    slice_sample(prop = 1, replace = TRUE) %>%
    mutate(id_groups_tmp = row_number()) %>%
    unnest(cols = c("data"))

  return(resampled_df)
}

#' Setup parallel backend for bootstrap
#'
#' @importFrom doRNG registerDoRNG
#' @importFrom foreach registerDoSEQ
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @param use_parallel A logical value indicating whether to use parallel
#'   processing.
#' @param n_cores The number of cores to be used for parallel processing.
#' @param seed_value The seed value for the random number generator.
#' @noRd
.SetupParallelBackend <- function(use_parallel, n_cores, seed_value) {
  if (use_parallel) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    doRNG::registerDoRNG(seed_value)
    on.exit({
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    })
  } else {
    # When use_parallel is FALSE, use sequential processing
    foreach::registerDoSEQ()
  }
  invisible(use_parallel)
}
