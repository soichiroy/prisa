#' Run bootstrap
#'
#' @importFrom foreach foreach %do% %dopar%
#' @import dplyr
#' @importFrom doRNG registerDoRNG
#' @importFrom foreach registerDoSEQ
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
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
  n_main_estimates <- point_estimate$n_main_estimates
  n_proxy_estimates <- point_estimate$n_proxy_estimates

  # Register parallel backend
  foreach::registerDoSEQ()
  if (use_parallel) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    doRNG::registerDoRNG(seed_value)
    on.exit({
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    })
  }

  var_cluster <- options$cluster_var_name
  if (isFALSE(options$debug_mode) && !is.null(var_cluster)) {
    var_cluster <- NULL
    message(
      "Clustering is currently disabled. Use debug_mode = TRUE to enable."
    )
  }

  # Get all objects from global environment to export to workers
  # This ensures user-defined functions (main_model, proxy_model) have access
  # to their dependencies
  global_objects <- ls(envir = .GlobalEnv, all.names = TRUE)
  loaded_packages <- loadedNamespaces()

  # Bootstrap for the labeled data to estimate the variance-covariance matrix
  # of the main estimators based on the labeled data.
  boot_estimate_labeled <- foreach(
    i = seq_len(n_boot),
    .combine = rbind,
    .packages = loaded_packages,
    .inorder = FALSE,
    .export = c(
      ".GetPointEstimates",
      ".ResampleDataFrame",
      global_objects
    )
  ) %dopar%
    {
      data_labeled_resampled <- list(
        dat_labeled = .ResampleDataFrame(data_list$dat_labeled, var_cluster)
      )
      estimates <- .GetPointEstimates(
        main_model,
        proxy_model,
        data_labeled_resampled,
        args_main_model,
        args_proxy_model,
        use_label_only = TRUE
      )

      c(estimates$tau_ell, estimates$delta_ell)
    }

  # unlist main and proxy estimates and compute the covariance
  vcov_labeled <- cov(boot_estimate_labeled, use = "complete.obs")

  # Exit the function if use_full is FALSE (skip the bootstrap for the entire
  # data).
  if (isFALSE(use_full)) {
    # Return vcov estimate
    boot_estimates <- list(vcov_labeled = vcov_labeled)
    return(.ProcessCovarianceEstimates(
      boot_estimates,
      n_main_estimates = n_main_estimates,
      use_full = use_full
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
    .packages = loaded_packages,
    .export = c(
      ".ResampleDataFrame",
      ".FitModel",
      global_objects
    )
  ) %dopar%
    {
      data_main_resampled <- .ResampleDataFrame(data_main, var_cluster)
      .FitModel(proxy_model, data_main_resampled, args_proxy_model)
    }

  # VCOV(delta_ell - delta_full)
  idx_delta_ell <- (n_main_estimates + 1):(n_main_estimates + n_proxy_estimates)
  delta_ell <- boot_estimate_labeled[, idx_delta_ell, drop = FALSE]
  boot_estimate_diff <- delta_ell - boot_estimate_full
  vcov_full <- cov(boot_estimate_diff, use = "complete.obs")

  # Cov(tau_main_ell, delta_ell - delta_full)
  vcov_main_diff <- cov(
    boot_estimate_diff,
    boot_estimate_labeled[, 1:n_main_estimates, drop = FALSE]
  )

  # Return bootstrapped variance-covariance estimates
  boot_estimates <- list(
    vcov_labeled = vcov_labeled,
    vcov_full = vcov_full,
    vcov_main_diff = vcov_main_diff
  )
  return(.ProcessCovarianceEstimates(
    boot_estimates,
    n_main_estimates = n_main_estimates,
    use_full = use_full
  ))
}

#' Process the bootstrap-based variance-covariance estimates
#'
#' @param cov_estimates A list of bootstrap outputs.
#' @param n_main_estimates The number of main estimates.
#' @param use_full Boolean. Passed from the SetOptions function.
#' @noRd
.ProcessCovarianceEstimates <- function(
  cov_estimates,
  n_main_estimates,
  use_full
) {
  if (n_main_estimates < 1) {
    stop("n_main_estimates must be greater than 1.")
  }

  # Variance covariance matrix of the main (labeled only) estimator
  idx_main <- 1:n_main_estimates
  vcov_tau <- cov_estimates$vcov_labeled[idx_main, idx_main, drop = FALSE]

  if (
    use_full &&
      !is.null(cov_estimates$vcov_full) &&
      !is.null(cov_estimates$vcov_main_diff)
  ) {
    # Variance-covariance matrix of the proxy estimator.
    # When use_full is TRUE, this vcov is estimated based on the difference of
    # delta (proxy estimators) between the labeled and full data.
    vcov_delta <- cov_estimates$vcov_full

    # Covariance matrix between the main and proxy estimators
    cov_delta_tau <- cov_estimates$vcov_main_diff

    return(list(
      vcov_tau = vcov_tau,
      vcov_delta = vcov_delta,
      cov_delta_tau = cov_delta_tau
    ))
  }

  # Variance covariance matrix of the proxy estimator
  # When use_full is FALSE, this vcov is estimated based on the labeled data.
  vcov_delta <- cov_estimates$vcov_labeled[-idx_main, -idx_main, drop = FALSE]

  # Covariance matrix between the main and proxy estimators
  cov_delta_tau <- cov_estimates$vcov_labeled[-idx_main, idx_main, drop = FALSE]
  
  list(
    vcov_tau = vcov_tau,
    vcov_delta = vcov_delta,
    cov_delta_tau = cov_delta_tau
  )
}

#' Resample data frame that allows for cluster sampling
#'
#' @param df A data frame to be resampled.
#' @param cluster_var A string representing the name of the cluster variable.
#'  If NULL, the function will perform simple random sampling.``
#' @return A resampled data frame.
#' @noRd
#' @importFrom dplyr slice_sample across all_of ungroup group_by
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
