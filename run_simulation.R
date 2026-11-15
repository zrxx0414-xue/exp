
source("drm_selection.R")
source("data_generating.R")
source("benchmark_aic_bic.R")
source("adaptive_basis.R")

library(dplyr)

run_linear_simulation <- function(structure = "Additional main effects", truth = "M2", snr = 1,
                                  n_train = 200L, n_test = 200L, repetitions = 500L,
                                  basis = "adaptive", kernel_bw="adaptive", dim_basis="adaptive",
                                  predictor_distribution = "Normal",
                                  error_distribution = "Normal", test_seed = 4100L,
                                  training_seed = 42L, alpha = 0.05, maxit = 2000L) {
  
  gamma <- if (truth == "M1") 0 else 1
  complexity <- candidate_parameter_count(structure)

  set.seed(test_seed)
  test_data <- generate_dataset(n_test, structure, truth, predictor_distribution = "Normal",
                                error_distribution = NULL, gamma = gamma)

  set.seed(training_seed)
  results <- vector("list", repetitions * 3L)
  result_id <- 1L
  signal_variances <- error_sds <- numeric(repetitions)

  for (replicate in seq_len(repetitions)) {
    training_data <- generate_dataset(n_train, structure, truth,
                                      predictor_distribution = predictor_distribution,
                                      error_distribution = error_distribution, gamma = gamma,
                                      snr = snr)
    signal_variances[replicate] <- attr(training_data, "signal_variance")
    error_sds[replicate] <- attr(training_data, "error_sd")
    models <- fit_candidate_models(training_data, structure)
    prediction_1 <- as.numeric(predict(models$M1, newdata = test_data))
    prediction_2 <- as.numeric(predict(models$M2, newdata = test_data))

    drm <- safe_drm_model_selection(test_data$y, prediction_1, prediction_2, basis,
                                    complexity = complexity, alpha = alpha, maxit = maxit,
                                    kernel_bw = kernel_bw, dim_basis = dim_basis)
    results[[result_id]] <- transform(drm, replicate = replicate, method = "DRM")
    result_id <- result_id + 1L

    benchmarks <- benchmark_aic_bic(models)
    for (method in names(benchmarks)) {
      results[[result_id]] <- data.frame(
        decision = unname(benchmarks[[method]]), reason = method, rejected = NA,
        equality_wald = NA_real_, equality_p = NA_real_, wald_1 = NA_real_, wald_2 = NA_real_,
        basis = NA_character_, valid = TRUE, error_message = NA_character_,
        replicate = replicate, method = method, stringsAsFactors = FALSE
      )
      result_id <- result_id + 1L
    }
  }

  raw <- bind_rows(results) |>
    mutate(structure = structure, truth = truth, snr = snr, n_train = n_train, n_test = n_test,
           train_signal_variance = signal_variances[replicate],
           error_sd = error_sds[replicate],
           correct = decision == truth)

  summary <- raw |>
    group_by(method) |>
    summarise(valid_runs = sum(!is.na(decision)), correct_selections = sum(correct, na.rm = TRUE),
              correct_rate = mean(correct, na.rm = TRUE), .groups = "drop")

  failures <- raw |>
    filter(method == "DRM") |>
    mutate(outcome = case_when(
      !valid ~ "Invalid DRM fit", correct ~ "Correct",
      !rejected ~ "Failed to reject equality",
      rejected ~ "Wrong Wald direction"
    )) |>
    count(outcome, name = "runs") |>
    mutate(proportion = runs / sum(runs))

  error_summary <- raw |>
    filter(method == "DRM", !valid) |>
    count(failure_stage, error_message, name = "runs")

  list(settings = list(structure = structure, truth = truth, snr = snr, gamma = gamma,
                       snr_definition = "Full signal sample variance / error sample variance, conditional on training design",
                       test_response = "Noiseless true signal",
                       n_train = n_train, n_test = n_test, repetitions = repetitions,
                       basis = basis, kernel_bw= kernel_bw, dim_basis=dim_basis, 
                       predictor_distribution = predictor_distribution,
                       error_distribution = error_distribution, test_seed = test_seed,
                       training_seed = training_seed),
       raw = raw, summary = summary, drm_outcomes = failures,
       error_summary = error_summary)
}

