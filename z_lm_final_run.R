# Run from this script's directory. SNR changes error size; regression coefficients stay fixed.

source("run_simulation.R")

result <- run_linear_simulation(structure = "Additional interaction term", truth = "M2", snr = 0.5, n_train = 200, 
    n_test = 200, repetitions = 500, basis = "adaptive", kernel_bw = "adaptive", dim_basis = 1, predictor_distribution = "Normal", 
    error_distribution = "Normal", test_seed = 4100, training_seed = 42)

result$summary

result$drm_outcomes

result$raw


result$runtime
result$error_summary
