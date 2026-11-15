#setwd("your directory")
source("run_simulation.R")

result <- run_linear_simulation(
  structure = "Additional interaction term",  #Additional main effects/Additional quadratic term/Additional interaction term
  truth = "M2",    #"M1"/"M2"
  snr = 0.5,
  n_train = 200,
  n_test = 200,
  repetitions = 500,
  basis = "adaptive",
  kernel_bw = "adaptive",
  dim_basis = 1,  # use 1, 2, or "adaptive"
  predictor_distribution = "Normal",   #"Normal"/"Uniform"/"Heavy-tailed t(3)"/"Bimodal"/"Right-skewed"
  error_distribution = "Normal",      #"Normal"/"Uniform"/"Heavy-tailed t(3)"/"Bimodal"/"Right-skewed"
  test_seed = 4100,
  training_seed = 42
)

result$summary   #"correct slection rate"
result$drm_outcomes    #detail for each replication: correct/Failed to reject equality/Wrong Wald direction/Invalid DRM fit
result$raw      #complete results: equalty test wald/p vlue; individual wald value; aic/bic; etc.
