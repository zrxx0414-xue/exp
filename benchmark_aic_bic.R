# AIC and BIC benchmark selection

benchmark_values <- function(models) {
  stopifnot(all(c("M1", "M2") %in% names(models)))

  data.frame(
    model = c("M1", "M2"),
    AIC = c(AIC(models$M1), AIC(models$M2)),
    BIC = c(BIC(models$M1), BIC(models$M2)),
    stringsAsFactors = FALSE
  )
}

benchmark_aic_bic <- function(models) {
  values <- benchmark_values(models)
  c(AIC = values$model[which.min(values$AIC)], BIC = values$model[which.min(values$BIC)])
}

# Backward-compatible name used by the existing experiment scripts.
benchmark_decisions <- function(models, n_train = NULL, parameter_count = NULL) {
  benchmark_aic_bic(models)
}
