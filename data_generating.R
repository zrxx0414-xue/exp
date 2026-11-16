# Linear regression simulation and model selection.

structure_levels <- c("Additional main effects", "Additional quadratic term", "Additional interaction term")

truth_levels <- c("M1", "M2")

standardised_random <- function(n, distribution = "Normal") {
    values <- switch(distribution, Normal = rnorm(n), Uniform = runif(n, -sqrt(3), sqrt(3)), `Heavy-tailed t(3)` = rt(n, 
        df = 3)/sqrt(3), Bimodal = {
        component <- sample(c(-1, 1), n, replace = TRUE)
        (component + rnorm(n, sd = 0.5))/sqrt(1.25)
    }, `Right-skewed` = rexp(n) - 1, stop("Unknown distribution: ", distribution))
    as.numeric(scale(values))
}

generate_predictors <- function(n, distribution = "Normal") {
    data.frame(x1 = standardised_random(n, distribution), x2 = standardised_random(n, distribution), x3 = standardised_random(n, 
        distribution))
}

shared_signal <- function(data, structure) {
    switch(structure, `Additional main effects` = 1 + data$x1, `Additional quadratic term` = 1 + data$x1, `Additional interaction term` = 1 + 
        data$x1 + data$x2, stop("Unknown model structure: ", structure))
}

additional_component <- function(data, structure) {
    switch(structure, `Additional main effects` = (data$x2 + data$x3)/sqrt(2), `Additional quadratic term` = (data$x1^2 - 
        1)/sqrt(2), `Additional interaction term` = data$x1 * data$x2, stop("Unknown model structure: ", structure))
}

signal_components <- function(data, structure) {
    list(data = data, structure = structure, shared = shared_signal(data, structure), additional = additional_component(data, 
        structure))
}

generate_response <- function(components, gamma, error_distribution = NULL, error_sd = 1) {
    response <- components$shared + gamma * components$additional
    if (is.null(error_distribution)) {
        return(response)
    }
    response + error_sd * standardised_random(length(response), error_distribution)
}

generate_dataset <- function(n, structure, truth, predictor_distribution = "Normal", error_distribution = NULL, 
    gamma = 1, snr = 1) {
    predictors <- generate_predictors(n, predictor_distribution)
    components <- signal_components(predictors, structure)
    effective_gamma <- if (truth == "M1") 
        0
    else gamma
    mu <- components$shared + effective_gamma * components$additional
    signal_variance <- var(mu)
    if (!is.null(error_distribution) && (length(snr) != 1L || is.na(snr) || snr <= 0)) {
        stop("snr must be positive (Inf gives zero noise).")
    }
    error_sd <- if (is.null(error_distribution)) 
        0
    else sqrt(signal_variance/snr)
    predictors$y <- generate_response(components, effective_gamma, error_distribution, error_sd = error_sd)
    attr(predictors, "signal_variance") <- signal_variance
    attr(predictors, "error_sd") <- error_sd
    predictors
}

fit_candidate_models <- function(data, structure) {
    switch(structure, `Additional main effects` = list(M1 = lm(y ~ x1, data = data), M2 = lm(y ~ x1 + x2 + x3, 
        data = data)), `Additional quadratic term` = list(M1 = lm(y ~ x1, data = data), M2 = lm(y ~ x1 + I(x1^2), 
        data = data)), `Additional interaction term` = list(M1 = lm(y ~ x1 + x2, data = data), M2 = lm(y ~ x1 * 
        x2, data = data)), stop("Unknown model structure: ", structure))
}

candidate_parameter_count <- function(structure) {
    switch(structure, `Additional main effects` = c(M1 = 2L, M2 = 4L), `Additional quadratic term` = c(M1 = 2L, 
        M2 = 3L), `Additional interaction term` = c(M1 = 3L, M2 = 4L), stop("Unknown model structure: ", structure))
}

candidate_parameter_counts <- candidate_parameter_count

fit_candidate_predictions <- function(x_train, y_train, x_test, structure) {
    training_data <- transform(x_train, y = y_train)
    models <- fit_candidate_models(training_data, structure)
    list(M1 = as.numeric(predict(models$M1, newdata = x_test)), M2 = as.numeric(predict(models$M2, newdata = x_test)), 
        models = models)
}

generate_normal_predictors <- function(n) generate_predictors(n, "Normal")

generate_signal <- function(data, structure, gamma) {
    shared_signal(data, structure) + gamma * additional_component(data, structure)
}

