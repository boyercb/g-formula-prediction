# Load helper functions ---------------------------------------------------

get_data <- function(path) {
  if(!is.null(CSV_ROOT_DIR)) {
    paste0(CSV_ROOT_DIR, path)
  } else {
    stop("Must specify location of CSV directory!")
  }
}

specd <- function(x, k) trimws(format(round(x, k), nsmall=k))

trim_glm <- function(fit) {
  fit$y <- c()
  fit$model <- c()
  
  fit$residuals <- c()
  fit$fitted.values <- c()
  fit$effects <- c()
  fit$qr$qr <- c()
  fit$linear.predictors <- c()
  fit$weights <- c()
  fit$prior.weights <- c()
  fit$data <- c()
  
  fit$family$variance <- c()
  fit$family$dev.resids <- c()
  fit$family$aic <- c()
  fit$family$validmu <- c()
  fit$family$simulate <- c()
  attr(fit$terms, ".Environment") <- c()
  attr(fit$formula, ".Environment") <- c()
  
  return(fit)
}

trim_truncreg <- function(fit) {
  fit$gradientObs <- c()
  fit$fitted.values <- c()
  fit$y <- c()
  
  return(fit)
}

trim_multinom <- function(fit) {
  fit$fitted.values <- c()
  fit$residuals <- c()
  fit$weights <- c()
  fit$y <- c()
  
  return(fit)
}

add_rmse <- function(fit) {
  return (sqrt(mean((fit$y - stats::fitted(fit))^2)))
}

add_stderr <- function(fit) {
  if (any(class(fit) == 'multinom')) {
    return(summary(fit)$coefficients)
  } else {
    return (stats::coefficients(summary(fit))[, "Std. Error"])
  }
}

add_vcov <- function(fit) {
  return(stats::vcov(fit))
}