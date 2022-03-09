# Load helper functions ---------------------------------------------------

get_data <- function(path) {
  if(!is.null(CSV_ROOT_DIR)) {
    paste0(CSV_ROOT_DIR, path)
  } else {
    stop("Must specify location of CSV directory!")
  }
}

specd <- function(x, k) trimws(format(round(x, k), nsmall=k))

print_ci <- function(est, lwr, upr, k = 2) paste0(specd(est, k), " (", specd(lwr, k), ", ", specd(upr, k), ")")

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
  return(sqrt(mean((fit$y - stats::fitted(fit))^2)))
}

add_grf_rmse <- function(rf) {
  if (rf$type %in% c("binomial grf", "survival")) {
    return(sqrt(mean((as.numeric(rf$Y.orig) - 1 - rf$predictions[, 2])^2)))
  } else {
    return(sqrt(mean((rf$Y.orig - rf$predictions)^2)))
  }
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


winsorize <- function(x, pct) {
  upr <- quantile(x, pct + (1 - pct) / 2, na.rm = TRUE)
  lwr <- quantile(x, (1 - pct) / 2, na.rm = TRUE)
  
  x[x > upr] <- upr
  x[x < lwr] <- lwr
  
  return(x)
}

draw_bootstrap_samples <- function(data, id, time) {
  ids <- unique(data[[id]])
  samples <- sample(x = ids, size = length(ids), replace = TRUE)
  
  rows <- sapply(samples, function(x) which(data[[id]] == x), simplify = FALSE)
  times <- sapply(rows, function(x) length(x))
  
  boot <- data[unlist(rows), ]
  boot$id <- rep(1:length(samples), times = times)
  
  return(boot)
}

zero_out_lags <- function(data, lags) {
  if (lags == 0) {
    data[, c(paste0("lag1_", covs_tv), paste0("lag2_", covs_tv))] <- 0
  }
  else if (lags == 1) {
    data[, paste0("lag2_", covs_tv)] <- 0
  }
  return(data)
}
