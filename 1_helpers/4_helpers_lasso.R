select_covariates <-
  function(outcome,
           covariates,
           data,
           glmnet = FALSE,
           lambda = "lambda.1se", 
           contvars = NULL,
           degree = NULL,
           knots = NULL,
           logit = FALSE,
           X.dependent = FALSE) {
    data <- na.omit(data[, c(outcome, covariates)])
  if (!is.null(contvars)) {
    if (!is.null(degree)) {
      poly_terms <- sapply(1:degree, function(x) {
        paste0("poly(", covariates[contvars], ", degree = ", degree, ")[,", x, "]") 
      })
      f <- as.formula(
        paste0("~ (", 
               paste0(
                 paste(poly_terms, collapse = " + "),
                 " + ",
                 paste(covariates[-contvars], collapse = " + ")
              ),
              ")^2 - 1"
            )
        )
    } else if (!is.null(knots)) {
      sn_terms <- sapply(1:(knots-1), function(x) {
        paste0("ns(", covariates[contvars], ", knots = ", knots, ")[,", x, "]") 
      })
      f <- as.formula(
        paste0("~ (", 
               paste0(
                 paste(sn_terms, collapse = " + "),
                 " + ",
                 paste(covariates[-contvars], collapse = " + ")
               ),
               ")^2 - 1"
        )
      )
    }
  } else {
    if (covariates %in% c("time")) {
      f <- as.formula(
        paste0("~ (", 
               paste0(
                 "as.factor(time) + ",
                 paste(covariates[-c("time")], collapse = " + ")
               ),
               ")^2 - 1"
        )
      )
    } else {
      f <- ~.^2 - 1
    }
  }
  x <- model.matrix(f, data = data[, covariates])
  y <- data[[outcome]]
  
  if (glmnet) {
    if (logit) {
      fit <- cv.glmnet(x, y, alpha = 1, nfolds = 10, family = "binomial")
    } else {
      fit <- cv.glmnet(x, y, alpha = 1, nfolds = 10)
    }
    inds <- which(abs(coef(fit, s = fit[[lambda]])) > 0)
    selected_covariates <- row.names(coef(fit, s = fit[[lambda]]))
    selected_covariates <- selected_covariates[inds]
    #selected_covariates <- coef(fit, s = fit$lambda.min)
  } else {
    if (logit) {
      fit <- rlassologit(x, y, penalty = list(X.dependent.lambda = X.dependent))
    } else {
      fit <- rlasso(x, y, penalty = list(X.dependent.lambda = X.dependent))
    }
    
    selected_covariates <- names(coef(fit))[abs(coef(fit)) > 0]
  }
 
  
  ret <- list(
    rlasso = fit,
    selected_covariates = selected_covariates[-1]
  )
  
  return(ret)
}

if (FALSE) {
  test <- select_covariates(
    outcome = "sbp",
    covariates = c(
      "sex",
      "age",
      "lag1_smk",
      "lag1_bmi",
      "lag1_dm",
      "lag1_sbp",
      "lag1_hdl",
      "lag1_tc",
      "lag1_hrx",
      "lag1_liprx",
      "smk",
      "bmi",
      "dm",
      "hrx",
      "liprx",
      "tc",
      "hdl",
      "time"
    ),
    data = analytic_long,
    contvars = c(2, 4, 6, 7, 8, 12, 16, 17, 18),
    degree = 3
  )
}