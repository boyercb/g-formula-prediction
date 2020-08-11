#' Fit GLM on Covariate
#'
#' This internal function fits a generalized linear model (GLM) for a single covariate using the observed data.
#'
#' @param formula 
#' @param family 
#' @param link 
#' @param data 
#' @param control 
#'
#' @return            Fitted model for the covariate.
#' @keywords internal

fit_glm <- function(formula, family, link = NULL, data, control = NULL) {
  if (is.null(link)) {
    famtext <- paste(family, "()", sep = "")
  } else {
    famtext <- paste(family, "(link = ", link, ")", sep = "")
  }
  
  # Fit GLM for covariate using user-specified formula
  if (!is.null(control)){
    fit <- stats::glm(
      formula, 
      family = eval(parse(text = famtext)),
      data = data,
      y = TRUE,
      control = control
    )
  } else {
    fit <- stats::glm(
      formula, 
      family = eval(parse(text = famtext)),
      data = data,
      y = TRUE
    )
  }
  
  fit$rmse <- add_rmse(fit)
  fit$stderrs <- add_stderr(fit)
  fit$vcov <- add_vcov(fit)
  fit$outcome <- all.vars(update(formula, . ~ 1))
  
  fit$type <- family

  return(fit)
}

#' Fit Multinomial Model on Covariate
#'
#' This internal function fits a multinomial regression model for a categorical covariate using the observed data.
#'
#' @param formula 
#' @param data 
#' @param control 
#'
#' @return            Fitted model for the covariate at index \eqn{j}.
#' @keywords internal

fit_multinomial <- function(formula, data, control = NULL){
  if (!is.null(control)) {
    args <- c(list(formula = formula, data = data), trace = FALSE, control)
    fit <- do.call(nnet::multinom, args = args)
  } else {
    fit <- nnet::multinom(formula = formula,
                          data = data,
                          trace = FALSE)
  }
  
  fit$stderr <- add_stderr(fit)
  fit$vcov <- add_vcov(fit)
  fit$outcome <- all.vars(update(formula, . ~ 1))
  fit$type <- 'categorical'
  
  return(fit)
}

#' Fit Zero-Inflated Normal Model on Covariate
#'
#' This internal function models a zero-inflated normal distribution through the combined
#' use of a generalized linear model (GLM) fit on a zero vs. non-zero indicator
#' and a GLM fit on all non-zero values.
#'
#' @param formula 
#' @param family 
#' @param link 
#' @param data 
#' @param control 
#'
#' @return            Fitted model for the covariate at index \eqn{j}.
#' @keywords internal
#' @import data.table

fit_zeroinfl_normal <- function(formula, link = NULL, data, control = NULL) {
  if (is.na(link)) {
    famtext <- "gaussian()"
  } else {
    famtext <- paste("gaussian(link = ", link, ")", sep = "")
  }
  
  cov.name <- all.vars(update(formula, . ~ 1))
  
  data[, paste("I_", cov.name, sep = "")] <- data[[cov.name]]
  data[[paste("I_", cov.name, sep = "")]][data[[cov.name]] != 0] <- 1
  
  # Take log to ensure that no negative values are predicted
  data[, paste("log_", cov.name, sep = "")] <- 0
  data[data[[cov.name]] != 0][, paste("log_", cov.name, sep = "")] <-
    log(data[data[[cov.name]] != 0][[cov.name]])
  
  # Fit binomial model on indicator of whether covariate is 0 or not
  if (!is.null(control)) {
    fit0 <- stats::glm(stats::as.formula(paste("I_", formula, sep = "")), family = stats::binomial(),
                       control = control[[1]], data = data)
  } else {
    fit0 <- stats::glm(stats::as.formula(paste("I_", formula, sep = "")), family = stats::binomial(),
                       data = data)
  }
  
  
  # Fit Gaussian model on data points for which covariate does not equal 0
  if (!is.null(control)) {
    fit <- stats::glm(stats::as.formula(paste("log_", formula, sep = "")),
                       family = eval(parse(text = famtext)),
                       control = ccontrol[[2]],
                       data = data[data[[cov.name]] != 0])
  } else {
    fit <- stats::glm(stats::as.formula(paste("log_", formula, sep = "")),
                       family = eval(parse(text = famtext)),
                       data = data[data[[cov.name]] != 0])
  }
  
  fit0$rmse <- add_rmse(fit0)
  fit$rmse <- add_rmse(fit)
  
  fit0$stderr <- add_stderr(fit0)
  fit$stderr <- add_stderr(fit)
  
  fit0$vcov <- add_vcov(fit0)
  fit$vcov <- add_vcov(fit)
  
  fit0$outcome <- all.vars(update(formula, . ~ 1))
  fit$outcome <- all.vars(update(formula, . ~ 1))
  
  fit0$type <- 'zero-inflated normal'
  fit$type <- 'zero-inflated normal'
  
  return(list(fit0, fit))
}

#' Fit Bounded Normal Model on Covariate
#'
#' This internal function models a covariate using a "bounded normal" distribution
#' by first standardizing the covariate values to the range [0, 1], noninclusive,
#' then fitting a generalized linear model (GLM) under the Gaussian family function.
#'
#' @param formula 
#' @param family 
#' @param link 
#' @param data 
#' @param control 
#'
#' @return            Fitted model for the covariate at index \eqn{j}.
#' @keywords internal
#' @import data.table

fit_bounded_continuous <- function(formula, link = NULL, data, control = NULL) {
  
  cov.name <- all.vars(update(formula, . ~ 1))
  
  data[, paste("norm_", cov.name, sep = "")] <-
    (data[[cov.name]] - min(data[[cov.name]])) /
    (max(data[[cov.name]]) - min(data[[cov.name]]))
  
  if (!is.null(link)) {
    if (!is.null(control)) {
      fit <- stats::glm(
        stats::as.formula(paste("norm_", formula, sep = "")),
        family = stats::gaussian(link = link), 
        data = data, 
        y = TRUE,
        control = control
      )
    } else {
      fit <- stats::glm(
        stats::as.formula(paste("norm_", formula, sep = "")),
        family = stats::gaussian(link = link), 
        data = data, 
        y = TRUE
      )
    }
  } else {
    if (!!is.null(control)) {
      fit <- stats::glm(
        stats::as.formula(paste("norm_", formula, sep = "")),
        family = stats::gaussian(), 
        data = data, 
        y = TRUE,
        control = control
      )
    } else {
      fit <- stats::glm(
        stats::as.formula(paste("norm_", formula, sep = "")),
        family = stats::gaussian(), 
        data = data, 
        y = TRUE
      )
    }
  }
  
  fit$rmse <- add_rmse(fit)
  fit$stderr <- add_stderr(fit)
  fit$vcov <- add_vcov(fit)
  fit$outcome <- cov.name
  fit$type <- 'bounded normal'
  return(fit)
}

#' Fit Truncated Normal Model on Covariate
#'
#' This internal function models a covariate using a normal distribution truncated on one
#' side at a user-specified cutoff.
#'
#' @param formula 
#' @param family 
#' @param link 
#' @param data 
#' @param control 
#' @param direction 
#' @param point 
#'
#' @return            Fitted model for the covariate at index \eqn{j}.
#' @keywords internal

fit_trunc_normal <- function(formula, link = NULL, data, control = NULL, direction = NULL, point = NULL) {
  if (!is.null(control)) {
    args <- c(list(
      formula = formula,
      data = data,
      point = point,
      direction = direction,
      y = TRUE
    ),
    control)
    
    fit <- do.call(truncreg::truncreg, args = args)
  } else {
    fit <- truncreg::truncreg(
      formula, 
      data = data, 
      point = point,
      direction = direction, 
      y = TRUE
    )
  }
  
  fit$rmse <- add_rmse(fit)
  fit$stderr <- add_stderr(fit)
  fit$vcov <- add_vcov(fit)
  fit$outcome <- all.vars(update(formula, . ~ 1))
  fit$type <- "truncated normal"
  return(fit)
}

#' Fit Covariate Models
#'
#' This internal function fits a model for each covariate using the observed data.
#'
#' @param models 
#' @param data 
#'
#' @return                A list of fitted models, one for each covariate in \code{covnames}.
#' @keywords internal
#' @import data.table

fit_X_models <- function(models, data, time) {
  # Create list of models for covariates
  cov.names <- names(models)
  subdata <- data[data[[time]] > 0, ]
  
  X.fits <- lapply(
    models,
    function(model) {
      if (!is.null(model$restrict)) {
          subdata <- subset(subdata, eval(parse(text = model$restrict$subset)))
      }
      
      fit <- switch(
        model$family,
        'binomial' = fit_glm(
          formula = model$formula, 
          family = 'binomial', 
          link = model$link, 
          data = subdata, 
          control = model$control
          ),
        
        'normal' = fit_glm(
          formula = model$formula, 
          family = 'gaussian',
          link = model$link, 
          data = subdata, 
          control = model$control
          ),
        
        'categorical' = fit_multinomial(
          formula = model$formula, 
          data = subdata, 
          control = model$control
          ),
       
         'zero-inflated normal' = fit_zeroinfl_normal(
          formula = model$formula, 
          link = model$link, 
          data = subdata, 
          control = model$control
          ),
        
        'bounded normal' = fit_bounded_continuous(
          formula = model$formula, 
          link = model$link, 
          data = subdata, 
          control = model$control
          ),
        
        'truncated normal' = fit_trunc_normal(
          formula = model$formula,
          link = model$link,
          data = subdata,
          control = model$control,
          point = model$point,
          direction = model$direction
        )
      )
      
      if (model$family == 'normal') {
        fit$type <- 'normal'
      }
      
      fit$restrict <- model$restrict 
      
      return(fit)
    })
  
  names(X.fits) <- cov.names
  return(X.fits)
}

#' Fit Outcome Model
#'
#' This internal function fits a generalized linear model (GLM) for the outcome variable using the observed data.
#'
#' @param data           A data frame, list or environment (or object coercible by
#'                       as.data.frame to a data frame) containing the variables in the model.
#'                       If not found in data, the variables are taken from environment(formula), 
#'                       typically the environment from which glm is called.
#' @param EOF            Only use end of follow up time point when fitting outcome model
#' @param model          Model statement for the outcome variable.
#' @param time           Time variable.
#'
#' @return               Fitted model for the outcome variable.
#' @keywords internal
#' @import data.table

fit_Y_model <- function(model, data, EOF = FALSE, time = NULL) {
  family <- switch(
    model$family,
    "normal" = stats::gaussian(link = model$link),
    "binomial" = stats::binomial(link = model$link)
  )
  
  outcome <- all.vars(update(model$formula, . ~ 1))
  
  if (EOF) {
    if (is.null(time)) {
      stop("Must specify time variable for EOF.")
    }
    data <- data[data[[time]] == max(unique(data[[time]]))]
  }
  
  if (!is.null(model$subset)) {
    
    if (family == "normal") {
      data[, paste("norm_", outcome, sep = "")] <-
        (data[[outcome]] - min(data[[outcome]]))/(max(data[[outcome]]) - min(data[[outcome]]))
      
      fitY <- stats::glm(
        stats::as.formula(paste("norm_", model, sep = "")), 
        family = family,
        data = subset(data, eval(parse(text = model$subset))),
        y = TRUE
      )
    } else {
      fitY <- stats::glm(
        formula = model$formula,
        family = family,
        data = data,
        y = TRUE
      )
    }
  } else { # Case where there are no restrictions on outcome variable
    # Fit GLM for outcome variable using user-specified model and entire dataset
    fitY <- stats::glm(
      formula = model$formula,
      family = family,
      data = data,
      y = TRUE
    )
  }
  
  fitY$rmse <- add_rmse(fitY)
  fitY$stderr <- add_stderr(fitY)
  fitY$vcov <- add_vcov(fitY)
  fitY$outcome <- outcome
  
  if (EOF) {
    if (family == "normal") {
      fitY$type <- 'continuous'
    } else {
      fitY$type <- 'binomial'
    }
  } else {
    fitY$type <- 'survival'
  }
  
  return(fitY)
}

#' Fit Competing Event Model
#'
#' This internal function fits a generalized linear model (GLM) for the competing event variable using the observed data.
#'
#' @param data                   A data frame, list or environment (or object coercible by
#'                               as.data.frame to a data frame) containing the variables in the model.
#'                               If not found in data, the variables are taken from environment(formula), 
#'                               typically the environment from which glm is called.
#' @param model                  Model statement for competing event variable.
#'
#' @return                       Fitted model for the competing event variable.
#' @keywords internal

fit_D_model <- function(model, data) {
  if (!model$family == "binomial") {
    stop("Currently only binomial models for competing events are allowed")
  }
  
  if (!is.null(model$subset)) {
    # Check for restrictions on compevent event variable modeling
    # Fit GLM assuming binomial distribution
    # Restrict data used for modeling to rows where condition is true
    fitD <- stats::glm(
      formula = model$formula,
      family = stats::binomial(link = model$link),
      data = subset(data, eval(parse(text = model$subset))),
      y = TRUE
    )
  } else {
    # Case where there are no restrictions on competing event variable
    #Fit GLM assuming binomial distribution
    fitD <- stats::glm(
      formula = model$formula,
      family = stats::binomial(link = model$link),
      data = data,
      y = TRUE
    )
  }
  
  fitD$rmse <- add_rmse(fitD)
  fitD$stderr <- add_stderr(fitD)
  fitD$vcov <- add_vcov(fitD)
  fitD$outcome <- all.vars(update(model$formula, . ~ 1))
  fitD$type <- 'survival'
  return(fitD)
}

