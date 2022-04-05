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

fit_X_models <- function(models, data, time = "time", t0 = 0, gam = FALSE) {
  # Create list of models for covariates
  cov.names <- names(models)
  subdata <- data[data[[time]] > t0, ]
  
  X.fits <- lapply(
    models,
    function(model) {
      if (!is.null(model$subset)) {
        subdata <- subset(subdata, eval(parse(text = model$subset)))
      }
      if (!is.null(model$restrict)) {
          subdata <- subset(subdata, eval(parse(text = model$restrict$subset)))
      }
      
      if (gam) {
        fit <- switch(
          model$family,
          'binomial' = mgcv::gam(
            formula = model$formula, 
            family = binomial(link = model$link), 
            data = subdata
          ),
          
          'normal' = mgcv::gam(
            formula = model$formula, 
            family = gaussian,
            data = subdata
          )
        )
        fit$rmse <- add_rmse(fit)
        fit$stderrs <- add_stderr(fit)
        fit$vcov <- add_vcov(fit)
        fit$outcome <- all.vars(update(model$formula, . ~ 1))
        if (model$family == 'normal') {
          fit$type <- 'normal'
        } else if (model$family == 'binomial') {
          fit$type <- 'binomial'
        }
          
      } else {
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
      }
      
      
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

fit_Y_model <- function(model, data, EOF = FALSE, time = NULL, gam = FALSE) {
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
  
  if (gam) {
    fit_function <- mgcv::gam
  } else {
    fit_function <- stats::glm
  }
  
  if (!is.null(model$subset)) {
    
    if (model$family == "normal") {
      data[, paste("norm_", outcome, sep = "")] <-
        (data[[outcome]] - min(data[[outcome]]))/(max(data[[outcome]]) - min(data[[outcome]]))
      
      fitY <- fit_function(
        stats::as.formula(paste("norm_", model, sep = "")), 
        family = family,
        data = subset(data, eval(parse(text = model$subset)))
      )
    } else {
      fitY <- fit_function(
        formula = model$formula,
        family = family,
        data = subset(data, eval(parse(text = model$subset)))
      )
    }
  } else { # Case where there are no restrictions on outcome variable
    # Fit GLM for outcome variable using user-specified model and entire dataset
    fitY <- fit_function(
      formula = model$formula,
      family = family,
      data = data
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

fit_D_model <- function(model, data, gam = FALSE) {
  if (!model$family == "binomial") {
    stop("Currently only binomial models for competing events are allowed")
  }
  
  if (gam) {
    fit_function <- mgcv::gam
  } else {
    fit_function <- stats::glm
  }
  
  if (!is.null(model$subset)) {
    # Check for restrictions on compevent event variable modeling
    # Fit GLM assuming binomial distribution
    # Restrict data used for modeling to rows where condition is true
    fitD <- fit_function(
      formula = model$formula,
      family = stats::binomial(link = model$link),
      data = subset(data, eval(parse(text = model$subset)))
    )
  } else {
    # Case where there are no restrictions on competing event variable
    #Fit GLM assuming binomial distribution
    fitD <- fit_function(
      formula = model$formula,
      family = stats::binomial(link = model$link),
      data = data
    )
  }
  
  fitD$rmse <- add_rmse(fitD)
  fitD$stderr <- add_stderr(fitD)
  fitD$vcov <- add_vcov(fitD)
  fitD$outcome <- all.vars(update(model$formula, . ~ 1))
  fitD$type <- 'survival'
  return(fitD)
}

define_models <- function(spec, outcome, compevent, time = TRUE, logvars = FALSE, sexints = FALSE) {
  if (time) {
    tvar <- "as.factor(time)"
  } else {
    tvar <- NULL
  }
  
  X.spec <- spec[!names(spec) %in% c(outcome, compevent)]
  Y.spec <- spec[[outcome]]
  D.spec <- spec[[compevent]]
  
  if (logvars) {
    # savenames <- names(X.spec)
    # X.spec <- lapply(X.spec, function(x) {
    #   str_replace_all(
    #     x,
    #     "((?:lag[12]_)?age)|((?:lag[12]_)?tc)|((?:lag[12]_)?hdl)|((?:lag[12]_)?ldl)|((?:lag[12]_)?sbp)",
    #     "log(\\1\\2\\3\\4\\5)"
    #   )
    # })
    # names(X.spec) <- savenames
    
    savenames <- names(Y.spec)
    Y.spec <- str_replace_all(
      Y.spec,
      #"((?:lag[12]_)?age)|((?:lag[12]_)?tc)|((?:lag[12]_)?hdl)|((?:lag[12]_)?ldl)|((?:lag[12]_)?sbp)",
      "((?<!lag[12]_)age)|((?<!lag[12]_)tc)|((?<!lag[12]_)hdl)|((?<!lag[12]_)ldl)|((?<!lag[12]_)sbp)",
      "log(\\1\\2\\3\\4\\5)"
    )
    names(Y.spec) <- savenames
    
    savenames <- names(D.spec)
    D.spec <- str_replace_all(
      D.spec,
      #"((?:lag[12]_)?age)|((?:lag[12]_)?tc)|((?:lag[12]_)?hdl)|((?:lag[12]_)?ldl)|((?:lag[12]_)?sbp)",
      "((?<!lag[12]_)age)|((?<!lag[12]_)tc)|((?<!lag[12]_)hdl)|((?<!lag[12]_)ldl)|((?<!lag[12]_)sbp)",
      "log(\\1\\2\\3\\4\\5)"
    )
    names(D.spec) <- savenames
  }

  X.models <- lapply(
    X = names(X.spec), 
    FUN = function(name) {
      if (sexints) {
        X.formula <- reformulate(
          c(str_remove(X.spec[[name]], "sex"), tvar),
          name
        )
        X.formula <- update(X.formula, . ~ sex * .)
        
      } else {
        X.formula <- reformulate(
          c(X.spec[[name]], tvar),
          name
        )
      }
      switch(
        name,
        # simple model for age
        "age" = list(
          formula = X.formula,
          link = "identity",
          family = "normal"
        ),
        
        # logit model for probability of smoking
        "smk" = list(
          formula = X.formula, 
          link = "logit",
          family = "binomial"
        ),
        
        # linear model of BMI
        "bmi" = list(
          formula = X.formula, 
          link = "identity",
          family = "normal"
        ),
        
        # logit model for diabetes (failure)
        "dm" = list(
          formula = X.formula, 
          link = "logit",
          family = "binomial",
          restrict = list(
            subset = 'lag1_dm == 0',
            otherwise = 1
          )
        ),
        
        # logit model for hypertension meds
        "hrx" = list(
          formula = X.formula, 
          link = "logit",
          family = "binomial"
        ),
        
        # logit model for lipids meds
        "liprx" = list(
          formula = X.formula, 
          link = "logit",
          family = "binomial"
        ),
        
        # linear model for total cholesterol
        "tc" = list(
          formula = X.formula, 
          link = "log",
          family = "normal"
        ),
        
        # linear model for total cholesterol
        "ldl" = list(
          formula = X.formula, 
          link = "log",
          family = "normal"
        ),
        
        # linear model for HDL cholesterol
        "hdl" = list(
          formula = X.formula, 
          link = "log",
          family = "normal"
        ),
        
        # linear model for systolic blood pressure
        "sbp" = list(
          formula = X.formula, 
          link = "log",
          family = "normal"
        )
      )
    })

  names(X.models) <- names(X.spec)
  
  
  if (sexints) {
    Y.formula <- reformulate(
      c(str_remove(Y.spec, "sex"), tvar),
      outcome
    )
    Y.formula <- update(Y.formula, . ~ sex * .)
    
    D.formula <- reformulate(
      c(str_remove(D.spec, "sex"), tvar),
      compevent
    )
    D.formula <- update(D.formula, . ~ sex * .)
  } else {
    
    Y.formula <- reformulate(
      c(Y.spec, tvar),
      outcome
    )
    
    D.formula <- reformulate(
      c(D.spec, tvar),
      compevent
    )
  }
  
  # pooled logistic model of CVD events
  Y.model <- list(
    formula = Y.formula, 
    link = "logit",
    family = "binomial"
  )
  
  # pooled logistic model of death due to other causes
  D.model <- list(
    formula = D.formula, 
    link = "logit",
    family = "binomial"
  )
  
  ret <- list(
    "X" = X.models,
    "Y" = Y.model,
    "D" = D.model
  )
  
  return(ret)
}



#' Fit Covariate Models (ICE)
#'
#' This internal function fits a model for each covariate using the observed data.
#'
#' @param models 
#' @param data 
#'
#' @return                A list of fitted models, one for each covariate in \code{covnames}.
#' @keywords internal
#' @import data.table

fit_X_models_ice <- function(models, data, time = "time", t0 = 0, k) {
  
  # here we can skip first time step for covariate models
  X.fits <- lapply(2:k, function(i) {
    mod <- lapply(models, function(model) {
      if (i != k) {
        f <- model$formula
        pat <- paste0("^lag[", paste0(i:(k-1), collapse = ""), "]")
        
        if (any(str_detect(attr(terms(f), "term.labels"), pat))) {
          model$formula <- update(
            f, 
            drop.terms(
              terms(f),
              str_which(attr(terms(f), "term.labels"), pat),
              keep.response = TRUE
            )
          )
        }
      } 
      model$subset = paste0("time == ", i - 1)
      
      return(model)
    })
    
    fit_X_models(
      model = mod, 
      data = data, 
      time = time,
      t0 = t0
    )
  })
    
    return(X.fits)
}

#' Fit Outcome Model (ICE)
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

fit_Y_model_ice <- function(model, data, EOF = FALSE, time = NULL, k) {

  fitY <- lapply(1:k, function(i) {
    if (i != k) {
      f <- model$formula
      pat <- paste0("^lag[", paste0(i:(k-1), collapse = ""), "]")
      
      if (any(str_detect(attr(terms(f), "term.labels"), pat))) {
        model$formula <- update(
          f, 
          drop.terms(
            terms(f),
            str_which(attr(terms(f), "term.labels"), pat),
            keep.response = TRUE
          )
        )
      }
    } 
    
    model$subset <- paste0("time == ", i - 1)
    
    fit_Y_model(
      model = model, 
      data = data, 
      EOF = EOF, 
      time = time
    )
  })
  
  return(fitY)
}

#' Fit Competing Event Model (ICE)
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

fit_D_model_ice <- function(model, data, time, k) {

  fitD <- lapply(1:k, function(i) {
    if (i != k) {
      f <- model$formula
      pat <- paste0("^lag[", paste0(i:(k-1), collapse = ""), "]")
      
      if (any(str_detect(attr(terms(f), "term.labels"), pat))) {
        model$formula <- update(
          f, 
          drop.terms(
            terms(f),
            str_which(attr(terms(f), "term.labels"), pat),
            keep.response = TRUE
          )
        )
      }
    } 
    
    model$subset <- paste0("time == ", i - 1)
    
    fit_D_model(
      model = model, 
      data = data
    )
  })
  
  return(fitD)
}
