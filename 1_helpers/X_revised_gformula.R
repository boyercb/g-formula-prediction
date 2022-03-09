#' Title
#'
#' @param Y.fit 
#' @param X.fit 
#' @param D.fit 
#' @param data 
#' @param id 
#' @param time 
#' @param base.covs 
#' @param mc.start 
#' @param mc.stop 
#' @param mc.samples 
#' @param mc.sims 
#' @param boot 
#' @param boot.reps 
#' @param interventions 
#' @param restrictions 
#' @param ci.method 
#' @param ci.level 
#' @param num.threads 
#' @param verbose 
#' @param seed 
#'
#' @return
#' @export
#'
#' @examples
gformula_mc2 <- function(Y.fit,
                        X.fit,
                        D.fit = NULL,
                        data,
                        id = "id",
                        time = "time",
                        treatment = NULL,
                        base.covs = NULL,
                        tv.covs = NULL,
                        hist.vars = NULL,
                        hist.fun = NULL,
                        ice = FALSE,
                        mc.start = 0,
                        mc.stop = NULL,
                        mc.samples = NULL,
                        mc.sims = 1,
                        bound.sims = TRUE,
                        boot = FALSE,
                        boot.reps = 500,
                        intervention = NULL,
                        ci.method = NULL,
                        ci.level = 95,
                        custom.pred.fun = stats::predict(),
                        num.threads = NULL,
                        verbose = FALSE,
                        merge = TRUE,
                        pool = FALSE,
                        last_only = FALSE,
                        seed = runif(1, 0, .Machine$integer.max)) {
  
  set.seed(seed)
  
  # local variables
  unique.ids <- unique(data[[id]])
  unique.times <- unique(data[[time]])
  n.unique.ids <- length(unique.ids)
  n.unique.times <- length(unique.times)
  
  obs <- nrow(data)
  orig <- data
  orig$order <- 1:obs
  
  if (ice) {
    Y.fit.ice <- Y.fit
    X.fit.ice <- X.fit
    
    Y.fit <- Y.fit[[1]]
    X.fit <- X.fit[[1]]
  }
  
  # collect variable names
  covs <- names(X.fit)
  outcome <- Y.fit$outcome
  
  if (!is.null(D.fit)) {
    if (ice) {
      D.fit.ice <- D.fit
      D.fit <- D.fit[[1]]
    } 
    comp.risk <- D.fit$outcome
  } 
  
  if (is.null(mc.stop)) {
    mc.stop <- n.unique.times - 1
  }
  
  # not specified set samples equal to number of individuals in data set
  if (is.null(mc.samples)) {
    mc.samples <- n.unique.ids
  }
  
  # make sure data is sorted by id and time
  data <- data[order(data[[id]], data[[time]]), ]
  
  # Determine ranges of observed covariates and outcome
  X.range <- lapply(seq_along(covs), function(x) range(data[[covs[x]]], na.rm = TRUE))
  Y.range <- range(data[[outcome]], na.rm = TRUE)
  
  if (!is.null(D.fit)) {
    D.range <- range(data[[comp.risk]])
  }
  
  # If the number of desired simulation units is greater than the number of
  # unique units sample with replacement
  if (mc.samples > n.unique.ids) {
    samples <- sample(unique.ids, mc.samples, replace = TRUE)
    data <- do.call(rbind, lapply(samples, function(x) data[data[[id]] == x, ]))
  }
  
  # If number of simulations PER unit is greater than 1 
  if (mc.sims > 1) {
    data <- data[rep(1:obs, each = mc.sims), ]
    data$sim <- rep(1:mc.sims, obs)
    data$group <- rep(1:obs, each = mc.sims)
    
    # make sure data is sorted by id and time
    data <- data[order(data[[id]], data[['sim']], data[[time]]), ]
    
    #n.unique.ids <- n.unique.ids * mc.sims
  } else {
    data$sim <- 1
  }
  
  X.rmses <- lapply(X.fit, get, x = "rmse")
  X.types <- lapply(X.fit, get, x = "type")
  
  X.restrictions <- lapply(X.fit, function (object) {
    if (exists(object, x = "restrict")) {
      get(object, x = "restrict")
    } else
      NULL
  })
  
  for (t in mc.start:mc.stop) {
    
    if (t == mc.start) {
      # Set simulated covariate values at time t = 0 equal to observed covariate values
      sims <- data[data[[time]], c(id, 'sim', 'group', time, covs, base.covs, tv.covs, hist.vars)]
      rows <- which(sims[[time]] == t)

      n.start <- nrow(sims[sims[[time]] == t, ])
      n.start.unique <- length(unique(sims[[id]][sims[[time]] == t]))
    } else {
      # Set initial simulated values at time t to simulated values at time t - 1, to be
      # updated later
      rows <- which(sims[[time]] == t)
      
      sims[rows, ] <- sims[sims[[time]] == t - 1, ]
      
      sims[rows, time] <- t
      
      sims[rows, hist.vars] <- update_histories(sims[rows, ], hist.vars, hist.fun)
      
      if (!is.null(tv.covs)) {
        sims[rows, tv.covs] <- data[data[[time]] == t, tv.covs]
      }
      
      if (ice) {
        X.fit <- X.fit.ice[[t]]
      } 
      print("b")
      for (i in seq_along(X.fit)) {
        
        sims[rows, covs[i]] <-
          switch(
            X.types[[i]],
            'binomial' = draw_binomial(
              N = n.start, 
              size = 1, 
              p.fit = X.fit[[i]], 
              data = sims[rows, ]
            ),
            'normal' = draw_normal(
              N = n.start, 
              mu.fit = X.fit[[i]],
              sd.hat = X.rmses[[i]], 
              data = sims[rows, ]
            ),
            'categorical' = draw_multinomial(
              p.fit = X.fit[[i]],
              data = sims[rows, ]
            ),
            'zero-inflated normal' = draw_zero_inflated_normal(
              N = n.start,
              size = 1, 
              p.fit = X.fit[[i]][[1]],
              mu.fit = X.fit[[i]][[2]],
              sd.hat = X.rmses[[i]],
              data = sims[rows, ]
            ),
            'bounded normal' = draw_bounded_normal(
              N = n.start, 
              mu.fit = X.fit[[i]],
              sd.hat = X.rmses[[i]], 
              range = c(min(sims[rows, covs[i]]), max(sims[rows, covs[i]])),
              data = sims[rows, ]
            ),
            'truncated normal' = draw_truncated_normal(
              N = n.start, 
              mu.fit = X.fit[[i]],
              sd.hat = X.rmses[[i]], 
              direction = X.fit[[i]]$direction,
              point = X.fit[[i]]$point,
              data = sims[rows, ]
            ),
            'binomial grf' = draw_grf_binomial(
              N = n.start, 
              size = 1, 
              p.fit = X.fit[[i]], 
              data = sims[rows, X.fit[[i]]$covs]
            ),
            'normal grf' = draw_grf_normal(
              N = n.start, 
              mu.fit = X.fit[[i]],
              sd.hat = X.rmses[[i]], 
              data = sims[rows, X.fit[[i]]$covs]
            ),
            'custom' = draw_custom(
              N = n.start,
              mu.fit = X.fit[[i]],
              sd.hat = X.rmses[[i]], 
              data = sims[rows, ],
              pred.fun = custom.pred.fun
            )
          )
        
        print("c")
        if (bound.sims) {
          # Set simulated covariate values outside the observed range to the observed min / max
          if (X.types[[i]] %in% c('normal', 'bounded normal', 'truncated normal')) {
            
            if (min(sims[rows, covs[i]]) < X.range[[i]][1]) {
              sims[rows, covs[i]][sims[rows, covs[i]] < X.range[[i]][1]] <- X.range[[i]][1]
            }
            if (max(sims[rows, covs[i]]) > X.range[[i]][2]) {
              sims[rows, covs[i]][sims[rows, covs[i]] > X.range[[i]][2]] <- X.range[[i]][2]
            }
            
          } else if (X.types[[i]] == 'zero-inflated normal') {
            
            if (min(sims[rows, covs[i]][sims[rows, covs[i]] > 0]) < X.range[[i]][1]) {
              sims[rows, covs[i]][sims[rows, covs[i]] < X.range[[i]][1] & sims[rows, covs[i]] > 0] <- X.range[[i]][1]
            }
            if (max(sims[rows, covs[i]]) > X.range[[i]][2]) {
              sims[rows, covs[i]][sims[rows, covs[i]] > X.range[[i]][2]] <- X.range[[i]][2]
            }
            
          }
        }
        # restrictions
        if (!is.null(X.restrictions[[i]])) {
          subset_rows <- with(sims[rows, ], !eval(parse(text = X.restrictions[[i]]$subset)))
          sims[rows, covs[i]][subset_rows] <- X.restrictions[[i]]$otherwise
        }
        
      }
      
      # interventions
      if (!is.null(intervention) & !is.null(treatment)) {
        sims[rows, treatment] <- intervention(sims[rows, ], sims[rows, treatment])
      }
    }
    
    if (ice) {
      Y.fit <- Y.fit.ice[[t + 1]]
      if (!is.null(D.fit)) {
        D.fit <- D.fit.ice[[t + 1]]
      }
    }
    
    if (t == mc.start) {
      if (Y.fit$type == 'survival') {
        sims[, c('Pd', 'D', 'Py', 'Y', 'prodd0', 'prodp0', 'prodp1', 'poprisk')] <- NA
      } else if (Y.fit$type == 'continous') {
        sims[['Ey']] <- NA
      } else if (Y.fit$type == 'binomial') {
        sims[['Py']] <- NA
      }
    }
      
    # Predict outcome value at time t using parametric models
    if (Y.fit$type == 'survival') {
      if (any(class(Y.fit) == "grf")) {
        if (is.factor(Y.fit$Y.orig)) {
          sims$Py[rows] <- stats::predict(Y.fit, type = 'response', newdata = sims[rows, Y.fit$covs])$predictions[,2]
        } else {
          sims$Py[rows] <- stats::predict(Y.fit, type = 'response', newdata = sims[rows, Y.fit$covs])$predictions
        }
      } else {
        sims$Py[rows] <- stats::predict(Y.fit, type = 'response', newdata = sims[rows, ])
      }
      
    } else if (Y.fit$type == 'continuous') {
      if (t < (n.unique.times - 1)) {
        sims$Ey[rows] <- NA
      } else if (t == (n.unique.times - 1)) {
        sims$Ey[rows] <- stats::predict(Y.fit, type = 'response', newdata = sims[rows, ])
      }
    } else if (Y.fit$type == 'binomial') {
      if (t < (n.unique.times - 1)) {
        sims$Py[rows] <- NA
      } else if (t == (n.unique.times - 1)) {
        sims$Py[rows] <- stats::predict(Y.fit, type = 'response', newdata = sims[rows, ])
      }
    }
    
    # Simulate outcome variable
    if (Y.fit$type == 'survival') {
      sims$Y[rows] <- stats::rbinom(n.start, 1, sims$Py[rows])
    }
    
    if (bound.sims) {
      # Set simulated outcome values outside the observed range to the observed min / max
      if (min(sims$Y[rows]) < Y.range[1]) {
        sims$Y[rows][sims$Y[rows] < Y.range[1]] <- Y.range[1]
      }
      if (max(sims$Y[rows]) > Y.range[2]) {
        sims$Y[rows][sims$Y[rows] > Y.range[2]] <- Y.range[2]
      }
    }
    if (Y.fit$type == 'survival') {
      if (!is.null(D.fit)) {
        # Predict competing event probabilities
        if (any(class(D.fit) == "grf")) {
          if (is.factor(D.fit$Y.orig)) {
            sims$Pd[rows] <- stats::predict(D.fit, type = 'response', newdata = sims[rpws, D.fit$covs])$predictions[,2]
          } else {
            sims$Pd[rows] <- stats::predict(D.fit, type = 'response', newdata = sims[rows, D.fit$covs])$predictions
          }
        } else {
          sims$Pd[rows] <- stats::predict(D.fit, type = 'response', newdata = sims[rows, ])
        }
        sims$D[rows] <- stats::rbinom(n.start, 1, sims$Pd)
        
        if (bound.sims) {
          
          # Set simulated competing event values outside the observed range to the observed min / max
          if (min(sims$D[rows]) < D.range[1]) {
            sims$D[rows][sims$D[rows] < D.range[1]] <- D.range[1]
          }
          if (max(sims$D[rows]) > D.range[2]) {
            sims$D[rows][sims$D[rows] > D.range[2]] <- D.range[2]
          }
        }
      } else {
        sims$D[rows] <- 0
        sims$Pd[rows] <- 0
      }
      
      if (t == mc.start) {
        # Calculate probability of Y rather than D at time t
        sims$prodp1[rows] <- sims$Py[rows] * (1 - sims$Pd[rows])
        
        # Calculate probability of survival or D = 1 
        sims$prodp0[rows] <- 1 - sims$Py[rows]
        
        sims$prodd0[rows] <- 1 - sims$Pd[rows]
        
        sims$poprisk[rows] <- sims$prodp1[rows]
        
      } else {
        
        # Calculate prodp1 as product of previous values
        sims$prodp1[rows] <- 
          sims$Py[rows] * (1 - sims$Pd[rows]) * sims$prodp0[rows][sims[[time]] == t - 1] * 
            sims$prodd0[rows][sims[[time]] == t - 1]

        sims$poprisk[rows] <- sims$poprisk[rows][sims[[time]] == t - 1] + sims$prodp1[rows]
          
        # Calculate probability of survival or D = 1 
        sims$prodp0[rows] <- (1 - sims$Py[rows]) * sims$prodp0[rows][sims[[time]] == t - 1]
        sims$prodd0[rows] <- (1 - sims$Pd[rows]) * sims$prodd0[rows][sims[[time]] == t - 1]
      }
      # If D occurs, flip a coin to determine which happens first Y or D
      Y.first <- stats::rbinom(length(sims$Y[rows][sims$D[rows] == 1]), 1, sims$prodp1[rows])
      
      # Update Y based on result of coin flip
      sims$Y[rows][sims$D[rows] == 1] <- ifelse(Y.first, sims$Y[rows][sims$D[rows] == 1], NA)
      
    }
    
    # Add simulated data for time t to aggregate simulated data over time
    #sims <- rbind(sims[sims[[time]] < t, ], sim)
    sims <- sims[order(sims[[id]], sims[['sim']], sims[[time]]), ]
  }
  
  if (Y.fit$type == 'survival') {
    #sims$poprisk <- stats::ave(sims$prodp1, sims[[id]], sims[['sim']], FUN = cumsum)
    sims$survival <- sims$prodp0
  }
  
  sims <- sims[order(sims[[id]], sims[[time]], sims[['sim']]), ]
  #predictions <- t(tapply(sims$poprisk, list(sims[[time]], sims[[id]]), FUN = mean))
  predictions <- rowsum(sims$poprisk, sims$group) / tabulate(sims$group)
  predictions <- data.frame(pred = predictions)
  #predictions[[id]] <- as.numeric(row.names(predictions))
  predictions[, c(id, time)] <- sims[!duplicated(sims[, c(id, time)]), c(id, time)]

  # predictions <-
  #   reshape(
  #     predictions,
  #     varying = list(1:(mc.stop - mc.start + 1)),
  #     v.names = "pred",
  #     idvar = id,
  #     timevar = time,
  #     direction = "long"
  #   )
  
  # predictions <- predictions[order(predictions[[id]], predictions[[time]]), ]
  # predictions[[time]] <- predictions[[time]] + (mc.start - 1)
  # row.names(predictions) <- NULL
  if (merge == TRUE) {
    predictions <- merge(orig, predictions, by = c(id, time), all.x = TRUE, sort = TRUE)
    predictions <- predictions[order(predictions[['order']]), 'pred']
    names(predictions) <- NULL
  } 
  
  if (pool == TRUE) {
    predictions <-  merge(sims, predictions, by = c(id, time), all.x = TRUE, sort = TRUE)
  }
  
  if (last_only == TRUE) {
    predictions <- predictions[predictions[[time]] == mc.stop, ]
  }
  
  return(predictions)
  
}
