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
gformula_mc <- function(Y.fit,
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
  
  # initialize sim dataset
  sim <- data[data[[time]] == mc.start, ]
  
  n.start <- nrow(sim)
  n.start.unique <- n.start
  
  if (mc.start == mc.stop) {
    mc.sims <- 1
  }
  
  # If number of simulations PER unit is greater than 1 
  if (mc.sims > 1) {
    # make copies
    sim <- sim[rep(1:n.start, each = mc.sims), ]
    sim$sim <- rep(1:mc.sims, n.start)
    sim$group <- rep(1:n.start, each = mc.sims)
    
    #n.unique.ids <- n.unique.ids * mc.sims
  } else {
    sim <- data[data[[time]] == mc.start, ]
    sim$sim <- 1
    sim$group <- rep(1:n.start, each = mc.sims)
  }
  
  # make sure data is sorted by id and time
  sim <- sim[order(sim[[id]], sim[['sim']], sim[[time]]), ]
  sim <- sim[, c(id, 'sim', 'group', time, covs, base.covs, tv.covs, hist.vars)]
  
  X.rmses <- lapply(X.fit, get, x = "rmse")
  X.types <- lapply(X.fit, get, x = "type")
  
  X.restrictions <- lapply(X.fit, function (object) {
    if (exists(object, x = "restrict")) {
      get(object, x = "restrict")
    } else
      NULL
  })
  
  for (t in mc.start:mc.stop) {
    ind <- which(mc.start:mc.stop %in% t)
    if (t == mc.start) {
      # Set simulated covariate values at time t = 0 equal to observed covariate values
      sims <- list(sim)
      n.start.sims <- nrow(sim)
    } else {
      # Set initial simulated values at time t to simulated values at time t - 1, to be
      # updated later
      sim <- sims[[ind - 1]]
      
      sim[[time]] <- t
      
      sim[, hist.vars] <- update_histories(sim, hist.vars, hist.fun)
      
      if (!is.null(tv.covs)) {
        warning("tv.covs not implemented")
        #sim[, tv.covs] <- sims[sims[[time]] == t, tv.covs]
      }
      
      if (ice) {
        X.fit <- X.fit.ice[[t]]
      } 
      
      for (i in seq_along(X.fit)) {
        sim[[covs[i]]] <-
          switch(
            X.types[[i]],
            'binomial' = draw_binomial(
              N = n.start.sims, 
              size = 1, 
              p.fit = X.fit[[i]], 
              data = sim
            ),
            'normal' = draw_normal(
              N = n.start.sims, 
              mu.fit = X.fit[[i]],
              sd.hat = X.rmses[[i]], 
              data = sim
            ),
            'categorical' = draw_multinomial(
              p.fit = X.fit[[i]],
              data = sim
            ),
            'zero-inflated normal' = draw_zero_inflated_normal(
              N = n.start.sims,
              size = 1, 
              p.fit = X.fit[[i]][[1]],
              mu.fit = X.fit[[i]][[2]],
              sd.hat = X.rmses[[i]],
              data = sim
            ),
            'bounded normal' = draw_bounded_normal(
              N = n.start.sims, 
              mu.fit = X.fit[[i]],
              sd.hat = X.rmses[[i]], 
              range = c(min(sim[[covs[i]]]), max(sim[[covs[i]]])),
              data = sim
            ),
            'truncated normal' = draw_truncated_normal(
              N = n.start.sims, 
              mu.fit = X.fit[[i]],
              sd.hat = X.rmses[[i]], 
              direction = X.fit[[i]]$direction,
              point = X.fit[[i]]$point,
              data = sim
            ),
            'binomial grf' = draw_grf_binomial(
              N = n.start.sims, 
              size = 1, 
              p.fit = X.fit[[i]], 
              data = sim[, X.fit[[i]]$covs]
            ),
            'normal grf' = draw_grf_normal(
              N = n.start.sims, 
              mu.fit = X.fit[[i]],
              sd.hat = X.rmses[[i]], 
              data = sim[, X.fit[[i]]$covs]
            ),
            'custom' = draw_custom(
              N = n.start.sims,
              mu.fit = X.fit[[i]],
              sd.hat = X.rmses[[i]], 
              data = sim,
              pred.fun = custom.pred.fun
            )
          )
        
        if (bound.sims) {
          # Set simulated covariate values outside the observed range to the observed min / max
          if (X.types[[i]] %in% c('normal', 'bounded normal', 'truncated normal')) {
            
            if (min(sim[[covs[i]]]) < X.range[[i]][1]) {
              sim[[covs[i]]][sim[[covs[i]]] < X.range[[i]][1]] <- X.range[[i]][1]
            }
            if (max(sim[[covs[i]]]) > X.range[[i]][2]) {
              sim[[covs[i]]][sim[[covs[i]]] > X.range[[i]][2]] <- X.range[[i]][2]
            }
            
          } else if (X.types[[i]] == 'zero-inflated normal') {
            
            if (min(sim[[covs[i]]][sim[[covs[i]]] > 0]) < X.range[[i]][1]) {
              sim[[covs[i]]][sim[[covs[i]]] < X.range[[i]][1] & sim[[covs[i]]] > 0] <- X.range[[i]][1]
            }
            if (max(sim[[covs[i]]]) > X.range[[i]][2]) {
              sim[[covs[i]]][sim[[covs[i]]] > X.range[[i]][2]] <- X.range[[i]][2]
            }
            
          }
        }
        # restrictions
        if (!is.null(X.restrictions[[i]])) {
          rows <- with(sim, !eval(parse(text = X.restrictions[[i]]$subset)))
          sim[[covs[i]]][rows] <- X.restrictions[[i]]$otherwise
        }
        
      }
      
      
      # interventions
      if (!is.null(intervention) & !is.null(treatment)) {
        sim[,treatment] <- intervention(sim, sim[,treatment])
      }
    }
    
    if (ice) {
      Y.fit <- Y.fit.ice[[t + 1]]
      if (!is.null(D.fit)) {
        D.fit <- D.fit.ice[[t + 1]]
      }
    }
    
    # Predict outcome value at time t using parametric models
    if (Y.fit$type == 'survival') {
      if (any(class(Y.fit) == "grf")) {
        if (is.factor(Y.fit$Y.orig)) {
          sim$Py <- stats::predict(Y.fit, type = 'response', newdata = sim[, Y.fit$covs])$predictions[,2]
        } else {
          sim$Py <- stats::predict(Y.fit, type = 'response', newdata = sim[, Y.fit$covs])$predictions
        }
      } else {
        sim$Py <- stats::predict(Y.fit, type = 'response', newdata = sim)
      }
      
    } else if (Y.fit$type == 'continuous') {
      if (t < (n.unique.times - 1)) {
        sim$Ey <- NA
      } else if (t == (n.unique.times - 1)) {
        sim$Ey <- stats::predict(Y.fit, type = 'response', newdata = sim)
      }
    } else if (Y.fit$type == 'binomial') {
      if (t < (n.unique.times - 1)) {
        sim$Py <- NA
      } else if (t == (n.unique.times - 1)) {
        sim$Py <- stats::predict(Y.fit, type = 'response', newdata = sim)
      }
    }

    # Simulate outcome variable
    if (Y.fit$type == 'survival') {
      sim$Y <- stats::rbinom(n.start.sims, 1, sim$Py)

    }
    if (bound.sims) {
      # Set simulated outcome values outside the observed range to the observed min / max
      if (min(sim$Y) < Y.range[1]) {
        sim$Y[sim$Y < Y.range[1]] <- Y.range[1]
      }
      if (max(sim$Y) > Y.range[2]) {
        sim$Y[sim$Y > Y.range[2]] <- Y.range[2]
      }
    }
    if (Y.fit$type == 'survival') {
      if (!is.null(D.fit)) {
        # Predict competing event probabilities
        if (any(class(D.fit) == "grf")) {
          if (is.factor(D.fit$Y.orig)) {
            sim$Pd <- stats::predict(D.fit, type = 'response', newdata = sim[, D.fit$covs])$predictions[,2]
          } else {
            sim$Pd <- stats::predict(D.fit, type = 'response', newdata = sim[, D.fit$covs])$predictions
          }
        } else {
          sim$Pd <- stats::predict(D.fit, type = 'response', newdata = sim)
        }
        sim$D <- stats::rbinom(n.start.sims, 1, sim$Pd)
        
        if (bound.sims) {
          
          # Set simulated competing event values outside the observed range to the observed min / max
          if (min(sim$D) < D.range[1]) {
            sim$D[sim$D < D.range[1]] <- D.range[1]
          }
          if (max(sim$D) > D.range[2]) {
            sim$D[sim$D > D.range[2]] <- D.range[2]
          }
        }
      } else {
        sim$D <- 0
        sim$Pd <- 0
      }
      
      if (t == mc.start) {
        # Calculate probability of Y rather than D at time t
        sim$prodp1 <- sim$Py * (1 - sim$Pd)
        
        # Calculate probability of survival or D = 1 
        sim$prodp0 <- 1 - sim$Py
        
        sim$prodd0 <- 1 - sim$Pd
        
        sim$poprisk <- sim$prodp1
    
      } else {
        
        # Calculate prodp1 as product of previous values
        sim$prodp1 <- 
          sim$Py * (1 - sim$Pd) * sims[[ind - 1]]$prodp0 * sims[[ind - 1]]$prodd0
        
        sim$poprisk <- sims[[ind - 1]]$poprisk + sim$prodp1
        
        # Calculate probability of survival or D = 1 
        sim$prodp0 <- (1 - sim$Py) * sims[[ind - 1]]$prodp0
        
        # Calculate probability of survival from D
        sim$prodd0 <- (1 - sim$Pd) * sims[[ind - 1]]$prodd0
        
      }
      # If D occurs, flip a coin to determine which happens first Y or D
      Y.first <- stats::rbinom(length(sim$Y[sim$D == 1]), 1, sim$prodp1)
      
      # Update Y based on result of coin flip
      sim$Y[sim$D == 1] <- ifelse(Y.first, sim$Y[sim$D == 1], NA)
      
    }
    
    if (t == mc.start) {
      if (Y.fit$type == 'survival') {
        sims[[ind]][, c('Pd', 'D', 'Py', 'Y', 'prodp0', 'prodp1', 'prodd0', 'poprisk')] <-
          sim[, c('Pd', 'D', 'Py', 'Y', 'prodp0', 'prodp1', 'prodd0', 'poprisk')]
      } else if (Y.fit$type == 'continous') {
        sims[[ind]][['Ey']] <- 
          sim[['Ey']]
      } else if (Y.fit$type == 'binomial') {
        sims[[ind]][['Py']] <- 
          sim[['Ey']]
      }
    } else {
    
      # Add simulated data for time t to aggregate simulated data over time
      sims[[ind]] <- sim
    }
    # sims <- sims[order(sims[[id]], sims[['sim']], sims[[time]]), ]
    
  }

  # sims <- dplyr::bind_rows(sims)
  # 
  # if (Y.fit$type == 'survival') {
  #   sims$survival <- sims$prodp0
  # }
  
  avgs <- 
    sapply(sims, function(x) rowsum(x$poprisk, x$group) / tabulate(x$group))
  ids <- sims[[1]][(1:n.start - 1) * mc.sims + 1, ][[id]]
  
  predictions <- data.frame(pred = as.vector(avgs))
  predictions[[id]] <- ids
  predictions[[time]] <- rep(mc.start:mc.stop, each = n.start.unique)
  #print(predictions)
  # sims$group <- rep(1:(n.start * (mc.stop - mc.start + 1)), each = mc.sims)
  # 
  # predictions <- rowsum(sims$poprisk, sims$group) / tabulate(sims$group)
  # predictions <- data.frame(pred = predictions)
  # predictions[, c(id, time)] <- sims[(1:max(sims$group) - 1) * mc.sims + 1, c(id, time)]
  # 
  
  predictions <- predictions[order(predictions[[id]], predictions[[time]]), ]
  row.names(predictions) <- NULL
  
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

draw_grf_binomial <- function(N, size, p.fit, data) {
  if (is.factor(p.fit$Y.orig)) {
    stats::rbinom(
      n = N, 
      size = size, 
      prob = stats::predict(p.fit, type = 'response', newdata = data)$predictions[,2]
    )
  } else {
    stats::rbinom(
      n = N, 
      size = size, 
      prob = stats::predict(p.fit, type = 'response', newdata = data)$predictions
    )
  }
}

draw_grf_normal <- function(N, mu.fit, sd.hat, data) {
  stats::rnorm(
    n = N, 
    mean = stats::predict(mu.fit, type = 'response', newdata = data)$predictions, 
    sd = sd.hat
  )
}


draw_binomial <- function(N, size, p.fit, data) {
  stats::rbinom(
    n = N, 
    size = size, 
    prob = stats::predict(p.fit, type = 'response', newdata = data)
  )
}

draw_normal <- function(N, mu.fit, sd.hat, data) {
  stats::rnorm(
    n = N, 
    mean = stats::predict(mu.fit, type = 'response', newdata = data), 
    sd = sd.hat
  )
}

draw_multinomial <- function(p.fit, data) {
  stats::predict(p.fit, type = 'class', newdata = data)
}

draw_zero_inflated_normal <- function(N, size, p.fit, mu.fit, sd.hat, data) {
  nonzero <- draw_binomial(
    N = N,
    size = size, 
    p.fit = p.fit,
    data = data
  )
  
  response <- draw_normal(
    N = N,
    mu.fit = mu.fit,
    sd.hat = sd.hat,
    data = data
  )
  
  nonzero * response
}

draw_bounded_normal <- function(N, mu.fit, sd.hat, range, data) {
  response <- draw_normal(
    N = N,
    mu.fit = mu.fit,
    sd.hat = sd.hat,
    data = data
  )
  
  (response * (range[2] - range[1])) + range[1]
}

draw_trunc_normal <- function(N, mu.fit, sd.hat, direction, point, data) {
  if (direction == 'left') {
    a <- point
    b <- Inf
  } else if (direction == 'right') {
    a <- -Inf
    b <- point
  }
  
  truncnorm::rtruncnorm(
    x = N,
    mean = predict(mu.fit, type = 'response', newdata = data),
    sd = sd.hat,
    a = a,
    b = b
  )
}

draw_custom <- function(N, mu.fit, sd.hat, data, pred.fun = stats::predict()) {
  stats::rnorm(
    n = N, 
    mean = stats::predict(mu.fit, type = 'response', newdata = data), 
    sd = rmse.hat
  )
}

update_histories <- function(data, hist.vars, hist.fun) {
  switch(
    hist.fun,
    "lag" = data[, gsub("^lag3_", "lag2_", gsub("^lag2_", "lag1_", gsub("^lag1_", "", hist.vars)))]
  )
}


