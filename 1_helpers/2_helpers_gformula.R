# take fitted models of joint distribution and simulate data for a given
# baseline covariate profile

predict.gformula <- function(object, obs_data, newdata = NULL, id, t0 = 0, covnames, covtypes, covparams,
                             covfits_custom = NA, covpredict_custom = NA,
                             histvars = NULL, histories = NA, basecovs = NA,
                             outcome_name, ymodel,
                             compevent_name = NULL, compevent_model = NA,
                             intvars = NULL, interventions = NULL,
                             int_times = NULL, int_descript = NULL, ref_int = 0, intcomp = NA,
                             visitprocess = NA, restrictions = NA,
                             yrestrictions = NA, compevent_restrictions = NA,
                             baselags = FALSE,
                             nsimul = NA, seed,
                             nsamples = 0, parallel = FALSE, ncores = NA,
                             ci_method = 'percentile', threads,
                             show_progress = TRUE, return_sims = FALSE, ...){
  if (object$time_points > 1){
    if (object$comprisk){
      fitD <- object$fits[[length(object$fits)]]
      fitY <- object$fits[[length(object$fits) - 1]]
      fitcov <- object$fits[1:(length(object$fits) - 2)]
    } else {
      fitY <- object$fits[[length(object$fits)]]
      fitcov <- object$fits[1:(length(object$fits) - 1)]
      fitD <- NA
    }
  } else {
    fitY <- unlist(object$fits)
    fitD <- NA
    fitcov <- NA
  }
  outcome_type <- "survival"
  max_visits <- NA
  
  lag_indicator <- lagavg_indicator <- cumavg_indicator <- c()
  lag_indicator <- gfoRmula:::update_lag_indicator(covparams$covmodels, lag_indicator)
  lagavg_indicator <- gfoRmula:::update_lagavg_indicator(covparams$covmodels, lagavg_indicator)
  cumavg_indicator <- gfoRmula:::update_cumavg_indicator(covparams$covmodels, cumavg_indicator)
  
  if (!missing(ymodel)){
    lag_indicator <- gfoRmula:::update_lag_indicator(ymodel, lag_indicator)
    lagavg_indicator <- gfoRmula:::update_lagavg_indicator(ymodel, lagavg_indicator)
    cumavg_indicator <- gfoRmula:::update_cumavg_indicator(ymodel, cumavg_indicator)
  }
  if (!(length(compevent_model) == 1 && is.na(compevent_model))){
    lag_indicator <- gfoRmula:::update_lag_indicator(compevent_model, lag_indicator)
    lagavg_indicator <- gfoRmula:::update_lagavg_indicator(compevent_model, lagavg_indicator)
    cumavg_indicator <- gfoRmula:::update_cumavg_indicator(compevent_model, cumavg_indicator)
  }
  histvals <- list(lag_indicator = lag_indicator, lagavg_indicator = lagavg_indicator,
                   cumavg_indicator = cumavg_indicator)
  
  comprisk <- object$comprisk
  
  if (!missing(threads)){
    data.table::setDTthreads(threads = threads)
  }
  else {
    threads <- data.table::getDTthreads()
  }
  
  min_time <- min(newdata[[object$time_name]])
  below_zero_indicator <- min_time < 0
  
  newdata <- data.table::copy(newdata)
  newdata$oldid <- newdata[[id]]

  obs_data <- data.table::copy(obs_data)
  
  sample_size <- length(unique(newdata[[id]]))
  
  time_points <- object$time_points

  
  for (i in seq_along(covnames)){
    if (covtypes[i] == 'absorbing'){
      restrictions <- c(restrictions[!is.na(restrictions)],
                        list(c(covnames[i], paste("lag1_", covnames[i], "==0", sep = ""),
                               carry_forward, 1)))
      covtypes[i] <- 'binary'
    }
  }
  
  # Create 1-indexed numerical IDs for observed datasets
  ids <- data.table::as.data.table(sort(unique(newdata[[id]])))
  ids[, 'newid' := seq_len(.N)]
  data.table::setkeyv(newdata, id)
  newdata <- newdata[J(ids), allow.cartesian = TRUE]
  newdata_geq_0 <- newdata[newdata[[object$time_name]] >= 0]
  obs_data_geq_0 <- obs_data[obs_data[[object$time_name]] >= 0]
  
  # Set default number of simulated individuals to equal number of individuals in
  # observed dataset
  if (is.na(nsimul)){
    nsimul <- length(unique(newdata$newid))
  }
  
  
  # Generate seeds for simulations and bootstrapping
  set.seed(seed)
  newseeds <- sample.int(2^30, size = nsamples + 1)
  subseed <- newseeds[1]
  bootseeds <- newseeds[2:(nsamples + 1)]
  
  # Determine ranges of observed covariates and outcome
  ranges <- lapply(seq_along(covnames), FUN = function(i){
    if (covtypes[i] == 'normal' || covtypes[i] == 'bounded normal' ||
        covtypes[i] == 'truncated normal') {
      range(obs_data_geq_0[[covnames[i]]])
    } else if (covtypes[i] == 'zero-inflated normal'){
      range(obs_data_geq_0[obs_data_geq_0[[covnames[i]]] > 0][[covnames[i]]])
    } else {
      NA
    }
  })
  
  yrange <- range(obs_data_geq_0[[outcome_name]])

  if (comprisk){
    compevent_range <- range(obs_data_geq_0[[compevent_name]])
  } else {
    compevent_range <- NA
  }
  
  newdata_noresample <- data.table::copy(newdata)
  len <- length(unique(newdata$newid))
  # If the number of user desired simulations differs from the number of individuals in
  # the observed dataset, sample the desired number of observed IDs with replacement
  # if (nsimul < len){
  #   ids <- data.table::as.data.table(sort(sample(unique(newdata$newid), nsimul, replace = TRUE)))
  #   colnames(ids) <- "newid"
  #   ids[, 'sid' := seq_len(.N)]
  #   newdata <- merge(ids, newdata, all.x = TRUE, by = "newid")
  #   newdata[, 'newid' := newdata$sid]
  #   newdata[, 'sid' := NULL]
  # } else if (nsimul > len){
    ids <- data.table::as.data.table(rep(unique(newdata$newid), each = nsimul))
    ids[, 'newid' := 1:(nsimul*len)]
    colnames(ids) <- c("newid", "sid")
    data.table::setkeyv(newdata, "newid")
    newdata <- newdata[J(ids), allow.cartesian = TRUE]
    newdata[, 'newid' := newdata$sid]
    newdata[, 'sid' := NULL]
  #}
  
  # Add natural course to list of interventions
  if (!is.null(interventions)){
    comb_interventions <- c(list(list(c(gfoRmula:::natural))), interventions)
    comb_intvars <- c(list('none'), intvars)
  } else {
    comb_interventions <- list(list(c(gfoRmula:::natural)))
    comb_intvars <- list('none')
  }
  
  if (is.null(int_times)){
    comb_int_times <- list()
    for (i in seq_along(comb_interventions)){
      comb_int_times[[i]] <- lapply(seq_along(comb_interventions[[i]]),
                                    FUN = function(i) {0:(time_points - 1)})
    }
  } else {
    comb_int_times <- c(list(list(0:(time_points - 1))), int_times)
  }
  
  if (parallel){
    cl <- gfoRmula:::prep_cluster(ncores = ncores, threads = threads, covtypes = covtypes)
    pools <- parallel::parLapply(cl, seq_along(comb_interventions), gformula.simulate,
                                 fitcov = fitcov, fitY = fitY, fitD = fitD,
                                 yrestrictions = yrestrictions,
                                 compevent_restrictions = compevent_restrictions,
                                 restrictions = restrictions,
                                 outcome_name = outcome_name, compevent_name = compevent_name,
                                 time_name = object$time_name,
                                 intvars = comb_intvars, interventions = comb_interventions,
                                 int_times = comb_int_times, histvars = histvars,
                                 histvals = histvals, histories = histories,
                                 covparams = covparams, covnames = covnames, covtypes = covtypes,
                                 covpredict_custom = covpredict_custom, basecovs = basecovs,
                                 comprisk = comprisk, ranges = ranges,
                                 yrange = yrange, compevent_range = compevent_range,
                                 outcome_type = outcome_type,
                                 subseed = subseed, t0 = t0, time_points = time_points,
                                 obs_data = newdata, parallel = parallel, max_visits = max_visits,
                                 baselags = baselags, below_zero_indicator = below_zero_indicator,
                                 min_time = min_time, show_progress = FALSE, ...)
    parallel::stopCluster(cl)
  } else {
    if (show_progress){
      pb <- progress::progress_bar$new(total = time_points,
                                       clear = FALSE,
                                       format = 'Simulation progress [:bar] :percent, Elapsed time :elapsed, Est. time remaining :eta')
    }
    pools <- lapply(seq_along(comb_interventions), FUN = function(i){
      gformula.simulate(fitcov = fitcov, fitY = fitY, fitD = fitD,
               yrestrictions = yrestrictions,
               compevent_restrictions = compevent_restrictions,
               restrictions = restrictions,
               outcome_name = outcome_name, compevent_name = compevent_name,
               time_name = object$time_name,
               intvars = comb_intvars[[i]], interventions = comb_interventions[[i]],
               int_times = comb_int_times[[i]], histvars = histvars, histvals = histvals,
               histories = histories, covparams = covparams,
               covnames = covnames, covtypes = covtypes,
               covpredict_custom = covpredict_custom, basecovs = basecovs, comprisk = comprisk,
               ranges = ranges, yrange = yrange, compevent_range = compevent_range,
               outcome_type = outcome_type,
               subseed = subseed, t0 = t0, time_points = time_points,
               obs_data = newdata, parallel = parallel, max_visits = max_visits,
               baselags = baselags, below_zero_indicator = below_zero_indicator,
               min_time = min_time, show_progress = TRUE, pb = pb,...)
    })
  }
  nat_pool <- pools[[1]] # Natural course data
  pools <- pools[-1] # List of intervention datasets
  
  # Calculate mean risk over all subjects at each time for natural course
  nat_result <- data.frame(id = unique(nat_pool[["oldid"]]))
  nat_result[,2:5] <- t(tapply(nat_pool$poprisk, list(nat_pool[[object$time_name]], nat_pool[["oldid"]]), FUN = mean))
  colnames(nat_result) <- c(id, "0", "1", "2", "3")
  
  # TODO: add functionality for multiple interventions
  
  if (return_sims == TRUE) {
    return(list(
      "Y.hat" = nat_result,
      "sims" = nat_pool
    ))
  } else {
    return(list(
      "Y.hat" = nat_result
    ))
  }
}

gformula.simulate <- function(o, fitcov, fitY, fitD,
                     yrestrictions, compevent_restrictions, restrictions,
                     outcome_name, compevent_name, time_name,
                     intvars, interventions, int_times, histvars, histvals, histories,
                     comprisk, ranges, yrange, compevent_range,
                     outcome_type, subseed, obs_data, time_points, t0, parallel,
                     covnames, covtypes, covparams, covpredict_custom,
                     basecovs, max_visits, baselags, below_zero_indicator,
                     min_time, show_progress, pb, ...){
  set.seed(subseed)
  
  # Mechanism of passing intervention variable and intervention is different for parallel
  # and non-parallel versions
  if (parallel){
    intvar <- intvars[[o]]
    intervention <- interventions[[o]]
    int_time <- int_times[[o]]
  } else {
    intvar <- intvars
    intervention <- interventions
    int_time <- int_times
  }
  
  if (!is.null(fitcov)){
    rmses <- lapply(seq_along(fitcov), FUN = gfoRmula:::rmse_calculate, fits = fitcov, covnames = covnames,
                    covtypes = covtypes, obs_data = obs_data, outcome_name = outcome_name,
                    time_name = time_name, restrictions = restrictions,
                    yrestrictions = yrestrictions, compevent_restrictions = compevent_restrictions)
  }
  
  # Initialize
  ids_unique <- unique(obs_data$newid)
  data_len <- length(ids_unique)
  restrict_ids <- rep(0, data_len)
  restrict_counts <- rep(list(rep(0, data_len)), length(restrictions))
  if (!is.na(restrictions[[1]][[1]])){
    restrict_covs <- lapply(restrictions, FUN = function(restriction){restriction[[1]]})
  }
  
  # Create histories_int and histvars_int, which are the necessary histories to create after the intervention
  if (!(length(intvar) == 1 && intvar == 'none')) {
    intvar_vec <- unique(unlist(intvar))
    histvars_int <- histories_int <- rep(list(NA), length(histvars))
    for (l in seq_along(histvars)){
      histvars_temp <- histvars[[l]][histvars[[l]] %in% intvar_vec]
      if (length(histvars_temp) > 0){
        histvars_int[[l]] <- histvars_temp
        histories_int[[l]] <- histories[[l]]
      }
    }
    histvars_int <- histvars_int[!is.na(histvars_int)]
    histories_int <- histories_int[!is.na(histories_int)]
  }
  
  for (t in ((1:time_points) - 1)){
    if (show_progress){
      pb$tick()
    }
    if (t <= t0){
      # Set simulated covariate values at time t = 0 equal to observed covariate values
      if (!is.na(basecovs[[1]])){
        pool <- obs_data[obs_data[[time_name]] <= t, ][, .SD, .SDcols = c(covnames, basecovs, time_name, 'oldid')]
      } else {
        pool <- obs_data[obs_data[[time_name]] <= t, ][, .SD, .SDcols = c(covnames, time_name, 'oldid')]
      }
      data.table::set(pool, j = 'id', value = obs_data[obs_data[[time_name]] <= t]$newid)
      if (!is.na(basecovs[[1]])){
        data.table::setcolorder(pool, c('id', time_name, covnames, basecovs, 'oldid'))
      } else {
        data.table::setcolorder(pool, c('id', time_name, covnames, 'oldid'))
      }
      newdf <- pool[pool[[time_name]] == 0]
      # Update datatable with specified treatment regime / intervention for this
      # simulation
      gfoRmula:::intfunc(newdf, pool = pool, intervention, intvar, unlist(int_time), time_name, t)
      if (ncol(newdf) > ncol(pool)){
        pool <- rbind(pool[pool[[time_name]] < t], newdf, fill = TRUE)
        pool <- pool[order(id, get(time_name))]
      } else {
        pool[pool[[time_name]] == t] <- newdf
      }
      # Update datatable with new covariates that are functions of history of existing
      # covariates
      gfoRmula:::make_histories(pool = pool, histvars = histvars, histvals = histvals,
                     histories = histories, time_name = time_name, t = t, id = 'id',
                     max_visits = max_visits, baselags = baselags, below_zero_indicator = below_zero_indicator)
      newdf <- pool[pool[[time_name]] == t]
      # Generate outcome probabilities
      if (outcome_type == 'survival'){
        data.table::set(newdf, j = 'Py', value = stats::predict(fitY, type = 'response', newdata = newdf))
      } else if (outcome_type == 'continuous_eof'){
        if (t < (time_points - 1)){
          data.table::set(newdf, j = 'Ey', value = as.double(NA))
        } else if (t == (time_points - 1)){
          data.table::set(newdf, j = 'Ey', value = stats::predict(fitY, type = 'response', newdata = newdf))
        }
      } else if (outcome_type == 'binary_eof'){
        if (t < (time_points - 1)){
          data.table::set(newdf, j = 'Py', value = as.double(NA))
        } else if (t == (time_points - 1)){
          data.table::set(newdf, j = 'Py', value = stats::predict(fitY, type = 'response', newdata = newdf))
        }
      }
      if (!is.na(yrestrictions[[1]][[1]])){ # Check if there are restrictions on outcome
        # variable simulation
        for (yrestriction in yrestrictions){
          # Set non-modeled outcome variable values equal to user-specified value
          if (outcome_type == 'survival'){
            newdf[!eval(parse(text = paste("newdf$", yrestriction[1]))),
                  "Py" := as.double(yrestriction[2])]
          }
        }
      }
      # Simulate outcome variable
      if (outcome_type == 'survival'){
        data.table::set(newdf, j = 'Y', value = stats::rbinom(data_len, 1, newdf$Py))
      }
      # Set simulated outcome values outside the observed range to the observed min / max
      if (length(newdf[newdf$Y < yrange[1]]$Y) != 0){
        data.table::set(newdf[newdf$Y < yrange[1]], j = 'Y', value = yrange[1])
      }
      if (length(newdf[newdf$Y > yrange[2]]$Y) != 0){
        data.table::set(newdf[newdf$Y > yrange[2]], j = 'Y', value = yrange[2])
      }
      if (outcome_type == 'survival')
      {
        if (comprisk){
          # Predict competing event probabilities
          data.table::set(newdf, j = 'Pd', value = stats::predict(fitD, type = 'response', newdata = newdf))
          if (!is.na(compevent_restrictions[[1]][[1]])){ # Check if there are restrictions
            # on competing event variable simulation
            for (compevent_restriction in compevent_restrictions){
              # Set non-modeled competing event values equal to user-specified value
              newdf[!eval(parse(text = compevent_restriction[1])),
                    "Pd" := as.double(compevent_restriction[2])]
            }
          }
          # Simulate competing event variable
          data.table::set(newdf, j = 'D', value = stats::rbinom(data_len, 1, newdf$Pd))
          # Set simulated competing event values outside the observed range to the observed
          # min / max
          if (length(newdf[newdf$D < compevent_range[1]]$D) != 0){
            data.table::set(newdf[newdf$D < compevent_range[1]], j = 'D', value = compevent_range[1])
          }
          if (length(newdf[newdf$D > compevent_range[2]]$D) != 0){
            data.table::set(newdf[newdf$D > compevent_range[2]], j = 'D', value = compevent_range[2])
          }
          # Calculate probability of death by main event rather than competing event at
          # time t
          data.table::set(newdf, j = 'prodp1', value = newdf$Py * (1 - newdf$Pd))
        } else {
          data.table::set(newdf, j = 'D', value = 0)
          # Calculate probability of death by main event without competing event
          if (outcome_type == 'survival'){
            data.table::set(newdf, j = 'prodp1', value = newdf$Py)
          }
        }
        data.table::set(newdf[newdf$D == 1], j = 'Y', value = NA)
        data.table::set(newdf, j = 'prodp0', value = 1 - newdf$Py)
      }
      # If competing event occurs, outcome cannot also occur because
      # both presumably lead to death
      # Calculate probability of survival or death by competing event at time t
      pool <- rbind(pool[pool[[time_name]] < t], newdf, fill = TRUE)
      pool <- pool[order(id, get(time_name))]
      col_types <- sapply(pool, class)
    } else {
      # Set initial simulated values at time t to simulated values at time t - 1, to be
      # updated later
      newdf <- pool[pool[[time_name]] == t - 1]
      data.table::set(newdf, j = time_name, value = rep(t, data_len))
      if ('categorical time' %in% covtypes){
        time_name_f <- paste(time_name, "_f", sep = "")
        newdf[, (time_name_f) :=
                obs_data[get(time_name) == t, get(time_name_f)][1]]
      }
      pool <- rbind(newdf, pool)
      gfoRmula:::make_histories(pool = pool, histvars = histvars, histvals = histvals, histories = histories,
                     time_name = time_name, t = t, id = 'id', max_visits = max_visits,
                     baselags = baselags, below_zero_indicator = below_zero_indicator)
      newdf <- pool[pool[[time_name]] == t]
      for (i in seq_along(covnames)){
        cast <- get(paste0('as.',unname(col_types[covnames[i]])))
        if (covtypes[i] == 'binary'){
          data.table::set(newdf, j = covnames[i],
              value = cast(gfoRmula:::predict_binomial(data_len, 1, stats::predict(fitcov[[i]], type = 'response',
                                                                        newdata = newdf))))
        } else if (covtypes[i] == 'normal'){
          data.table::set(newdf, j = covnames[i],
              value = cast(gfoRmula:::predict_normal(data_len, stats::predict(fitcov[[i]], type = 'response',
                                                                   newdata = newdf),
                                          est_sd = rmses[[i]])))
        } else if (covtypes[i] == 'categorical'){
          data.table::set(newdf, j = covnames[i],
              value = cast(stats::predict(fitcov[[i]], type = 'class', newdata = newdf)))
        } else if (covtypes[i] == 'zero-inflated normal'){
          data.table::set(newdf, j = paste("I_", covnames[i], sep = ""),
              value = gfoRmula:::predict_binomial(data_len, 1, stats::predict(fitcov[[i]][[1]], type = 'response',
                                                                   newdata = newdf)))
          # Where indicator is nonzero, simulate from Gaussian model
          data.table::set(newdf, j = covnames[i],
              value = cast(exp(gfoRmula:::predict_normal(data_len,
                                              stats::predict(fitcov[[i]][[2]], type = 'response',
                                                             newdata = newdf),
                                              est_sd =  rmses[[i]]))))
          data.table::set(newdf, j = covnames[i],
              value = cast(newdf[[paste("I_", covnames[i], sep = "")]] * newdf[[covnames[i]]]))
          # Remove indicator
          data.table::set(newdf, j = paste("I_", covnames[i], sep = ""), value = NULL)
        } else if (covtypes[i] == 'bounded normal'){
          if (!is.na(restrictions[[1]][[1]])){
            restrictnames <- lapply(seq_along(restrictions), FUN = function(r){
              restrictions[[r]][[1]]})
            # Create list of conditions where covariates are modeled
            conditions <- lapply(seq_along(restrictions), FUN = function(r){
              restrictions[[r]][[2]]
            })
            if (covnames[i] %in% restrictnames){
              j <- which(restrictnames %in% covnames[i])
              condition <- ""
              if (length(j) > 1){
                for (k in j){
                  condition_var <- sub("<.*$", "", restrictions[[k]][[2]])
                  condition_var <- sub(">.*$", "", condition_var)
                  condition_var <- sub("=.*$", "", condition_var)
                  condition_var <- sub("!.*$", "", condition_var)
                  condition_var <- sub("%in%.*$", "", condition_var)
                  condition_var <- sub(" ", "", condition_var)
                  if (condition_var %in% names(obs_data)){
                    if (condition[1] == ""){
                      condition <- conditions[j]
                    } else {
                      condition <- paste(condition, conditions[j], sep = "||")
                    }
                  }
                }
              } else {
                condition_var <- sub("<.*$", "", restrictions[[j]][[2]])
                condition_var <- sub(">.*$", "", condition_var)
                condition_var <- sub("=.*$", "", condition_var)
                condition_var <- sub("!.*$", "", condition_var)
                condition_var <- sub("%in%.*$", "", condition_var)
                condition_var <- sub(" ", "", condition_var)
                if (condition_var %in% names(obs_data)){
                  condition <- conditions[j]
                }
              }
              if (condition[1] != ""){
                sub_obs_data <- subdata.table::set(obs_data[obs_data[[time_name]] >= 0], eval(parse(text = condition)))
              } else {
                sub_obs_data <- obs_data[obs_data[[time_name]] >= 0]
              }
            } else {
              sub_obs_data <- obs_data[obs_data[[time_name]] >= 0]
            }
          } else {
            sub_obs_data <- obs_data[obs_data[[time_name]] >= 0]
          }
          data.table::set(newdf, j = paste("norm_", covnames[i], sep = ""),
              value = gfoRmula:::predict_normal(data_len, stats::predict(fitcov[[i]], type = 'response',
                                                              newdata = newdf), est_sd = rmses[[i]]))
          data.table::set(newdf, j = covnames[i],
              value = cast((newdf[[paste("norm_", covnames[i], sep = "")]] *
                              (max(sub_obs_data[[covnames[i]]]) - min(sub_obs_data[[covnames[i]]]))) +
                             min(sub_obs_data[[covnames[i]]])))
          data.table::set(newdf, j = paste("norm_", covnames[i], sep = ""), value = NULL)
        } else if (covtypes[i] == 'truncated normal'){
          if (covparams$direction[i] == 'left'){
            data.table::set(newdf, j = covnames[i],
                value = cast(gfoRmula:::predict_trunc_normal(data_len, stats::predict(fitcov[[i]], type = 'response',
                                                                           newdata = newdf),
                                                  est_sd = rmses[[i]], a = covparams$point[i], b = Inf)))
          } else if (covparams$direction[i] == 'right'){
            data.table::set(newdf, j = covnames[i],
                value = cast(gfoRmula:::predict_trunc_normal(data_len, stats::predict(fitcov[[i]], type = 'response',
                                                                           newdata = newdf),
                                                  est_sd = rmses[[i]], a = - Inf, b = covparams$point[i])))
          }
        } else if (covtypes[i] == 'custom'){
          if (!is.na(restrictions[[1]][[1]])){
            restrictnames <- lapply(seq_along(restrictions), FUN = function(r){
              restrictions[[r]][[1]]})
            # Create list of conditions where covariates are modeled
            conditions <- lapply(seq_along(restrictions), FUN = function(r){
              restrictions[[r]][[2]]
            })
            if (covnames[i] %in% restrictnames){
              j <- which(restrictnames %in% covnames[i])
              condition <- ""
              if (length(j) > 1){
                for (k in j){
                  condition_var <- sub("<.*$", "", restrictions[[k]][[2]])
                  condition_var <- sub(">.*$", "", condition_var)
                  condition_var <- sub("=.*$", "", condition_var)
                  condition_var <- sub("!.*$", "", condition_var)
                  condition_var <- sub("%in%.*$", "", condition_var)
                  condition_var <- sub(" ", "", condition_var)
                  if (condition_var %in% names(obs_data)){
                    if (condition[1] == ""){
                      condition <- conditions[j]
                    } else {
                      condition <- paste(condition, conditions[j], sep = "||")
                    }
                  }
                }
              } else {
                condition_var <- sub("<.*$", "", restrictions[[j]][[2]])
                condition_var <- sub(">.*$", "", condition_var)
                condition_var <- sub("=.*$", "", condition_var)
                condition_var <- sub("!.*$", "", condition_var)
                condition_var <- sub("%in%.*$", "", condition_var)
                condition_var <- sub(" ", "", condition_var)
                if (condition_var %in% names(obs_data)){
                  condition <- conditions[j]
                }
              }
            }
          }
          data.table::set(newdf, j = covnames[i],
              value = cast(covpredict_custom[[i]](obs_data, newdf, fitcov[[i]],time_name, t,
                                                  condition, covnames[i], ...)))
        }
        if (covtypes[i] == 'normal' || covtypes[i] == 'bounded normal' ||
            covtypes[i] == 'truncated normal'){
          
          if (length(newdf[newdf[[covnames[i]]] < ranges[[i]][1]][[covnames[i]]]) != 0){
            newdf[newdf[[covnames[i]]] < ranges[[i]][1], (covnames[i]) := cast(ranges[[i]][1])]
          }
          if (length(newdf[newdf[[covnames[i]]] > ranges[[i]][2]][[covnames[i]]]) != 0){
            newdf[newdf[[covnames[i]]] > ranges[[i]][2], (covnames[i]) := cast(ranges[[i]][2])]
          }
        } else if (covtypes[i] == 'zero-inflated normal') {
          if (length(newdf[newdf[[covnames[i]]] < ranges[[i]][1] & newdf[[covnames[i]]] > 0][[covnames[i]]]) != 0){
            newdf[newdf[[covnames[i]]] < ranges[[i]][1] & newdf[[covnames[i]]] > 0, (covnames[i]) := cast(ranges[[i]][1])]
          }
          if (length(newdf[newdf[[covnames[i]]] > ranges[[i]][2]][[covnames[i]]]) != 0){
            newdf[newdf[[covnames[i]]] > ranges[[i]][2], (covnames[i]) := cast(ranges[[i]][2])]
          }
        }
        # Check if there are restrictions on covariate simulation
        if (!is.na(restrictions[[1]][[1]])){
          lapply(seq_along(restrictions), FUN = function(r){
            if (restrictions[[r]][[1]] == covnames[i]){
              restrict_ids <- newdf[!eval(parse(text = restrictions[[r]][[2]]))]$id
              if (length(restrict_ids) != 0){
                restrictions[[r]][[3]](newdf, pool[pool[[time_name]] < t & pool[[time_name]] >= 0], restrictions[[r]], time_name, t)
              }
            }
          })
        }
        pool[pool[[time_name]] == t] <- newdf
        if (covnames[i] %in% unlist(histvars)){
          ind <- unlist(lapply(histvars, FUN = function(x) {
            covnames[i] %in% x
          }))
          gfoRmula:::make_histories(pool = pool, histvars = rep(list(covnames[i]), sum(ind)),
                         histvals = histvals, histories = histories[ind],
                         time_name = time_name, t = t, id = 'id', max_visits = max_visits,
                         baselags = baselags, below_zero_indicator = below_zero_indicator)
          newdf <- pool[pool[[time_name]] == t]
        }
      }
      # Update datatable with specified treatment regime / intervention for this
      # simulation
      newdf <- pool[pool[[time_name]] == t]
      gfoRmula:::intfunc(newdf, pool, intervention, intvar, unlist(int_time), time_name, t)
      # Update datatable with new covariates that are functions of history of existing
      # covariates
      pool[pool[[time_name]] == t] <- newdf
      
      if (!(length(intvar) == 1 && intvar == 'none') && length(histvars_int) > 0){
        gfoRmula:::make_histories(pool = pool, histvars = histvars_int, histvals = histvals,
                       histories = histories_int, time_name = time_name, t = t, id = 'id',
                       max_visits = max_visits, baselags = baselags, below_zero_indicator = below_zero_indicator)
      }
      
      newdf <- pool[pool[[time_name]] == t]
      
      # Predict outcome probabilities
      if (outcome_type == 'survival'){
        if (comprisk){
          # Predict competing event probabilities
          data.table::set(newdf, j = 'Pd', value = stats::predict(fitD, type = 'response', newdata = newdf))
          if (!is.na(compevent_restrictions[[1]][[1]])){ # Check if there are restrictions
            # on competing event variable
            # simulation
            for (compevent_restriction in compevent_restrictions){
              # Set non-modeled competing event values equal to user-specified value
              newdf[!eval(parse(text = compevent_restriction[1])), "Pd" := as.double(compevent_restriction[2])]
            }
          }
          # Simulate competing event variable
          data.table::set(newdf, j = 'D', value = stats::rbinom(data_len, 1, newdf$Pd))
          # Set simulated competing event values outside the observed range to the observed
          # min / max
          if (length(newdf[newdf$D < compevent_range[1]]$D) != 0){
            newdf[newdf$D < compevent_range[1], "D" := compevent_range[1]]
          }
          if (length(newdf[newdf$D > compevent_range[2]]$D) != 0){
            newdf[newdf$D > compevent_range[2], "D" := compevent_range[2]]
          }
        } else {
          data.table::set(newdf, j = 'D', value = 0)
        }
        data.table::set(newdf, j = 'Py', value = stats::predict(fitY, type = 'response', newdata = newdf))
      } else if (outcome_type == 'continuous_eof'){
        if (t < (time_points - 1)){
          data.table::set(newdf, j = 'Ey', value = as.double(NA))
        } else if (t == (time_points - 1)){
          data.table::set(newdf, j = 'Ey', value = stats::predict(fitY, type = 'response', newdata = newdf))
        }
      } else if (outcome_type == 'binary_eof'){
        if (t < (time_points - 1)){
          data.table::set(newdf, j = 'Py', value = as.double(NA))
        } else if (t == (time_points - 1)){
          data.table::set(newdf, j = 'Py', value = stats::predict(fitY, type = 'response', newdata = newdf))
        }
      }
      if (!is.na(yrestrictions[[1]][[1]])){ # Check if there are restrictions on outcome
        # variable simulation
        for (yrestriction in yrestrictions){
          # Set non-modeled outcome variable values equal to user-specified value
          if (outcome_type == 'survival'){
            newdf[!eval(parse(text = paste("newdf$", yrestriction[1]))), "Py" := as.double(yrestriction[2])]
          } else if (outcome_type == 'continuous_eof' && t == (time_points - 1)){
            newdf[!eval(parse(text = paste("newdf$", yrestriction[1]))), "Ey" := as.double(yrestriction[2])]
          } else if (outcome_type == 'binary_eof' && t == (time_points - 1)){
            newdf[!eval(parse(text = paste("newdf$", yrestriction[1]))), "Py" := as.double(yrestriction[2])]
          }
        }
      }
      # Calculate probability of survival or death from competing event (if any) at time t
      if (outcome_type == 'survival'){
        data.table::set(newdf, j = 'prodp0', value = 1 - newdf$Py)
      }
      # Simulate outcome variable
      if (outcome_type == 'survival'){
        data.table::set(newdf, j = 'Y', value = stats::rbinom(data_len, 1, newdf$Py))
      }
      # Set simulated outcome values outside the observed range to the observed min / max
      if (length(newdf[newdf$Y < yrange[1]]$Y) != 0){
        newdf[newdf$Y < yrange[1], 'Y' := yrange[1]]
      }
      if (length(newdf[newdf$Y > yrange[2]]$Y) != 0){
        newdf[newdf$Y > yrange[2], 'Y' := yrange[2]]
      }
      if (outcome_type == 'survival')
        newdf[newdf$D == 1, 'Y' := NA]
      # If competing event occurs, outcome cannot also occur because
      # both presumably lead to death
      # Calculate probability of death from main event at time t
      if (comprisk){
        data.table::set(newdf, j = 'prodp1',
            value = newdf$Py * tapply(pool[pool[[time_name]] < t & pool[[time_name]] >= 0]$prodp0,
                                      pool[pool[[time_name]] < t & pool[[time_name]] >= 0]$id, FUN = prod) *
              tapply(1 - pool[pool[[time_name]] < t & pool[[time_name]] >= 0]$Pd,
                     pool[pool[[time_name]] < t & pool[[time_name]] >= 0]$id,
                     FUN = prod) * (1 - newdf$Pd))
      } else if (outcome_type == 'survival'){
        data.table::set(newdf, j = 'prodp1', value = newdf$Py * tapply(pool[pool[[time_name]] < t & pool[[time_name]] >= 0]$prodp0,
                                                           pool[pool[[time_name]] < t & pool[[time_name]] >= 0]$id,
                                                           FUN = prod))
      }
      # Add simulated data for time t to aggregate simulated data over time
      pool[pool[[time_name]] == t] <- newdf
    }
  }
  colnames(pool)[colnames(pool) == time_name] <- 't0'
  data.table::setorder(pool, id, t0)
  colnames(pool)[colnames(pool) == 't0'] <- time_name
  # Calculate probabiity of death from main event at or before time t for each individual
  # at each time point
  pool <- pool[pool[[time_name]] >= 0]
  if (outcome_type == 'survival'){
    pool[, 'poprisk' := stats::ave(pool$prodp1, by = pool$id, FUN = cumsum)]
    pool[, 'survival' := stats::ave(pool$prodp0, by = pool$id, FUN = cumprod)]
  }
  pool2 <- data.table::copy(pool)

  return (pool2)
}

if (FALSE) {
  fit_ldl <- 
    gformula_survival(
      obs_data = drop_na(analytic_long),
      id = "pid",
      time_name = "time",
      time_points = 4,
      covnames = covs_dvs[covs_dvs != "liprx"], 
      covtypes = covtypes,
      covparams = covparams,
      outcome_name = "event_chd",
      ymodel = ymodel,
      compevent_name = "event_dth",
      compevent_model = compevent_model,
      restrictions = restrictions,
      basecovs = covs_fixed[!covs_fixed %in% covs_refs],
      histvars = list(covs_dvs[covs_dvs != "liprx"]),
      histories = c(lagged),
      nsimul = GFORM_SIM,
      seed = 1234,
      show_progress = TRUE,
      sim_data_b = TRUE,
      model_fits = TRUE
    )
  
  newdata <- fit_ldl$sim_data$`Natural course`[1:400]

  predict.gformula(
    object = fit_ldl,
    obs_data = data.table::as.data.table(drop_na(analytic_long)),
    newdata = data.table::as.data.table(newdata),
    id = "id",
    covnames = covs_dvs[covs_dvs != "liprx"], 
    covtypes = covtypes,
    covparams = covparams,
    outcome_name = "event_chd",
    ymodel = ymodel,
    compevent_name = "event_dth",
    compevent_model = compevent_model,
    restrictions = restrictions,
    basecovs = covs_fixed[!covs_fixed %in% covs_refs],
    histvars = list(covs_dvs[covs_dvs != "liprx"]),
    histories = c(lagged),
    nsamples = 0,
    nsimul = GFORM_SIM,
    seed = 1234,
    show_progress = TRUE,
    # parallel = TRUE,
    # ncores = parallel::detectCores() - 1,
    model_fits = TRUE
  )
}