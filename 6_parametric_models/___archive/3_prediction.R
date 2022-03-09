# Use parametric g-formula to predict risks -------------------------------



# option 1  ---------------------------------------------------------------

# start from exam X, iteratively update using data from subsequent exams until
# you reach the final exam, the number of exams each prediction is conditional
# on is determined by the lag structure of the parametric models, e.g.
#
# parametric models w/ 2 lags:
#   iter 0 - conditional on exam 0 alone project 0 through 5
#   iter 1 - conditional on exams 0 to 1 project 1 through 5
#   iter 2 - conditional on exams 0 to 2 project 2 through 5
#   iter 3 - conditional on exams 1 to 3 project 3 through 5
#   iter 4 - conditional on exams 2 to 4 project 4 through 5
#   iter 5 - conditional on exams 3 to 5 project 5


run_params <- list(
  name = list(
    "no lags",
    "one lag",
    "two lags"
  ),
  models = parametric_models_uncens_noD,
  data = list(filter(offspring_long, time < 5 & !pid %in% censored_ids)),
  stop = 4,
  sims = 100
)

plan(multisession, workers = 14)

with_progress({
  p <- progressor(steps = length(0:run_params$stop) * length(run_params$name))
  
  gformula_predictions_opt1 <- 
    future_pmap(
      run_params,
      function(name, models, data, stop, sims) {
        updated_preds <- 
          map(
            set_names(0:stop, paste0("pred_", 0:stop)),
            function (x) {
              pred <- gformula_mc(
                Y.fit = models[["Y"]],
                X.fit = models[["X"]],
                D.fit = models[["D"]],
                data = data,
                id = "pid",
                time = "time",
                base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
                hist.vars = c(paste0("lag1_", covs_tv), paste0("lag2_", covs_tv)),
                hist.fun = "lag",
                mc.sims = sims,
                mc.start = x,
                mc.stop = stop,
                merge = FALSE
              )
              p()
              return(pred)
            }
          )
        
        updated_preds <- reduce(updated_preds, left_join, by = c("pid", "time"))
        names(updated_preds) <- c("pid", "time", paste0("pred_", 0:stop))
        
        updated_preds <- left_join(
          observed_events,
          updated_preds,
          by = c("pid", "time")
        )
        
        return(updated_preds)
      },
      .options = furrr_options(
        seed = TRUE
      )
    )
})

plan(sequential)

names(gformula_predictions_opt1) <- run_params$name


# option 2 ----------------------------------------------------------------

# always start from same exam (e.g. exam 2) and iteratively add previous exam
# information to obtain updated predictions -- here the advantage is that we are
# comparing predictions over the same time horizon and observed cases/deaths,
# etc.

# parametric models w/ 2 lags:
#   iter 0 - conditional on exam 2 alone project 2 through 5
#   iter 1 - conditional on exams 1 to 2 project 2 through 5

run_params <- list(
  name = list(
    "no lags",
    "one lag",
    "two lags"
  ),
  models = parametric_models_1to5_noD,
  data = list(filter(offspring_long, time == 2 & !pid %in% censored_ids)),
  start = 2,
  stop = 4,
  sims = 100
)

plan(multisession, workers = 14)

with_progress({
  p <- progressor(steps = length(0:2) * length(run_params$name))
  
  gformula_predictions_opt2 <- 
    future_pmap(
      run_params,
      function(name, models, data, start, stop, sims) {
        updated_preds <- 
          map(
            set_names(0:2, paste0("pred_", 0:2)),
            function (x) {
              pred <- gformula_mc(
                Y.fit = models[["Y"]],
                X.fit = models[["X"]],
                D.fit = models[["D"]],
                data = zero_out_lags(data, x),
                id = "pid",
                time = "time",
                base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
                hist.vars = c(paste0("lag1_", covs_tv), paste0("lag2_", covs_tv)),
                hist.fun = "lag",
                mc.sims = sims,
                mc.start = start,
                mc.stop = stop,
                merge = FALSE
              )
              p()
              return(pred)
            }
          )
        
        updated_preds <- reduce(updated_preds, left_join, by = c("pid", "time"))
        names(updated_preds) <- c("pid", "time", paste0("pred_", 0:2))
        
        updated_preds <- left_join(
          observed_events,
          updated_preds,
          by = c("pid", "time")
        )
        
        return(updated_preds)
      },
      .options = furrr_options(
        seed = TRUE
      )
    )
})

plan(sequential)

names(gformula_predictions_opt2) <- run_params$name
