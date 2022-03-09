# Use parametric g-formula to predict risks -------------------------------


# experiment 1 ------------------------------------------------------------

# always start from same exam (e.g. exam 2) and iteratively add previous exam
# information to obtain updated predictions -- here the advantage is that we are
# comparing predictions over the same time horizon and observed cases/deaths,
# etc.

# parametric models w/ 2 lags:
#   iter 0 - conditional on exam 2 alone project 2 through 5
#   iter 1 - conditional on exams 1 to 2 project 2 through 5

run_params <- list(
  name = c(
    "no lags, net risk",
    "one lag, net risk",
    "two lags, net risk",
    "no lags, ASCVD-specific risk",
    "one lag, ASCVD-specific risk",
    "two lags, ASCVD-specific risk"
  ),
  models = c(
    models_ex1_noD,
    models_ex1
  ),
  lags = rep(0:2, each = 2),
  start = 0,
  stop = 2,
  sims = 100
)

with_progress({
  p <- progressor(steps = length(run_params$name))
  
  gformula_predictions_ex1 <- 
    pmap(
      run_params,
      function(name, models, lags, start, stop, sims, compevent) {
        
        pred <- gformula_mc(
          Y.fit = models[["Y"]],
          X.fit = models[["X"]],
          D.fit = models[["D"]],
          data = zero_out_lags(offspring_long, 2-lags),
          id = "pid",
          time = "time",
          base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
          hist.vars = c("lag1_age", "lag2_age", covs_lag1_ex1, covs_lag2_ex1),
          hist.fun = "lag",
          mc.sims = sims,
          mc.start = start,
          mc.stop = stop,
          merge = FALSE
        )
        p()
        
        pred <- left_join(
          observed_events,
          pred,
          by = c("pid", "time")
        )
        
        return(pred)
      }
      # .options = furrr_options(
      #   seed = TRUE
      # )
    )
})

names(gformula_predictions_ex1) <- run_params$name

# experiment 2 ------------------------------------------------------------

# always start from same exam (e.g. exam 2) and iteratively add previous exam
# information to obtain updated predictions -- here the advantage is that we are
# comparing predictions over the same time horizon and observed cases/deaths,
# etc.

# parametric models w/ 2 lags:
#   iter 0 - conditional on exam 2 alone project 2 through 5
#   iter 1 - conditional on exams 1 to 2 project 2 through 5


run_params <- list(
  name = c(
    "bmi, net risk",
    "bmi + ldl, net risk",
    "bmi + ldl + liprx, net risk",
    "bmi, ASCVD-specific risk",
    "bmi + ldl, ASCVD-specific risk",
    "bmi + ldl + liprx, ASCVD-specific risk"
  ),
  models = c(
    models_ex2_noD,
    models_ex2
  ),
  id = rep(1:3, 2),
  start = 0,
  stop = 2,
  sims = 100
)


with_progress({
  p <- progressor(steps = length(run_params$name))
  
  gformula_predictions_ex2 <- 
    pmap(
      run_params,
      function(name, models, data, id, start, stop, sims, compevent) {
        pred <- gformula_mc(
          Y.fit = models[["Y"]],
          X.fit = models[["X"]],
          D.fit = models[["D"]],
          data = data,
          id = "pid",
          time = "time",
          base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
          hist.vars = c("lag1_age", "lag2_age", covs_lag1_ex2[[id]], covs_lag2_ex2[[id]]),
          hist.fun = "lag",
          mc.sims = sims,
          mc.start = start,
          mc.stop = stop,
          merge = FALSE
        )
        p()
        
        pred <- left_join(
          observed_events,
          pred,
          by = c("pid", "time")
        )
        
        return(pred)
      },
      data = offspring_long
      # .options = furrr_options(
      #   seed = TRUE
      # )
    )
})

names(gformula_predictions_ex2) <- run_params$name


# experiment 3 ------------------------------------------------------------

# always start from same exam (e.g. exam 2) and iteratively add previous exam
# information to obtain updated predictions -- here the advantage is that we are
# comparing predictions over the same time horizon and observed cases/deaths,
# etc.

# parametric models w/ 2 lags:
#   iter 0 - conditional on exam 2 alone project 2 through 5
#   iter 1 - conditional on exams 1 to 2 project 2 through 5


run_params <- list(
  name = c(
    "gam, net risk",
    "gam, ASCVD-specific risk"
  ),
  models = c(
    models_ex3_noD,
    models_ex3
  ),
  id = 3,
  start = 0,
  stop = 2,
  sims = 100
)
with_progress({
  p <- progressor(steps = length(run_params$name))
  
  gformula_predictions_ex3 <- 
    pmap(
      run_params,
      function(name, models, data, id, start, stop, sims, compevent) {
        pred <- gformula_mc(
          Y.fit = models[["Y"]],
          X.fit = models[["X"]],
          D.fit = models[["D"]],
          data = data,
          id = "pid",
          time = "time",
          base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
          hist.vars = c("lag1_age", "lag2_age", covs_lag1_ex2[[id]], covs_lag2_ex2[[id]]),
          hist.fun = "lag",
          mc.sims = sims,
          mc.start = start,
          mc.stop = stop,
          merge = FALSE
        )
        p()
        
        pred <- left_join(
          observed_events,
          pred,
          by = c("pid", "time")
        )
        
        return(pred)
      },
      data = offspring_long
      # .options = furrr_options(
      #   seed = TRUE
      # )
    )
})

names(gformula_predictions_ex3) <- run_params$name


# conventional models -----------------------------------------------------

# always start from same exam (e.g. exam 2) and iteratively add previous exam
# information to obtain updated predictions -- here the advantage is that we are
# comparing predictions over the same time horizon and observed cases/deaths,
# etc.

# parametric models w/ 2 lags:
#   iter 0 - conditional on exam 2 alone project 2 through 5
#   iter 1 - conditional on exams 1 to 2 project 2 through 5


run_params <- list(
  name = c(
    "logistic, net risk",
    "coxph, net risk",
    "logistic, ASCVD-specific risk",
    "coxph, ASCVD-specific risk"
  ),
  models = list(
    conventional_models[[1]],
    conventional_models[[2]],
    conventional_models[[1]],
    conventional_models[[2]]
  ),
  data = list(
    offspring_wide,
    offspring_coxph,
    offspring_wide,
    offspring_coxph
  )
)

with_progress({
  p <- progressor(steps = length(run_params$name))
  
  conventional_predictions <- 
    pmap(
      run_params,
      function(name, models, data, sims, compevent) {
        if (str_detect(name, "coxph")) {
          data$pred <- 1 - predict(models, newdata = mutate(data, enddate = 3652), type = "survival")
        } else {
          data$pred <- predict(models, newdata = data, type = "response")
        }
        
        pred <- left_join(
          observed_events_wide,
          select(data, pid, pred),
          by = "pid"
        )
        
        return(pred)
      }
      # .options = furrr_options(
      #   seed = TRUE
      # )
    )
})

names(conventional_predictions) <- run_params$name
