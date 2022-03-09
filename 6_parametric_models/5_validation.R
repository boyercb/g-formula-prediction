# Use parametric g-formula to predict risks -------------------------------

BOOTS <- 500

# option 2 ----------------------------------------------------------------

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
  def = c(definitions_ex1, definitions_ex1),
  data = list(offspring_long),
  lags = rep(0:2, each = 2),
  start = 0,
  stop = 2,
  sims = 100
)


plan(multisession, workers = 10)

with_progress({
  p <- progressor(steps = length(run_params$name) * BOOTS)
  
  boot_stats_ex1 <- future_map_dfr(
    1:BOOTS,
    function(bsim, run_params) {
      validation_preds_ex1 <- 
        pmap(
          run_params,
          function(name, def, data, start, stop, lags, sims) {
            boot <- draw_bootstrap_samples(data, "pid", "time")
            oob <- data[!data[['pid']] %in% boot[['pid']], ]
            
            if (str_detect(name, "ASCVD")) {
              defD <- def[["D"]]
            } else {
              defD <- NULL
            }
            
            models <- list(
              "X" = fit_X_models(
                models = def[["X"]],
                data = boot,
                time = "time"
              ),
              "Y" = fit_Y_model(
                model = def[["Y"]],
                data = boot
              ),
              "D" = if (!is.null(defD)) {
                fit_D_model(
                  model = defD,
                  data = boot
                )
              } else {
                NULL
              }
            )
            
            pred_boot <- gformula_mc(
              Y.fit = models[["Y"]],
              X.fit = models[["X"]],
              D.fit = models[["D"]],
              data = zero_out_lags(boot, 2-lags),
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
            
            pred_test <- gformula_mc(
              Y.fit = models[["Y"]],
              X.fit = models[["X"]],
              D.fit = models[["D"]],
              data = zero_out_lags(data, 2-lags),
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
            
            pred_boot <- left_join(
              observed_events,
              pred_boot,
              by = c("pid", "time")
            )
            
            pred_test <- left_join(
              observed_events,
              pred_test,
              by = c("pid", "time")
            )
            
            return(list(
              boot = pred_boot,
              test = pred_test
            ))
            
          }
        )
      
      names(validation_preds_ex1) <- run_params$name
      
      stats <- map2_dfr(validation_preds_ex1,
                        names(validation_preds_ex1),
                        function(x, name, time = 2) {
                          if (str_detect(name, "ASCVD")) {
                            outcome <- "event_ascvd_zero"
                          } else {
                            outcome <- "event_ascvd"
                          }
                          
                          boot_stats <- 
                            val.prob(p = x$boot[x$boot[["time"]] == time,][["pred"]],
                                     y = x$boot[x$boot[["time"]] == time,][[outcome]],
                                     pl = FALSE)
                          
                          test_stats <- 
                            val.prob(p = x$test[x$test[["time"]] == time,][["pred"]],
                                     y = x$test[x$test[["time"]] == time,][[outcome]],
                                     pl = FALSE)
                          
                          bind_rows(boot_stats, test_stats, .id = "test")
                        }, .id = "model")
      
      return(stats)
    },
    run_params = run_params,
    .id = "bsim",
    .options = furrr_options(
      seed = TRUE
    )
  )
  
})

plan(sequential)

# option 2 ----------------------------------------------------------------

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
  def = c(definitions_ex2, definitions_ex2),
  data = list(offspring_long),
  id = rep(1:3, 2),
  start = 0,
  stop = 2,
  sims = 100
)

plan(multisession, workers = 10)

with_progress({
  p <- progressor(steps = length(run_params$name) * BOOTS)
  
  boot_stats_ex2 <- future_map_dfr(
    1:BOOTS,
    function(bsim, run_params) {
      validation_preds_ex2 <- 
        pmap(
          run_params,
          function(name, def, data, id, start, stop, sims) {
            boot <- draw_bootstrap_samples(data, "pid", "time")
            oob <- data[!data[['pid']] %in% boot[['pid']], ]
            
            if (str_detect(name, "ASCVD")) {
              defD <- def[["D"]]
            } else {
              defD <- NULL
            }
            
            models <- list(
              "X" = fit_X_models(
                models = def[["X"]],
                data = boot,
                time = "time"
              ),
              "Y" = fit_Y_model(
                model = def[["Y"]],
                data = boot
              ),
              "D" = if (!is.null(defD)) {
                fit_D_model(
                  model = defD,
                  data = boot
                )
              } else {
                NULL
              }
            )
            
            pred_boot <- gformula_mc(
              Y.fit = models[["Y"]],
              X.fit = models[["X"]],
              D.fit = models[["D"]],
              data = boot,
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
            
            pred_test <- gformula_mc(
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
            
            pred_boot <- left_join(
              observed_events,
              pred_boot,
              by = c("pid", "time")
            )
            
            pred_test <- left_join(
              observed_events,
              pred_test,
              by = c("pid", "time")
            )
            
            return(list(
              boot = pred_boot,
              test = pred_test
            ))
            
          }
        )
      
      names(validation_preds_ex2) <- run_params$name
      
      stats <- map2_dfr(validation_preds_ex2,
                        names(validation_preds_ex2),
                        function(x, name, time = 2) {
                          if (str_detect(name, "ASCVD")) {
                            outcome <- "event_ascvd_zero"
                          } else {
                            outcome <- "event_ascvd"
                          }
                          
                          boot_stats <- 
                            val.prob(p = x$boot[x$boot[["time"]] == time,][["pred"]],
                                     y = x$boot[x$boot[["time"]] == time,][[outcome]],
                                     pl = FALSE)
                          
                          test_stats <- 
                            val.prob(p = x$test[x$test[["time"]] == time,][["pred"]],
                                     y = x$test[x$test[["time"]] == time,][[outcome]],
                                     pl = FALSE)
                          
                          bind_rows(boot_stats, test_stats, .id = "test")
                          
                        }, .id = "model")
      
      return(stats)
    },
    run_params = run_params,
    .id = "bsim",
    .options = furrr_options(
      seed = TRUE
    )
  )
  
})

plan(sequential)

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
    "gam, net risk",
    "gam, ASCVD-specific risk"
  ),
  def = c(definitions_ex3, definitions_ex3),
  data = list(offspring_long),
  id = 3,
  start = 0,
  stop = 2,
  sims = 100
)

plan(multisession, workers = 10)

with_progress({
  p <- progressor(steps = length(run_params$name) * BOOTS)
  
  boot_stats_ex3 <- future_map_dfr(
    1:BOOTS,
    function(bsim, run_params) {
      validation_preds_ex3 <- 
        pmap(
          run_params,
          function(name, def, data, id, start, stop, sims) {
            boot <- draw_bootstrap_samples(data, "pid", "time")
            oob <- data[!data[['pid']] %in% boot[['pid']], ]
            
            if (str_detect(name, "ASCVD")) {
              defD <- def[["D"]]
            } else {
              defD <- NULL
            }
            
            models <- list(
              "X" = fit_X_models(
                models = def[["X"]],
                data = boot,
                time = "time"
              ),
              "Y" = fit_Y_model(
                model = def[["Y"]],
                data = boot
              ),
              "D" = if (!is.null(defD)) {
                fit_D_model(
                  model = defD,
                  data = boot
                )
              } else {
                NULL
              }
            )
            
            pred_boot <- gformula_mc(
              Y.fit = models[["Y"]],
              X.fit = models[["X"]],
              D.fit = models[["D"]],
              data = boot,
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

            pred_test <- gformula_mc(
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
            
            pred_boot <- left_join(
              observed_events,
              pred_boot,
              by = c("pid", "time")
            )
            
            pred_test <- left_join(
              observed_events,
              pred_test,
              by = c("pid", "time")
            )
            
            return(list(
              boot = pred_boot,
              test = pred_test
            ))
            
          }
        )
      
      names(validation_preds_ex3) <- run_params$name
      
      stats <- map2_dfr(validation_preds_ex3,
                        names(validation_preds_ex3),
                        function(x, name, time = 2) {
                          if (str_detect(name, "ASCVD")) {
                            outcome <- "event_ascvd_zero"
                          } else {
                            outcome <- "event_ascvd"
                          }
                          
                          boot_stats <- 
                            val.prob(p = x$boot[x$boot[["time"]] == time,][["pred"]],
                                     y = x$boot[x$boot[["time"]] == time,][[outcome]],
                                     pl = FALSE)
                          
                          test_stats <- 
                            val.prob(p = x$test[x$test[["time"]] == time,][["pred"]],
                                     y = x$test[x$test[["time"]] == time,][[outcome]],
                                     pl = FALSE)
                          
                          bind_rows(boot_stats, test_stats, .id = "test")
                          
                        }, .id = "model")
      
      return(stats)
    },
    run_params = run_params,
    .id = "bsim",
    .options = furrr_options(
      seed = TRUE
    )
  )
  
})

plan(sequential)

# option 2 ----------------------------------------------------------------

# always start from same exam (e.g. exam 2) and iteratively add previous exam
# information to obtain updated predictions -- here the advantage is that we are
# comparing predictions over the same time horizon and observed cases/deaths,
# etc.

# parametric models w/ 2 lags:
#   iter 0 - conditional on exam 2 alone project 2 through 5
#   iter 1 - conditional on exams 1 to 2 project 2 through 5

run_params <- list(
  name = list("PCE, net risk",
              "PCE, ASCVD-specific risk",
              "PCE (calibrated), net risk",
              "PCE (calibrated), ASCVD-specific risk",
              "FRS, net risk",
              "FRS, ASCVD-specific risk",
              "FRS (calibrated), net risk",
              "FRS (calibrated), ASCVD-specific risk"),
  def = "TEST",
  data = list(
    observed_events_wide
  ),
  id = 3,
  start = 0,
  stop = 2,
  sims = 100
)

# plan(multisession, workers = 10)

with_progress({
  p <- progressor(steps = length(run_params$name) * BOOTS)
  
  boot_stats_accaha <- future_map_dfr(
    1:BOOTS,
    function(bsim, run_params) {
      
      stats <- 
        pmap_dfr(
          run_params,
          function(name, def, data, id, start, stop, sims, time = 2) {
            boot <- draw_bootstrap_samples(data, "pid", "time")
            oob <- data[!data[['pid']] %in% boot[['pid']], ]
            
            outcome <-
              ifelse(str_detect(name, "ASCVD"),
                     "event_ascvd_zero_2",
                     "event_ascvd_2")
            variable <- 
              ifelse(str_detect(name, "PCE"),
                     "ascvd_10yr_accaha_risk",
                     "ascvd_10yr_frs_risk")
            
            calibrated <- 
              ifelse(str_detect(name, "calibrated"),
                     "",
                     "_updated")
            
            variable <- paste0(variable, calibrated)
            
            boot_stats <- 
              val.prob(p = boot[[paste0(variable, "_0")]],
                       y = boot[[outcome]],
                       pl = FALSE)
            
            test_stats <- 
              val.prob(p = data[[paste0(variable, "_0")]],
                       y = data[[outcome]],
                       pl = FALSE)
            
            bind_rows(boot_stats, test_stats, .id = "test")
            
            
          }, .id = "model")
      p()
      
      stats$model <- case_when(
        stats$model == 1 ~ "PCE, net risk",
        stats$model == 2 ~ "PCE, ASCVD-specific risk",
        stats$model == 3 ~ "PCE (calibrated), net risk",
        stats$model == 4 ~ "PCE (calibrated), ASCVD-specific risk",
        stats$model == 5 ~ "FRS, net risk",
        stats$model == 6 ~ "FRS, ASCVD-specific risk",
        stats$model == 7 ~ "FRS (calibrated), net risk",
        stats$model == 8 ~ "FRS (calibrated), ASCVD-specific risk"
      )
      return(stats)

    },
    run_params = run_params,
    .id = "bsim",
    .options = furrr_options(
      seed = TRUE
    )
  )
  
})

# plan(sequential)



# conventional ------------------------------------------------------------

run_params <- list(
  name = list(
    "logistic, net risk",
    "coxph, net risk",
    "logistic, ASCVD-specific risk",
    "coxph, ASCVD-specific risk"
  ),
  data = list(
    offspring_wide,
    offspring_coxph,
    offspring_wide,
    offspring_coxph
  )
)

# plan(multisession, workers = 10)

with_progress({
  p <- progressor(steps = length(run_params$name) * BOOTS)
  
  boot_stats_conventional <- future_map_dfr(
    1:BOOTS,
    function(bsim, run_params) {
      validation_preds_conventional <- 
        pmap(
          run_params,
          function(name, def, data) {
            boot <- draw_bootstrap_samples(data, "pid", "time")
            oob <- data[!data[['pid']] %in% boot[['pid']], ]
            
            if (str_detect(name, "coxph")) {
              fit <- coxph(
                Surv(enddate, event_ascvd) ~ strata(sex) + log(age_0) + I(log(age_0)^2) + smk_0 + dm_0 + 
                  hrx_0 + log(tc_0) + log(hdl_0) + log(sbp_0) + log(age_0):log(tc_0) +
                  log(age_0):log(hdl_0) + log(age_0):smk_0 + log(age_0):log(sbp_0) + log(age_0):hrx_0 + 
                  log(age_0):hrx_0:log(sbp_0) + hrx_0:log(sbp_0),
                data = boot
              )
              boot$pred <- 1 - predict(fit, newdata = mutate(boot, enddate = 3652), type = "survival")
              data$pred <- 1 - predict(fit, newdata = mutate(data, enddate = 3652), type = "survival")
            } else {
              fit <- glm(
                formula = event_ascvd_2 ~ sex * (
                  poly(log(age_0), 2) + smk_0 + dm_0 + hrx_0 +
                    log(tc_0) + log(hdl_0) + log(sbp_0) + log(age_0):log(tc_0) +
                    log(age_0):log(hdl_0) + log(age_0):smk_0 + log(age_0):log(sbp_0) + log(age_0):hrx_0 +
                    log(age_0):hrx_0:log(sbp_0) + hrx_0:log(sbp_0)
                ), 
                family = binomial(link = 'logit'),
                data = boot
              )
              
              boot$pred <- predict(fit, newdata = boot, type = "response")
              data$pred <- predict(fit, newdata = data, type = "response")
            }
            
            
            pred_boot <- left_join(
              select(boot, pid, pred),
              observed_events_wide,
              by = c("pid")
            )
            
            pred_test <- left_join(
              select(data, pid, pred),
              observed_events_wide,
              by = c("pid")
            )
            
            return(list(
              boot = pred_boot,
              test = pred_test
            ))
            
          }
        )
      
      names(validation_preds_conventional) <- run_params$name
      
      stats <- map2_dfr(validation_preds_conventional,
                        names(validation_preds_conventional),
                        function(x, name, time = 2) {
                          outcome <-
                            ifelse(str_detect(name, "ASCVD"),
                                   "event_ascvd_zero_2",
                                   "event_ascvd_2")
                          
                          boot_stats <- 
                            val.prob(p = x$boot[["pred"]],
                                     y = x$boot[[outcome]],
                                     pl = FALSE)
                          
                          test_stats <- 
                            val.prob(p = x$test[["pred"]],
                                     y = x$test[[outcome]],
                                     pl = FALSE)
                          
                          bind_rows(boot_stats, test_stats, .id = "test")
                          
                        }, .id = "model")
      
      return(stats)
    },
    run_params = run_params,
    .id = "bsim",
    .options = furrr_options(
      seed = TRUE
    )
  )
  
})

# plan(sequential)

# write_rds(boot_stats_ex1, "../../3_data_encrypted/rds/boot_stats_ex1.rds")
# write_rds(boot_stats_ex2, "../../3_data_encrypted/rds/boot_stats_ex2.rds")
# write_rds(boot_stats_ex3, "../../3_data_encrypted/rds/boot_stats_ex3.rds")
# write_rds(boot_stats_accaha, "../../3_data_encrypted/rds/boot_stats_accaha.rds")
# write_rds(boot_stats_conventional, "../../3_data_encrypted/rds/boot_stats_conventional.rds")



