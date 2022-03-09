# define simulation function ----------------------------------------------

run_simulation <- function(N, k, alpha, beta, gamma, delta, eta, sigma, mc.sims, pb = NULL) {
 
  train <- datagen(N, k, alpha, beta, gamma, delta, eta, sigma)
  test <- datagen(N, k, alpha, beta, gamma, delta, eta, sigma)
  
  # update progress bar if exists
  if (!is.null(pb)) {
    i <- getTxtProgressBar(pb)
    setTxtProgressBar(pb, ifelse(is.na(i), 1, i + 1))
  }
  
  # fit models for g-formula
  Y.fit_g <- 
    fit_Y_model(
      model = list(
        formula = Y ~ X + W + lag1_X + lag1_W,
        link = "logit",
        family = "binomial"
      ),
      data = filter(train$long, D != 1)
    )
  
  X.fit_g <-
    fit_X_models(
        models = list(
        "X" = list(
          formula = X ~ lag1_X + lag1_W, 
          family = "normal"
        ),
        "W" = list(
          formula = W ~ X + lag1_X + lag1_W, 
          link = "logit",
          family = "binomial"
        )
      ),
      data = train$long,
      time = "time"
    )
  
  if (!is.null(delta)) {
    D.fit_g <-
      fit_Y_model(
        model = list(
          formula = D ~ X + W,
          link = "logit",
          family = "binomial"
        ),
        data = train$long
      )
  } else {
    D.fit_g <- NULL
  }
  
  stats <-
    map(0:(k - 2), function(x) {

      # run g-formula
      Y.hat_g <-
        gformula_mc(
          Y.fit_g,
          X.fit_g,
          D.fit_g,
          data = filter(test$long, time >= x & !is.na(X)),
          id = "id",
          time = "time",
          mc.start = x,
          mc.stop = k - 1,
          hist.vars = c("lag1_X", "lag1_W"),
          hist.fun = "lag",
          mc.sims = mc.sims,
          merge = FALSE
        )
      
      # fit conventional models
      Y.fit_c_base <-
        glm(
          formula = reformulate(paste0(c("X", "W"), 0), paste0("Y", k - 1)),
          data = train$wide,
          family = binomial(link = "logit")
        )
      
      Y.fit_c_cond <-
        glm(
          formula = reformulate(paste0(c("X", "W"), 0:x), paste0("Y", k - 1)),
          data = train$wide,
          family = binomial(link = "logit")
        )
      
      # get predictions
      Y.fit_c_base <-
        predict(Y.fit_c_base, type = 'response', newdata = test$wide)
      
      Y.fit_c_cond <-
        predict(Y.fit_c_cond, type = 'response', newdata = test$wide)
      
      test$wide$pred_c_base <- Y.fit_c_base
      test$wide$pred_c_cond <- Y.fit_c_cond
      
      # calculate calibration and validation stats
      preds <-
        left_join(Y.hat_g,
                  filter(test$long, time >= x & !is.na(X)),
                  by = c("id", "time")) %>%
        rename(pred_g = pred)
      
      preds <-
        left_join(select(test$wide, id, pred_c_base, pred_c_cond, paste0("Y", k - 1)),
                  filter(preds, time == k - 1),
                  by = "id")
      
      stats_g <-
        rms::val.prob(preds$pred_g, preds[, paste0("Y", k - 1)], pl = FALSE)
      stats_c_base <-
        rms::val.prob(preds$pred_c_base, preds[, paste0("Y", k - 1)], pl = FALSE)
      stats_c_cond <-
        rms::val.prob(preds$pred_c_cond, preds[, paste0("Y", k - 1)], pl = FALSE)
      
      stats <-
        bind_rows(stats_g,
                  stats_c_base,
                  stats_c_cond,
                  .id = "id") %>%
        mutate(model = case_when(id %in% c(1) ~ "g-formula",
                                 id %in% c(2) ~ "conventional (baseline-only)",
                                 id %in% c(3) ~ "conventional (updated)")
               )
      
      return(stats)
    })
  
  return(bind_rows(stats, .id = "num"))
}


# run simulation ----------------------------------------------------------

if (rerun_simulation) {

  SIMS <- 500
  
  sim_params <- 
    expand.grid(
      N = rep(2500, SIMS)
    )
  
  # plan for parallel multiprocess
  future::plan(future::multiprocess, workers = 4)
  
  # run simulations for full parameter grid
  tictoc::tic()
  
  
  risk_factor_results_noDC <-
    furrr::future_pmap(
      as.list(sim_params),
      run_simulation,
      # These below are constant across sims:
      k = 5,
      alpha = c(0, 2, -1.5),
      beta = c(log(0.1), log(1.5), log(0.5), log(1.5)),
      gamma = NULL,
      delta = NULL,
      eta = c(log(0.1), -log(1.5), log(1.5), -log(1), log(1)),
      sigma = 1,
      mc.sims = 50,
      .progress = TRUE
    )
  tictoc::toc()
  
  # plan for parallel multiprocess
  future::plan(future::multiprocess, workers = 4)
  
  # run simulations for full parameter grid
  tictoc::tic()
  
  
  risk_factor_results <-
    furrr::future_pmap(
      as.list(sim_params),
      run_simulation,
      # These below are constant across sims:
      k = 5,
      alpha = c(0, 1.2, -1.5),
      beta = c(log(0.1), log(1.5), log(1.2), log(1.2)),
      gamma = c(log(0.1), log(1.2), log(1.2)),
      delta = c(log(0.03), log(1.2), log(0.85)),
      eta = c(log(0.05), log(1.2), log(0.5), log(1.1), log(0.85)),
      sigma = 1,
      mc.sims = 50,
      .progress = TRUE
    )
  tictoc::toc()
  
  write_rds(risk_factor_results, get_data("risk_factor_results"))
  
} else {
  
  risk_factor_results <- read_rds(get_data("risk_factor_results"))
}

