# define simulation function ----------------------------------------------

run_simulation <- function(N, k, alpha, beta, gamma, sigma, pb = NULL) {
 
  train <- datagen(N, k, alpha, beta, gamma, sigma)
  test <- datagen(N, k, alpha, beta, gamma, sigma)
  
  # update progress bar if exists
  if (!is.null(pb)) {
    i <- getTxtProgressBar(pb)
    setTxtProgressBar(pb, ifelse(is.na(i), 1, i + 1))
  }
  
  # fit models for g-formula
  Y.fit_g <- 
    fit_Y_model(
      model = list(
        formula = Y ~ X + W,
        link = "logit",
        family = "binomial"
      ),
      data = train$long
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
  
  stats <-
    map(0:(k - 2), function(x) {
      
      # run g-formula
      Y.hat_g <-
        gformula_mc(
          Y.fit_g,
          X.fit_g,
          data = filter(test$long, !is.na(X) | time > x),
          id = "id",
          time = "time",
          mc.start = x,
          mc.stop = k,
          hist.vars = c("lag1_X", "lag1_W"),
          hist.fun = "lag",
          mc.sims = 500
        )
      
      # fit conventional models
      Y.fit_c_base <-
        glm(
          formula = reformulate(paste0(c("X", "W"), 0), "Y3"),
          data = train$wide,
          family = binomial(link = "logit")
        )
      
      Y.fit_c_cond <-
        glm(
          formula = reformulate(paste0(c("X", "W"), 0:x), "Y3"),
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
        cbind(filter(test$long, !is.na(X) | time > x), 'pred_g' = Y.hat_g)
      
      preds <-
        left_join(select(test$wide, id, pred_c_base, pred_c_cond, Y3),
                  filter(preds, time == 3),
                  by = "id")
      
      stats_g <- rms::val.prob(preds$pred_g, preds$Y3, pl = FALSE)
      stats_c_base <- rms::val.prob(preds$pred_c_base, preds$Y3, pl = FALSE)
      stats_c_cond <- rms::val.prob(preds$pred_c_cond, preds$Y3, pl = FALSE)
      
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
      N = rep(5000, SIMS)
    )
  
  pb <- txtProgressBar(max = nrow(sim_params), initial = NA, style = 3)
  
  risk_factor_results <-
    pmap(
      as.list(sim_params),
      run_simulation,
      # These below are constant across sims:
      k = 4,
      alpha = c(0, 2, -2),
      beta = c(log(0.25), log(3), log(3), log(2)),
      gamma = c(log(0.25), log(2), log(2)),
      sigma = 1,
      pb = pb
    )
  
  close(pb)
  
  write_rds(risk_factor_results, get_data("risk_factor_results"))
  
} else {
  
  risk_factor_results <- read_rds(get_data("risk_factor_results"))
}

