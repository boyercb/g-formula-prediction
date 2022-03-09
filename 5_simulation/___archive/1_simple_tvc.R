# define simulation function ----------------------------------------------

run_simulation <- function(N,
                           a0, a1, a2, a3, a4,
                           b0, b1, b2, b3, 
                           c0, c1, c2, c3, c4,
                           sigma, pb) {
  
  # draw from data generation process
  train <- datagen(
    N = N, 
    a0 = a0,
    a1 = a1,
    a2 = a2,
    a3 = a3,
    a4 = a4,
    b0 = b0,
    b1 = b1,
    b2 = b2,
    b3 = b3,
    c0 = c0,
    c1 = c1,
    c2 = c2,
    c3 = c3,
    c4 = c4, 
    sigma = sigma,
    output = "both"
  )
  
  # draw from data generation process
  test <- datagen(
    N = N, 
    a0 = a0,
    a1 = a1,
    a2 = a2,
    a3 = a3,
    a4 = a4,
    b0 = b0,
    b1 = b1,
    b2 = b2,
    b3 = b3,
    c0 = c0,
    c1 = c1,
    c2 = c2,
    c3 = c3,
    c4 = c4, 
    sigma = sigma,
    output = "both"
  )
  
  # draw from data generation process
  test_no_treat <- datagen(
    N = N, 
    a0 = a0,
    a1 = a1,
    a2 = a2,
    a3 = a3,
    a4 = a4,
    b0 = -Inf,
    b1 = 0,
    b2 = 0,
    b3 = 0,
    c0 = c0,
    c1 = c1,
    c2 = c2,
    c3 = c3,
    c4 = c4, 
    sigma = sigma,
    output = "both"
  )

  # draw from data generation process
  test_const_treat <- datagen(
    N = N, 
    a0 = a0,
    a1 = a1,
    a2 = a2,
    a3 = a3,
    a4 = a4,
    b0 = 0,
    b1 = 0,
    b2 = 0,
    b3 = Inf,
    c0 = c0,
    c1 = c1,
    c2 = c2,
    c3 = c3,
    c4 = c4, 
    sigma = sigma,
    output = "both"
  )
  
  # draw from data generation process
  test_emod <- datagen(
    N = N, 
    a0 = a0,
    a1 = a1,
    a2 = a2,
    a3 = a3,
    a4 = a4,
    b0 = b0,
    b1 = b1,
    b2 = b2,
    b3 = b3,
    c0 = c0,
    c1 = c1,
    c2 = c2,
    c3 = c3,
    c4 = c4, 
    c5 = log(0.25),
    sigma = sigma,
    output = "both"
  )
  
  # update progress bar if exists
  if (!is.null(pb)) {
    i <- getTxtProgressBar(pb)
    setTxtProgressBar(pb, ifelse(is.na(i), 1, i + 1))
  }
  
  # define models 
  X.models <-
    list(
      "L" = list(
        formula = L ~ lag1_L + lag1_A, 
        family = "normal"
      ),
      "A" = list(
        formula = A ~ L + lag1_L + lag1_A, 
        link = "logit",
        family = "binomial"
      )
    )
  
  Y.model <- 
    list(
      formula = Y ~ L + A + lag1_L + lag1_A,
      link = "logit",
      family = "binomial"
    )

  # fit models for g-formula
  Y.fit <- 
    fit_Y_model(
      Y.model,
      data = train$long
    )
  
  X.fit <-
    fit_X_models(
      X.models,
      data = train$long,
      time = "time"
    )
  
  # fit conventional model
  Y.fit_c <-
    glm(Y1 ~ A0 + L0, data = train$wide, family = binomial(link = "logit"))
  
  # fit conventional model w/out A0
  Y.fit_c_noA <-
    glm(Y1 ~ L0, data = train$wide, family = binomial(link = "logit"))
  
  # run g-formula
  Y.hat_g <-
    gformula_mc(
      Y.fit,
      X.fit,
      data = test$long,
      id = "id",
      time = "time",
      hist.vars = c("lag1_A", "lag1_L"),
      hist.fun = "lag",
      mc.sims = 50
    )
  
  Y.hat_g_no_treat <-
    gformula_mc(
      Y.fit,
      X.fit,
      data = test_no_treat$long,
      id = "id",
      time = "time",
      treatment = "A",
      intervention = function(data, x) { 0 },
      hist.vars = c("lag1_A", "lag1_L"),
      hist.fun = "lag",
      mc.sims = 50
    )
  
  Y.hat_g_const_treat <-
    gformula_mc(
      Y.fit,
      X.fit,
      data = test_const_treat$long,
      id = "id",
      time = "time",
      treatment = "A",
      intervention = function(data, x) { data[["lag1_A"]] },
      hist.vars = c("lag1_A", "lag1_L"),
      hist.fun = "lag",
      mc.sims = 50
    )
  
  Y.hat_g_emod <-
    gformula_mc(
      Y.fit,
      X.fit,
      data = test_emod$long,
      id = "id",
      time = "time",
      hist.vars = c("lag1_A", "lag1_L"),
      hist.fun = "lag",
      mc.sims = 50
    )
  
  Y.hat_c <-
    predict(Y.fit_c, type = 'response', newdata = test$wide)

  Y.hat_c_no_treat <-
    predict(Y.fit_c, type = 'response', newdata = test_no_treat$wide)
  
  Y.hat_c_const_treat <-
    predict(Y.fit_c, type = 'response', newdata = test_const_treat$wide)
  
  Y.hat_c_emod <-
    predict(Y.fit_c, type = 'response', newdata = test_emod$wide)
  
  Y.hat_c_noA <-
    predict(Y.fit_c_noA, type = 'response', newdata = test$wide)
  
  Y.hat_c_noA_no_treat <-
    predict(Y.fit_c_noA, type = 'response', newdata = test_no_treat$wide)
  
  Y.hat_c_noA_const_treat <-
    predict(Y.fit_c_noA, type = 'response', newdata = test_const_treat$wide)
  
  Y.hat_c_noA_emod <-
    predict(Y.fit_c_noA, type = 'response', newdata = test_emod$wide)
  
  test$wide$pred_c <- Y.hat_c
  test_no_treat$wide$pred_c <- Y.hat_c_no_treat
  test_const_treat$wide$pred_c <- Y.hat_c_const_treat
  test_emod$wide$pred_c <- Y.hat_c_emod
  
  test$wide$pred_c_noA <- Y.hat_c_noA
  test_no_treat$wide$pred_c_noA <- Y.hat_c_noA_no_treat
  test_const_treat$wide$pred_c_noA <- Y.hat_c_noA_const_treat
  test_emod$wide$pred_c_noA <- Y.hat_c_noA_emod
  
  # calculate calibration and validation stats
  preds <- 
    cbind(test$long, 'pred_g' = Y.hat_g) 
  
  preds <-
    left_join(select(test$wide, id, pred_c, pred_c_noA, Y1),
              filter(preds, time == 1),
              by = "id") 
  
  stats_g <- rms::val.prob(preds$pred_g, preds$Y1, pl = FALSE)
  stats_c <- rms::val.prob(preds$pred_c, preds$Y1, pl = FALSE)
  stats_c_noA <- rms::val.prob(preds$pred_c_noA, preds$Y1, pl = FALSE)
  
  preds_no_treat <- 
    cbind(test_no_treat$long, 'pred_g' = Y.hat_g_no_treat) 
  
  preds_no_treat <-
    left_join(select(test_no_treat$wide, id, pred_c, pred_c_noA, Y1),
              filter(preds_no_treat, time == 1),
              by = "id")
  
  stats_g_no_treat <- rms::val.prob(preds_no_treat$pred_g, preds_no_treat$Y1, pl = FALSE)
  stats_c_no_treat <- rms::val.prob(preds_no_treat$pred_c, preds_no_treat$Y1, pl = FALSE)
  stats_c_noA_no_treat <- rms::val.prob(preds_no_treat$pred_c_noA, preds_no_treat$Y1, pl = FALSE)
  
  preds_const_treat <- 
    cbind(test_const_treat$long, 'pred_g' = Y.hat_g_const_treat) 
  
  preds_const_treat <-
    left_join(select(test_const_treat$wide, id, pred_c, pred_c_noA, Y1),
              filter(preds_const_treat, time == 1),
              by = "id")
  
  stats_g_const_treat <- rms::val.prob(preds_const_treat$pred_g, preds_const_treat$Y1, pl = FALSE)
  stats_c_const_treat <- rms::val.prob(preds_const_treat$pred_c, preds_const_treat$Y1, pl = FALSE)
  stats_c_noA_const_treat <- rms::val.prob(preds_const_treat$pred_c_noA, preds_const_treat$Y1, pl = FALSE)
  
  preds_emod <- 
    cbind(test_emod$long, 'pred_g' = Y.hat_g_emod) 
  
  preds_emod <-
    left_join(select(test_emod$wide, id, pred_c, pred_c_noA, Y1),
              filter(preds_emod, time == 1),
              by = "id")
  
  stats_g_emod <- rms::val.prob(preds_emod$pred_g, preds_emod$Y1, pl = FALSE)
  stats_c_emod <- rms::val.prob(preds_emod$pred_c, preds_emod$Y1, pl = FALSE)
  stats_c_noA_emod <- rms::val.prob(preds_emod$pred_c_noA, preds_emod$Y1, pl = FALSE)
  
  
  stats <-
    bind_rows(
      stats_g,
      stats_c, 
      stats_c_noA,
      stats_g_no_treat, 
      stats_c_no_treat, 
      stats_c_noA_no_treat,
      stats_g_const_treat, 
      stats_c_const_treat, 
      stats_c_noA_const_treat,
      stats_g_emod, 
      stats_c_emod, 
      stats_c_noA_emod,
      .id = "id"
    ) %>%
    mutate(
      model = case_when(
        id %in% c(1, 4, 7, 10) ~ "g-formula",
        id %in% c(2, 5, 8, 11) ~ "conventional",
        id %in% c(3, 6, 9, 12) ~ "conventional (no A0)"
      ),
      test_data = case_when(
        id %in% 1:3 ~ "natural course",
        id %in% 4:6 ~ "no treatment",
        id %in% 7:9 ~ "constant treatment",
        id %in% 10:12 ~ "effect modification"
      )
    )
  
  stats <- mutate(
    stats,
    a2 = a2
  )

  return(stats)
}


# run simulation ----------------------------------------------------------

if (rerun_simulation) {
  
SIMS <- 500

sim_params <- 
  expand.grid(
    N = rep(5000, SIMS),
    a2 = c(-3, -2.5, -2, -1.5, -1, -0.5, 0)
  )

pb <- txtProgressBar(max = nrow(sim_params), initial = NA, style = 3)

simple_tvc_results <-
  pmap(
    as.list(sim_params),
    run_simulation,
    # These below are constant across sims:
    a0 = 0,
    a1 = 1,
    a3 = 0,
    a4 = 0,
    b0 = log(0.25),
    b1 = log(2),
    b2 = 0,
    b3 = log(2),
    c0 = log(0.25),
    c1 = log(1.5),
    c2 = log(0.5),
    c3 = log(1.5),
    c4 = log(0.5),
    sigma = 1,
    pb = pb
  )

close(pb)

write_rds(simple_tvc_results, get_data("simple_tvc_results"))

} else {
  
  simple_tvc_results <- read_rds(get_data("simple_tvc_results"))
}



