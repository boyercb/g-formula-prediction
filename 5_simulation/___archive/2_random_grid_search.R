grid_search <- function(a0, a1, a2, a3, a4 = 0,
                        c0, c1, c2, c3 = 0,
                        d0, d1, d2, pb) {
  
  sim <- datagen(
    N = 2500, 
    a0 = a0,
    a1 = a1,
    a2 = a2,
    a3 = a3,
    a4 = a4,
    c0 = c0,
    c1 = c1,
    c2 = c2,
    c3 = c3,
    d0 = d0,
    d1 = d1,
    d2 = d2
  )
  
  i <- getTxtProgressBar(pb)
  setTxtProgressBar(pb, ifelse(is.na(i), 1, i + 1))
  
  # create long dataset
  sim_long <-
    sim %>% 
    pivot_longer(
      -id,
      names_to = c("variable", "time"),
      names_pattern = "([LAY])([0-9])"
    ) %>%
    pivot_wider(names_from = "variable") %>%
    mutate(time = as.numeric(time))
  
  # add lagged variables for L and A
  sim_long <- mutate(
    sim_long, 
    across(c("L", "A"), lag, default = 0, .names = "lag1_{col}")
  )
  
  # 70/30 split into training and validation sets
  train_rows <- sample(1:nrow(sim), size = floor(0.7 * nrow(sim)))
  
  train <- sim[train_rows, ]
  test <- sim[-train_rows, ] 
  
  train_long <- left_join(train, sim_long, by = "id")
  test_long <- left_join(test, sim_long, by = "id")
  
  
  # define models 
  X.models <-
    list(
      "L" = list(
        formula = L ~ lag1_L + lag1_A, 
        link = "logit",
        family = "binomial"
      ),
      "A" = list(
        formula = A ~ L + lag1_L + lag1_A, 
        link = "logit",
        family = "binomial"
      )
    )
  
  D.model <- 
    list(
      formula = D ~ L + A + as.factor(time),
      link = "logit",
      family = "binomial"
    )
  
  Y.model <- 
    list(
      formula = Y ~ L + A + as.factor(time),
      link = "logit",
      family = "binomial"
    )
  
  
  # fit models 
  Y.fit <- 
    fit_Y_model(
      Y.model,
      data = train_long
    )
  
  X.fit <- 
    fit_X_models(
      X.models,
      data = train_long,
      time = "time"
    )
  
  # run g-formula 
  Y.hat_gformula <-
    gformula_mc(
      Y.fit,
      X.fit,
      #D.fit,
      data = test_long,
      id = "id",
      time = "time",
      hist.vars = c("lag1_A", "lag1_L"),
      hist.fun = "lag",
      mc.sims = 50
    )
  
  
  # run conventional model
  Y.fit_conventional <- 
    glm(Y1 ~ A0 + L0, data = train, family = binomial(link = "logit"))
  
  Y.hat_conventional <-
    predict(Y.fit_conventional, type = 'response', newdata = test)
  
  test$pred_c <- Y.hat_conventional
  
  # calculate calibration and validation stats
  test_w_pred <-
    cbind(test_long, 'pred_g' = Y.hat_gformula)
  
  test_w_pred <- 
    left_join(select(test, id, pred_c), filter(test_w_pred, time == 1), by = "id")
  
  stats_g <- rms::val.prob(test_w_pred$pred_g, test_w_pred$Y1, pl = FALSE)
  stats_c <- rms::val.prob(test_w_pred$pred_c, test_w_pred$Y1, pl = FALSE)
  
  stats <- bind_rows(stats_g, stats_c,.id = "id") %>% 
    mutate(calibration = ifelse(id == 1, "g", "c"))
  
  stats <- mutate(
    stats,
    a0 = a0,
    a1 = a1,
    a2 = a2,
    a3 = a3,
    a4 = a4,
    c0 = c0,
    c1 = c1,
    c2 = c2,
    c3 = c3,
    d0 = d0,
    d1 = d1,
    d2 = d2
  )
  
  return(stats)
}


# run grid search ---------------------------------------------------------

parameter_grid <-
  expand.grid(
    a0 = runif(4,-2, 2),
    a1 = runif(4,-2, 2),
    a2 = runif(4,-2, 2),
    a3 = runif(4,-2, 2),
    #a4 = runif(4,23,23),
    c0 = runif(4,-2, 2),
    c1 = runif(4,-2, 2),
    c2 = runif(4,-2, 2),
    #c3 = runif(4,23,23),
    d0 = runif(4,-2, 2),
    d1 = runif(4,-2, 2),
    d2 = runif(4,-2, 2)
  ) %>%
  slice_sample(prop = 0.01)

pb <- txtProgressBar(max = nrow(parameter_grid), initial = NA, style = 3)

stats <- 
  pmap(
    c(as.list(parameter_grid)),
    grid_search,
    pb = pb
  ) %>% bind_rows()
