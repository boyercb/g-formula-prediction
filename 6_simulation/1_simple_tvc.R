# define simulation function ----------------------------------------------

run_simulation <- function(N,
                           a0, a1, a2, a3, a4,
                           b0, b1, b2, b3, 
                           c0, c1, c2, c3, c4,
                           sigma, pb) {
  
  # draw from data generation process
  sim <- datagen(
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
    type = "wide"
  )

  # update progress bar if exists
  if (!is.null(pb)) {
    i <- getTxtProgressBar(pb)
    setTxtProgressBar(pb, ifelse(is.na(i), 1, i + 1))
  }
  
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
  sim_long <-
    sim_long %>%
    group_by(id) %>%
    mutate(across(c("L", "A"), lag, default = 0, .names = "lag1_{col}"))
  
  # 70/30 split into training and validation sets
  train_rows <- sample(1:N, size = floor(0.7 * N))
  
  train <- sim[train_rows, ]
  test <- sim[-train_rows, ] 
  
  train_long <- left_join(train, sim_long, by = "id")
  test_long <- left_join(test, sim_long, by = "id")
  
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
      data = test_long,
      id = "id",
      time = "time",
      hist.vars = c("lag1_A", "lag1_L"),
      hist.fun = "lag",
      mc.sims = 100
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
    b0 = b0,
    b1 = b1,
    b2 = b2,
    b3 = b3,
    c0 = c0,
    c1 = c1,
    c2 = c2,
    c3 = c3,
    c4 = c4,
    sigma = sigma
  )

  return(stats)
}


# run simulation ----------------------------------------------------------

SIMS <- 1000
pb <- txtProgressBar(max = SIMS, initial = NA, style = 3)

results <-
  rerun(
    SIMS,
    run_simulation(
      N = 15000,
      a0 = 0,
      a1 = 1,
      a2 = -3,
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
  )



