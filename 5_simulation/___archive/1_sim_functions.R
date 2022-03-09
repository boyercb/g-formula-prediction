# function to run monte carlo simulations
run_simulation <- 
  function(n,
           k,
           a,
           b,
           c = NULL,
           d = NULL,
           e,
           treatment = NULL,
           shift = NULL,
           sigma = NULL,
           s0 = NULL,
           n.sims = 1000,
           mc.sims = 100,
           X.models,
           Y.model,
           D.model = NULL,
           base.model,
           updated.model
           ) {
    
    
    p <- progressor(steps = n.sims)
    
    future_map_dfr(
      1:n.sims,
      function(sim) {
        # draw training sample from SEM 
        train <- generate_sem_data(
          n = n,
          k = k,
          a = a,
          b = b,
          c = c,
          d = d,
          e = e,
          sigma = sigma,
          s0 = s0,
          long = TRUE
        ) 
        
        # draw testing sample from SEM 
        test <- generate_sem_data(
          n = n,
          k = k,
          a = a,
          b = b,
          c = c,
          d = d,
          e = e,
          treatment = treatment,
          shift = shift,
          sigma = sigma,
          s0 = s0,
          long = TRUE
        ) 
        
        # draw estimates from sample
        estimates <- draw_estimates(
          train = train,
          test = test,
          k = k,
          treatment = treatment,
          shift = shift,
          mc.sims = mc.sims,
          X.models = X.models,
          Y.model = Y.model,
          D.model = D.model,
          base.model,
          updated.model
        )
        
        # update progress bar
        p()
        return(estimates)
      },
      .id = "sim",
      .options = furrr_options(
        seed = TRUE
      )
    )
  }

# helper function which draws estimates
draw_estimates <- function(
  train,
  test,
  k,
  treatment = NULL,
  shift = NULL,
  mc.sims,
  X.models,
  Y.model,
  D.model = NULL,
  base.model,
  updated.model
) {
  
  # fit the models for joint distribution
  X.fit <- fit_X_models(X.models, data = train$long, time = "time")
  Y.fit <- fit_Y_model(Y.model, data = train$long)
  
  if (!is.null(D.model)) {
    D.fit <- fit_D_model(D.model, data = train$long)
  } else {
    D.fit <- NULL
  }
  
  # fit conventional (base) model
  base.fit <- fit_Y_model(base.model, data = train$wide)
  
  # fit conventional (updated) models
  updated.fit <- lapply(updated.model, function (x) {
    fit_Y_model(x, data = train$wide)
  })
  
  # iteratively condition on time steps between 0 and k - 1
  estimates_gformula <- map_dfr(
    0:(k-1),
    function(x) {
      # run monte carlo sim to get g-formula probability estimates
      prob_gformula <- gformula_mc(
        Y.fit = Y.fit,
        X.fit = X.fit,
        D.fit = D.fit,
        data = subset(test$long, !is.na(X)),
        id = "id",
        time = "time",
        hist.vars = names(test$long)[grepl("lag", names(test$long))],
        hist.fun = "lag",
        treatment = treatment,
        intervention = shift,
        mc.sims = mc.sims,
        mc.start = x,
        mc.stop = k-1,
        merge = FALSE,
        last_only = TRUE
      )
      
      prob_gformula <- left_join(prob_gformula, test$wide, by = "id")
      
      # estimate relevant summary statistics
      rms::val.prob(
        p = prob_gformula[["pred"]],
        y = prob_gformula[[paste0("Y", k-1)]], 
        pl = FALSE
      )
    },
    .id = "cond"
  )
  
  # get conventional probability estimates
  prob_base <- predict(base.fit, newdata = test$wide, type = "response")
  
  # estimate relevant summary statistics
  estimates_base <- rms::val.prob(
    p = prob_base, 
    y = test$wide[[paste0("Y", k-1)]], 
    pl = FALSE
  )
  
  estimates_updated <- map2_dfr(
    0:(k-1),
    updated.fit,
    function(x, fit) {
      # get (updated) conventional probability estimates
      prob_updated <- predict(fit, newdata = test$wide, type = "response")
      
      # estimate relevant summary statistics
      rms::val.prob(
        p = prob_updated,
        y = test$wide[[paste0("Y", k-1)]], 
        pl = FALSE
      )
    },
    .id = "cond"
  )
  
  estimates <-
    bind_rows(estimates_gformula, estimates_base, estimates_updated, .id = "id")
  
  # create model variable
  estimates$model <- case_when(
    estimates$id == 1 ~ "g-formula",
    estimates$id == 2 ~ "conventional (baseline-only)",
    estimates$id == 3 ~ "conventional (updated)"
  )
  
  return(estimates)
}


# tests
if (FALSE) {
  # testing the draw_estimates function
  train <- generate_sem_data(
    n = 2000,
    k = 4,
    #          Int,  X, lag1_X, lag2_X, lag3_X, W, lag1_W, lag2_W, lag3_W
    a = c(       0,  0,      1,      0,      0, 0,   -0.5,      0,      0),
    b = c(       0,  1,    0.5,      0,      0, 0,      1,      0,      0),
    e = c(log(0.1),  1,   0.75,   0.50,   0.25,-1,  -0.75,  -0.50,  -0.25),
    sigma = 0.75,
    s0 = 1,
    long = TRUE
  )
  
  test <- generate_sem_data(
    n = 2000,
    k = 4,
    #          Int,  X, lag1_X, lag2_X, lag3_X, W, lag1_W, lag2_W, lag3_W
    a = c(       0,  0,      1,      0,      0, 0,   -0.5,      0,      0),
    b = c(       0,  1,    0.5,      0,      0, 0,      1,      0,      0),
    e = c(log(0.1),  1,   0.75,   0.50,   0.25,-1,  -0.75,  -0.50,  -0.25),
    sigma = 0.75,
    s0 = 1,
    long = TRUE
  )
  
  X.models <- 
    list(
      "X" = list(
        formula = X ~ lag1_X + lag1_W,
        family = "normal"
      ),
      "W" = list(
        formula = W ~ X + lag1_X + lag1_W, 
        link = "logit",
        family = "binomial"
      )
    )
  
  Y.model <- list(
    formula = Y ~ (X + W) * as.factor(time) + lag1_X + lag1_W + lag2_X + lag2_W + lag3_X + lag3_W,
    link = "logit",
    family = "binomial"
  )
  
  BASE.MODEL <- list(
    formula = Y3 ~ X0 + W0,
    link = "logit",
    family = "binomial"
  )
  
  UPDATED.MODEL <- list(
    "base only" = list(
      formula = Y3 ~ X0 + W0,
      link = "logit",
      family = "binomial"
    ),
    "time 0 to 1" = list(
      formula = Y3 ~ X0 + W0 + X1 + W1,
      link = "logit",
      family = "binomial"
    ),
    "time 0 to 2" = list(
      formula = Y3 ~ X0 + W0 + X1 + W1 + X2 + W2,
      link = "logit",
      family = "binomial"
    ),
    "time 0 to 3" = list(
      formula = Y3 ~ X0 + W0 + X1 + W1 + X2 + W2 + X3 + W3,
      link = "logit",
      family = "binomial"
    )
  )
  
  # X.fit <- fit_X_models(X.models, data = sample$long, time = "time")
  # Y.fit <- fit_Y_model(Y.model, data = sample$long)
  # alt.fit <- fit_Y_model(alt.model, data = sample$wide)
  df <- draw_estimates(
    train = train,
    test = test,
    k = 4,
    mc.sims = 100,
    X.models = X.models,
    Y.model = Y.model,
    base.model = BASE.MODEL,
    updated.model = UPDATED.MODEL
  )
  
  # prob <- gformula_mc(
  #   Y.fit = Y.fit,
  #   X.fit = X.fit,
  #   data = subset(sample$long, !is.na(X)),
  #   id = "id",
  #   time = "time",
  #   hist.vars = names(sample$long)[grepl("lag", names(sample$long))],
  #   hist.fun = "lag",
  #   mc.sims = 100,
  #   mc.start = 0,
  #   mc.stop = 3,
  #   merge = FALSE,
  #   last_only = TRUE
  # )
  # 
  # prob <- left_join(prob, sample$wide, by = "id") 
  # 
  # rms::val.prob(
  #   p = prob[["pred"]],
  #   y = prob[["Y3"]],
  #   pl = TRUE
  # )
  # 
  # # get conventional probability estimates
  # prob_alt <- predict(alt.fit, type = "response")
  # 
  # # estimate relevant summary statistics
  # estimates_alt <-  rms::val.prob(
  #   p = prob_alt, 
  #   y = sample$wide[["Y3"]], 
  #   pl = TRUE,
  #   riskdist = "predicted"
  # )
  # 
  # ggplot(prob, aes(x = pred, y = Y1)) + 
  #   geom_smooth() +
  #   geom_abline()
  # 
  # lm(Y1 ~ pred, data = prob)
  # lm(Y1 ~ pred, data = cbind(sample$wide, "pred" = prob_alt))
  # 
  # df <- sample$long %>%
  #   filter(!is.na(X)) %>%
  #   mutate(
  #     pred = predict(Y.fit, newdata = ., type = "response")
  #   ) %>%
  #   group_by(id) %>%
  #   mutate(
  #     cumpred = cumsum(pred * cumprod(lag(1 - pred, default = 1)))
  #   ) %>%
  #   ungroup()
  # 
  # df <- left_join(df, sample$wide, by = "id")
  # 
  # rms::val.prob(
  #   p = filter(df, time == 0) %>% pull(pred),
  #   y = filter(df, time == 0) %>% pull(Y0),
  #   pl = TRUE
  # )
  # 
  # rms::val.prob(
  #   p = filter(df, time == 1) %>% pull(pred),
  #   y = filter(df, time == 1) %>% pull(Y1),
  #   pl = TRUE
  # )
  # 
  # rms::val.prob(
  #   p = filter(df, time == 2) %>% pull(pred),
  #   y = filter(df, time == 2) %>% pull(Y2),
  #   pl = TRUE
  # )
  # 
  # rms::val.prob(
  #   p = filter(df, time == 3) %>% pull(pred),
  #   y = filter(df, time == 3) %>% pull(Y3),
  #   pl = TRUE
  # )
}