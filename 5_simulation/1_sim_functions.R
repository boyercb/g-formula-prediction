
# create simulation function ----------------------------------------------

run_simulation <- function(scenario, nsims) {
  
  # get parameter values from scenario object
  alpha <- scenario$alpha
  beta <- scenario$beta
  gamma <- scenario$gamma
  delta <- scenario$delta
  epsilon <- scenario$epsilon

  error <-  scenario$error
  rho <- scenario$rho
  sigma <- scenario$sigma
  
  sigma[lower.tri(sigma)] <- rho
  sigma[upper.tri(sigma)] <- rho
  
  mc.sims <- scenario$mc.sims
  sizes <- scenario$sizes
  p <- scenario$p
  k <- scenario$k
  
  treatment <- scenario$treatment
  shift <- scenario$shift
  
  specification <- scenario$specification
  
  has_net_risk <- any(str_detect(specification, "elimination risk"))
  
  if (has_net_risk) {
    # define simulation parameters
    sims <- expand_grid(
      sim = 1:nsims,
      n = sizes,
      specification = specification
    )
  } else {
    # define simulation parameters
    sims <- expand_grid(
      sim = 1:nsims,
      n = sizes,
      specification = "one"
    )
  }
 
  
  # create training and validation datasets
  sims <- 
    sims %>%
    mutate(
      train = map(
        n,
        function(n) {
          generate_sem_data(
            n = n,
            k = k,
            p = p,
            a = alpha,
            b = beta,
            c = gamma,
            d = delta,
            e = epsilon,
            error = error,
            sigma = sigma,
            long = TRUE
          )
        }),
      test = map(
        n,
        function(n) {
          generate_sem_data(
            n = max(sizes),
            k = k,
            p = p,
            a = alpha,
            b = beta,
            c = gamma,
            d = delta,
            e = epsilon,
            error = error,
            sigma = sigma,
            long = TRUE
          )
        })
    ) 
  
  # expand grid to include conditioning points
  sims <-
    sims %>%
    expand_grid(time = 0:(k - 1)) 
  
  # add true value of estimand
  sims <-
    sims %>%
    mutate(
      truth = pmap(
        list(
          test,
          time,
          specification
        ),
        function(test, time, spec) {
          calculate_sem_truth(
            n = max(sizes),
            k = k,
            p = p,
            a = alpha,
            b = beta,
            c = gamma,
            d = delta,
            e = epsilon,
            L0 = test$wide[, calc_position(
              sapply(paste0("_", 0:time), function(x) paste0("L", 1:p, x)) %>% 
                as.character(),
              test$wide)
            ],
            A0 = test$wide[, calc_position(paste0("A", 0:time), test$wide)],
            time = time,
            treatment = if (has_net_risk) {
              if (spec == "elimination risk") {
                "D"
              } else {
                treatment
              }
            } else {
              treatment
            },
            shift = if (has_net_risk) {
              if (spec == "elimination risk") {
                function(data, treatment, time) { 0 }
              } else {
                shift
              }
            } else {
              shift
            }
          )
        })
    )
  
  if (has_net_risk) {
    # expand grid to include estimators 
    sims <-
      sims %>%
      expand_grid(
        estimator = c(
          "conventional",
          "conventional (updated)",
          "gformula (pooled)",
          "gformula (unpooled)"
        )
      ) 
  } else {
    # expand grid to include estimators 
    sims <-
      sims %>%
      select(-specification) %>%
      expand_grid(
        estimator = c(
          "conventional",
          "conventional (updated)",
          "gformula (pooled)",
          "gformula (unpooled)"
        ),
        specification = specification
      ) 
  }
  
  # defined model and prediction functions for each estimator
  sims <-
    sims %>%
    mutate(
      model = map(estimator, function(x) {
        switch(
          x,
          "conventional" = function(data, time, spec) {
            switch(
              as.character(time),
              "0" = glm(
                formula = create_formula("conventional", "Y3", spec, time),
                family = binomial(link = "logit"),
                data = as.data.frame(data$wide)
              ),
              "1" = glm(
                formula = create_formula("conventional", "Y3", spec, time),
                family = binomial(link = "logit"),
                data = as.data.frame(data$wide)
              ),
              "2" = glm(
                formula = create_formula("conventional", "Y3", spec, time),
                family = binomial(link = "logit"),
                data = as.data.frame(data$wide)
              ),
              "3" = glm(
                formula = create_formula("conventional", "Y3", spec, time),
                family = binomial(link = "logit"),
                data = as.data.frame(data$wide)
              )
            ) 
          },
          "conventional (updated)" = function(data, time, spec) {
            switch(
              as.character(time),
              "0" = glm(
                formula = create_formula("conventional (updated)", "Y3", spec, time),
                family = binomial(link = "logit"),
                data = as.data.frame(data$wide)
              ),
              "1" = glm(
                formula = create_formula("conventional (updated)", "Y3", spec, time),
                family = binomial(link = "logit"),
                data = as.data.frame(data$wide)
              ),
              "2" = glm(
                formula = create_formula("conventional (updated)", "Y3", spec, time),
                family = binomial(link = "logit"),
                data = as.data.frame(data$wide)
              ),
              "3" = glm(
                formula = create_formula("conventional (updated)", "Y3", spec, time),
                family = binomial(link = "logit"),
                data = as.data.frame(data$wide)
              )
            )
          },
          "gformula (pooled)" = function(data, time, spec) {
            list(
              X = fit_X_models(
                models = list(
                  "L1" = list(
                    formula = create_formula("gformula (pooled)", "L1", spec, time, long = TRUE),
                    family = "normal"
                  ),
                  "L2" = list(
                    formula = create_formula("gformula (pooled)", "L2", spec, time, long = TRUE),
                    family = "normal"
                  ),
                  "L5" = list(
                    formula = create_formula("gformula (pooled)", "L5", spec, time, long = TRUE),
                    family = "normal"
                  ), 
                  "A" = list(
                    formula = create_formula("gformula (pooled)", "A", spec, time, long = TRUE),
                    link = "logit",
                    family = "binomial"
                  )
                ),
                data = as.data.frame(data$long)
              ),
              Y = fit_Y_model(
                model = list(
                  formula = create_formula("gformula (pooled)", "Y", spec, time, long = TRUE),
                  link = "logit",
                  family = "binomial"
                ),
                data = as.data.frame(data$long)
              ),
              D = if (spec == "outcome-specific risk") {
                fit_D_model(
                  model = list(
                    formula = create_formula("gformula (pooled)", "D", spec, time, long = TRUE),
                    link = "logit",
                    family = "binomial"
                  ),
                  data = as.data.frame(data$long)
                )
              } else {
                NULL
              }
            )
          },
          "gformula (unpooled)" = function(data, time, spec) {
            list(
              X = fit_X_models_ice(
                models = list(
                  "L1" = list(
                    formula = create_formula("gformula (unpooled)", "L1", spec, time, long = TRUE),
                    family = "normal"
                  ),
                  "L2" = list(
                    formula = create_formula("gformula (unpooled)", "L2", spec, time, long = TRUE),
                    family = "normal"
                  ),
                  "L5" = list(
                    formula = create_formula("gformula (unpooled)", "L5", spec, time, long = TRUE),
                    family = "normal"
                  ), 
                  "A" = list(
                    formula = create_formula("gformula (unpooled)", "A", spec, time, long = TRUE),
                    link = "logit",
                    family = "binomial"
                  )
                ),
                data = as.data.frame(data$long),
                k = k
              ),
              Y = fit_Y_model_ice(
                model = list(
                  formula = create_formula("gformula (unpooled)", "Y", spec, time, long = TRUE),
                  link = "logit",
                  family = "binomial"
                ),
                data = as.data.frame(data$long),
                k = k
              ),
              D = if (spec == "outcome-specific risk") {
                fit_D_model_ice(
                  model = list(
                    formula = create_formula("gformula (unpooled)", "D", spec, time, long = TRUE),
                    link = "logit",
                    family = "binomial"
                  ),
                  data = as.data.frame(data$long),
                  k = k
                )
              } else {
                NULL
              }
            )
          })
      }),
      p_function = map(estimator, function(x) {
        switch(
          x,
          "conventional" = function(fit, data, time) {
            predict(fit, as.data.frame(data$wide), type = "response")
          },
          "conventional (updated)" = function(fit, data, time) {
            predict(fit, as.data.frame(data$wide), type = "response")
          },
          "gformula (pooled)" = function(fit, data, time) {
            p <- gformula_mc(
              Y.fit = fit[["Y"]],
              X.fit = fit[["X"]],
              D.fit = fit[["D"]],
              data = subset(as.data.frame(data$long), !is.na(L1)),
              id = "id",
              time = "time",
              base.covs = c("L3", "L4"),
              hist.vars = c(
                paste0("lag1_", paste0("L", c(1, 2, 5))), 
                paste0("lag2_", paste0("L", c(1, 2, 5))),
                paste0("lag3_", paste0("L", c(1, 2, 5))),
                paste0("lag", 1:3, "_A")
              ),
              hist.fun = "lag",
              treatment = treatment,
              intervention = shift,
              mc.sims = mc.sims,
              mc.start = time,
              mc.stop = k - 1,
              merge = FALSE,
              last_only = TRUE
            )
            
            p <- left_join(as.data.frame(data$wide), p, by = "id")
            p[["pred"]]
          },
          "gformula (unpooled)" = function(fit, data, time) {
            p <- gformula_mc(
              Y.fit = fit[["Y"]],
              X.fit = fit[["X"]],
              D.fit = fit[["D"]],
              data = subset(as.data.frame(data$long), !is.na(L1)),
              id = "id",
              time = "time",
              ice = TRUE,
              base.covs = c("L3", "L4"),
              hist.vars = c(
                paste0("lag1_", paste0("L", c(1, 2, 5))), 
                paste0("lag2_", paste0("L", c(1, 2, 5))),
                paste0("lag3_", paste0("L", c(1, 2, 5))),
                paste0("lag", 1:3, "_A")
              ),
              hist.fun = "lag",
              treatment = treatment,
              intervention = shift,
              mc.sims = mc.sims,
              mc.start = time,
              mc.stop = k - 1,
              merge = FALSE,
              last_only = TRUE
            )
            
            p <- left_join(as.data.frame(data$wide), p, by = "id")
            p[["pred"]]
          }
        )
      })
    )
  
  
  # run simulations ---------------------------------------------------------
  
  tic = Sys.time()
  with_progress({
    p <- progressor(steps = nrow(sims))
    
    # predict on training and test data 
    sims <-
      sims %>%
      mutate(
        p_test = future_pmap(
          list(model, estimator, specification, time, p_function, train, test),
          get_predictions,
          .options = furrr_options(
            seed = TRUE
          ),
          p = p
        ),
        p_truth = map(
          truth,
          function(truth) tibble(p = truth[, calc_position(paste0("Y_cum",k-1), truth)])
        )
      )
  })
  Sys.time() - tic
  # plan(sequential)
  
  
  # calculate stats ---------------------------------------------------------
  
  sims <-
    sims %>%
    mutate(    
      bias = map2(
        p_test,
        p_truth,
        function(test, truth) mean(abs(test$p - truth$p), na.rm = TRUE)
      ),
      rmspe = map2(
        p_test,
        p_truth,
        function(test, truth) sqrt(mean((test$p - truth$p)^2, na.rm = TRUE))
      ),
      # test_stats = map(
      #   p_test,
      #   function(x) {
      #     # add small tolerance in case of perfect predictions
      #     p <- if_else(x$p == 0, x$p + 0.0001, if_else(x$p == 1, x$p - 0.0001, x$p))
      #     val.prob(p, x$y, pl = FALSE)
      #   }
      # )
    )
  
  return(sims)
}



# helper functions --------------------------------------------------------


create_formula <- function(estimator, variable, specification, time, long = FALSE) {
  
  if (long) {
    if (specification == "incorrect covariate") {
      f <- switch(
        variable,
        "L1" = L1 ~ 
          I(exp(lag1_L1/2)) + I(exp(lag2_L1/2)) + I(exp(lag3_L1/2)) + 
          I(lag1_L2/(1 + exp(lag1_L1)) + 10) + I(lag2_L2/(1 + exp(lag2_L1)) + 10) + I(lag3_L2/(1 + exp(lag3_L1)) + 10) + 
          lag1_A + lag2_A + lag3_A,
        "L2" = L2 ~ 
          I(exp(L1/2)) + I(exp(lag1_L1/2)) + I(exp(lag2_L1/2)) + I(exp(lag3_L1/2)) + 
          I(lag1_L2/(1 + exp(lag1_L1)) + 10) + I(lag2_L2/(1 + exp(lag2_L1)) + 10) + I(lag3_L2/(1 + exp(lag3_L1)) + 10) + 
          lag1_A + lag2_A + lag3_A,
        "L5" = L5 ~ I(((lag1_L1 * lag1_L5)/25 + 0.6)^3), 
        "A" = A ~ 
          I(exp(L1/2)) + I(exp(lag1_L1/2)) + I(exp(lag2_L1/2)) + I(exp(lag3_L1/2)) + 
          I(L2/(1 + exp(L1)) + 10) + I(lag1_L2/(1 + exp(lag1_L1)) + 10) + I(lag2_L2/(1 + exp(lag2_L1)) + 10) + I(lag3_L2/(1 + exp(lag3_L1)) + 10) + 
          L3 + 
          lag1_A + lag2_A + lag3_A,
        "Y" =  Y ~ 
          L1 + lag1_L1 + lag2_L1 + lag3_L1 + 
          L2 + lag1_L2 + lag2_L2 + lag3_L2 + 
          L3 + 
          L4 + 
          L5 + lag1_L5 + lag2_L5 + lag3_L5 + 
          A + lag1_A + lag2_A + lag3_A,
        "D" = D ~ 
          L1 + lag1_L1 + lag2_L1 + lag3_L1 + 
          L2 + lag1_L2 + lag2_L2 + lag3_L2 + 
          L3 + 
          L4 + 
          L5 + lag1_L5 + lag2_L5 + lag3_L5 + 
          A + lag1_A + lag2_A + lag3_A,
      )
    } else if (specification == "incorrect outcome") {
      f <- switch(
        variable,
        "L1" = L1 ~ 
          lag1_L1 + lag2_L1 + lag3_L1 + 
          lag1_L2 + lag2_L2 + lag3_L2 + 
          lag1_A + lag2_A + lag3_A,
        "L2" = L2 ~ 
          L1 + lag1_L1 + lag2_L1 + lag3_L1 + 
          lag1_L2 + lag2_L2 + lag3_L2 + 
          lag1_A + lag2_A + lag3_A,
        "L5" = L5 ~ lag1_L5,
        "A" = A ~ 
          L1 + lag1_L1 + lag2_L1 + lag3_L1 + 
          L2 + lag1_L2 + lag2_L2 + lag3_L2 + 
          L3 + 
          lag1_A + lag2_A + lag3_A,
        "Y" =  Y ~ 
          I(exp(L1/2)) + I(exp(lag1_L1/2)) + I(exp(lag2_L1/2)) + I(exp(lag3_L1/2)) + 
          I(L2/(1 + exp(L1)) + 10) + I(lag1_L2/(1 + exp(lag1_L1)) + 10) + I(lag2_L2/(1 + exp(lag2_L1)) + 10) + I(lag3_L2/(1 + exp(lag3_L1)) + 10) + 
          L3 + 
          L4 + 
          I(((L1 * L5)/25 + 0.6)^3) + I(((lag1_L1 * lag1_L5)/25 + 0.6)^3) + I(((lag2_L1 * lag2_L5)/25 + 0.6)^3) + I(((lag3_L1 * lag3_L5)/25 + 0.6)^3) + 
          A + lag1_A + lag2_A + lag3_A,
        "D" = D ~ 
          I(exp(L1/2)) + I(exp(lag1_L1/2)) + I(exp(lag2_L1/2)) + I(exp(lag3_L1/2)) + 
          I(L2/(1 + exp(L1)) + 10) + I(lag1_L2/(1 + exp(lag1_L1)) + 10) + I(lag2_L2/(1 + exp(lag2_L1)) + 10) + I(lag3_L2/(1 + exp(lag3_L1)) + 10) + 
          L3 + 
          L4 + 
          I(((L1 * L5)/25 + 0.6)^3) + I(((lag1_L1 * lag1_L5)/25 + 0.6)^3) + I(((lag2_L1 * lag2_L5)/25 + 0.6)^3) + I(((lag3_L1 * lag3_L5)/25 + 0.6)^3) + 
          A + lag1_A + lag2_A + lag3_A,
      )
    } else {
      f <- switch(
        variable,
        "L1" = L1 ~ 
          lag1_L1 + lag2_L1 + lag3_L1 + 
          lag1_L2 + lag2_L2 + lag3_L2 + 
          lag1_A + lag2_A + lag3_A,
        "L2" = L2 ~ 
          L1 + lag1_L1 + lag2_L1 + lag3_L1 + 
          lag1_L2 + lag2_L2 + lag3_L2 + 
          lag1_A + lag2_A + lag3_A,
        "L5" = L5 ~ lag1_L5, 
        "A" = A ~ 
          L1 + lag1_L1 + lag2_L1 + lag3_L1 + 
          L2 + lag1_L2 + lag2_L2 + lag3_L2 + 
          L3 + 
          lag1_A + lag2_A + lag3_A,
        "Y" =  Y ~ 
          L1 + lag1_L1 + lag2_L1 + lag3_L1 + 
          L2 + lag1_L2 + lag2_L2 + lag3_L2 + 
          L3 + 
          L4 + 
          L5 + lag1_L5 + lag2_L5 + lag3_L5 + 
          A + lag1_A + lag2_A + lag3_A,
        "D" = D ~ 
          L1 + lag1_L1 + lag2_L1 + lag3_L1 + 
          L2 + lag1_L2 + lag2_L2 + lag3_L2 + 
          L3 + 
          L4 + 
          L5 + lag1_L5 + lag2_L5 + lag3_L5 + 
          A + lag1_A + lag2_A + lag3_A,
      )
    }
    
  } else {
    if (estimator == "conventional") {
      if (specification == "incorrect outcome") {
        f <- switch(
          as.character(time),
          "0" = Y3 ~ A0 + I(exp(L1_0/2)) + I(L2_0/(1 + L1_0) + 10) + L3_0 + L4_0 + I(((L1_0 * L5_0)/25 + 0.6)^3),
          "1" = Y3 ~ A1 + I(exp(L1_1/2)) + I(L2_1/(1 + L1_1) + 10) + L3_0 + L4_0 + I(((L1_1 * L5_1)/25 + 0.6)^3),
          "2" = Y3 ~ A2 + I(exp(L1_2/2)) + I(L2_2/(1 + L1_2) + 10) + L3_0 + L4_0 + I(((L1_2 * L5_2)/25 + 0.6)^3),
          "3" = Y3 ~ A3 + I(exp(L1_3/2)) + I(L2_3/(1 + L1_3) + 10) + L3_0 + L4_0 + I(((L1_3 * L5_3)/25 + 0.6)^3)
        )
      } else {
        f <- switch(
          as.character(time),
          "0" = Y3 ~ A0 + L1_0 + L2_0 + L3_0 + L4_0 + L5_0,
          "1" = Y3 ~ A1 + L1_1 + L2_1 + L3_0 + L4_0 + L5_1,
          "2" = Y3 ~ A2 + L1_2 + L2_2 + L3_0 + L4_0 + L5_2,
          "3" = Y3 ~ A3 + L1_3 + L2_3 + L3_0 + L4_0 + L5_3
        )
      }
    } else {
      if (specification == "incorrect outcome") {
        f <- switch(
          as.character(time),
          "0" = Y3 ~ A0 + I(exp(L1_0/2)) + I(L2_0/(1 + L1_0) + 10) + L3_0 + L4_0 + I(((L1_0 * L5_0)/25 + 0.6)^3),
          "1" =  Y3 ~ 
            A0 + I(exp(L1_0/2)) + I(L2_0/(1 + L1_0) + 10) + L3_0 + L4_0 + I(((L1_0 * L5_0)/25 + 0.6)^3) + 
            A1 + I(exp(L1_1/2)) + I(L2_1/(1 + L1_1) + 10) + I(((L1_1 * L5_1)/25 + 0.6)^3),
          "2" = Y3 ~ 
            A0 + I(exp(L1_0/2)) + I(L2_0/(1 + L1_0) + 10) + L3_0 + L4_0 + I(((L1_0 * L5_0)/25 + 0.6)^3) + 
            A1 + I(exp(L1_1/2)) + I(L2_1/(1 + L1_1) + 10) + I(((L1_1 * L5_1)/25 + 0.6)^3) +
            A2 + I(exp(L1_2/2)) + I(L2_2/(1 + L1_2) + 10) + I(((L1_2 * L5_2)/25 + 0.6)^3),
          "3" = Y3 ~ 
            A0 + I(exp(L1_0/2)) + I(L2_0/(1 + L1_0) + 10) + L3_0 + L4_0 + I(((L1_0 * L5_0)/25 + 0.6)^3) + 
            A1 + I(exp(L1_1/2)) + I(L2_1/(1 + L1_1) + 10) + I(((L1_1 * L5_1)/25 + 0.6)^3) +
            A2 + I(exp(L1_2/2)) + I(L2_2/(1 + L1_2) + 10) + I(((L1_2 * L5_2)/25 + 0.6)^3) + 
            A3 + I(exp(L1_3/2)) + I(L2_3/(1 + L1_3) + 10) + I(((L1_3 * L5_3)/25 + 0.6)^3)
        )
      } else {
        f <- switch(
          as.character(time),
          "0" = Y3 ~ A0 + L1_0 + L2_0 + L3_0 + L4_0 + L5_0,
          "1" =  Y3 ~ 
            A0 + L1_0 + L2_0 + L3_0 + L4_0 + L5_0 + 
            A1 + L1_1 + L2_1 + L5_1,
          "2" = Y3 ~ 
            A0 + L1_0 + L2_0 + L3_0 + L4_0 + L5_0 + 
            A1 + L1_1 + L2_1 + L5_1 + 
            A2 + L1_2 + L2_2 + L5_2,
          "3" = Y3 ~ 
            A0 + L1_0 + L2_0 + L3_0 + L4_0 + L5_0 + 
            A1 + L1_1 + L2_1 + L5_1 + 
            A2 + L1_2 + L2_2 + L5_2 + 
            A3 + L1_3 + L2_3 + L5_3
        )
      }
    }
  }
  return(f)
  
}


fit_models <- function(model, estimator, specification, time, train) {
  
  fit <- model(train, time, specification)
  
  return(fit)
}

get_predictions <- function(model, estimator, specification, time, p_function, train, test, p) {
  
  fit <- model(train, time, specification)
  
  pr <- p_function(fit, test, time)
  y <- test$wide[, calc_position("Y3", test$wide)]
  p()
  tibble(p = pr, y)
}
