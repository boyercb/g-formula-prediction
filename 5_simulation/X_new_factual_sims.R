
# set up fixed model parameter values -------------------------------------

# scenario 1:
# L1 and L2 are time-varying risk factors that determine whether treatment A is
# administered and there is treatment-confounder feedback. L5 is a time-varying
# risk factor that is unaffected by treatment. L3 and L4 are fixed baseline risk
# factors with L3 affecting whether future treatment is administered. Lagged
# effects, where present, are exponentially decreasing.
alpha <- matrix(
  # Int,               L1,               L2,            L3,            L4,              L5,                 A
  c(  0, c(0, 2^(-(1:3))), c(0, 2^(-(1:3))),     rep(0, 4),     rep(0, 4),       rep(0, 4), c(0, -1.5^(-(1:3))),  # L1
      0,       2^(-(1:4)), c(0, 2^(-(1:3))),     rep(0, 4),     rep(0, 4),       rep(0, 4), c(0, -1.5^(-(1:3))),  # L2
      0,        rep(0, 4),        rep(0, 4), c(0, 1, 0, 0),     rep(0, 4),       rep(0, 4),         rep(0, 4),  # L3
      0,        rep(0, 4),        rep(0, 4),     rep(0, 4), c(0, 1, 0, 0),       rep(0, 4),         rep(0, 4),  # L4
      0,        rep(0, 4),        rep(0, 4),     rep(0, 4),     rep(0, 4), c(0, 1.1, 0, 0),         rep(0, 4)), # L5
  nrow = 5,
  ncol = 25,
  byrow = TRUE
)

#                 Int,         L1,         L2,               L3,              L4,         L5,                A
beta <-    c(log(0.1), 2^(-(1:4)), 3^(-(1:4)), -c(0, 0.5, 0, 0),       rep(0, 4),  rep(0, 4),  c(0, 2^(-(1:3))))
epsilon <- c(log(0.05), 2^(-(1:4)), 3^(-(1:4)),  c(0.5, 0, 0, 0), c(0.5, 0, 0, 0), 1^(-(1:4)),       -3^(-(1:4)))

#            L1   L2   L3   L4   L5
error <-  c(0.2, 0.2,   0,   0, 0.2)

rho <- 0.1
sigma <- diag(rep(1, 5))
sigma[lower.tri(sigma)] <- rho
sigma[upper.tri(sigma)] <- rho

SIMS <- 500
GSIMS <- 100
SIZES <- c(500, 1500, 3000)
p <- 5
k <- 4


# initialize simulation grid ----------------------------------------------

# define simulation parameters
factual_sims_df <- expand_grid(
  sim = 1:SIMS,
  n = SIZES
)

# create training and validation datasets
factual_sims_df <- 
  factual_sims_df %>%
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
          n = max(SIZES),
          k = k,
          p = p,
          a = alpha,
          b = beta,
          e = epsilon,
          error = error,
          sigma = sigma,
          long = TRUE
        )
      })
  ) 

# expand grid to include conditioning points
factual_sims_df <-
  factual_sims_df %>%
  expand_grid(time = 0:(k - 1)) 

# add true value of estimand
factual_sims_df <-
  factual_sims_df %>%
  mutate(
    truth = map2(
      test,
      time,
      function(test, time) {
        calculate_sem_truth(
          n = max(SIZES),
          k = k,
          p = p,
          a = alpha,
          b = beta,
          e = epsilon,
          L0 = test$wide[, calc_position(
            sapply(paste0("_", 0:time), function(x) paste0("L", 1:5, x)) %>% 
              as.character(),
            test$wide)
            ],
          A0 = test$wide[, calc_position(paste0("A", 0:time), test$wide)],
          time = time
        )
      })
  )

# expand grid to include estimators 
factual_sims_df <-
  factual_sims_df %>%
  expand_grid(
    estimator = c(
      "conventional",
      "conventional (updated)",
      "gformula (pooled)",
      "gformula (unpooled)"
    ),
    specification = c(
      "correct",
      "incorrect covariate",
      "incorrect outcome"
    )
  ) 

# defined model and prediction functions for each estimator
factual_sims_df <-
  factual_sims_df %>%
  mutate(
    model = map(estimator, function(x) {
      switch(
        x,
        "conventional" = function(data, time, spec) {
          if (spec == "incorrect outcome") {
            data$wide <- misspecify(as.data.frame(data$wide))
          } 
          switch(
            as.character(time),
            "0" = glm(
              formula = Y3 ~ A0 + L1_0 + L2_0 + L3_0 + L4_0 + L5_0,
              family = binomial(link = "logit"),
              data = as.data.frame(data$wide)
            ),
            "1" = glm(
              formula = Y3 ~ 
                A1 + L1_1 + L2_1 + L3_0 + L4_0 + L5_1,
              family = binomial(link = "logit"),
              data = as.data.frame(data$wide)
            ),
            "2" = glm(
              formula = Y3 ~ 
                A2 + L1_2 + L2_2 + L3_0 + L4_0 + L5_2,
              family = binomial(link = "logit"),
              data = as.data.frame(data$wide)
            ),
            "3" = glm(
              formula = Y3 ~ 
                A3 + L1_3 + L2_3 + L3_0 + L4_0 + L5_3,
              family = binomial(link = "logit"),
              data = as.data.frame(data$wide)
            )
          ) 
        },
        "conventional (updated)" = function(data, time, spec) {
          if (spec == "incorrect outcome") {
            data$wide <- misspecify(as.data.frame(data$wide))
          } 
          switch(
            as.character(time),
            "0" = glm(
              formula = Y3 ~ A0 + L1_0 + L2_0 + L3_0 + L4_0 + L5_0,
              family = binomial(link = "logit"),
              data = as.data.frame(data$wide)
            ),
            "1" = glm(
              formula = Y3 ~ 
                A0 + L1_0 + L2_0 + L3_0 + L4_0 + L5_0 + 
                A1 + L1_1 + L2_1 + L5_1,
              family = binomial(link = "logit"),
              data = as.data.frame(data$wide)
            ),
            "2" = glm(
              formula = Y3 ~ 
                A0 + L1_0 + L2_0 + L3_0 + L4_0 + L5_0 + 
                A1 + L1_1 + L2_1 + L5_1 + 
                A2 + L1_2 + L2_2 + L5_2,
              family = binomial(link = "logit"),
              data = as.data.frame(data$wide)
            ),
            "3" = glm(
              formula = Y3 ~ 
                A0 + L1_0 + L2_0 + L3_0 + L4_0 + L5_0 + 
                A1 + L1_1 + L2_1 + L5_1 + 
                A2 + L1_2 + L2_2 + L5_2 + 
                A3 + L1_3 + L2_3 + L5_3,
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
                  formula = L1 ~ 
                    lag1_L1 + lag2_L1 + lag3_L1 + 
                    lag1_L2 + lag2_L2 + lag3_L2 + 
                    lag1_A + lag2_A + lag3_A,
                  family = "normal"
                ),
                "L2" = list(
                  formula = L2 ~ 
                    L1 + lag1_L1 + lag2_L1 + lag3_L1 + 
                    lag1_L2 + lag2_L2 + lag3_L2 + 
                    lag1_A + lag2_A + lag3_A,
                  family = "normal"
                ),
                "L5" = list(
                  formula = L5 ~ lag1_L5,
                  family = "normal"
                ), 
                "A" = list(
                  formula = A ~ 
                    L1 + lag1_L1 + lag2_L1 + lag3_L1 + 
                    L2 + lag1_L2 + lag2_L2 + lag3_L2 + 
                    lag1_A + lag2_A + lag3_A,
                  link = "logit",
                  family = "binomial"
                )
              ),
              data = if (spec == "incorrect covariate") {
                misspecify(as.data.frame(data$long), long = TRUE)
              } else {
                as.data.frame(data$long)
              }
            ),
            Y = fit_Y_model(
              model = list(
                formula = Y ~ 
                  L1 + lag1_L1 + lag2_L1 + lag3_L1 + 
                  L2 + lag1_L2 + lag2_L2 + lag3_L2 + 
                  L3 + 
                  L4 + 
                  L5 + lag1_L5 + lag2_L5 + lag3_L5 + 
                  A + lag1_A + lag2_A + lag3_A,
                link = "logit",
                family = "binomial"
              ),
              data = if (spec == "incorrect outcome") {
                misspecify(as.data.frame(data$long), long = TRUE)
              } else {
                as.data.frame(data$long)
              }
            )
          )
        },
        "gformula (unpooled)" = function(data, time, spec) {
          list(
            X = fit_X_models_ice(
              models = list(
                "L1" = list(
                  formula = L1 ~ 
                    lag1_L1 + lag2_L1 + lag3_L1 + 
                    lag1_L2 + lag2_L2 + lag3_L2 + 
                    lag1_A + lag2_A + lag3_A,
                  family = "normal"
                ),
                "L2" = list(
                  formula = L2 ~ 
                    L1 + lag1_L1 + lag2_L1 + lag3_L1 + 
                    lag1_L2 + lag2_L2 + lag3_L2 + 
                    lag1_A + lag2_A + lag3_A,
                  family = "normal"
                ),
                "L5" = list(
                  formula = L5 ~ lag1_L5,
                  family = "normal"
                ), 
                "A" = list(
                  formula = A ~ 
                    L1 + lag1_L1 + lag2_L1 + lag3_L1 + 
                    L2 + lag1_L2 + lag2_L2 + lag3_L2 + 
                    lag1_A + lag2_A + lag3_A,
                  link = "logit",
                  family = "binomial"
                )
              ),
              data = if (spec == "incorrect covariate") {
                misspecify(as.data.frame(data$long), long = TRUE)
              } else {
                as.data.frame(data$long)
              },
              k = k
            ),
            Y = fit_Y_model_ice(
              model = list(
                formula = Y ~ 
                  L1 + lag1_L1 + lag2_L1 + lag3_L1 + 
                  L2 + lag1_L2 + lag2_L2 + lag3_L2 + 
                  L3 + 
                  L4 + 
                  L5 + lag1_L5 + lag2_L5 + lag3_L5 + 
                  A + lag1_A + lag2_A + lag3_A,
                link = "logit",
                family = "binomial"
              ),
              data = if (spec == "incorrect outcome") {
                misspecify(as.data.frame(data$long), long = TRUE)
              } else {
                as.data.frame(data$long)
              },
              k = k
            )
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
            mc.sims = GSIMS,
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
            mc.sims = GSIMS,
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

# plan(multicore, workers = 12)
tic = Sys.time()
with_progress({
  p <- progressor(steps = nrow(factual_sims_df))

  # predict on training and test data 
  factual_sims_df <-
    factual_sims_df %>%
    mutate(
      # fit = future_pmap(
      #   list(model, estimator, specification, time, train),
      #   fit_models,
      #   .options = furrr_options(
      #     seed = TRUE
      #   )
      # ),
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
        function(truth) tibble(p = truth[, k])
      )
    )
})
Sys.time() - tic
# plan(sequential)


# calculate stats ---------------------------------------------------------

factual_sims_df <-
  factual_sims_df %>%
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
    # train_stats = map(
    #   p_train,
    #   ~val.prob(.x$p, .x$y, pl = FALSE)
    # ),
    # test_stats = map(
    #   p_test,
    #   ~val.prob(.x$p + 0.01 * I(.x$p == 0) - 0.01 * I(.x$p == 1), .x$y, pl = FALSE)
    # )
  )


# TODO: 
# 1. cut back the generate truth function to just output p
# 2. make sure there's no memory leak
# 3. store data as matrices instead of tibbles for speed/memory?
# future isn't working
misspecify <- function(data, long = FALSE) {
  
  # wide_df <- as.data.frame(data$wide)
  # long_df <- as.data.frame(data$long)
  
  for (i in 1:(k - 1)) {
    if (long) {
      data[data[["time"]] == i, "L1"] <- exp(data[data[["time"]] == i, "L1"] / 2)
      data[data[["time"]] == i, "L2"] <- data[data[["time"]] == i, "L2"] / (1 + exp(data[data[["time"]] == i, "L1"])) + 10
      data[data[["time"]] == i, "L5"] <- (data[data[["time"]] == i, "L1"] * data[data[["time"]] == i, "L5"] / 25 + 0.6) ^ 3
      
      data[data[["time"]] == i, "lag1_L1"] <- exp(data[data[["time"]] == i, "lag1_L1"] / 2)
      data[data[["time"]] == i, "lag1_L2"] <- data[data[["time"]] == i, "lag1_L2"] / (1 + exp(data[data[["time"]] == i, "lag1_L1"])) + 10
      data[data[["time"]] == i, "lag1_L5"] <- (data[data[["time"]] == i, "lag1_L1"] * data[data[["time"]] == i, "lag1_L5"] / 25 + 0.6) ^ 3
      
      data[data[["time"]] == i, "lag2_L1"] <- exp(data[data[["time"]] == i, "lag2_L1"] / 2)
      data[data[["time"]] == i, "lag2_L2"] <- data[data[["time"]] == i, "lag2_L2"] / (1 + exp(data[data[["time"]] == i, "lag2_L1"])) + 10
      data[data[["time"]] == i, "lag2_L5"] <- (data[data[["time"]] == i, "lag2_L1"] * data[data[["time"]] == i, "lag2_L5"] / 25 + 0.6) ^ 3
      
      data[data[["time"]] == i, "lag3_L1"] <- exp(data[data[["time"]] == i, "lag3_L1"] / 2)
      data[data[["time"]] == i, "lag3_L2"] <- data[data[["time"]] == i, "lag3_L2"] / (1 + exp(data[data[["time"]] == i, "lag3_L1"])) + 10
      data[data[["time"]] == i, "lag3_L5"] <- (data[data[["time"]] == i, "lag3_L1"] * data[data[["time"]] == i, "lag3_L5"] / 25 + 0.6) ^ 3
    } else {
      data[[paste0("L1_", i)]] <- exp(data[[paste0("L1_", i)]] / 2)
      data[[paste0("L2_", i)]] <- data[[paste0("L2_", i)]] / (1 + exp(data[[paste0("L1_", i)]])) + 10
      data[[paste0("L5_", i)]] <- (data[[paste0("L1_", i)]] * data[[paste0("L5_", i)]] / 25 + 0.6) ^ 3
    }
    
  }
  
  return(data)
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


