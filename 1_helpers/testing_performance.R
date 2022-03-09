profvis::profvis(
  gformula_mc(
    Y.fit = fits[["Y"]],
    X.fit = fits[["X"]],
    data = subset(as.data.frame(factual_sims_df$test[[3]]$long), !is.na(L1)),
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
    mc.start = 0,
    mc.stop = 3,
    merge = FALSE,
    last_only = TRUE
  )
)

fits <- factual_sims_df$model[[3]](factual_sims_df$train[[3]], factual_sims_df$time[[3]])

p1 <- gformula_mc(
  Y.fit = fits[["Y"]],
  X.fit = fits[["X"]],
  data = subset(as.data.frame(factual_sims_df$test[[9]]$long), !is.na(L1)),
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
  mc.start = 0,
  mc.stop = 3,
  merge = FALSE,
  seed = 12345
)

p2 <- gformula_mc2(
  Y.fit = factual_sims_df$fit[[9]][["Y"]],
  X.fit = factual_sims_df$fit[[9]][["X"]],
  data = factual_sims_df$test[[9]]$long,
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
  mc.start = 0,
  mc.stop = 3,
  merge = FALSE,
  pool = TRUE,
  seed = 12345,
  bound.sims = FALSE
)

head(p1)
head(p2)


fit_X_models(
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
  data = as.data.frame(factual_sims_df$train[[1]]$long)
)


fit_X_models_ice(
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
  data = as.data.frame(factual_sims_df$train[[1]]$long),
  k = 4
)

fit_Y_model_ice(
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
  data = as.data.frame(factual_sims_df$train[[1]]$long),
  k = 4
)
