# Fit parametric models ---------------------------------------------------
censored_ids <- observed_events_wide$pid[is.na(observed_events_wide$event_chd_zero_5)]
Y.fit_lag0 <- 
  fit_Y_model(
    Y.model_lag0,
    data = filter(analytic_long, time < 3 & !pid %in% censored_ids)
  )

D.fit_lag0 <- 
  fit_D_model(
    D.model_lag0,
    data = filter(analytic_long, time < 3 & !pid %in% censored_ids)
  )

Y.fit_lag1 <- 
  fit_Y_model(
    Y.model_lag1,
    data = filter(analytic_long, time < 3 & !pid %in% censored_ids)
  )

X.fit_lag1 <- 
  fit_X_models(
    X.models_lag1,
    data = filter(analytic_long, time < 3 & !pid %in% censored_ids),
    time = "time"
  )

D.fit_lag1 <- 
  fit_D_model(
    D.model_lag1,
    data = filter(analytic_long, time < 3 & !pid %in% censored_ids)
  )

Y.fit_lag2 <- 
  fit_Y_model(
    Y.model_lag2,
    data = filter(analytic_long, time < 3 & !pid %in% censored_ids)
  )

X.fit_lag2 <- 
  fit_X_models(
    X.models_lag2,
    data = filter(analytic_long, time < 3 & !pid %in% censored_ids),
    time = "time"
  )

D.fit <- 
  fit_D_model(
    D.model,
    data = filter(analytic_long, time < 3 & !pid %in% censored_ids)
  )


Y.fit <- 
  fit_Y_model(
    Y.model,
    data = filter(analytic_long, time < 3 & !pid %in% censored_ids)
  )

X.fit <- 
  fit_X_models(
    X.models,
    data = filter(analytic_long, time < 3 & !pid %in% censored_ids),
    time = "time"
  )

D.fit_lag2 <- 
  fit_D_model(
    D.model_lag2,
    data = filter(analytic_long, time < 3 & !pid %in% censored_ids)
  )

X.fit_lasso <- 
  fit_X_models(
    X.lasso_lag1,
    data = filter(analytic_long, time <=2 & !pid %in% censored_ids),
    time = "time"
  )

Y.fit_lasso <- 
  fit_Y_model(
    Y.lasso_lag1,
    data = filter(analytic_long, time <=2 & !pid %in% censored_ids)
  )

D.fit_lasso <- 
  fit_D_model(
    D.lasso_lag1,
    data = filter(analytic_long, time <=2 & !pid %in% censored_ids)
  )


X.fit_lasso <- 
  fit_X_models(
    X.lasso_lag2,
    data = analytic_long,
    time = "time"
  )

Y.fit_lasso <- 
  fit_Y_model(
    Y.lasso_lag2,
    data = analytic_long
  )

D.fit_lasso <- 
  fit_D_model(
    D.lasso_lag2,
    data = analytic_long
  )

X.fit_ints <- 
  fit_X_models(
    X.ints,
    data = filter(analytic_long, time < 3 & !pid %in% censored_ids),
    time = "time"
  )

Y.fit_ints <- 
  fit_Y_model(
    Y.ints,
    data = filter(analytic_long, time < 3 & !pid %in% censored_ids)
  )

D.fit_ints <- 
  fit_D_model(
    D.ints,
    data = filter(analytic_long, time < 3 & !pid %in% censored_ids)
  )


# Fit random forest models ------------------------------------------------

X.rf_lag1 <-
  map2(
    .x = list(
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "smk", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "smk", "bmi", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "smk", "bmi", "dm", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "smk", "bmi", "dm", "hrx", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "smk", "bmi", "dm", "hrx", "liprx", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "smk", "bmi", "dm", "hrx", "liprx", "tc", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "smk", "bmi", "dm", "hrx", "liprx", "tc", "hdl", "time")
      ),
    .y = covs_tv[covs_tv != "age"],
    .f = function(x, y) {
      if (y == "dm") {
        df <- subset(analytic_long, time > 0 & lag1_dm == 0)
      } else {
        df <- subset(analytic_long, time > 0)
      }
      
      rf <- 
        regression_forest(
          Y = df[[y]],
          X = df[, x],
          tune.parameters = "all",
        )
      
      rf$rmse <- add_grf_rmse(rf)
      rf$covs <- x
      
      if (y == "dm") {
        rf$restrict <- 
          list(
            subset = 'lag1_dm == 0',
            otherwise = 1
          )
      }
      
      if (y %in% c("smk", "dm", "hrx", "liprx")) {
        rf$type <- "binomial grf"
      } else {
        rf$type <- "normal grf"
      }
      return(rf)
  })

X.rf_lag1 <- c(
  X.fit_lag1[1],
  X.rf_lag1
)
names(X.rf_lag1) <- covs_tv

Y.rf_lag1 <-
  regression_forest(
    Y = analytic_long$event_chd[!is.na(analytic_long$event_chd)],
    X = analytic_long[!is.na(analytic_long$event_chd), c(
      "sex",
      "age",
      "lag1_smk",
      "lag1_bmi",
      "lag1_dm",
      "lag1_sbp",
      "lag1_hdl",
      "lag1_tc",
      "lag1_hrx",
      "lag1_liprx",
      "smk",
      "bmi",
      "dm",
      "hrx",
      "liprx",
      "tc",
      "hdl",
      "sbp",
      "time"
    )],
    tune.parameters = "all"
  )

Y.rf_lag1$outcome <- "event_chd"
Y.rf_lag1$covs <- c(
  "sex",
  "age",
  "lag1_smk",
  "lag1_bmi",
  "lag1_dm",
  "lag1_sbp",
  "lag1_hdl",
  "lag1_tc",
  "lag1_hrx",
  "lag1_liprx",
  "smk",
  "bmi",
  "dm",
  "hrx",
  "liprx",
  "tc",
  "hdl",
  "sbp",
  "time"
)
Y.rf_lag1$type <- "survival"
Y.rf_lag1$rmse <- add_grf_rmse(Y.rf_lag1)

D.rf_lag1 <-
  regression_forest(
    Y = analytic_long$event_dth,
    X = analytic_long[, c(
      "sex",
      "age",
      "lag1_smk",
      "lag1_bmi",
      "lag1_dm",
      "lag1_sbp",
      "lag1_hdl",
      "lag1_tc",
      "lag1_hrx",
      "lag1_liprx",
      "smk",
      "bmi",
      "dm",
      "hrx",
      "liprx",
      "tc",
      "hdl",
      "sbp",
      "time"
    )],
    tune.parameters = "all"
  )

D.rf_lag1$outcome <- "event_dth"
D.rf_lag1$type <- "survival"
D.rf_lag1$covs <- c(
  "sex",
  "age",
  "lag1_smk",
  "lag1_bmi",
  "lag1_dm",
  "lag1_sbp",
  "lag1_hdl",
  "lag1_tc",
  "lag1_hrx",
  "lag1_liprx",
  "smk",
  "bmi",
  "dm",
  "hrx",
  "liprx",
  "tc",
  "hdl",
  "sbp",
  "time"
)
D.rf_lag1$rmse <- add_grf_rmse(D.rf_lag1)


X.rf_lag2 <-
  map2(
    .x = list(
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "lag2_smk", "lag2_bmi", "lag2_dm", "lag2_sbp", "lag2_hdl", "lag2_tc", "lag2_hrx", "lag2_liprx", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "lag2_smk", "lag2_bmi", "lag2_dm", "lag2_sbp", "lag2_hdl", "lag2_tc", "lag2_hrx", "lag2_liprx", "smk", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "lag2_smk", "lag2_bmi", "lag2_sbp", "lag2_hdl", "lag2_tc", "lag2_hrx", "lag2_liprx", "smk", "bmi", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "lag2_smk", "lag2_bmi", "lag2_dm", "lag2_sbp", "lag2_hdl", "lag2_tc", "lag2_hrx", "lag2_liprx", "smk", "bmi", "dm", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "lag2_smk", "lag2_bmi", "lag2_dm", "lag2_sbp", "lag2_hdl", "lag2_tc", "lag2_hrx", "lag2_liprx", "smk", "bmi", "dm", "hrx", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "lag2_smk", "lag2_bmi", "lag2_dm", "lag2_sbp", "lag2_hdl", "lag2_tc", "lag2_hrx", "lag2_liprx", "smk", "bmi", "dm", "hrx", "liprx", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "lag2_smk", "lag2_bmi", "lag2_dm", "lag2_sbp", "lag2_hdl", "lag2_tc", "lag2_hrx", "lag2_liprx", "smk", "bmi", "dm", "hrx", "liprx", "tc", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "lag2_smk", "lag2_bmi", "lag2_dm", "lag2_sbp", "lag2_hdl", "lag2_tc", "lag2_hrx", "lag2_liprx", "smk", "bmi", "dm", "hrx", "liprx", "tc", "hdl", "time")
    ),
    .y = covs_tv[covs_tv != "age"],
    .f = function(x, y) {
      if (y == "dm") {
        df <- subset(analytic_long, time > 0 & lag1_dm == 0)
      } else {
        df <- subset(analytic_long, time > 0)
      }
      
      rf <- 
        regression_forest(
          Y = df[[y]],
          X = df[, x],
          tune.parameters = "all",
        )
      
      rf$rmse <- add_grf_rmse(rf)
      rf$covs <- x
      
      if (y == "dm") {
        rf$restrict <- 
          list(
            subset = 'lag1_dm == 0',
            otherwise = 1
          )
      }
      
      if (y %in% c("smk", "dm", "hrx", "liprx")) {
        rf$type <- "binomial grf"
      } else {
        rf$type <- "normal grf"
      }
      return(rf)
    })

X.rf_lag2 <- c(
  X.fit_lag2[1],
  X.rf_lag2
)
names(X.rf_lag2) <- covs_tv

Y.rf_lag2 <-
  regression_forest(
    Y = analytic_long$event_chd[!is.na(analytic_long$event_chd)],
    X = analytic_long[!is.na(analytic_long$event_chd), c(
      "sex",
      "age",
      "lag1_smk",
      "lag1_bmi",
      "lag1_dm",
      "lag1_sbp",
      "lag1_hdl",
      "lag1_tc",
      "lag1_hrx",
      "lag1_liprx",
      "lag2_smk",
      "lag2_bmi",
      "lag2_dm",
      "lag2_sbp",
      "lag2_hdl",
      "lag2_tc",
      "lag2_hrx",
      "lag2_liprx",
      "smk",
      "bmi",
      "dm",
      "hrx",
      "liprx",
      "tc",
      "hdl",
      "sbp",
      "time"
    )],
    tune.parameters = "all"
  )

Y.rf_lag2$outcome <- "event_chd"
Y.rf_lag2$covs <- c(
  "sex",
  "age",
  "lag1_smk",
  "lag1_bmi",
  "lag1_dm",
  "lag1_sbp",
  "lag1_hdl",
  "lag1_tc",
  "lag1_hrx",
  "lag1_liprx",
  "lag2_smk",
  "lag2_bmi",
  "lag2_dm",
  "lag2_sbp",
  "lag2_hdl",
  "lag2_tc",
  "lag2_hrx",
  "lag2_liprx",
  "smk",
  "bmi",
  "dm",
  "hrx",
  "liprx",
  "tc",
  "hdl",
  "sbp",
  "time"
)
Y.rf_lag2$type <- "survival"
Y.rf_lag2$rmse <- add_grf_rmse(Y.rf_lag2)

D.rf_lag2 <-
  regression_forest(
    Y = analytic_long$event_dth,
    X = analytic_long[, c(
      "sex",
      "age",
      "lag1_smk",
      "lag1_bmi",
      "lag1_dm",
      "lag1_sbp",
      "lag1_hdl",
      "lag1_tc",
      "lag1_hrx",
      "lag1_liprx",
      "lag2_smk",
      "lag2_bmi",
      "lag2_dm",
      "lag2_sbp",
      "lag2_hdl",
      "lag2_tc",
      "lag2_hrx",
      "lag2_liprx",
      "smk",
      "bmi",
      "dm",
      "hrx",
      "liprx",
      "tc",
      "hdl",
      "sbp",
      "time"
    )],
    tune.parameters = "all"
  )

D.rf_lag2$outcome <- "event_dth"
D.rf_lag2$type <- "survival"
D.rf_lag2$covs <- c(
  "sex",
  "age",
  "lag1_smk",
  "lag1_bmi",
  "lag1_dm",
  "lag1_sbp",
  "lag1_hdl",
  "lag1_tc",
  "lag1_hrx",
  "lag1_liprx",
  "lag2_smk",
  "lag2_bmi",
  "lag2_dm",
  "lag2_sbp",
  "lag2_hdl",
  "lag2_tc",
  "lag2_hrx",
  "lag2_liprx",
  "smk",
  "bmi",
  "dm",
  "hrx",
  "liprx",
  "tc",
  "hdl",
  "sbp",
  "time"
)
D.rf_lag2$rmse <- add_grf_rmse(D.rf_lag2)


#  probability forests ----------------------------------------------------

X.pf_lag1 <-
  map2(
    .x = list(
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "smk", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "smk", "bmi", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "smk", "bmi", "dm", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "smk", "bmi", "dm", "hrx", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "smk", "bmi", "dm", "hrx", "liprx", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "smk", "bmi", "dm", "hrx", "liprx", "tc", "time"),
      c("sex", "age", "lag1_smk", "lag1_bmi", "lag1_dm", "lag1_sbp", "lag1_hdl", "lag1_tc", "lag1_hrx", "lag1_liprx", "smk", "bmi", "dm", "hrx", "liprx", "tc", "hdl", "time")
    ),
    .y = covs_tv[covs_tv != "age"],
    .f = function(x, y) {
      if (y == "dm") {
        df <- subset(analytic_long, time > 0 & lag1_dm == 0)
      } else {
        df <- subset(analytic_long, time > 0)
      }
      
      if (y %in% c("smk", "dm", "hrx", "liprx")) {
        
        rf <- 
          probability_forest(
            Y = factor(df[[y]]),
            X = df[, x]
          )
        rf$type <- "binomial grf"
        
        rf$rmse <- add_grf_rmse(rf)
        rf$covs <- x
        
        if (y == "dm") {
          rf$restrict <- 
            list(
              subset = 'lag1_dm == 0',
              otherwise = 1
            )
        }
        
      } else {
        rf <- 
          regression_forest(
            Y = df[[y]],
            X = df[, x],
            tune.parameters = "all",
          )
        rf$type <- "normal grf"
        
        rf$rmse <- add_grf_rmse(rf)
        rf$covs <- x
        
      }
      return(rf)
    })

X.pf_lag1 <- c(
  X.fit_lag1[1],
  X.pf_lag1
)
names(X.pf_lag1) <- covs_tv

Y.pf_lag1 <-
  probability_forest(
    Y = factor(analytic_long$event_chd[!is.na(analytic_long$event_chd)]),
    X = analytic_long[!is.na(analytic_long$event_chd), c(
      "sex",
      "age",
      "lag1_smk",
      "lag1_bmi",
      "lag1_dm",
      "lag1_sbp",
      "lag1_hdl",
      "lag1_tc",
      "lag1_hrx",
      "lag1_liprx",
      "smk",
      "bmi",
      "dm",
      "hrx",
      "liprx",
      "tc",
      "hdl",
      "sbp",
      "time"
    )], 
    num.trees = 20000
  )

Y.pf_lag1$outcome <- "event_chd"
Y.pf_lag1$covs <- c(
  "sex",
  "age",
  "lag1_smk",
  "lag1_bmi",
  "lag1_dm",
  "lag1_sbp",
  "lag1_hdl",
  "lag1_tc",
  "lag1_hrx",
  "lag1_liprx",
  "smk",
  "bmi",
  "dm",
  "hrx",
  "liprx",
  "tc",
  "hdl",
  "sbp",
  "time"
)
Y.pf_lag1$type <- "survival"
Y.pf_lag1$rmse <- add_grf_rmse(Y.pf_lag1)

D.pf_lag1 <-
  probability_forest(
    Y = factor(analytic_long$event_dth),
    X = analytic_long[, c(
      "sex",
      "age",
      "lag1_smk",
      "lag1_bmi",
      "lag1_dm",
      "lag1_sbp",
      "lag1_hdl",
      "lag1_tc",
      "lag1_hrx",
      "lag1_liprx",
      "smk",
      "bmi",
      "dm",
      "hrx",
      "liprx",
      "tc",
      "hdl",
      "sbp",
      "time"
    )]
  )

D.pf_lag1$outcome <- "event_dth"
D.pf_lag1$type <- "survival"
D.pf_lag1$covs <- c(
  "sex",
  "age",
  "lag1_smk",
  "lag1_bmi",
  "lag1_dm",
  "lag1_sbp",
  "lag1_hdl",
  "lag1_tc",
  "lag1_hrx",
  "lag1_liprx",
  "smk",
  "bmi",
  "dm",
  "hrx",
  "liprx",
  "tc",
  "hdl",
  "sbp",
  "time"
)
D.pf_lag1$rmse <- add_grf_rmse(D.pf_lag1)


# conventional.fit_0_2 <- 
#   glm(
#     formula = event_chd_3 ~ sex + age0 + educ_1 + educ_2 + educ_3 +
#       marital_1 + marital_2 + eversmk + pre_dpd + pre_bmi + pre_dm + pre_sbp +
#       pre_cpd + pre_ldl + pre_hrx + pre_liprx + cpd_0 + dpd_0 + bmi_0 +
#       dm_0 + sbp_0 + ldl_0 + hrx_0 + liprx_0,
#     family = binomial(link = 'logit'),
#     data = analytic_wide
#   )
# 
# conventional.fit_0_6 <- 
#   glm(
#     formula = event_chd_6 ~ sex + age0 + educ_1 + educ_2 + educ_3 +
#       marital_1 + marital_2 + eversmk + pre_dpd + pre_bmi + pre_dm + pre_sbp +
#       pre_cpd + pre_ldl + pre_hrx + pre_liprx + cpd_0 + dpd_0 + bmi_0 +
#       dm_0 + sbp_0 + ldl_0 + hrx_0 + liprx_0,
#     family = binomial(link = 'logit'),
#     data = analytic_wide
#   )


# Create output tables ----------------------------------------------------

sink("9_results/tables/cov_models_lag1.tex")
texreg(
  l = X.fit_lag1,
  booktabs = TRUE,
  use.packages = FALSE,
  table = FALSE
) %>% print()
sink()

sink("9_results/tables/cov_models_lag2.tex")
texreg(
  l = X.fit_lag2,
  booktabs = TRUE,
  use.packages = FALSE,
  table = FALSE
) %>% print()
sink()

sink("9_results/tables/out_models.tex")
texreg(
  l = list(
    "\\shortstack{(1) \\\\ event\\_chd}" = Y.fit_lag0,
    "\\shortstack{(2) \\\\ event\\_chd}" = Y.fit_lag1,
    "\\shortstack{(3) \\\\ event\\_chd}" = Y.fit_lag2,
    "\\shortstack{(4) \\\\ event\\_dth}" = D.fit_lag0, 
    "\\shortstack{(5) \\\\ event\\_dth}" = D.fit_lag1,
    "\\shortstack{(6) \\\\ event\\_dth}" = D.fit_lag2),
  booktabs = TRUE,
  use.packages = FALSE,
  table = FALSE
) %>% print()
sink()
