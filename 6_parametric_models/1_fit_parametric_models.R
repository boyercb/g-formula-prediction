
# Fit models for experiment 1 ---------------------------------------------

models_ex1 <- lapply(
  definitions_ex1, 
  function (def) {
    list(
      "X" = fit_X_models(
        models = def[["X"]],
        data = offspring_long,
        time = "time"
      ),
      "Y" = fit_Y_model(
        model = def[["Y"]],
        data = offspring_long
      ),
      "D" = fit_D_model(
        model = def[["D"]],
        data = offspring_long
      )
    )
  }) 
  
names(models_ex1) <- names(definitions_ex1)

# drop censored observations

models_ex1_uncens <- lapply(
  definitions_ex1, 
  function (def) {
    list(
      "X" = fit_X_models(
        models = def[["X"]],
        data = filter(offspring_long, !pid %in% censored_ids),
        time = "time"
      ),
      "Y" = fit_Y_model(
        model = def[["Y"]],
        data = filter(offspring_long, !pid %in% censored_ids)
      ),
      "D" = fit_D_model(
        model = def[["D"]],
        data = filter(offspring_long, !pid %in% censored_ids)
      )
    )
  }) 

names(models_ex1_uncens) <- names(definitions_ex1)

# drop D model

models_ex1_noD <- lapply(
  models_ex1,
  function(fit) {
    list(
      "X" = fit[["X"]],
      "Y" = fit[["Y"]],
      "D" = NULL
    )
  })

names(models_ex1_noD) <- names(definitions_ex1)


# Fit models for experiment 2 ---------------------------------------------

models_ex2 <- lapply(
  definitions_ex2, 
  function (def) {
    list(
      "X" = fit_X_models(
        models = def[["X"]],
        data = offspring_long,
        time = "time"
      ),
      "Y" = fit_Y_model(
        model = def[["Y"]],
        data = offspring_long
      ),
      "D" = fit_D_model(
        model = def[["D"]],
        data = offspring_long
      )
    )
  }) 

names(models_ex2) <- names(definitions_ex2)

# drop censored observations

models_ex2_uncens <- lapply(
  definitions_ex2, 
  function (def) {
    list(
      "X" = fit_X_models(
        models = def[["X"]],
        data = filter(offspring_long, !pid %in% censored_ids),
        time = "time"
      ),
      "Y" = fit_Y_model(
        model = def[["Y"]],
        data = filter(offspring_long, !pid %in% censored_ids)
      ),
      "D" = fit_D_model(
        model = def[["D"]],
        data = filter(offspring_long, !pid %in% censored_ids)
      )
    )
  }) 

names(models_ex2_uncens) <- names(definitions_ex2)

# drop D model

models_ex2_noD <- lapply(
  models_ex2,
  function(fit) {
    list(
      "X" = fit[["X"]],
      "Y" = fit[["Y"]],
      "D" = NULL
    )
  })

names(models_ex2_noD) <- names(definitions_ex2)


# Fit models for experiment 3 ---------------------------------------------

models_ex3 <- lapply(
  definitions_ex3, 
  function (def) {
    list(
      "X" = fit_X_models(
        models = def[["X"]],
        data = offspring_long,
        time = "time",
        gam = TRUE
      ),
      "Y" = fit_Y_model(
        model = def[["Y"]],
        data = offspring_long,
        gam = TRUE
      ),
      "D" = fit_D_model(
        model = def[["D"]],
        data = offspring_long,
        gam = TRUE
      )
    )
  }) 

names(models_ex3) <- names(definitions_ex3)
models_ex3$gam$X$smk$type <- "binomial"
models_ex3$gam$X$dm$type <- "binomial"
models_ex3$gam$X$hrx$type <- "binomial"
models_ex3$gam$X$liprx$type <- "binomial"

# drop censored observations

models_ex3_uncens <- lapply(
  definitions_ex3, 
  function (def) {
    list(
      "X" = fit_X_models(
        models = def[["X"]],
        data = filter(offspring_long, !pid %in% censored_ids),
        time = "time",
        gam = TRUE
      ),
      "Y" = fit_Y_model(
        model = def[["Y"]],
        data = filter(offspring_long, !pid %in% censored_ids),
        gam = TRUE
      ),
      "D" = fit_D_model(
        model = def[["D"]],
        data = filter(offspring_long, !pid %in% censored_ids),
        gam = TRUE
      )
    )
  }) 

names(models_ex3_uncens) <- names(definitions_ex3)
models_ex3_uncens$gam$X$smk$type <- "binomial"
models_ex3_uncens$gam$X$dm$type <- "binomial"
models_ex3_uncens$gam$X$hrx$type <- "binomial"
models_ex3_uncens$gam$X$liprx$type <- "binomial"

# drop D model

models_ex3_noD <- lapply(
  models_ex3,
  function(fit) {
    list(
      "X" = fit[["X"]],
      "Y" = fit[["Y"]],
      "D" = NULL
    )
  })

names(models_ex3_noD) <- names(definitions_ex3)


# Fit comparison models ---------------------------------------------------


conventional_models <- list(
  logistic =
    glm(
      formula = event_ascvd_2 ~ sex * (
        poly(log(age_0), 2) + smk_0 + dm_0 + hrx_0 +
          log(tc_0) + log(hdl_0) + log(sbp_0) + log(age_0):log(tc_0) +
          log(age_0):log(hdl_0) + log(age_0):smk_0 + log(age_0):log(sbp_0) + log(age_0):hrx_0 +
          log(age_0):hrx_0:log(sbp_0) + hrx_0:log(sbp_0)
      ), 
      family = binomial(link = 'logit'),
      data = offspring_wide
    ),

  coxph =
    coxph(
      Surv(enddate, event_ascvd) ~ strata(sex) + log(age_0) + I(log(age_0)^2) + smk_0 + dm_0 + 
        hrx_0 + log(tc_0) + log(hdl_0) + log(sbp_0) + log(age_0):log(tc_0) +
        log(age_0):log(hdl_0) + log(age_0):smk_0 + log(age_0):log(sbp_0) + log(age_0):hrx_0 + 
        log(age_0):hrx_0:log(sbp_0) + hrx_0:log(sbp_0),
      data = offspring_coxph
    )
)

conventional_models_uncens <- list(
  logistic =
    glm(
      formula = event_ascvd_1 ~ sex * (
        poly(log(age_0), 2) + smk_0 + dm_0 + hrx_0 +
          log(tc_0) + log(hdl_0) + log(sbp_0) + log(age_0):log(tc_0) +
          log(age_0):log(hdl_0) + log(age_0):smk_0 + log(age_0):log(sbp_0) + log(age_0):hrx_0 +
          log(age_0):hrx_0:log(sbp_0) + hrx_0:log(sbp_0)
      ), 
      family = binomial(link = 'logit'),
      data = filter(offspring_wide, !pid %in% censored_ids)
    ),

  coxph =
    coxph(
      Surv(enddate, event_ascvd) ~ strata(sex) + log(age_0) + I(log(age_0)^2) + smk_0 + dm_0 + 
        hrx_0 + log(tc_0) + log(hdl_0) + log(sbp_0) + log(age_0):log(tc_0) +
        log(age_0):log(hdl_0) + log(age_0):smk_0 + log(age_0):log(sbp_0) + log(age_0):hrx_0 +
        log(age_0):hrx_0:log(sbp_0) + hrx_0:log(sbp_0),
      data = filter(offspring_coxph,  !pid %in% censored_ids)
    )
)


# Create output tables ----------------------------------------------------

sink("9_results/tables/cov_models_lag0.tex")
texreg(
  l = models_ex1_uncens[["lag0"]][["X"]],
  booktabs = TRUE,
  use.packages = FALSE,
  table = FALSE
) %>% print()
sink()

sink("9_results/tables/cov_models_lag1.tex")
texreg(
  l = models_ex1_uncens[["lag1"]][["X"]],
  booktabs = TRUE,
  use.packages = FALSE,
  table = FALSE
) %>% print()
sink()

sink("9_results/tables/cov_models_lag2.tex")
texreg(
  l = models_ex1_uncens[["lag2"]][["X"]],
  booktabs = TRUE,
  use.packages = FALSE,
  table = FALSE
) %>% print()
sink()

sink("9_results/tables/out_models.tex")
texreg(
  l = list(
    "\\shortstack{(1) \\\\ event\\_chd}" = models_ex1_uncens[["lag0"]][["Y"]],
    "\\shortstack{(2) \\\\ event\\_chd}" = models_ex1_uncens[["lag1"]][["Y"]],
    "\\shortstack{(3) \\\\ event\\_chd}" = models_ex1_uncens[["lag2"]][["Y"]],
    "\\shortstack{(4) \\\\ event\\_dth}" = models_ex1_uncens[["lag0"]][["D"]], 
    "\\shortstack{(5) \\\\ event\\_dth}" = models_ex1_uncens[["lag1"]][["D"]],
    "\\shortstack{(6) \\\\ event\\_dth}" = models_ex1_uncens[["lag2"]][["D"]]),
  booktabs = TRUE,
  use.packages = FALSE,
  table = FALSE
) %>% print()
sink()
