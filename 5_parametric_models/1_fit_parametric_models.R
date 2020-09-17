
# Fit parametric models ---------------------------------------------------

Y.fit <- 
  fit_Y_model(
    Y.model,
    data = analytic_long
  )

X.fit <- 
  fit_X_models(
    X.models,
    data = analytic_long,
    time = "time"
  )

D.fit <- 
  fit_D_model(
    D.model,
    data = analytic_long
  )

conventional.fit_0_2 <- 
  glm(
    formula = event_chd_3 ~ sex + age0 + educ_1 + educ_2 + educ_3 +
      marital_1 + marital_2 + eversmk + pre_dpd + pre_bmi + pre_dm + pre_sbp +
      pre_cpd + pre_ldl + pre_hrx + pre_liprx + cpd_0 + dpd_0 + bmi_0 +
      dm_0 + sbp_0 + ldl_0 + hrx_0 + liprx_0,
    family = binomial(link = 'logit'),
    data = analytic_wide
  )

conventional.fit_0_6 <- 
  glm(
    formula = event_chd_6 ~ sex + age0 + educ_1 + educ_2 + educ_3 +
      marital_1 + marital_2 + eversmk + pre_dpd + pre_bmi + pre_dm + pre_sbp +
      pre_cpd + pre_ldl + pre_hrx + pre_liprx + cpd_0 + dpd_0 + bmi_0 +
      dm_0 + sbp_0 + ldl_0 + hrx_0 + liprx_0,
    family = binomial(link = 'logit'),
    data = analytic_wide
  )


# Create output tables ----------------------------------------------------

sink("9_results/tables/cov_models.tex")
texreg(
  l = X.fit,
  booktabs = TRUE,
  use.packages = FALSE,
  table = FALSE
) %>% print()
sink()

sink("9_results/tables/out_models.tex")
texreg(
  l = list("event_chd" = Y.fit, "event_dth" = D.fit),
  booktabs = TRUE,
  use.packages = FALSE,
  table = FALSE
) %>% print()
sink()