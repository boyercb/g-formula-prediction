
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