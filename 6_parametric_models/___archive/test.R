stan_data_mc <- list(
  N = nrow(analytic_wide),
  K = 6,
  N_new = nrow(analytic_wide) * 6,
  v_start = v,
  x1_start = analytic_long$smk,
  x2_start = (analytic_long$cpd - means$cpd) / sds$cpd,
  x3_start = analytic_long$drk,
  x4_start = (analytic_long$dpd - means$dpd) / sds$dpd,
  x5_start = (analytic_long$bmi - means$bmi) / sds$bmi,
  x6_start = analytic_long$dm,
  x7_start = (analytic_long$sbp - means$sbp) / sds$sbp,
  x8_start = (analytic_long$ldl - means$ldl) / sds$ldl,
  x9_start = analytic_long$hrx,
  x10_start = analytic_long$liprx,
  return_survival = FALSE,
  comprisk = TRUE
)

mod_mc <- cmdstan_model("5_parametric_models/gformula_rng.stan", compile = FALSE)
mod_mc$compile(cpp_options = list(STAN_THREADS = TRUE, STAN_OPENCL = TRUE))
mod_mc$exe_file()

fit_mc <-
  mod_mc$generate_quantities(
    fitted_params = fit,
    data = stan_data_mc,
    parallel_chains = 4,
    seed = 58798235
  )

preds <- as_draws_df(fit_mc$draws())
preds <- select(preds, starts_with("pred"))

y.hat <- tibble(
  pid = rep(analytic_wide$pid, each = 6),
  time = rep(0:5, nrow(analytic_wide)),
  pred = map_dbl(preds, median)
)

map(0:5,  function (y) {
  p.name <- "pred"
  y.name <- paste0("event_chd_", y + 1)
  
  df <- left_join(filter(y.hat, time == y),
                  select(observed_events, "pid", y.name))
  
  rms::val.prob(df[[p.name]], df[[y.name]])
})
