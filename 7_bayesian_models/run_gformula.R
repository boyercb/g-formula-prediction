library(cmdstanr)
library(posterior)
library(bayesplot)

set.seed(950285)

bayes_ids <- analytic_wide$pid
bayes_sample <- filter(analytic_long, pid %in% bayes_ids)

covs_cont <- c("age", "bmi", "tc", "hdl", "sbp")

means <- bayes_sample %>% 
  #filter(time == 0) %>%
  summarise(across(covs_cont, mean)) 

sds <- bayes_sample %>% 
  #filter(time == 0) %>%
  summarise(across(covs_cont, sd))

v_start <- bayes_sample %>% 
  filter(time == 0) %>%
  select(all_of(covs_fixed[!covs_fixed %in% covs_refs]), pid) 

v <- left_join(select(bayes_sample, pid, time), v_start) %>%
     select(-pid, -time)

v_start <- select(v_start, -pid)


# fit parametric joint model in stan --------------------------------------

stan_data <- list(
  N = length(bayes_ids),
  K = 6,
  N_obs = nrow(bayes_sample),
  K_obs = count(bayes_sample, pid) %>% pull(n),
  N_y = sum(!is.na(bayes_sample$event_chd)),
  y_not_missing = which(!is.na(bayes_sample$event_chd)),
  v = v,
  x1 = (bayes_sample$age - means$age) / sds$age,
  x2 = bayes_sample$smk,
  x3 = (bayes_sample$bmi - means$bmi) / sds$bmi,
  x4 = bayes_sample$dm,
  x5 = bayes_sample$hrx,
  x6 = bayes_sample$liprx,
  x7 = (bayes_sample$tc - means$tc) / sds$tc,
  x8 = (bayes_sample$hdl - means$hdl) / sds$hdl,
  x9 = (bayes_sample$sbp - means$sbp) / sds$sbp,
  y = filter(bayes_sample, !is.na(event_chd)) %>% pull(event_chd),
  d = bayes_sample$event_dth,
  t0 = as.numeric(bayes_sample$time == 0),
  t1 = as.numeric(bayes_sample$time == 1),
  t2 = as.numeric(bayes_sample$time == 2),
  t3 = as.numeric(bayes_sample$time == 3),
  t4 = as.numeric(bayes_sample$time == 4),
  t5 = as.numeric(bayes_sample$time == 5)
)

mod <- cmdstan_model("7_bayesian_models/gformula_bayes.stan", compile = FALSE)
mod$compile(cpp_options = list(STAN_THREADS = TRUE, STAN_OPENCL = TRUE))
mod$exe_file()

if (rerun_sample_from_posterior) {
  fit <-
    mod$sample(
      data = stan_data,
      iter_warmup = 500,
      iter_sampling = 2000, 
      parallel_chains = 4,
      seed = 58798235
    )
  
  fit$save_object(file = "../../3_data_encrypted/rds/gformula_bayes_fit_w_d.RDS")
}

fit$cmdstan_diagnose()
fit$cmdstan_summary()

bayes_gformula_draws <- as_draws_df(fit$draws())

# make predictions based on posterior samples -----------------------------

stan_data_mc <- list(
  N = nrow(analytic_wide),
  K = 6,
  N_new = nrow(analytic_wide) * 6,
  v_start = v_start,
  x1_start = (analytic_wide$age_0 - means$age) / sds$age,
  x2_start = analytic_wide$smk_0,
  x3_start = (analytic_wide$bmi_0 - means$bmi) / sds$bmi,
  x4_start = analytic_wide$dm_0,
  x5_start = analytic_wide$hrx_0,
  x6_start = analytic_wide$liprx_0,
  x7_start = (analytic_wide$tc_0 - means$tc) / sds$tc,
  x8_start = (analytic_wide$hdl_0 - means$hdl) / sds$hdl,
  x9_start = (analytic_wide$sbp_0 - means$sbp) / sds$sbp,
  return_survival = FALSE,
  comprisk = TRUE
)

mod_mc <- cmdstan_model("7_bayesian_models/gformula_rng.stan", compile = FALSE)
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

#stanfit <- rstan::read_stan_csv(fit$output_files())



