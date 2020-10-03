library(cmdstanr)
library(posterior)
library(bayesplot)

covs_cont <- c("cpd", "dpd", "bmi", "sbp", "ldl")
covs_fixed_cont <- c("age0", "pre_dpd", "pre_cpd", "pre_bmi", "pre_sbp", "pre_ldl")

means <- analytic_long %>% 
  filter(time == 0) %>%
  summarise(across(covs_cont, mean)) 

sds <- analytic_long %>% 
  filter(time == 0) %>%
  summarise(across(covs_cont, sd))

v_start <- analytic_long %>% 
  filter(time == 0) %>%
  select(all_of(covs_fixed[!covs_fixed %in% covs_refs]), pid) %>%
  mutate(across(covs_fixed_cont, function(x) (x - mean(x)) / sd(x)))

v <- left_join(select(analytic_long, pid, time), v_start) %>%
     select(-pid, -time)

v_start <- select(v_start, -pid)


# fit parametric joint model in stan --------------------------------------

stan_data <- list(
  N = nrow(analytic_wide),
  K = 6,
  N_obs = nrow(analytic_long),
  K_obs = count(analytic_long, pid) %>% pull(n),
  N_new = nrow(analytic_wide) * 6,
  v = v,
  x1 = analytic_long$smk,
  x2 = (analytic_long$cpd - means$cpd) / sds$cpd,
  x3 = analytic_long$drk,
  x4 = (analytic_long$dpd - means$dpd) / sds$dpd,
  x5 = (analytic_long$bmi - means$bmi) / sds$bmi,
  x6 = analytic_long$dm,
  x7 = (analytic_long$sbp - means$sbp) / sds$sbp,
  x8 = (analytic_long$ldl - means$ldl) / sds$ldl,
  x9 = analytic_long$hrx,
  x10 = analytic_long$liprx,
  y = analytic_long$event_chd,
  t1 = as.numeric(analytic_long$time == 1),
  t2 = as.numeric(analytic_long$time == 2),
  t3 = as.numeric(analytic_long$time == 3),
  t4 = as.numeric(analytic_long$time == 4),
  t5 = as.numeric(analytic_long$time == 5)
)

mod <- cmdstan_model("5_parametric_models/gformula_bayes.stan", compile = FALSE)
mod$compile(cpp_options = list(STAN_THREADS = TRUE, STAN_OPENCL = TRUE))
mod$exe_file()

if (rerun_sample_from_posterior) {
  fit <-
    mod$sample(
      data = stan_data,
      iter_warmup = 1000,
      iter_sampling = 2000, 
      parallel_chains = 4,
      seed = 58798235
    )
  
  fit$save_object(file = get_data("gformula_bayes_fit.RDS"))
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
  x1_start = analytic_wide$smk_0,
  x2_start = (analytic_wide$cpd_0 - means$cpd) / sds$cpd,
  x3_start = analytic_wide$drk_0,
  x4_start = (analytic_wide$dpd_0 - means$dpd) / sds$dpd,
  x5_start = (analytic_wide$bmi_0 - means$bmi) / sds$bmi,
  x6_start = analytic_wide$dm_0,
  x7_start = (analytic_wide$sbp_0 - means$sbp) / sds$sbp,
  x8_start = (analytic_wide$ldl_0 - means$ldl) / sds$ldl,
  x9_start = analytic_wide$hrx_0,
  x10_start = analytic_wide$liprx_0
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

#stanfit <- rstan::read_stan_csv(fit$output_files())

