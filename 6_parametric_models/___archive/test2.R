# simulating full

stan_data_ind <- list(
  N = 1,
  K = 6,
  N_new = 6,
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

mod_ind <- cmdstan_model("5_parametric_models/gformula_rng_ind.stan", compile = FALSE)
mod_ind$compile(cpp_options = list(STAN_THREADS = TRUE, STAN_OPENCL = TRUE))
mod_ind$exe_file()

fit_ind <-
  mod_ind$generate_quantities(
    fitted_params = fit,
    data = stan_data_ind,
    parallel_chains = 4,
    seed = 58798235
  )

labs <- c(
  "p_x1_new" = "Smoker (0/1)",
  "x2_new" = "Cigarettes per day",
  "p_x3_new" = "Drinker (0/1)",
  "x4_new" = "Drinks per day",
  "x5_new" = "BMI",
  "p_x6_new" = "Diabetes (0/1)",
  "x7_new" = "SBP",
  "x8_new" = "LDL", 
  "p_x9_new" = "Anti-hypertensive meds (0/1)",
  "p_x10_new" = "Anti-lipid meds (0/1)"
)
preds <- as_draws_df(fit_ind$draws())

preds %>%
  select(-`.chain`, -`.iteration`, -`.draw`) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("var", "time"),
    names_pattern = "(.*)\\[(.)\\]"
  ) %>%
  filter(str_detect(var, "^(p_x[0-9]+)|^(x[24578]+)")) %>%
  group_by(var, time) %>%
  filter(time > 1) %>%
  summarise(
    pred_05 = quantile(value,  prob = 0.05),
    pred_10 = quantile(value,  prob = 0.10),
    pred_25 = quantile(value,  prob = 0.25),
    pred_50 = median(value),
    pred_75 = quantile(value,  prob = 0.75),
    pred_90 = quantile(value,  prob = 0.90),
    pred_95 = quantile(value,  prob = 0.95)
  ) %>%
  ungroup() %>%
  mutate(var = fct_relevel(
    var,
    "p_x1_new",
    "x2_new",
    "p_x3_new",
    "x4_new",
    "x5_new",
    "p_x6_new",
    "x7_new",
    "x8_new", 
    "p_x9_new",
    "p_x10_new")) %>%
  ggplot(., aes(x = time, y = pred_50, group = var)) + 
  facet_wrap(~var, scales = "free_y", nrow = 3, labeller = labeller(var = labs)) +
  geom_ribbon(aes(ymin = pred_05, ymax = pred_95), fill = "#d1e1ec", alpha = 0.7) +
  geom_ribbon(aes(ymin = pred_10, ymax = pred_90), fill = "#6497b1", alpha = 0.7) +
  geom_ribbon(aes(ymin = pred_25, ymax = pred_75), fill = "#011f4b", alpha = 0.7) +
  geom_line(linetype = "dashed") +
  labs(
    x = "exam",
    y = ""
  ) + 
  slides_theme()



preds %>%
  select(starts_with("p_new")) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("var", "time"),
    names_pattern = "(.*)\\[(.)\\]"
  ) %>%
  group_by(var, time) %>%
  summarise(
    pred_05 = quantile(value,  prob = 0.05),
    pred_10 = quantile(value,  prob = 0.10),
    pred_25 = quantile(value,  prob = 0.25),
    pred_50 = median(value),
    pred_75 = quantile(value,  prob = 0.75),
    pred_90 = quantile(value,  prob = 0.90),
    pred_95 = quantile(value,  prob = 0.95)
  ) %>%
  ggplot(., aes(x = time, y = pred_50, group = var)) + 
  geom_ribbon(aes(ymin = pred_05, ymax = pred_95), fill = "#DCBCBC") +
  geom_ribbon(aes(ymin = pred_10, ymax = pred_90), fill = "#B97C7C") +
  geom_ribbon(aes(ymin = pred_25, ymax = pred_75), fill = "#7C0000") +
  geom_line(linetype = "dashed") +
  labs(
    x = "exam",
    y = "cumulative risk of CHD\n"
  ) + 
  slides_theme()

preds %>%
  select(starts_with("p_new")) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("var", "time"),
    names_pattern = "(.*)\\[(.)\\]"
  ) %>%
  ggplot(., aes(x = time, y = value, group = var)) + 
  stat_gradientinterval(show_point = FALSE,
                        show_interval = FALSE, fill = "skyblue") +
  #geom_bin2d(binwidth = c(0.001, 0.001))+
  #stat_density_2d(aes(fill = ..density..), geom = "polygon", contour = FALSE) +
  stat_summary(fun = mean) +
  geom_smooth(
    linetype = "dashed",
    method = "lm",
    formula = y ~ factor(x),
    se = F,
    color = "black"
  ) +
  scale_fill_distiller(direction = 1, palette = "Reds") +
  labs(
    x = "exam",
    y = "cumulative risk of CHD\n"
  ) + 
  slides_theme() +
  theme(legend.position = "none")