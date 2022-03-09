bind_rows(
  boot_stats_accaha,
  boot_stats_conventional,
  boot_stats_ex1,
  boot_stats_ex2,
  boot_stats_ex3
) %>%
  filter(test == 1) %>%
  separate(model, c("model", "compevent"), ",") %>%
  mutate(
    compevent = if_else(
      compevent == " ASCVD-specific", 
      " ASCVD-specific risk",
      compevent
    ),
    model = fct_rev(factor(
      model,
      levels = c(
        "PCE",
        "PCE (calibrated)",
        "FRS",
        "FRS (calibrated)",
        "logistic",
        "coxph",
        "no lags",
        "one lag",
        "two lags",
        "bmi",
        "bmi + ldl",
        "bmi + ldl + liprx",
        "gam"
      )
    )),
    Brier = 1 - (Brier / 0.079),
    compevent = factor(trimws(compevent), levels = c("net risk",
                                                     "ASCVD-specific risk"))
  ) %>%
  filter(model != "logistic") %>%
  group_by(model, compevent) %>%
  summarise(across(
    .cols = c(-bsim,-test),
    .fns = list(
      mean = ~ mean(.x, na.rm = TRUE),
      # lwr = function(x) mean - 1.96 * sd(x, na.rm = TRUE),
      # upr = function(x) mean + 1.96 * sd(x, na.rm = TRUE)
      lwr = ~ quantile(.x, 0.025, na.rm = TRUE),
      upr = ~ quantile(.x, 0.975, na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  select(model,
         compevent,
         starts_with("C (ROC)"),
         starts_with("R2"),
         starts_with("Brier")) %>%
  pivot_longer(cols = c(
    starts_with("C (ROC)"),
    starts_with("R2"),
    starts_with("Brier")
  )) %>%
  separate(name, c("name", "type"), sep = "_") %>%
  pivot_wider(names_from = type, values_from = value) %>%
  ggplot(., aes(x = mean, y = model)) +
  facet_grid(compevent ~ name, scales = 'free_x') +
  geom_pointrange(aes(xmin = lwr, xmax = upr)) +
  geom_text(aes(label = specd(mean, 2)), nudge_y = 0.4, size = 3) +
  labs(x = "") +
  ggpubr::theme_pubclean()


bind_rows(
  boot_stats_accaha,
  boot_stats_conventional,
  boot_stats_ex1,
  boot_stats_ex2,
  boot_stats_ex3
) %>%
  filter(test == 1) %>%
  separate(model, c("model", "compevent"), ",") %>%
  mutate(
    compevent = if_else(
      compevent == " ASCVD-specific", 
      " ASCVD-specific risk",
      compevent
    ),
    model = fct_rev(factor(
      model,
      levels = c(
        "PCE",
        "PCE (calibrated)",
        "FRS",
        "FRS (calibrated)",
        "logistic",
        "coxph",
        "no lags",
        "one lag",
        "two lags",
        "bmi",
        "bmi + ldl",
        "bmi + ldl + liprx",
        "gam"
      )
    )),
    Brier = 1 - (Brier / 0.079),
    compevent = factor(trimws(compevent), levels = c("net risk",
                                                     "ASCVD-specific risk"))
  ) %>%
  filter(model != "logistic") %>%
  group_by(model, compevent) %>%
  summarise(across(
    .cols = c(-bsim,-test),
    .fns = list(
      mean = ~ mean(.x, na.rm = TRUE),
      # lwr = function(x) mean - 1.96 * sd(x, na.rm = TRUE),
      # upr = function(x) mean + 1.96 * sd(x, na.rm = TRUE)
      lwr = ~ quantile(.x, 0.025, na.rm = TRUE),
      upr = ~ quantile(.x, 0.975, na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  select(model,
         compevent,
         starts_with("Slope"),
         starts_with("Intercept")) %>%
  pivot_longer(cols = c(
    starts_with("Slope"),
    starts_with("Intercept")
  )) %>%
  separate(name, c("name", "type"), sep = "_") %>%
  pivot_wider(names_from = type, values_from = value) %>%
  ggplot(., aes(x = mean, y = model)) +
  facet_grid(compevent ~ name, scales = 'free_x') +
  geom_pointrange(aes(xmin = lwr, xmax = upr)) +
  geom_text(aes(label = specd(mean, 2)), nudge_y = 0.4, size = 3) +
  labs(x = "") +
  geom_vline(aes(xintercept = rep(c(1, 0), 24)), linetype = "dashed") +
  ggpubr::theme_pubclean()




df <- gformula_predictions_ex2$`bmi + ldl + liprx, ASCVD-specific risk` %>%
    group_by(pid) %>%
    mutate(
      ascvd_10yr_accaha_risk_lag2 = lag(ascvd_10yr_accaha_risk, 2),
      ascvd_10yr_frs_risk_lag2 = lag(ascvd_10yr_frs_risk, 2),
      ascvd_10yr_frs_simple_risk_lag2 = lag(ascvd_10yr_frs_simple_risk, 2)
    ) %>%
    ungroup() %>%
    filter(time == 2)
plot(df$pred, df$ascvd_10yr_accaha_risk_lag2)


ggplot(df, aes(x = ascvd_10yr_accaha_risk_lag2, y = pred)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_vline(xintercept = 0.075, linetype = "dashed") +
  geom_hline(yintercept = 0.075, linetype = "dashed") +
  geom_point(shape = 1, alpha = 0.25) + 
  geom_rug(alpha = 0.1, position = "jitter") +
  geom_smooth(se = FALSE) +
  #coord_fixed() +
  labs(
    x = "\nACC/AHA PCE 2013 10-year cumulative risk",
    y = "G-formula 10-year cumulative risk\n"
  ) +
  ggpubr::theme_pubr()


df %>%
  pivot_longer(cols = c(ascvd_10yr_accaha_risk_lag2, pred)) %>%
ggplot(., aes(x = value, fill = name)) +
  geom_histogram(bins = 50, position = "identity", alpha = 0.4) + 
  scale_fill_brewer(
    name = "",
    palette = "Set1",
    labels = c(
      "ACC/AHA PCE 2013",
      "G-formula"
    )
  ) +
  labs(
    x = "\npredicted risk",
    y = "count\n"
  ) +
  geom_vline(xintercept = 0.075, linetype = "dashed") +
  ggpubr::theme_pubr() +
  theme(legend.position = c(0.3, 0.6))
