# create plots ------------------------------------------------------------

risk_factor_results_noDC <- bind_rows(risk_factor_results_noDC, .id = "iter") 
risk_factor_results_noDC$num <- as.numeric(risk_factor_results_noDC$num) - 1
plot_df <-
  risk_factor_results_noDC %>%
  group_by(model, num) %>%
  summarise(
    across(Dxy:`S:p`, mean, .names = "{col}_mean"),
    across(Dxy:`S:p`, sd, .names = "{col}_se")
  )

discrimination_plot <-
  plot_df %>%
  select(`C (ROC)_mean`,
         `C (ROC)_se`,
         Brier_mean,
         Brier_se,
         model,
         num) %>%
  group_by(model, num) %>%
  pivot_longer(c("C (ROC)_mean", "C (ROC)_se", "Brier_mean", "Brier_se")) %>%
  separate(name, c("stat", "type"), sep = "_") %>%
  pivot_wider(names_from = "type") %>%
  mutate(stat = factor(stat, labels = c("Brier\nscore", "C-statistic\n(AUC)")))

calibration_plot <-
  plot_df %>%
  select(Slope_mean,
         Slope_se,
         Intercept_mean,
         Intercept_se,
         model,
         num) %>%
  group_by(model, num) %>%
  pivot_longer(c("Slope_mean", "Slope_se", "Intercept_mean", "Intercept_se")) %>%
  separate(name, c("stat", "type"), sep = "_") %>%
  pivot_wider(names_from = "type") %>%
  mutate(stat = factor(stat, labels = c("Intercept", "Slope")))


pdf("9_results/figures/risk_factor_discrimination_plot.pdf", width = 11, height = 8)
ggplot(discrimination_plot, aes(x = factor(num), y = mean, color = model, fill = model, group = model)) +
  facet_wrap(~ fct_rev(stat), scales = "free_y", as.table = FALSE) +
  geom_ribbon(aes(
    ymin = mean - 1.96 * se,
    ymax = mean + 1.96 * se
  ), alpha = 0.4) +
  geom_point() +
  geom_line() +
  labs(
    x = "Time (k)",
    y = ""
  ) +
  coord_cartesian(expand = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  slides_theme() +
  theme(
    text = element_text(size = 14),
    strip.text.y = element_text(angle = 0)
  )
dev.off()


pdf("9_results/figures/risk_factor_calibration_plot.pdf", width = 11, height = 8)
ggplot(calibration_plot, aes(x = factor(num), y = mean, color = model, fill = model, group = model)) +
  facet_wrap(~ stat, scales = "free_y", as.table = FALSE) +
  geom_ribbon(aes(
    ymin = mean - 1.96 * se,
    ymax = mean + 1.96 * se
  ), alpha = 0.4) +
  geom_point() +
  geom_line() +
  labs(
    x = "Time (k)",
    y = ""
  ) +
  coord_cartesian(expand = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  slides_theme() +
  theme(
    text = element_text(size = 14),
    strip.text.y = element_text(angle = 0)
  )
dev.off()
