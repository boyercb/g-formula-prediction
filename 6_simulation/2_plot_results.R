# create plots ------------------------------------------------------------

simple_tvc_results <- bind_rows(simple_tvc_results, .id = "iter") 

plot_df <-
  simple_tvc_results %>%
  group_by(model, test_data, a2) %>%
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
         test_data,
         a2) %>%
  group_by(model, test_data, a2) %>%
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
         test_data,
         a2) %>%
  group_by(model, test_data, a2) %>%
  pivot_longer(c("Slope_mean", "Slope_se", "Intercept_mean", "Intercept_se")) %>%
  separate(name, c("stat", "type"), sep = "_") %>%
  pivot_wider(names_from = "type") %>%
  mutate(stat = factor(stat, labels = c("Intercept", "Slope")))


pdf("9_results/figures/simulation_discrimination_plot.pdf", width = 11, height = 8)
ggplot(discrimination_plot, aes(x = a2, y = mean, color = model, fill = model)) +
  facet_grid(stat ~ test_data, scales = "free_y", as.table = FALSE) +
  geom_ribbon(aes(
    ymin = mean - 1.96 * se,
    ymax = mean + 1.96 * se
  ), alpha = 0.4) +
  geom_point() +
  geom_line() +
  labs(
    x = bquote(alpha[2]),
    y = ""
  ) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  latex_theme() +
  theme(
    text = element_text(size = 14),
    strip.text.y = element_text(angle = 0)
  )
dev.off()


pdf("9_results/figures/simulation_calibration_plot.pdf", width = 11, height = 8)
ggplot(calibration_plot, aes(x = a2, y = mean, color = model, fill = model)) +
  facet_grid(stat ~ test_data, scales = "free_y", as.table = FALSE) +
  geom_ribbon(aes(
    ymin = mean - 1.96 * se,
    ymax = mean + 1.96 * se
  ), alpha = 0.4) +
  geom_point() +
  geom_line() +
  labs(
    x = bquote(alpha[2]),
    y = ""
  ) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  latex_theme() +
  theme(
    text = element_text(size = 14),
    strip.text.y = element_text(angle = 0)
  )
dev.off()
