# create plots ------------------------------------------------------------

copy_baseline_model <- 
  factual_sim_results %>%
  expand(scenario, model, cond, sim) %>%
  filter(model == "conventional (baseline-only)" & cond > 1) %>%
  left_join(
    filter(
      factual_sim_results,
      model == "conventional (baseline-only)" & cond == 1
    ) %>% 
    select(-cond),
  by = c("scenario", "model", "sim"))

plot_df <-
  bind_rows(
    factual_sim_results,
    copy_baseline_model
  ) %>%
  mutate(
    scenario = paste0("scenario: ", scenario)
    ) %>%
  group_by(scenario, model, cond) %>%
  summarise(
    across(Dxy:`S:p`, mean, .names = "{col}_mean"),
    across(Dxy:`S:p`, sd, .names = "{col}_sd"),
    across(Dxy:`S:p`, ~quantile(.x, 0.9, na.rm = TRUE), .names = "{col}_q90"),
    across(Dxy:`S:p`, ~quantile(.x, 0.1, na.rm = TRUE), .names = "{col}_q10")
  ) %>%
  filter(model != "conventional (updated)")

cstat_plot <-
  plot_df %>%
  select(`C (ROC)_mean`,
         `C (ROC)_sd`,
         `C (ROC)_q10`,
         `C (ROC)_q90`,
         scenario,
         model,
         cond) %>%
  group_by(scenario, model, cond) %>%
  pivot_longer(c("C (ROC)_mean", "C (ROC)_sd", "C (ROC)_q10", "C (ROC)_q90")) %>%
  separate(name, c("stat", "type"), sep = "_") %>%
  pivot_wider(names_from = "type") %>%
  mutate(stat = factor(stat, labels = c("C-statistic\n(AUC)")))

brier_plot <-
  plot_df %>%
  select(`Brier_mean`,
         `Brier_sd`,
         `Brier_q10`,
         `Brier_q90`,
         scenario,
         model,
         cond) %>%
  group_by(scenario, model, cond) %>%
  pivot_longer(c("Brier_mean", "Brier_sd", "Brier_q10", "Brier_q90")) %>%
  separate(name, c("stat", "type"), sep = "_") %>%
  pivot_wider(names_from = "type") 

somers_d_plot <-
  plot_df %>%
  select(`Dxy_mean`,
         `Dxy_sd`,
         `Dxy_q10`,
         `Dxy_q90`,
         scenario,
         model,
         cond) %>%
  group_by(scenario, model, cond) %>%
  pivot_longer(c("Dxy_mean", "Dxy_sd", "Dxy_q10", "Dxy_q90")) %>%
  separate(name, c("stat", "type"), sep = "_") %>%
  pivot_wider(names_from = "type") 

R2_plot <-
  plot_df %>%
  select(`R2_mean`,
         `R2_sd`,
         `R2_q10`,
         `R2_q90`,
         scenario,
         model,
         cond) %>%
  group_by(scenario, model, cond) %>%
  pivot_longer(c("R2_mean", "R2_sd", "R2_q10", "R2_q90")) %>%
  separate(name, c("stat", "type"), sep = "_") %>%
  pivot_wider(names_from = "type") 

calibration_plot <-
  plot_df %>%
  select(Slope_mean,
         Slope_sd,
         Slope_q10,
         Slope_q90,
         scenario,
         model,
         cond) %>%
  group_by(scenario, model, cond) %>%
  pivot_longer(c("Slope_mean", "Slope_sd", "Slope_q10", "Slope_q90")) %>%
  separate(name, c("stat", "type"), sep = "_") %>%
  pivot_wider(names_from = "type") 

calibration_int_plot <-
  plot_df %>%
  select(Intercept_mean,
         Intercept_sd,
         Intercept_q10,
         Intercept_q90,
         scenario,
         model,
         cond) %>%
  group_by(scenario, model, cond) %>%
  pivot_longer(c("Intercept_mean", "Intercept_sd", "Intercept_q10", "Intercept_q90")) %>%
  separate(name, c("stat", "type"), sep = "_") %>%
  pivot_wider(names_from = "type") 

calibration_avg_plot <-
  plot_df %>%
  select(Eavg_mean,
         Eavg_sd,
         Eavg_q10,
         Eavg_q90,
         scenario,
         model,
         cond) %>%
  group_by(scenario, model, cond) %>%
  pivot_longer(c("Eavg_mean", "Eavg_sd", "Eavg_q10", "Eavg_q90")) %>%
  separate(name, c("stat", "type"), sep = "_") %>%
  pivot_wider(names_from = "type") 


pdf("9_results/figures/factual_sims_cstat_plot.pdf", width = 16, height = 8)
ggplot(cstat_plot, aes(x = factor(cond), y = mean, fill = model, group = model)) +
  facet_wrap(~ scenario, scales = "free_y") +
  geom_ribbon(aes(
    ymin = q10,
    ymax = q90
  ), alpha = 0.35) +
  geom_point(aes(color = model), size = 2) +
  geom_line(aes(color = model)) +
  labs(
    x = "Conditioning on k time points",
    y = "C-statistic\n(AUC)"
  ) +
  coord_cartesian(expand = FALSE) +
  scale_fill_brewer(name = "", type = "qual", palette = "Set2") +
  scale_color_brewer(name = "", type = "qual", palette = "Set2") +
  slides_theme() +
  theme(
    text = element_text(size = 14),
    strip.text.y = element_text(angle = 0)
  )
dev.off()

pdf("9_results/figures/factual_sims_brier_plot.pdf", width = 16, height = 8)
ggplot(brier_plot, aes(x = factor(cond), y = mean, fill = model, group = model)) +
  facet_wrap(~ scenario, scales = "free_y") +
  geom_ribbon(aes(
    ymin = q10,
    ymax = q90
  ), alpha = 0.35) +
  geom_point(aes(color = model), size = 2) +
  geom_line(aes(color = model)) +
  labs(
    x = "Conditioning on k time points",
    y = "Brier score"
  ) +
  coord_cartesian(expand = FALSE) +
  scale_fill_brewer(name = "", type = "qual", palette = "Set2") +
  scale_color_brewer(name = "", type = "qual", palette = "Set2") +
  slides_theme() +
  theme(
    text = element_text(size = 14),
    strip.text.y = element_text(angle = 0)
  )
dev.off()


pdf("9_results/figures/factual_sims_calibration_plot.pdf", width = 16, height = 8)
ggplot(calibration_avg_plot, aes(x = factor(cond), y = mean, fill = model, group = model)) +
  facet_wrap(~ scenario) +
  geom_ribbon(aes(
    ymin = q10,
    ymax = q90
  ), alpha = 0.35) +
  geom_point(aes(color = model), size = 2) +
  geom_line(aes(color = model)) +
  labs(
    x = "Conditioning on k time points",
    y = "Average absolute\ncalibration error"
  ) +
  coord_cartesian(expand = FALSE) +
  scale_fill_brewer(name = "", type = "qual", palette = "Set2") +
  scale_color_brewer(name = "", type = "qual", palette = "Set2") +
  slides_theme() +
  theme(
    text = element_text(size = 14),
    strip.text.y = element_text(angle = 0)
  )
dev.off()


pdf("9_results/figures/factual_sims_somers_d_plot.pdf", width = 16, height = 8)
ggplot(somers_d_plot, aes(x = factor(cond), y = mean, fill = model, group = model)) +
  facet_wrap(~ scenario, scales = "free_y") +
  geom_ribbon(aes(
    ymin = q10,
    ymax = q90
  ), alpha = 0.35) +
  geom_point(aes(color = model), size = 2) +
  geom_line(aes(color = model)) +
  labs(
    x = "Conditioning on k time points",
    y = "C-statistic\n(AUC)"
  ) +
  coord_cartesian(expand = FALSE) +
  scale_fill_brewer(name = "", type = "qual", palette = "Set2") +
  scale_color_brewer(name = "", type = "qual", palette = "Set2") +
  slides_theme() +
  theme(
    text = element_text(size = 14),
    strip.text.y = element_text(angle = 0)
  )
dev.off()

pdf("9_results/figures/factual_sims_R2_plot.pdf", width = 16, height = 8)
ggplot(R2_plot, aes(x = factor(cond), y = mean, fill = model, group = model)) +
  facet_wrap(~ scenario, scales = "free_y") +
  geom_ribbon(aes(
    ymin = q10,
    ymax = q90
  ), alpha = 0.35) +
  geom_point(aes(color = model), size = 2) +
  geom_line(aes(color = model)) +
  labs(
    x = "Conditioning on k time points",
    y = "C-statistic\n(AUC)"
  ) +
  coord_cartesian(expand = FALSE) +
  scale_fill_brewer(name = "", type = "qual", palette = "Set2") +
  scale_color_brewer(name = "", type = "qual", palette = "Set2") +
  slides_theme() +
  theme(
    text = element_text(size = 14),
    strip.text.y = element_text(angle = 0)
  )
dev.off()

pdf("9_results/figures/factual_sims_calibration_slope.pdf", width = 16, height = 8)
ggplot(calibration_plot, aes(x = factor(cond), y = mean, fill = model, group = model)) +
  facet_wrap(~ scenario) +
  geom_ribbon(aes(
    ymin = q10,
    ymax = q90
  ), alpha = 0.35) +
  geom_point(aes(color = model), size = 2) +
  geom_line(aes(color = model)) +
  labs(
    x = "Conditioning on k time points",
    y = "Average absolute\ncalibration error"
  ) +
  coord_cartesian(expand = FALSE) +
  scale_fill_brewer(name = "", type = "qual", palette = "Set2") +
  scale_color_brewer(name = "", type = "qual", palette = "Set2") +
  slides_theme() +
  theme(
    text = element_text(size = 14),
    strip.text.y = element_text(angle = 0)
  )
dev.off()

pdf("9_results/figures/factual_sims_calibration_intercept.pdf", width = 16, height = 8)
ggplot(calibration_int_plot, aes(x = factor(cond), y = mean, fill = model, group = model)) +
  facet_wrap(~ scenario) +
  geom_ribbon(aes(
    ymin = q10,
    ymax = q90
  ), alpha = 0.35) +
  geom_point(aes(color = model), size = 2) +
  geom_line(aes(color = model)) +
  labs(
    x = "Conditioning on k time points",
    y = "Average absolute\ncalibration error"
  ) +
  coord_cartesian(expand = FALSE) +
  scale_fill_brewer(name = "", type = "qual", palette = "Set2") +
  scale_color_brewer(name = "", type = "qual", palette = "Set2") +
  slides_theme() +
  theme(
    text = element_text(size = 14),
    strip.text.y = element_text(angle = 0)
  )
dev.off()