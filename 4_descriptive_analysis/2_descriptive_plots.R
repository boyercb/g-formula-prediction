# bivariate functional relationships

# Outcome: CHD events -----------------------------------------------------

# with raw data
pdf("9_results/figures/event_chd_functional_forms_raw.pdf", width = 7, height = 4)
offspring_long %>%
  select(all_of(vars_cont), all_of(dvs)) %>%
  pivot_longer( 
    cols = all_of(vars_cont), 
    names_to = "variable",
    values_to = "value"
  ) %>%
  ggplot(., aes(x = value, y = event_ascvd)) +
  facet_wrap(~variable, scales = "free") +
  geom_jitter(
    aes(y = 1.1 * event_ascvd + (1 - event_ascvd) * -0.1),
    height = 0.05,
    alpha = 0.2,
    shape = 16
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    y = "Probability of CHD event",
    x = ""
  ) +
  slides_theme()
dev.off()

# zoomed in on functional relationship
pdf("9_results/figures/event_chd_functional_forms.pdf", width = 7, height = 4)
offspring_long %>%
  filter(time == 0) %>%
  select(all_of(vars_cont), all_of(dvs)) %>%
  pivot_longer( 
    cols = all_of(vars_cont), 
    names_to = "variable",
    values_to = "value"
  ) %>%
  ggplot(., aes(x = value, y = event_ascvd)) +
  facet_wrap(~variable, scales = "free") +
  geom_smooth(method = "loess") +
  labs(
    y = "Probability of CHD event",     
    x = ""
  ) +
  slides_theme()
dev.off()

# comparison on and not on anti-hypertensive treatment
pdf("9_results/figures/event_chd_functional_forms_bp_treament.pdf", width = 7, height = 4)
offspring_long %>%
  select(c("hdl", "tc", "sbp", "hrx", "liprx"), all_of(dvs)) %>%
  pivot_longer( 
    cols = c("hdl", "tc", "sbp"), 
    names_to = "variable",
    values_to = "value"
  ) %>%
  ggplot(., aes(x = value, y = event_chd)) +
  facet_wrap(hrx~variable, scales = "free", labeller = label_both) +
  geom_smooth(method = "loess") +
  labs(
    y = "Probability of CHD event",    
    x = ""
  ) +
  slides_theme()
dev.off()

# comparison on and not on lipid lowering medications
pdf("9_results/figures/event_chd_functional_forms_lipid_treament.pdf", width = 7, height = 4)
offspring_long %>%
  select(c("hdl", "tc", "sbp", "hrx", "liprx"), all_of(dvs)) %>%
  pivot_longer( 
    cols = c("hdl", "tc", "sbp"), 
    names_to = "variable",
    values_to = "value"
  ) %>%
  ggplot(., aes(x = log(value), y = event_chd)) +
  facet_wrap(liprx~variable, scales = "free", labeller = label_both) +
  geom_smooth(method = "loess") +
  labs(
    y = "Probability of CHD event",     
    x = ""
  ) +
  slides_theme()
dev.off()


# Outcome: Non-CHD death  -------------------------------------------------

# with raw data
pdf("9_results/figures/event_dth_functional_forms_raw.pdf", width = 7, height = 4)
offspring_long %>%
  select(all_of(vars_cont), all_of(dvs)) %>%
  pivot_longer( 
    cols = all_of(vars_cont), 
    names_to = "variable",
    values_to = "value"
  ) %>%
  ggplot(., aes(x = value, y = event_dth)) +
  facet_wrap(~variable, scales = "free") +
  geom_jitter(
    aes(y = 1.1 * event_dth + (1 - event_dth) * -0.1),
    height = 0.05,
    alpha = 0.2,
    shape = 16
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    y = "Probability of death from non-CHD cause",     
    x = ""
  ) +
  slides_theme()
dev.off()

# zoomed in on functional relationship
pdf("9_results/figures/event_dth_functional_forms.pdf", width = 7, height = 4)
offspring_long %>%
  select(all_of(vars_cont), all_of(dvs)) %>%
  pivot_longer( 
    cols = all_of(vars_cont), 
    names_to = "variable",
    values_to = "value"
  ) %>%
  ggplot(., aes(x = value, y = event_dth)) +
  facet_wrap(~variable, scales = "free") +
  geom_smooth(method = "loess") +
  labs(
    y = "Probability of death from non-CHD cause",     
    x = ""
  ) +
  slides_theme()
dev.off()