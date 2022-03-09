# fit random forest model -------------------------------------------------

rf_nc_lag1 <- 
  gformula_mc(
    Y.fit = Y.rf_lag1,
    X.fit = X.rf_lag1,
    D.fit = D.rf_lag1,
    data = analytic_long,
    id = "pid",
    time = "time",
    base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
    hist.vars = paste0("lag1_", covs_tv),
    hist.fun = "lag",
    mc.sims = 10,
    mc.start = 0,
    mc.stop = 5,
    merge = FALSE,
    pool = TRUE
  )

rf_nc_lag2 <- 
  gformula_mc(
    Y.fit = Y.rf_lag2,
    X.fit = X.rf_lag2,
    D.fit = D.rf_lag2,
    data = analytic_long,
    id = "pid",
    time = "time",
    base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
    hist.vars = c(paste0("lag1_", covs_tv), paste0("lag2_", covs_tv)),
    hist.fun = "lag",
    mc.sims = 10,
    mc.start = 0,
    mc.stop = 5,
    merge = FALSE,
    pool = TRUE
  )

# clean up results
rf_nc_lag1 <- 
  rf_nc_lag1 %>%
    rename(event_chd = pred,
           event_chd_alt = Py,
           event_dth = Pd) %>%
    group_by(pid, time) %>%
    summarise(
      across(
        c(all_of(covs_tv), event_chd, event_chd_alt, event_dth),
        mean, 
        .names = "rf_{col}"
      ),
      .groups = "drop_last") %>%
    mutate(
      rf_event_dth = cumsum(rf_event_dth),
      rf_event_chd_alt = cumsum(rf_event_chd_alt)
    )

# clean up results
rf_nc_lag2 <- 
  rf_nc_lag2 %>%
  rename(event_chd = pred,
         event_chd_alt = Py,
         event_dth = Pd) %>%
  group_by(pid, time) %>%
  summarise(
    across(
      c(all_of(covs_tv), event_chd, event_chd_alt, event_dth),
      mean, 
      .names = "rf_{col}"
    ),
    .groups = "drop_last") %>%
  mutate(
    rf_event_dth = cumsum(rf_event_dth),
    rf_event_chd_alt = cumsum(rf_event_chd_alt)
  )

observed_events_covs <-
  left_join(mutate(observed_events, time = as.numeric(time)),
            select(analytic_long, all_of(covs_tv), pid, time),
            by = c("pid", "time"))


# Bootstrap differences in natural course ---------------------------------


BOOT_SMPS <- 500

pb <- txtProgressBar(max = BOOT_SMPS, initial = NA, style = 3)

nc_boots <- 
  map(
    .x = 1:BOOT_SMPS, 
    .f = function(x, pb) {
      setTxtProgressBar(pb, x)
      
      ids <- unique(analytic_long$pid)
      boot_ids <- sample(ids, size = length(ids), replace = TRUE)
      
      boot_rows <-
        sapply(boot_ids, function(x)
          which(analytic_long$pid == x))
      
      boot_rows_obs <-
        sapply(boot_ids, function(x)
          which(observed_events_covs$pid == x),
          simplify = F)
      
      boot_rows_rf <-
        sapply(boot_ids, function(x)
          which(rf_nc_lag1$pid == x),
          simplify = F)
      
      boot_times <-
        sapply(boot_ids, function(x)
          length(which(analytic_long$pid == x)))
      
      boot_df <- analytic_long[unlist(boot_rows), ]
      boot_df_obs <- observed_events_covs[unlist(boot_rows_obs), ]
      boot_df_rf_lag1 <- rf_nc_lag1[unlist(boot_rows_rf),]
      
      
      boot_df$id <- rep(1:length(boot_ids), times = boot_times)
      boot_df_rf$id <- rep(1:length(boot_ids), each = 6)
      boot_df_obs$id <- rep(1:length(boot_ids), each = 6)
      
      Y.fit_boot <- 
        fit_Y_model(
          Y.model,
          data = filter(boot_df, time < 6)
        )
      
      X.fit_boot <- 
        fit_X_models(
          X.models,
          data = filter(boot_df, time < 6),
          time = "time"
        )
      
      D.fit_boot <- 
        fit_D_model(
          D.model,
          data = filter(boot_df, time < 6)
        )
      
      
      natural_course <- 
        gformula_mc(
          Y.fit = Y.fit_boot,
          X.fit = X.fit_boot,
          D.fit = D.fit_boot,
          data = boot_df,
          id = "id",
          time = "time",
          base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
          hist.vars = paste0("lag1_", covs_tv),
          hist.fun = "lag",
          mc.sims = 10,
          mc.start = 0,
          mc.stop = 5,
          merge = FALSE,
          pool = TRUE
        )
      
      natural_course <-
        natural_course %>%
        rename(event_chd = pred,
               event_chd_alt = Py,
               event_dth = Pd) %>%
        group_by(id, time) %>%
        summarise(
          across(
            c(all_of(covs_tv), event_chd, event_chd_alt, event_dth),
            mean, 
            .names = "nc_{col}"
          ),
          .groups = "drop_last") %>%
        mutate(
          nc_event_dth = cumsum(nc_event_dth),
          nc_event_chd_alt = cumsum(nc_event_chd_alt)
        )

      observed <-
        boot_df_obs %>%
        select(id, time, all_of(covs_tv), event_chd, event_chd_zero, event_dth) %>%
        mutate(event_chd_alt = event_chd_zero)
      
      ev <- 
        left_join(natural_course, observed, by = c("id", "time")) %>%
        left_join(boot_df_rf, by = c("id", "time")) %>%
        pivot_longer(-c("id", "time")) %>%
        mutate(
          type = case_when(
            str_detect(name, "^nc_") ~ "natural course (parametric)", 
            str_detect(name, "^rf_") ~ "natural course (random forest)", 
            TRUE ~ "observed"
          ),
          name = str_remove(name, "^nc_"),
          name = str_remove(name, "^rf_")
        )
      
      eval_mv <- 
        ev %>%
        group_by(time, name, type) %>%
        summarise(
          mean_value = mean(value, na.rm = TRUE),
          .groups = "drop"
        ) 
      
      eval_md <- 
        ev %>%
        group_by(id, time, name) %>%
        summarise(
          nc_diff = first(value) - last(value), 
          rf_diff = nth(value, 2) - last(value), 
          .groups = "drop"
        ) %>%
        group_by(time, name) %>%
        summarise(
          mean_nc_diff = mean(nc_diff, na.rm = TRUE),
          mean_rf_diff = mean(rf_diff, na.rm = TRUE),
          .groups = "drop"
        )
      
      left_join(eval_mv, eval_md, by = c("name", "time"))
    }, pb = pb)

close(pb)


# Make plots --------------------------------------------------------------

# plot 1: mean covariate values
bind_rows(nc_boots) %>%
  group_by(time, name, type) %>%
  summarise(
    value = mean(mean_value, na.rm = TRUE),
    lwr = quantile(mean_value, 0.025),
    upr = quantile(mean_value, 0.975),
    .groups = "drop"
  ) %>%
  filter(!name %in% c(
    "event_chd",
    "event_chd_alt",
    "event_dth",
    "event_chd_zero",
    "pid"
  )) %>%
  mutate( 
    exam = time + 4,
    name = case_when(
      name == "age" ~ "Age",
      name == "bmi" ~ "BMI",
      name == "dm" ~ "Diabetes (%)",
      name == "hdl" ~ "HDL",
      name == "hrx" ~ "Anti-hypertensive drugs (%)",
      name == "liprx" ~ "Lipid-lowering drugs (%)",
      name == "sbp" ~ "SBP",
      name == "smk" ~ "Smoking (%)",
      name == "tc" ~ "Total cholesterol"
    )
  ) %>%
  ggplot(., aes(
    x = exam,
    y = value,
    group = type,
    color = type,
    fill = type,
    linetype = type
  )) +
  facet_wrap(~name, scales = "free_y") + 
  #geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3, color = NA) +
  geom_line() +
  #geom_point(shape = 16) +  
  labs(
    x = "Exam",
    y = ""
  ) +
  scale_color_manual(
    name = "",
    values = c(RColorBrewer::brewer.pal(3, name = "Set1")[1:2], "black")
    ) +
  scale_fill_manual(
    name = "",
    values = c(RColorBrewer::brewer.pal(3, name = "Set1")[1:2], "black")
  ) +
  scale_linetype_manual(
    name = "",
    values = c(1, 1, 2)
  ) +
  theme_bw() +
  theme(legend.position = "bottom")


# plot 2: outcome means
bind_rows(nc_boots) %>%
  group_by(time, name, type) %>%
  summarise(
    value = mean(mean_value, na.rm = TRUE),
    lwr = quantile(mean_value, 0.025),
    upr = quantile(mean_value, 0.975),
    .groups = "drop"
  ) %>%
  filter(name %in% c(
    "event_chd",
    "event_chd_alt",
    "event_dth"
  )) %>%
  mutate(
    exam = time + 4,
    name = case_when(
      name == "event_chd_alt"  & type == "observed" ~ "CHD (comp. risk)",
      name == "event_chd_alt" ~ "CHD",
      name == "event_chd"  & type == "observed" ~ "CHD",
      name == "event_chd" ~ "CHD (comp. risk)",
      name == "event_dth" ~ "non-CHD death"
    )
  ) %>%
  ggplot(., aes(
    x = exam,
    y = value,
    group = type,
    color = type,
    fill = type,
    linetype = type
  )) +
  facet_wrap(~name, scales = "free_y") + 
  #geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3, color = NA) +
  geom_line() +
  #geom_point(shape = 16) +  
  labs(
    x = "Exam",
    y = "Cumulative risk (%)"
  ) +
  scale_color_manual(
    name = "",
    values = c(RColorBrewer::brewer.pal(3, name = "Set1")[1:2], "black")
  ) +
  scale_fill_manual(
    name = "",
    values = c(RColorBrewer::brewer.pal(3, name = "Set1")[1:2], "black")
  ) +
  scale_linetype_manual(
    name = "",
    values = c(1, 1, 2)
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# plot 3: covariate differences
bind_rows(nc_boots) %>%
  filter(type == "natural course (parametric)") %>%
  group_by(time, name) %>%
  summarise(
    nc_value = mean(mean_nc_diff, na.rm = TRUE),
    nc_lwr = quantile(mean_nc_diff, 0.025),
    nc_upr = quantile(mean_nc_diff, 0.975),
    rf_value = mean(mean_rf_diff, na.rm = TRUE),
    rf_lwr = quantile(mean_rf_diff, 0.025),
    rf_upr = quantile(mean_rf_diff, 0.975),
    .groups = "drop"
  ) %>%
  filter(!name %in% c(
    "event_chd",
    "event_chd_alt",
    "event_dth",
    "event_chd_zero",
    "pid"
  )) %>%
  mutate(
    exam = time + 4,
    name = case_when(
      name == "age" ~ "Age",
      name == "bmi" ~ "BMI",
      name == "dm" ~ "Diabetes (%)",
      name == "hdl" ~ "HDL",
      name == "hrx" ~ "Anti-hypertensive drugs (%)",
      name == "liprx" ~ "Lipid-lowering drugs (%)",
      name == "sbp" ~ "SBP",
      name == "smk" ~ "Smoking (%)",
      name == "tc" ~ "Total cholesterol"
    )
  ) %>%
  pivot_longer(
    cols = -c("time", "name", "exam"),
    names_to = c("type", "variable"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  pivot_wider(
    names_from = "variable",
    values_from = "value"
  ) %>%
  mutate(
    type = case_when(
      type == "rf" ~ "natural course (random forest)",
      type == "nc" ~ "natural course (parametric)"
    )
  ) %>%
  ggplot(.,
         aes(
           x = exam,
           y = value,
           group = type,
           color = type,
           fill = type,
         )) +
  facet_wrap(~name, scales = "free_y") + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              alpha = 0.5,
              color = NA) +
  geom_line() +
  #geom_point() +
  scale_color_manual(
    name = "",
    values = RColorBrewer::brewer.pal(3, name = "Set1")[1:2]
  ) +
  scale_fill_manual(
    name = "",
    values = RColorBrewer::brewer.pal(3, name = "Set1")[1:2]
  ) +
  labs(
    x = "Exam",
    y = "Observed - Natural Course"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")