# fit random forest model -------------------------------------------------

observed_events_covs <-
  left_join(mutate(observed_events, time = as.numeric(time)),
            select(offspring_long, all_of(covs_tv), pid, time),
            by = c("pid", "time"))


# Bootstrap differences in natural course ---------------------------------

BOOT_SMPS <- 500

plan(multisession, workers = 14)

with_progress({
  p <- progressor(steps = BOOT_SMPS)
  
  nc_boots <- 
    future_map_dfr(
      .x = 1:BOOT_SMPS, 
      .f = function(x) {
        
        # draw bootstrapped samples
        boot <- draw_bootstrap_samples(filter(offspring_long, !pid %in% censored_ids), "pid", "time")
        
        # match to observed data
        boot_obs <- left_join(select(boot, "pid", "time", "id"), observed_events_covs, by = c("pid", "time"))
        
        # fit models on bootstrap sample
        Y.fit_boot <- 
          fit_Y_model(
            parametric_model_definitions[["lag2"]][["Y"]],
            data = filter(boot, time < 6)
          )
        
        X.fit_boot <- 
          fit_X_models(
            parametric_model_definitions[["lag2"]][["X"]],
            data = filter(boot, time < 6),
            time = "time"
          )
        
        D.fit_boot <- 
          fit_D_model(
            parametric_model_definitions[["lag2"]][["D"]],
            data = filter(boot, time < 6)
          )
        
        # estimate risk under the natural course
        natural_course <- 
          gformula_mc(
            Y.fit = Y.fit_boot,
            X.fit = X.fit_boot,
            D.fit = D.fit_boot,
            data = boot,
            id = "id",
            time = "time",
            base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
            hist.vars = c(paste0("lag1_", covs_tv), paste0("lag2_", covs_tv)),
            hist.fun = "lag",
            mc.sims = 50,
            mc.start = 0,
            mc.stop = 5,
            merge = FALSE,
            pool = TRUE
          )
        
        # average predictions
        natural_course <-
          natural_course %>%
          rename(event_ascvd = pred,
                 event_ascvd_alt = Py,
                 event_dth_ascvd = Pd) %>%
          group_by(id, time) %>%
          summarise(
            across(
              c(all_of(covs_tv), event_ascvd, event_ascvd_alt, event_dth_ascvd),
              mean
            ),
            .groups = "drop_last") %>%
          mutate(
            event_dth_ascvd = cumsum(event_dth_ascvd),
            event_ascvd_alt = cumsum(event_ascvd_alt)
          ) %>%
          group_by(time) %>%
          summarise(
            across(
              c(all_of(covs_tv), event_ascvd, event_ascvd_alt, event_dth_ascvd),
              mean, 
              .names = "nc_{col}"
            ),
            .groups = "drop_last")
  
        observed <-
          boot_obs %>%
          select(id, time, all_of(covs_tv), event_ascvd, event_ascvd_zero, event_dth_ascvd) %>%
          mutate(event_ascvd_alt = event_ascvd_zero) %>% 
          group_by(time) %>%
          summarise(
            across(
              c(all_of(covs_tv), event_ascvd, event_ascvd_alt, event_dth_ascvd),
              mean, 
              .names = "obs_{col}",
              na.rm = TRUE
            ),
            .groups = "drop"
          ) %>%
          mutate(
            obs_event_ascvd = 1 - cumprod(1 - obs_event_ascvd),
            obs_event_ascvd_alt = 1 - cumprod(1 - obs_event_ascvd_alt),
            obs_event_dth_ascvd = 1 - cumprod(1 - obs_event_dth_ascvd),
          )
        
        ev <- 
          left_join(natural_course, observed, by = c("time")) %>%
          pivot_longer(-c("time"), values_to = "mean_value") %>%
          mutate(
            type = case_when(
              str_detect(name, "^nc_") ~ "natural course (parametric)", 
              TRUE ~ "observed (nonparametric)"
            ),
            name = str_remove(name, "^(nc_)|(obs_)")
          )
        
        
        eval_md <- 
          ev %>%
          group_by(time, name) %>%
          summarise(
            nc_diff = first(mean_value) - last(mean_value), 
            .groups = "drop"
          ) 
        
        p()
        
        return(left_join(ev, eval_md, by = c("name", "time")))
      },     
      .id = "bsim",
      .options = furrr_options(
        seed = TRUE
      ))
})

plan(sequential)

# Make plots --------------------------------------------------------------

# plot 1: mean covariate values
pdf("9_results/figures/diag_covs.pdf", width = 8, height = 5)
nc_boots %>%
  group_by(time, name, type) %>%
  summarise(
    value = mean(mean_value, na.rm = TRUE),
    lwr = quantile(mean_value, 0.025),
    upr = quantile(mean_value, 0.975),
    .groups = "drop"
  ) %>%
  filter(!name %in% c(
    "event_ascvd",
    "event_ascvd_alt",
    "event_dth_ascvd",
    "event_ascvd_zero",
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
dev.off()


# plot 2: outcome means
pdf("9_results/figures/diag_outcomes.pdf", width = 8, height = 5)
nc_boots %>%
  group_by(time, name, type) %>%
  summarise(
    value = mean(mean_value, na.rm = TRUE),
    lwr = quantile(mean_value, 0.025),
    upr = quantile(mean_value, 0.975),
    .groups = "drop"
  ) %>%
  filter(name %in% c(
    "event_ascvd",
    "event_ascvd_alt",
    "event_dth_ascvd"
  )) %>%
  mutate(
    exam = time + 4,
    name = case_when(
      name == "event_ascvd_alt"  & type == "observed" ~ "CHD (comp. risk)",
      name == "event_ascvd_alt" ~ "CHD",
      name == "event_ascvd"  & type == "observed" ~ "CHD",
      name == "event_ascvd" ~ "CHD (comp. risk)",
      name == "event_dth_ascvd" ~ "non-CHD death"
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
dev.off()


# plot 3: covariate differences
pdf("9_results/figures/diag_covs_diff.pdf", width = 8, height = 5)
nc_boots %>%
  filter(type == "natural course (parametric)") %>%
  group_by(time, name) %>%
  summarise(
    nc_value = mean(nc_diff, na.rm = TRUE),
    nc_lwr = quantile(nc_diff, 0.025),
    nc_upr = quantile(nc_diff, 0.975),
    .groups = "drop"
  ) %>%
  filter(!name %in% c(
    "event_ascvd",
    "event_ascvd_alt",
    "event_dth_ascvd",
    "event_ascvd_zero",
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
    y = "Natural Course - Observed"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()
