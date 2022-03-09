# Use parametric g-formula to predict risks -------------------------------



# option 1  ---------------------------------------------------------------

# start from exam X, iteratively update using data from subsequent exams until
# you reach the final exam, the number of exams each prediction is conditional
# on is determined by the lag structure of the parametric models, e.g.
#
# parametric models w/ 2 lags:
#   iter 0 - conditional on exam 0 alone project 0 through 5
#   iter 1 - conditional on exams 0 to 1 project 1 through 5
#   iter 2 - conditional on exams 0 to 2 project 2 through 5
#   iter 3 - conditional on exams 1 to 3 project 3 through 5
#   iter 4 - conditional on exams 2 to 4 project 4 through 5
#   iter 5 - conditional on exams 3 to 5 project 5

BOOTS <- 500

run_params <- list(
  name = list(
    "no lags",
    "one lag",
    "two lags"
  ),
  def = parametric_model_definitions,
  data = list(filter(offspring_long, time < 5 & !pid %in% censored_ids)),
  stop = 4,
  sims = 100
)

plan(multisession, workers = 14)

with_progress({
  p <- progressor(steps = length(0:run_params$stop) * length(run_params$name) * BOOTS)
  
  boot_stats_opt1 <- future_map_dfr(
    1:BOOTS,
    function(bsim, run_params) {
      validation_preds_opt1 <- 
        pmap(
          run_params,
          function(name, def, data, stop, sims) {
            boot <- draw_bootstrap_samples(data, "pid", "time")
            oob <- data[!data[['pid']] %in% boot[['pid']], ]
            
            models <- list(
              "X" = fit_X_models(
                models = def[["X"]],
                data = boot,
                time = "time"
              ),
              "Y" = fit_Y_model(
                model = def[["Y"]],
                data = boot
              ),
              "D" = NULL
            )
            
            updated_preds <- 
              map(
                set_names(0:stop, paste0("pred_", 0:stop)),
                function (x) {
                  pred <- gformula_mc(
                    Y.fit = models[["Y"]],
                    X.fit = models[["X"]],
                    D.fit = models[["D"]],
                    data = boot,
                    id = "id",
                    time = "time",
                    base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
                    hist.vars = c(paste0("lag1_", covs_tv), paste0("lag2_", covs_tv)),
                    hist.fun = "lag",
                    mc.sims = sims,
                    mc.start = x,
                    mc.stop = stop,
                    merge = FALSE
                  )
                  p()
                  return(pred)
                }
              )
            
            updated_preds <- reduce(updated_preds, left_join, by = c("id", "time"))
            names(updated_preds) <- c("id", "time", paste0("pred_", 0:stop))
            
            updated_preds <- left_join(
              updated_preds,
              select(filter(boot, time == 0), "id", "pid"),
              by = c("id")
            )
            
            updated_preds <- left_join(
              updated_preds,
              observed_events,
              by = c("pid", "time")
            )
            
            
            return(updated_preds)
          }
        )
      
      names(validation_preds_opt1) <- run_params$name
      
      stats <- map_dfr(validation_preds_opt1,
                   function(x, start = 0, stop = 4) {
                     
                     grid <- expand_grid(
                       time = start:stop,
                       pred = start:stop
                     )
                     
                     grid <- filter(grid, time >= pred)
                     
                     pmap_dfr(
                       as.list(grid),
                       function(time, pred) {
                         df <- val.prob(
                           p = x[x[["time"]] == time, ][[paste0("pred_", pred)]], 
                           y = x[x[["time"]] == time, ][["event_ascvd"]],
                           pl = FALSE
                         )
                         
                         df[['start']] <- pred
                         df[['stop']] <- time
                         
                         return(df)
                       }
                     )
                   }, .id = "model")
      
      return(stats)
    },
    run_params = run_params,
    .id = "bsim",
    .options = furrr_options(
      seed = TRUE
    )
  )
  
})

plan(sequential)




# option 2 ----------------------------------------------------------------

# always start from same exam (e.g. exam 2) and iteratively add previous exam
# information to obtain updated predictions -- here the advantage is that we are
# comparing predictions over the same time horizon and observed cases/deaths,
# etc.

# parametric models w/ 2 lags:
#   iter 0 - conditional on exam 2 alone project 2 through 5
#   iter 1 - conditional on exams 1 to 2 project 2 through 5

run_params <- list(
  name = list(
    "no lags",
    "one lag",
    "two lags"
  ),
  def = parametric_model_definitions,
  data = list(filter(offspring_long, time < 5 & !pid %in% censored_ids)),
  start = 2,
  stop = 4,
  sims = 100
)

plan(multisession, workers = 14)

with_progress({
  p <- progressor(steps = length(0:2) * length(run_params$name) * BOOTS)
  
  boot_stats_opt2 <- future_map_dfr(
    1:BOOTS,
    function(bsim, run_params) {
      validation_preds_opt2 <- 
        pmap(
          run_params,
          function(name, def, data, start, stop, sims) {
            boot <- draw_bootstrap_samples(data, "pid", "time")
            oob <- data[!data[['pid']] %in% boot[['pid']], ]
            
            models <- list(
              "X" = fit_X_models(
                models = def[["X"]],
                data = filter(boot, time > 1 & time < 5),
                time = "time",
                t0 = 2
              ),
              "Y" = fit_Y_model(
                model = def[["Y"]],
                data = filter(boot, time > 1 & time < 5)
              ),
              "D" = NULL
            )
            
            updated_preds <- 
              map(
                set_names(0:2, paste0("pred_", 0:2)),
                function (x) {
                  pred <- gformula_mc(
                    Y.fit = models[["Y"]],
                    X.fit = models[["X"]],
                    D.fit = models[["D"]],
                    data = zero_out_lags(filter(boot, time == 2), x),
                    id = "id",
                    time = "time",
                    base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
                    hist.vars = c(paste0("lag1_", covs_tv), paste0("lag2_", covs_tv)),
                    hist.fun = "lag",
                    mc.sims = sims,
                    mc.start = start,
                    mc.stop = stop,
                    merge = FALSE
                  )
                  p()
                  return(pred)
                }
              )
            
            updated_preds <- reduce(updated_preds, left_join, by = c("id", "time"))
            names(updated_preds) <- c("id", "time", paste0("pred_", 0:2))
            
            updated_preds <- left_join(
              updated_preds,
              select(filter(boot, time == 2), "id", "pid"),
              by = c("id")
            )
            
            updated_preds <- left_join(
              updated_preds,
              observed_events,
              by = c("pid", "time")
            )
            # these seem a bit nutty!
            #print(updated_preds)
                        
            return(updated_preds)
          }
        )
      
      names(validation_preds_opt2) <- run_params$name
      
      stats <- map_dfr(validation_preds_opt2,
                   function(x, start = 2, stop = 4, lags = 0:2) {
                     
                     grid <- expand_grid(
                       time = start:stop,
                       pred = lags
                     )
                     
                     pmap_dfr(
                       as.list(grid),
                       function(time, pred) {
                         df <- val.prob(
                           p = x[x[["time"]] == time, ][[paste0("pred_", pred)]], 
                           y = x[x[["time"]] == time, ][["event_ascvd"]],
                           pl = FALSE
                         )
                         
                         df[['start']] <- pred
                         df[['stop']] <- time
                         
                         return(df)
                       }
                     )
                   }, .id = "model")
      
      return(stats)
    },
    run_params = run_params,
    .id = "bsim",
    .options = furrr_options(
      seed = TRUE
    )
  )
  
})

plan(sequential)

boot_stats_opt1_pct <- 
    boot_stats_opt1 %>%
    group_by(model, start, stop) %>%
    summarise(
      across(
        .cols = -bsim,
        .fns = list(
          mean = ~mean(.x, na.rm = TRUE),
          lwr = ~quantile(.x, 0.025, na.rm = TRUE),
          upr = ~quantile(.x, 0.975, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )
    )

val_stats_opt1 <-
  left_join(stats_opt1,
            boot_stats_opt1_pct)

val_stats_opt1 %>%
  mutate(
    stat = print_ci(`C (ROC)`, `C (ROC)_lwr`, `C (ROC)_upr`, 3)
  ) %>%
  select(stat, model, start, stop) %>%
  pivot_wider(
    names_from = stop,
    values_from = stat
  ) %>%
  mutate(
    start = paste0("hspace{2mm}Exam ", start + 4)
  ) %>%
  select(-model) %>%
  kable(format = "latex", booktabs = TRUE, linesep = "")

val_stats_opt1 %>%
  mutate(
    stat = print_ci(`Slope`, `Slope_lwr`, `Slope_upr`, 3)
  ) %>%
  select(stat, model, start, stop) %>%
  pivot_wider(
    names_from = stop,
    values_from = stat
  ) %>%
  mutate(
    start = paste0("hspace{2mm}Exam ", start + 4)
  ) %>%
  select(-model) %>%
  kable(format = "latex", booktabs = TRUE, linesep = "")


val_stats_opt1 %>%
  mutate(
    stat = print_ci(`Intercept`, `Intercept_lwr`, `Intercept_upr`, 3)
  ) %>%
  select(stat, model, start, stop) %>%
  pivot_wider(
    names_from = stop,
    values_from = stat
  ) %>%
  mutate(
    start = paste0("hspace{2mm}Exam ", start + 4)
  ) %>%
  select(-model) %>%
  kable(format = "latex", booktabs = TRUE, linesep = "")



boot_stats_opt2_pct <- 
  boot_stats_opt2 %>%
  group_by(model, start, stop) %>%
  summarise(
    across(
      .cols = -bsim,
      .fns = list(
        lwr = ~quantile(.x, 0.025, na.rm = TRUE),
        upr = ~quantile(.x, 0.975, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  )


val_stats_opt2 <-
  left_join(stats_opt2,
            boot_stats_opt2_pct)

val_stats_opt2 %>%
  mutate(
    stat = print_ci(`C (ROC)`, `C (ROC)_lwr`, `C (ROC)_upr`, 3)
  ) %>%
  select(stat, model, start, stop) %>%
  pivot_wider(
    names_from = stop,
    values_from = stat
  ) %>%
  mutate(
    start = start + 1
  ) %>%
  kable(format = "latex", booktabs = TRUE, linesep = "")

val_stats_opt2 %>%
  mutate(
    stat = print_ci(Slope, `Slope_lwr`, `Slope_upr`, 3)
  ) %>%
  select(stat, model, start, stop) %>%
  pivot_wider(
    names_from = stop,
    values_from = stat
  ) %>%
  mutate(
    start = start + 1
  ) %>%
  kable(format = "latex", booktabs = TRUE, linesep = "")

val_stats_opt2 %>%
  mutate(
    stat = print_ci(Intercept, `Intercept_lwr`, `Intercept_upr`, 3)
  ) %>%
  select(stat, model, start, stop) %>%
  pivot_wider(
    names_from = stop,
    values_from = stat
  ) %>%
  mutate(
    start = start + 1
  ) %>%
  kable(format = "latex", booktabs = TRUE, linesep = "")

# TODO: net re-classification, etc. 
# TODO: re-calibrated ASCVD pooled cohort risks
# TODO: finish simulations 
# TODO: look at ARIC data


# what do we need input on?
# 1. the specification of the outcome and covariate models
# 2. how to handle time (duration)
# 3. should we stick with exam-based time intervals or try to incorporate calendar or study time


# what do we normally do? we take just one exam (often most recent) as baseline and then build a model for 10-year cumulative incidence
# but we're discarding all this additional information about each subjects full history
# an easy way to solve this might be to just build conditional models but then you have to have one for each set of possible exams (and what do you do if people skip etc)
# however, this doesn't make the natural course explicit
# segway to simple example of why natural course may be relevant 

# what happens if natural course changes (to both g-formula and )


left_join(
  offspring_long %>% filter(time == 2) %>% select(pid, ascvd_10yr_accaha_risk),
  observed_events_wide
) %>%
  select(pid, ascvd_10yr_accaha_risk, ascvd_10yr_accaha_risk_2)

val.prob(observed_events_wide$ascvd_10yr_accaha_risk_2, observed_events_wide$event_ascvd_4)
val.prob(observed_events_wide$ascvd_10yr_frs_risk_2, observed_events_wide$event_ascvd_4)


df <- 
  gformula_predictions_opt2$`two lags` %>% 
  select(pid, time, pred_2) %>%
  filter(time == 4) 

comb <- left_join(
  observed_events_wide,
  df
) %>%
  select(pid, event_ascvd_4, ascvd_10yr_accaha_risk_2, pred_2) %>%
  drop_na()

nribin(
  event = comb$event_ascvd_4,
  p.std = comb$ascvd_10yr_accaha_risk_2,
  p.new = comb$pred_2,
  cut = 0.075
)

Hmisc::improveProb(x1 = comb$ascvd_10yr_accaha_risk_2, x2 = comb$pred_2, y = comb$event_ascvd_4)

