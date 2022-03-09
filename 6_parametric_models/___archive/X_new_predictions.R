# Run g-computation algorithm on full time course -------------------------

run_params <- list(
  name = list(
    # "no lags",
    # "one lag",
    "two lags"
  ),
  models = list(parametric_models_uncens[["lag2"]]),
  comparison = list(conventional_models_uncens),
  data_long = list(filter(offspring_long, time < 5 & !pid %in% censored_ids)),
  data_wide = list(offspring_wide)
)


pb <- txtProgressBar(max = length(1:3), initial = NA, style = 3)

results <- 
  pmap(
    run_params,
    function(name, models, comparison, data_long, data_wide) {
      i <- getTxtProgressBar(pb)
      setTxtProgressBar(pb, ifelse(is.na(i), 1, i + 1))
      
      stop <- 4
      sims <- 100
      
      # iteratively conditioning on past and projecting forward
      Y.hats_cond_wD <- 
        map(
          set_names(0:stop, paste0("pred_", 0:stop)),
          function (x) {
            gformula_mc(
              Y.fit = models[["Y"]],
              X.fit = models[["X"]],
              D.fit = models[["D"]],
              data = data_long,
              id = "pid",
              time = "time",
              base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
              hist.vars = c(paste0("lag1_", covs_tv), paste0("lag2_", covs_tv)),
              hist.fun = "lag",
              mc.sims = sims,
              mc.start = x,
              mc.stop = stop,
              merge = FALSE
            )
          })
      
      Y.hats_cond_noD <- 
        map(
          set_names(0:stop, paste0("pred_", 0:stop)),
          function (x) {
            gformula_mc(
              Y.fit = models[["Y"]],
              X.fit = models[["X"]],
              data = data_long,
              id = "pid",
              time = "time",
              base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
              hist.vars = c(paste0("lag1_", covs_tv), paste0("lag2_", covs_tv)),
              hist.fun = "lag",
              mc.sims = sims,
              mc.start = x,
              mc.stop = stop,
              merge = FALSE
            )
          })
      
      Y.hats_cond_wD <- reduce(Y.hats_cond_wD, left_join, by = c("pid", "time"))
      names(Y.hats_cond_wD) <- c("pid", "time", paste0("pred", 0:stop))
      stats_cond_wD <-
        map_dfr(0:stop,  function (begin) {
          df <- map_dfr(begin:stop,  function (end) {
            df <- left_join(
              filter(observed_events, time == end),
              filter(Y.hats_cond_wD, time == end),
              by = c("pid", "time")
            )
            
            rms::val.prob(df[[paste0("pred", begin)]], df[["event_chd_zero"]], pl = FALSE)
          }, .id = "cond_inner")
          df$start <- begin
          df$stop <- as.numeric(df$cond_inner) - 1
          return(df)
        }, .id = "cond_outer")
      
      Y.hats_cond_noD <- reduce(Y.hats_cond_noD, left_join, by = c("pid", "time"))
      names(Y.hats_cond_noD) <- c("pid", "time", paste0("pred", 0:stop))
      stats_cond_noD <-
        map_dfr(0:stop,  function (begin) {
          df <- map_dfr(begin:stop,  function (end) {
            df <- left_join(
              filter(observed_events, time == end),
              filter(Y.hats_cond_noD, time == end),
              by = c("pid", "time")
            )
            
            rms::val.prob(df[[paste0("pred", begin)]], df[["event_chd_zero"]], pl = FALSE)
          }, .id = "cond_inner")
          df$start <- begin
          df$stop <- as.numeric(df$cond_inner) - 1
          return(df)
        }, .id = "cond_outer")
      
      stats <-
        bind_rows(stats_cond_wD, stats_cond_noD, .id = "id")
      
      stats$name <- name
      stats$comprisk <- ifelse(stats$id %in% c(1), 1, 0)

      return(stats)
    }
  )

close(pb)

results <- bind_rows(results)

test <- filter(offspring_wide, !pid %in% censored_ids)
test$pred <- predict(conventional_models_uncens[["conventional_12yr"]], newdata = test, type = "response")
test <- left_join(select(test, pid, pred), observed_events_wide)
val.prob(test$pred, test$event_chd_3)

# results %>%
#   filter(type == 1) %>%
#   select(time, name, comprisk, `C (ROC)`, Brier) %>%
#   pivot_longer(cols = c("C (ROC)", "Brier"), 
#                names_to = "var") %>%
#   ggplot(.,
#          aes(
#            x = as.numeric(time) * 4,
#            y = value,
#            group = name,
#            color = name
#          )) +
#   facet_grid(var ~ factor(comprisk, labels = c("no competing risks", "competing risks")), scales = "free_y") +
#   geom_point() +
#   geom_line() +
#   scale_color_brewer(name = "", type = "qual", palette = "Set1") +
#   theme_bw() +
#   labs(
#     x = "conditioning on exam",
#     y = "projecting to exam",
#     title = "C-statistic"
#   )
#   theme(
#     legend.position = "bottom"
#   )

results %>%
  select(start, stop, name, comprisk, `C (ROC)`, Brier) %>%
  mutate(stop = start + stop) %>%
  pivot_longer(cols = c("C (ROC)", "Brier"), 
               names_to = "var") %>%
  filter(var == "C (ROC)") %>%
  ggplot(.,
         aes(
           x = factor(start + 4),
           y = factor(stop + 4),
           fill = value
         )) +
  facet_grid(name ~ factor(comprisk, labels = c("no competing risks", "competing risks")), scales = "free_y") +
  geom_tile() +
  geom_text(aes(label = round(value, 3)), color = "black") +
  scale_fill_distiller(name = "", type = "seq", palette = "Spectral") +
  theme_bw() +
  labs(
    x = "conditioning on exam",
    y = "projecting to exam",
    title = "C-statistic"
  ) +
  theme(
    legend.position = "bottom"
  )

results %>%
  select(start, stop, name, comprisk, `C (ROC)`, Brier) %>%
  mutate(stop = start + stop) %>%
  pivot_longer(cols = c("C (ROC)", "Brier"), 
               names_to = "var") %>%
  filter(var == "Brier") %>%
  ggplot(.,
         aes(
           x = factor(start + 4),
           y = factor(stop + 4),
           fill = value
         )) +
  facet_grid(name ~ factor(comprisk, labels = c("no competing risks", "competing risks")), scales = "free_y") +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black") +
  scale_fill_distiller(name = "", type = "seq", palette = "Spectral") +
  theme_bw() +
  labs(
    x = "conditioning on exam",
    y = "projecting to exam",
    title = "Brier score"
  ) +
  theme(
    legend.position = "bottom"
  )

results %>%
  select(start, stop, name, comprisk, `Slope`, `Intercept`) %>%
  mutate(stop = start + stop) %>%
  pivot_longer(cols = c("Slope", "Intercept"), 
               names_to = "var") %>%
  filter(var == "Slope") %>%
  ggplot(.,
         aes(
           x = factor(start + 4),
           y = factor(stop + 4),
           fill = value
         )) +
  facet_grid(name ~ factor(comprisk, labels = c("no competing risks", "competing risks")), scales = "free_y") +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black") +
  scale_fill_distiller(name = "", type = "seq", palette = "Spectral") +
  theme_bw() +
  labs(
    x = "conditioning on exam",
    y = "projecting to exam",
    title = "Calibration slope"
  ) +
  theme(
    legend.position = "bottom"
  )

results %>%
  select(start, stop, name, comprisk, `Slope`, `Intercept`) %>%
  mutate(stop = start + stop) %>%
  pivot_longer(cols = c("Slope", "Intercept"), 
               names_to = "var") %>%
  filter(var == "Intercept") %>%
  ggplot(.,
         aes(
           x = factor(start + 4),
           y = factor(stop + 4),
           fill = value
         )) +
  facet_grid(name ~ factor(comprisk, labels = c("no competing risks", "competing risks")), scales = "free_y") +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black") +
  scale_fill_distiller(name = "", type = "seq", palette = "Spectral") +
  theme_bw() +
  labs(
    x = "conditioning on exam",
    y = "projecting to exam",
    title = "Calibration intercept"
  ) +
  theme(
    legend.position = "bottom"
  )

results %>%
  select(start, stop, name, comprisk, `U:p`) %>%
  mutate(stop = start + stop) %>%
  pivot_longer(cols = c("U:p"), 
               names_to = "var") %>%
  ggplot(.,
         aes(
           x = factor(start + 4),
           y = factor(stop + 4),
           fill = value
         )) +
  facet_grid(name ~ factor(comprisk, labels = c("no competing risks", "competing risks")), scales = "free_y") +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black") +
  scale_fill_distiller(name = "", type = "seq", palette = "Spectral") +
  theme_bw() +
  labs(
    x = "conditioning on exam",
    y = "projecting to exam",
    title = "Calibration intercept"
  ) +
  theme(
    legend.position = "bottom"
  )



