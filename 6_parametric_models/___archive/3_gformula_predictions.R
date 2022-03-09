# Run g-computation algorithm on full time course -------------------------

run_params <- list(
  name = list(
    "no lags",
    "one lag",
    "two lags"
  ),
  Y.model = list(
    Y.fit_lag0,
    Y.fit_lag1,
    Y.fit_lag2
  ),
  X.models = list(
    X.fit_lag1,
    X.fit_lag1,
    X.fit_lag2
  ),
  D.model = list(
    D.fit_lag0,
    D.fit_lag1,
    D.fit_lag2
  )
)

pb <- txtProgressBar(max = length(1:3), initial = NA, style = 3)

results <- 
  pmap(
  run_params,
  function(name, Y.model, X.models, D.model) {
    i <- getTxtProgressBar(pb)
    setTxtProgressBar(pb, ifelse(is.na(i), 1, i + 1))
    
    stop <- 5
    sims <- 50
    # using only baseline project forward 
    # with D model
    Y.hat_wD <-
      gformula_mc(
        Y.fit = Y.model,
        X.fit = X.models,
        D.fit = D.model,
        data = offspring_long,
        id = "pid",
        time = "time",
        base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
        hist.vars = c(paste0("lag1_", covs_tv), paste0("lag2_", covs_tv)),
        hist.fun = "lag",
        mc.sims = sims,
        mc.start = 0,
        mc.stop = stop,
        merge = FALSE
      )
    
    # without D model
    Y.hat_noD <-
      gformula_mc(
        Y.fit = Y.model,
        X.fit = X.models,
        data = offspring_long,
        id = "pid",
        time = "time",
        base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
        hist.vars = c(paste0("lag1_", covs_tv), paste0("lag2_", covs_tv)),
        hist.fun = "lag",
        mc.sims = sims,
        mc.start = 0,
        mc.stop = stop,
        merge = FALSE
      )
    
    stats_wD <-
      map_dfr(0:stop,  function (y) {
        df <- left_join(filter(observed_events, time == y),
                        filter(Y.hat_wD, time == y), 
                        by = c("pid", "time"))
        
        val.prob(df[["pred"]], df[["event_chd_zero"]], pl = FALSE)
      }, .id = "time")
    
    stats_noD <-
      map_dfr(0:stop,  function (y) {
        df <- left_join(filter(observed_events, time == y),
                        filter(Y.hat_noD, time == y), 
                        by = c("pid", "time"))
        
        val.prob(df[["pred"]], df[["event_chd"]], pl = FALSE)
      }, .id = "time")
    
    # iteratively conditioning on past and projecting forward
    Y.hats_cond_wD <- 
      map(
        set_names(0:stop, paste0("pred_", 0:stop)),
        function (x) {
          gformula_mc(
            Y.fit = Y.model,
            X.fit = X.models,
            D.fit = D.model,
            data = offspring_long,
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
            Y.fit = Y.model,
            X.fit = X.models,
            data = offspring_long,
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
          
          rms::val.prob(df[[paste0("pred", begin)]], df[["event_chd"]], pl = FALSE)
        }, .id = "cond_inner")
        df$start <- begin
        df$stop <- as.numeric(df$cond_inner) - 1
        return(df)
      }, .id = "cond_outer")
    
    stats <-
      bind_rows(stats_wD, stats_noD, stats_cond_wD, stats_cond_noD, .id = "id")
    
    stats$name <- name
    stats$comprisk <- ifelse(stats$id %in% c(1, 3), 1, 0)
    stats$type <- ifelse(stats$id %in% c(1, 2), 1, 2)
    
    return(stats)
  }
)

close(pb)

results <- bind_rows(results)

results %>%
  filter(type == 1) %>%
  select(time, name, comprisk, `C (ROC)`, Brier) %>%
  pivot_longer(cols = c("C (ROC)", "Brier"), 
               names_to = "var") %>%
ggplot(.,
       aes(
         x = as.numeric(time) * 4,
         y = value,
         group = name,
         color = name
       )) +
  facet_grid(var ~ factor(comprisk, labels = c("no competing risks", "competing risks")), scales = "free_y") +
  geom_point() +
  geom_line() +
  scale_color_brewer(name = "", type = "qual", palette = "Set1") +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

results %>%
  filter(type == 2) %>%
  select(start, stop, name, comprisk, `C (ROC)`, Brier) %>%
  mutate(stop = start + stop) %>%
  pivot_longer(cols = c("C (ROC)", "Brier"), 
               names_to = "var") %>%
  filter(var == "C (ROC)") %>%
  ggplot(.,
         aes(
           x = factor(start),
           y = factor(stop),
           fill = value
         )) +
  facet_grid(name ~ factor(comprisk, labels = c("no competing risks", "competing risks")), scales = "free_y") +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black") +
  scale_fill_distiller(name = "", type = "seq", palette = "Spectral") +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

results %>%
  filter(type == 2) %>%
  select(start, stop, name, comprisk, `C (ROC)`, Brier) %>%
  mutate(stop = start + stop) %>%
  pivot_longer(cols = c("C (ROC)", "Brier"), 
               names_to = "var") %>%
  filter(var == "Brier") %>%
  ggplot(.,
         aes(
           x = factor(start),
           y = factor(stop),
           fill = value
         )) +
  facet_grid(name ~ factor(comprisk, labels = c("no competing risks", "competing risks")), scales = "free_y") +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black") +
  scale_fill_distiller(name = "", type = "seq", palette = "Spectral") +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

# Y.hat_lag1 <-
#   gformula_mc(
#     Y.fit = Y.fit_lag1,
#     X.fit = X.fit_lag1,
#     D.fit = D.fit_lag1,
#     data = offspring_long,
#     id = "pid",
#     time = "time",
#     base.covs = covs_fixed[!covs_fixed %i n% c(covs_refs)],
#     hist.vars = paste0("lag1_", covs_tv),
#     hist.fun = "lag",
#     mc.sims = 100,
#     mc.start = 0,
#     mc.stop = 5,
#     merge = FALSE
#   )
# 
# Y.hat_lag2 <-
#   gformula_mc(
#     Y.fit = Y.fit_lag2,
#     X.fit = X.fit_lag2,
#     D.fit = D.fit_lag2,
#     data = offspring_long,
#     id = "pid",
#     time = "time",
#     base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
#     hist.vars = c(paste0("lag1_", covs_tv), paste0("lag2_", covs_tv)),
#     hist.fun = "lag",
#     mc.sims = 100,
#     mc.start = 0,
#     mc.stop = 5,
#     merge = FALSE
#   )
# 
# Y.hat_rf <- 
#   gformula_mc(
#     Y.fit = Y.rf,
#     X.fit = X.rf,
#     D.fit = D.rf,
#     data = offspring_long,
#     id = "pid",
#     time = "time",
#     base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
#     hist.vars = paste0("lag1_", covs_tv),
#     hist.fun = "lag",
#     mc.sims = 100,
#     mc.start = 0,
#     mc.stop = 5,
#     merge = FALSE
#   )
# 
# Y.hat_comprisk <-
#   gformula_mc(
#     Y.fit = Y.fit,
#     X.fit = X.fit,
#      data = offspring_long,
#     id = "pid",
#     time = "time",
#     base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
#     hist.vars = paste0("lag1_", covs_tv),
#     hist.fun = "lag",
#     mc.sims = 100,
#     mc.start = 0,
#     mc.stop = 5,
#     merge = FALSE
#   )
# 
# Y.hat_comprisk_rf <- 
#   gformula_mc(
#     Y.fit = Y.rf,
#     X.fit = X.rf,
#     data = offspring_long,
#     id = "pid",
#     time = "time",
#     base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
#     hist.vars = paste0("lag1_", covs_tv),
#     hist.fun = "lag",
#     mc.sims = 100,
#     mc.start = 0,
#     mc.stop = 5,
#     merge = FALSE
#   )
# 
# map(0:5,  function (y) {
#   df <- left_join(filter(observed_events, time == y),
#                   filter(Y.hat_lag2, time == y))
#   
#   rms::val.prob(df[["pred"]], df[["event_chd_zero"]], m = 250)
# })
# 
# map(0:5,  function (y) {
#   df <- left_join(filter(Y.hat_rf, time == y),
#                   filter(observed_events, time == y))
#   
#   rms::val.prob(df[["pred"]], df[["event_chd"]])
# })
# 
# map(0:5,  function (y) {
#   df <- left_join(filter(Y.hat_comprisk, time == y),
#                   filter(observed_events, time == y))
#   
#   rms::val.prob(df[["pred"]], df[["event_chd"]])
# })
# 
# map(0:5,  function (y) {
#   df <- left_join(filter(Y.hat_comprisk_rf, time == y),
#                   filter(observed_events, time == y))
#   
#   rms::val.prob(df[["pred"]], df[["event_chd"]])
# })
# 
# map(0:5,  function (y) {
#   df <- left_join(filter(Y.hat, time == y),
#                   filter(observed_events, time == y))
#   
#   rms::val.prob(df[["pred"]], df[["event_chd_zero"]])
# })
# 
# map(0:5,  function (y) {
#   df <- left_join(filter(Y.hat_rf, time == y),
#                   filter(observed_events, time == y))
#   
#   rms::val.prob(df[["pred"]], df[["event_chd_zero"]])
# })
# 
# map(0:5,  function (y) {
#   df <- left_join(filter(Y.hat_comprisk, time == y),
#                   filter(observed_events, time == y))
#   
#   rms::val.prob(df[["pred"]], df[["event_chd_zero"]])
# })
# 
# map(0:5,  function (y) {
#   df <- left_join(filter(Y.hat_comprisk_rf, time == y),
#                   filter(observed_events, time == y))
#   
#   rms::val.prob(df[["pred"]], df[["event_chd_zero"]])
# })
# 
# # Run g-computation iteratively conditioning on past ----------------------
# 
# stop <- 5
# pb <- txtProgressBar(max = length(0:stop), initial = NA, style = 3)
# 
# Y.hats_cond <- 
#   map(
#     set_names(0:stop, paste0("pred_", 0:stop)),
#     function (x) {
#       i <- getTxtProgressBar(pb)
#       setTxtProgressBar(pb, ifelse(is.na(i), 1, i + 1))
#       gformula_mc(
#         Y.fit = Y.fit_lag2,
#         X.fit = X.fit_lag2,
#         D.fit = D.fit_lag2,
#         data = offspring_long,
#         id = "pid",
#         time = "time",
#         base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
#         hist.vars = c(paste0("lag1_", covs_tv), paste0("lag2_", covs_tv)),
#         hist.fun = "lag",
#         mc.sims = 100,
#         mc.start = x,
#         mc.stop = stop,
#         merge = FALSE
#       )
#     })
# 
# close(pb)
# 
# Y.hats_cond <- reduce(Y.hats_cond, left_join, by = c("pid", "time"))
# names(Y.hats_cond) <- c("pid", "time", paste0("pred", 0:stop))
# map(0:stop,  function (y) {
#   df <- left_join(filter(observed_events, time == stop),
#                   filter(Y.hats_cond, time == stop))
#   
#   rms::val.prob(df[[paste0("pred", y)]], df[["event_chd"]], m = 250)
# })
# 
# 
# pb <- txtProgressBar(max = length(0:stop), initial = NA, style = 3)
# 
# Y.hats_cond_comp <- 
#   map(
#     set_names(0:stop, paste0("pred_", 0:stop)),
#     function (x) {
#       i <- getTxtProgressBar(pb)
#       setTxtProgressBar(pb, ifelse(is.na(i), 1, i + 1))
#       gformula_mc(
#         Y.fit = Y.fit,
#         X.fit = X.fit,
#         D.fit = D.fit,
#         data = offspring_long,
#         id = "pid",
#         time = "time",
#         base.covs = covs_fixed[!covs_fixed %in% covs_refs],
#         hist.vars = paste0("lag1_", covs_tv),
#         hist.fun = "lag",
#         mc.sims = 100,
#         mc.start = x,
#         mc.stop = stop,
#         merge = FALSE
#       )
#     })
# 
# close(pb)
# 
# Y.hats_cond_comp <- reduce(Y.hats_cond_comp, left_join, by = c("pid", "time"))
# names(Y.hats_cond_comp) <- c("pid", "time", paste0("pred", 0:stop))
# map(0:stop,  function (y) {
#   df <- left_join(filter(observed_events, time == stop),
#                   filter(Y.hats_cond_comp, time == stop))
#   
#   rms::val.prob(df[[paste0("pred", y)]], df[["event_chd_zero"]], m = 250)
# })
# 
# pb <- txtProgressBar(max = length(0:4), initial = NA, style = 3)
# 
# Y.hats_cond_rf <- 
#   map(
#     set_names(0:4, paste0("pred_", 0:4)),
#     function (x) {
#       i <- getTxtProgressBar(pb)
#       setTxtProgressBar(pb, ifelse(is.na(i), 1, i + 1))
#       gformula_mc(
#         Y.fit = Y.rf_lag1,
#         X.fit = X.rf_lag1,
#         D.fit = D.rf_lag1,
#         data = offspring_long,
#         id = "pid",
#         time = "time",
#         base.covs = covs_fixed[!covs_fixed %in% covs_refs],
#         hist.vars = paste0("lag1_", covs_tv),
#         hist.fun = "lag",
#         mc.sims = 100,
#         mc.start = x,
#         mc.stop = 5,
#         merge = FALSE
#       )
#     })
# 
# close(pb)
# 
# Y.hats_cond_rf <- reduce(Y.hats_cond_rf, left_join, by = c("pid", "time"))
# names(Y.hats_cond_rf) <- c("pid", "time", paste0("pred", 0:4))
# map(0:4,  function (y) {
#   df <- left_join(filter(observed_events, time == 5),
#                   filter(Y.hats_cond_rf, time == 5))
#   
#   rms::val.prob(df[[paste0("pred", y)]], df[["event_chd_zero"]], m = 500)
# })

# 
# # Run g-computation iteratively projecting x time steps  ------------------
# 
# pb <- txtProgressBar(max = length(0:3), initial = NA, style = 3)
# 
# Y.hats_proj <- 
#   map(
#     set_names(1:5, paste0("pred_", 0:4)),
#     function (x) {
#       i <- getTxtProgressBar(pb)
#       setTxtProgressBar(pb, ifelse(is.na(i), 1, i + 1))
#       gformula_mc(
#         Y.fit = Y.fit,
#         X.fit = X.fit,
#         D.fit = D.fit,
#         data = offspring_long,
#         id = "pid",
#         time = "time",
#         base.covs = covs_fixed[!covs_fixed %in% covs_refs],
#         hist.vars = paste0("lag1_", covs_tv),
#         hist.fun = "lag",
#         mc.sims = 500,
#         mc.start = 0,
#         mc.stop = x,
#         merge = FALSE
#       )
#     })
# 
# close(pb)
# 
# Y.hats_proj <- reduce(Y.hats_proj, right_join, by = c("pid", "time"))
# names(Y.hats_proj) <- c("pid", "time", paste0("pred", 1:5))
# map(1:5,  function (y) {
#   p.name <- paste0("pred", y)
#   y.name <- paste0("event_chd_", y + 1)
#   
#   df <- left_join(filter(Y.hats_proj, time == y),
#                   select(offspring_wide, "pid", y.name))
#   
#   rms::val.prob(df[[p.name]], df[[y.name]])
# })
# 
# 
# # Run g-computation algorithm ---------------------------------------------
# 
# pb <- txtProgressBar(max = length(0:3), initial = NA, style = 3)
# 
# Y.hats_ <- 
#   map_dfc(
#     set_names(0:3, paste0("pred_", 0:3, "_", 0:3 + 2)),
#     function (x) {
#       i <- getTxtProgressBar(pb)
#       setTxtProgressBar(pb, ifelse(is.na(i), 1, i + 1))
#       gformula_mc(
#         Y.fit = Y.fit,
#         X.fit = X.fit,
#         D.fit = D.fit,
#         data = offspring_long,
#         id = "pid",
#         time = "time",
#         base.covs = covs_fixed[!covs_fixed %in% covs_refs],
#         hist.vars = paste0("lag1_", covs_tv),
#         hist.fun = "lag",
#         mc.sims = 500,
#         mc.start = x,
#         mc.stop = x + 2
#       )
#   })
#   
# close(pb)
# 
# offspring_long <- cbind(offspring_long, Y.hats)
# 
# 
# # get estimates of in sample discrimination/calibration -------------------
# 
# validation_stats <-
#   map(list(
#     "0_2" = 0:2,
#     "1_3" = 1:3,
#     "2_4" = 2:4,
#     "3_5" = 3:5
#   ),
#   function(x) {
#     map(x,  function (y) {
#       p.name <-  paste0("pred_", min(x), "_", max(x))
#       y.name <- paste0("event_chd_", y + 1)
#       
#       df <- filter(offspring_long, time == y)
#       
#       df <- left_join(
#         select(df, "pid", p.name), 
#         select(offspring_wide, "pid", y.name)
#       )
#       
#       rms::val.prob(df[[p.name]], df[[y.name]])
#     })
#   })
#       


