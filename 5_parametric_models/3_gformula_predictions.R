# Run gcomputation algorithm ----------------------------------------------

Y.hats <- 
  map_dfc(
    set_names(0:3, paste0("pred_", 0:3, "_", 0:3 + 2)),
    function (x) {
      gformula_mc(
        Y.fit,
        X.fit,
        D.fit,
        data = analytic_long,
        id = "pid",
        time = "time",
        base.covs = covs_fixed[!covs_fixed %in% covs_refs],
        hist.vars = paste0("lag1_", covs_tv),
        hist.fun = "lag",
        mc.sims = 10,
        mc.start = x,
        mc.stop = x + 2
      )
  })
  
analytic_long <- cbind(analytic_long, Y.hats)


# # predict risk on hold-out sample -----------------------------------------

stats <-
  map(0:3,
      function (x) {
        df <- filter(analytic_long, between(time, x, x + 2))
        rms::val.prob(df[, paste0("pred_", x, "_", x + 2)], df$event_chd)
      })
# analytic_hat <-
#   left_join(analytic_long, Y.hat, by = c("pid", "time"))
# 
# stats <- map(0:5, function(t) {
#   df <- filter(analytic_hat, time == t)
#   rms::val.prob(df$pred, df$event_chd)
# })
# 
# aucs <- map(0:3, function(t) {
#   df <- filter(analytic_hat, time == t)
#   pROC::roc(df$event_chd, df$pred)
# })
# 
# map(aucs, pROC::ggroc)
# 
# 
# # plot some sample trajectories -------------------------------------------
# 
# samples <- newdata[1:6,]
# 
# trajectories <- predict.gformula(
#   object = fit_pred,
#   obs_data = data.table::as.data.table(drop_na(analytic_long)),
#   newdata = data.table::as.data.table(samples),
#   id = "pid",
#   t0 = 0,
#   covnames = covs_dvs, 
#   covtypes = covtypes,
#   covparams = covparams,
#   outcome_name = "event_chd",
#   ymodel = ymodel,
#   compevent_name = "event_dth",
#   compevent_model = compevent_model,
#   restrictions = restrictions,
#   basecovs = covs_fixed[!covs_fixed %in% covs_refs],
#   histvars = list(covs_dvs),
#   histories = c(lagged),
#   nsamples = 0,
#   nsimul = 1000,
#   seed = 4548076,
#   model_fits = TRUE,
#   return_sims = TRUE
# )
# 
# ggplot(trajectories$sims, aes(x = time, y = survival)) +
#   #facet_grid(~factor(id)) +
#   geom_line(aes(group = factor(id)), alpha = 0.04) +
#   stat_summary(geom = "line", fun = mean, color = "blue", size = 1.2) +
#   stat_summary(geom = "point", fun = mean, color = "blue") +
#   g_theme() +
#   labs(
#     x = "Follow up",
#     y = "Predicted survival"
#   ) +
#   coord_cartesian(expand = F)
