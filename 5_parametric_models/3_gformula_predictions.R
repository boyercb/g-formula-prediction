# Run g-computation algorithm ---------------------------------------------

pb <- txtProgressBar(max = length(0:3), initial = NA, style = 3)

Y.hats <- 
  map_dfc(
    set_names(0:3, paste0("pred_", 0:3, "_", 0:3 + 2)),
    function (x) {
      i <- getTxtProgressBar(pb)
      setTxtProgressBar(pb, ifelse(is.na(i), 1, i + 1))
      gformula_mc(
        Y.fit = Y.fit,
        X.fit = X.fit,
        D.fit = D.fit,
        data = analytic_long,
        id = "pid",
        time = "time",
        base.covs = covs_fixed[!covs_fixed %in% covs_refs],
        hist.vars = paste0("lag1_", covs_tv),
        hist.fun = "lag",
        mc.sims = 500,
        mc.start = x,
        mc.stop = x + 2
      )
  })
  
analytic_long <- cbind(analytic_long, Y.hats)


# get estimates of in sample discrimination/calibration -------------------

validation_stats <-
  map(list(
    "0_2" = 0:2,
    "1_3" = 1:3,
    "2_4" = 2:4,
    "3_5" = 3:5
  ),
  function(x) {
    map(x,  function (y) {
      p.name <-  paste0("pred_", min(x), "_", max(x))
      y.name <- paste0("event_chd_", y + 1)
      
      df <- filter(analytic_long, time == y)
      
      df <- left_join(
        select(df, "pid", p.name), 
        select(analytic_wide, "pid", y.name)
      )
      
      rms::val.prob(df[[p.name]], df[[y.name]])
    })
  })
      


