
stats_ex1 <- map2_dfr(gformula_predictions_ex1,
                      names(gformula_predictions_ex1),
                      function(x, name, time = 2) {
                        outcome <-
                          ifelse(str_detect(name, "ASCVD"),
                                 "event_ascvd_zero",
                                 "event_ascvd")
                        
                        val.prob(p = x[x[["time"]] == time, ][["pred"]],
                                 y = x[x[["time"]] == time, ][[outcome]],
                                 pl = FALSE)
                      }, .id = "model")


stats_ex2 <- map2_dfr(gformula_predictions_ex2,
                      names(gformula_predictions_ex2),
                      function(x, name, time = 2) {
                        outcome <-
                          ifelse(str_detect(name, "ASCVD"),
                                 "event_ascvd_zero",
                                 "event_ascvd")
                        
                        val.prob(p = x[x[["time"]] == time,][["pred"]],
                                 y = x[x[["time"]] == time,][[outcome]],
                                 pl = FALSE)
                      }, .id = "model")

stats_ex3 <- map2_dfr(gformula_predictions_ex3,
                      names(gformula_predictions_ex3),
                      function(x, name, time = 2) {
                        outcome <-
                          ifelse(str_detect(name, "ASCVD"),
                                 "event_ascvd_zero",
                                 "event_ascvd")
                        
                        val.prob(p = x[x[["time"]] == time,][["pred"]],
                                 y = x[x[["time"]] == time,][[outcome]],
                                 pl = FALSE)
                      }, .id = "model")


df <- 
  list(
    "PCE, net risk" = gformula_predictions_ex2$`bmi, net risk`,
    "PCE, ASCVD-specific risk" = gformula_predictions_ex2$`bmi, ASCVD-specific risk`,
    "PCE (calibrated), net risk" = gformula_predictions_ex2$`bmi, net risk`,
    "PCE (calibrated), ASCVD-specific risk" = gformula_predictions_ex2$`bmi, ASCVD-specific risk`,
    "FRS, net risk" = gformula_predictions_ex2$`bmi, net risk`,
    "FRS, ASCVD-specific risk" = gformula_predictions_ex2$`bmi, ASCVD-specific risk`,
    "FRS (calibrated), net risk" = gformula_predictions_ex2$`bmi, net risk`,
    "FRS (calibrated), ASCVD-specific risk" = gformula_predictions_ex2$`bmi, ASCVD-specific risk`
  )

stats_accaha <- map2_df(df, names(df),
                        function(x, name, time = 2) {
                          outcome <-
                            ifelse(str_detect(name, "ASCVD"),
                                   "event_ascvd_zero",
                                   "event_ascvd")
                          variable <- 
                            ifelse(str_detect(name, "PCE"),
                                   "ascvd_10yr_accaha_risk",
                                   "ascvd_10yr_frs_risk")
                          
                          calibrated <- 
                            ifelse(str_detect(name, "calibrated"),
                                   "",
                                   "_updated")
                          
                          variable <- paste0(variable, calibrated)
                          
                          x <- x %>%
                            group_by(pid) %>%
                            mutate(across(all_of(variable), ~ lag(.x, 2), .names = "{.col}_lag2")) %>%
                            ungroup()
                          
                          val.prob(p = x[x[["time"]] == time, ][[paste0(variable, "_lag2")]],
                                   y = x[x[["time"]] == time, ][[outcome]],
                                   pl = FALSE)
                        }, .id = "model")


stats_conventional <- map2_dfr(conventional_predictions,
                      names(conventional_predictions),
                      function(x, name, time = 2) {
                        outcome <-
                          ifelse(str_detect(name, "ASCVD"),
                                 "event_ascvd_zero_2",
                                 "event_ascvd_2")
                        
                        val.prob(p = x[["pred"]],
                                 y = x[[outcome]],
                                 pl = FALSE)
                      }, .id = "model")


# df <- gformula_predictions_ex2$bmi %>%
#   group_by(pid) %>%
#   mutate(
#     ascvd_10yr_accaha_risk_lag2 = lag(ascvd_10yr_accaha_risk, 2),
#     ascvd_10yr_frs_risk_lag2 = lag(ascvd_10yr_frs_risk, 2),
#     ascvd_10yr_frs_simple_risk_lag2 = lag(ascvd_10yr_frs_simple_risk, 2)
#   ) %>%
#   ungroup() %>%
#   filter(time == 2) 
# 
# stats_accaha <- val.prob(
#   df$ascvd_10yr_accaha_risk_lag2,
#   df$event_ascvd_zero
# )
# 
# 
# # CHANGE LAG BELOW TO ONE!!!!
# 
# df <- observed_events %>%
#   group_by(pid) %>%
#   mutate(
#     ascvd_10yr_accaha_risk_lag2 = lag(ascvd_10yr_accaha_risk, 2),
#     ascvd_10yr_frs_risk_lag2 = lag(ascvd_10yr_frs_risk, 2),
#     ascvd_10yr_frs_simple_risk_lag2 = lag(ascvd_10yr_frs_simple_risk, 2)
#   ) %>%
#   ungroup()
# 
# val.prob(
#   df$ascvd_10yr_accaha_risk_lag2[df$time == 2],
#   df$event_ascvd[df$time == 2]
# )
# 
# val.prob(
#   df$ascvd_10yr_frs_risk_lag2[df$time == 2],
#   df$event_chd_zero[df$time == 2]
# )
# 
# val.prob(
#   df$ascvd_10yr_frs_simple_risk_lag2[df$time == 4],
#   df$event_chd_zero[df$time == 4]
# )
# 
# 
# val.prob(
#   df$ascvd_10yr_frs_risk_lag2,
#   df$event_ascvd
# )
# 
# mgcv::gam(
#   event_ascvd ~ sex * age + smk + bmi + dm + hrx + liprx + hdl + ldl + sbp + 
#     lag1_smk + lag1_bmi + lag1_dm + lag1_hrx + lag1_liprx + lag1_hdl + lag1_ldl + lag1_sbp,
#   data = offspring_long,
#   family = binomial(link = "logit")
# )
