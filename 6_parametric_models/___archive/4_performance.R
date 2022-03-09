stats_opt1 <- map_dfr(gformula_predictions_opt1,
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


stats_opt2 <- map_dfr(gformula_predictions_opt2,
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


# CHANGE LAG BELOW TO ONE!!!!

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
#   df$ascvd_10yr_accaha_risk_lag2[df$time == 4],
#   df$event_chd_zero[df$time == 4]
# )
# 
# val.prob(
#   df$ascvd_10yr_frs_risk_lag2[df$time == 4],
#   df$event_chd_zero[df$time == 4]
# )
# 
# val.prob(
#   df$ascvd_10yr_frs_simple_risk_lag2[df$time == 4],
#   df$event_chd_zero[df$time == 4]
# )
