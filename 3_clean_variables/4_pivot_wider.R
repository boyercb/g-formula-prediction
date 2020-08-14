# pivot to wide format for lmtp package
analytic_wide <- 
  analytic_long %>%
  drop_na(.) %>%
  mutate(event_chd = replace(event_chd, event_dth == 1, NA)) %>%
  pivot_longer(
    cols = c(covs_tv, paste0("lag1_", covs_tv),  dvs),
    names_to = "variable"
  ) %>%
  pivot_wider(
    names_from = c("variable", "time"),
    values_from = "value"
  ) 

# define appropriate censoring indicator for lmtp
analytic_wide <- 
  mutate(
    analytic_wide,
    cens_0 = if_else(is.na(event_chd_0), 0, 1),
    cens_1 = if_else(is.na(event_chd_1), 0, 1),
    cens_2 = if_else(is.na(event_chd_2), 0, 1),
    cens_3 = if_else(is.na(event_chd_3), 0, 1),
    cens_4 = if_else(is.na(event_chd_4), 0, 1),
    cens_5 = if_else(is.na(event_chd_5), 0, 1),
  )

# code outcomes in appropriate way for lmtp (i.e. carry forward last obs.)
analytic_wide <- event_locf(analytic_wide, paste0("event_chd_", 0:5))

analytic_wide <- 
  mutate(
    analytic_wide,
    event_chd_6 = event_chd_5,
    event_chd_5 = event_chd_4,
    event_chd_4 = event_chd_3,
    event_chd_3 = event_chd_2,
    event_chd_2 = event_chd_1,
    event_chd_1 = event_chd_0,
    event_chd_0 = 0
  )
