# pivot to wide format for lmtp package
offspring_wide <- 
  offspring_long %>%
#  mutate(event_chd = replace(event_chd, event_dth == 1, NA)) %>%
  pivot_longer(
    cols = c(
      all_of(covs_tv),
      paste0("lag1_", covs_tv),
      paste0("lag2_", covs_tv),
      all_of(dvs),
      all_of(risk_scores),
      "date",
      "enddate"
    ),
    names_to = "variable"
  ) %>%
  pivot_wider(
    names_from = c("variable", "time"),
    values_from = "value"
  ) 

# define appropriate censoring indicator for lmtp
offspring_wide <- 
  mutate(
    offspring_wide,
    cens_0 = if_else(is.na(event_chd_0), 0, 1),
    cens_1 = if_else(is.na(event_chd_1), 0, 1),
    cens_2 = if_else(is.na(event_chd_2), 0, 1)
    # cens_3 = if_else(is.na(event_chd_3), 0, 1),
    # cens_4 = if_else(is.na(event_chd_4), 0, 1),
    # cens_5 = if_else(is.na(event_chd_5), 0, 1),
  )

# code outcomes in appropriate way for lmtp (i.e. carry forward last obs.)
offspring_wide <- event_locf(offspring_wide, paste0("event_chd_", 0:2))
offspring_wide <- event_locf(offspring_wide, paste0("event_dth_", 0:2))
offspring_wide <- event_locf(offspring_wide, paste0("event_ascvd_", 0:2))
offspring_wide <- event_locf(offspring_wide, paste0("event_dth_ascvd_", 0:2))

# compile observed events
observed_events_wide <- 
  offspring_wide %>%
  select(pid,
         starts_with("event_chd_"),
         starts_with("event_dth_"),
         starts_with("event_ascvd_"),
         starts_with("event_dth_ascvd_"),
         starts_with("cens_"),
         starts_with("ascvd_10yr_")
         ) %>%
  select(pid,
         ends_with("0"),
         ends_with("1"),
         ends_with("2"),
         ends_with("3"),
         ends_with("4"),
         ends_with("5"))

# fix coding of CHD event to be 0 to infinity if death occurs to be consistent
# with Fine and Gray
observed_events_wide <- 
  observed_events_wide %>% 
  rowwise() %>%
  mutate(
    died = max(c_across(matches("event_dth_[0-5]")), na.rm = T),
    died_ascvd = max(c_across(matches("event_dth_ascvd_[0-5]")), na.rm = T)
  ) %>% 
  ungroup() %>%
  mutate(
    event_chd_zero_0 = replace(event_chd_0, is.na(event_chd_0) & died == 1, 0),
    event_chd_zero_1 = replace(event_chd_1, is.na(event_chd_1) & died == 1, 0),
    event_chd_zero_2 = replace(event_chd_2, is.na(event_chd_2) & died == 1, 0),
    # event_chd_zero_3 = replace(event_chd_3, is.na(event_chd_3) & died == 1, 0),
    # event_chd_zero_4 = replace(event_chd_4, is.na(event_chd_4) & died == 1, 0),
    # event_chd_zero_5 = replace(event_chd_5, is.na(event_chd_5) & died == 1, 0),
    event_ascvd_zero_0 = replace(event_ascvd_0, is.na(event_ascvd_0) & died_ascvd == 1, 0),
    event_ascvd_zero_1 = replace(event_ascvd_1, is.na(event_ascvd_1) & died_ascvd == 1, 0),
    event_ascvd_zero_2 = replace(event_ascvd_2, is.na(event_ascvd_2) & died_ascvd == 1, 0),
    # event_ascvd_zero_3 = replace(event_ascvd_3, is.na(event_ascvd_3) & died_ascvd == 1, 0),
    # event_ascvd_zero_4 = replace(event_ascvd_4, is.na(event_ascvd_4) & died_ascvd == 1, 0),
    # event_ascvd_zero_5 = replace(event_ascvd_5, is.na(event_ascvd_5) & died_ascvd == 1, 0),
  )

observed_events <-
  observed_events_wide %>%
  pivot_longer(
    cols = -c("pid", "died", "died_ascvd"), 
    names_to = c("variable", "time"),
    names_pattern = c("(.*)_([0-6])$"),
    values_to = "values"
  ) %>%
  pivot_wider(
    names_from = "variable",
    values_from = "values"
  ) %>% 
  mutate(
    time = as.numeric(time)
  )

# identify all who are loss to follow up by the end of study
censored_ids <-
  observed_events_wide$pid[is.na(observed_events_wide$event_ascvd_zero_2)]



offspring_coxph <- left_join(
  offspring_long %>%
    group_by(pid) %>%
    filter(time == max(time)) %>%
    select(pid, time, enddate, event_ascvd) %>%
    ungroup(),
  offspring_wide %>%
    select(pid, sex, paste0(covs_tv, "_0")),
  by = "pid"
)


# offspring_wide <- 
#   mutate(
#     offspring_wide,
#     event_chd_6 = event_chd_5,
#     event_chd_5 = event_chd_4,
#     event_chd_4 = event_chd_3,
#     event_chd_3 = event_chd_2,
#     event_chd_2 = event_chd_1,
#     event_chd_1 = event_chd_0,
#     event_chd_0 = 0
#   )

