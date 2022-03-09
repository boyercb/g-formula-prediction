
# Create and clean fixed and time-varying covariates ----------------------

offspring_long <- offspring_long %>%
  mutate(
    # female sex (0/1)
    sex = sex - 1,
    
    # educational level (categorical)
    educ_1 = if_else(educ3 == 6, 1, 0),              # less than high school
    educ_2 = if_else(educ3 == 1 | educ3 == 2, 1, 0), # high school/some college
    educ_3 = if_else(educ3 == 3, 1, 0),              # bachelor's degree
    educ_4 = if_else(educ3 == 4 | educ3 == 5, 1, 0), # postgraduate education
    
    # marital status (categorical)
    marital_1 = if_else(marital == 1, 1, 0),            # single
    marital_2 = if_else(marital == 2, 1, 0),            # married
    marital_3 = if_else(marital %in% c(3, 4, 5), 1, 0), # widowed, divorced, or separated
    
    # ever smoked (0/1)
    eversmk = case_when(
      currsmk1 == 1 | currsmk2 == 1 | currsmk3 == 1 ~ 1, 
      currsmk1 == 1 & (is.na(currsmk2) | is.na(currsmk3)) ~ 1,
      currsmk2 == 1 & (is.na(currsmk1) | is.na(currsmk3)) ~ 1,
      currsmk3 == 1 & (is.na(currsmk1) | is.na(currsmk2)) ~ 1,
      currsmk1 == 0 & currsmk2 == 0 & currsmk3 == 0 ~ 0, 
      currsmk1 == 0 & (is.na(currsmk2) | is.na(currsmk3)) ~ 0,
      currsmk2 == 0 & (is.na(currsmk1) | is.na(currsmk3)) ~ 0,
      currsmk3 == 0 & (is.na(currsmk1) | is.na(currsmk2)) ~ 0,
      TRUE ~ NA_real_
    ),
    
    # baseline standard drinks per day (continuous)
    pre_dpd = (beer_week3 + wine_week3 + liquor_week3) / 7,
    
    # baseline standard drinks per day (categorical)
    pre_dpd_1 = if_else(pre_dpd == 0, 1, 0),               # 0 drinks per day
    pre_dpd_2 = if_else(pre_dpd > 0 & pre_dpd < 2, 1, 0),  # 0 to 2 drinks per day
    pre_dpd_3 = if_else(pre_dpd >= 2 & pre_dpd < 4, 1, 0), # 2 to 4 drinks per day
    pre_dpd_4 = if_else(pre_dpd >= 4, 1, 0),               # >4 drinks per day
    
    # baseline cigarettes per day (categorical)
    pre_cpd_1 = if_else(cpd3 == 0, 1, 0),             # None
    pre_cpd_2 = if_else(cpd3 == 1, 1, 0),             # 1 or fewer
    pre_cpd_3 = if_else(cpd3 > 1 & cpd3 < 5, 1, 0),   # 2 to 4 
    pre_cpd_4 = if_else(cpd3 >= 5 & cpd3 < 25, 1, 0), # 5 to 24
    pre_cpd_5 = if_else(cpd3 > 25, 1, 0),             # 25 or more
    
    # standard drinks per day (continuous)
    dpd = case_when(
      exam %in% c(6, 7) ~ (beer_week + white_wine_week + red_wine_week + liquor_week) / 7,
      exam %in% c(4, 5, 8, 9) ~ (beer_week + wine_week + liquor_week) / 7
    ),
    
    # currently drinks (0/1)
    drk = if_else(dpd > 0, 1, 0),
    
    # standard drinks per day (categorical)
    dpd_1 = if_else(dpd == 0, 1, 0),           # 0 drinks per day
    dpd_2 = if_else(dpd > 0 & dpd < 2, 1, 0),  # 0 to 2 drinks per day
    dpd_3 = if_else(dpd >= 2 & dpd < 4, 1, 0), # 2 to 4 drinks per day
    dpd_4 = if_else(dpd >= 4, 1, 0),           # >4 drinks per day
    
    # cigarettes per day (categorical)
    cpd_1 = if_else(cpd == 0, 1, 0),             # None
    cpd_2 = if_else(cpd == 1, 1, 0),             # 1 or fewer
    cpd_3 = if_else(cpd > 1 & cpd < 5, 1, 0),    # 2 to 4 
    cpd_4 = if_else(cpd >= 5 & cpd < 25, 1, 0),  # 5 to 24
    cpd_5 = if_else(cpd > 25, 1, 0),            # 25 or more
    
    # winsorize continuous measures
    sbp = winsorize(sbp, 0.995),
    tc = winsorize(tc, 0.995),
    hdl = winsorize(hdl, 0.995),
    bmi = winsorize(bmi, 0.995),
    calc_ldl = winsorize(calc_ldl, 0.995),
    
    # exam indicators (0/1)
    exam = as.numeric(exam),
    time = exam - 5,
    time_f = factor(time)

  ) %>%
  # rename variables to final covariate names
  rename(
    pre_age = age3,
    pre_cpd = cpd3,
    pre_bmi = bmi3,
    pre_sbp = sbp3,
    pre_ldl = calc_ldl3,
    pre_tc = tc3,
    pre_hdl = hdl3,
    pre_hrx = hrx3,
    pre_liprx = liprx3,
    pre_dm = curr_diab3,
    ldl = calc_ldl,
    dm = curr_diab,
    smk = currsmk
  ) %>%
  # drop intermediate variables
  select(
    -matches("wine"), 
    -starts_with("beer"), 
    -starts_with("liquor"),
    -matches(".*[a-z][1-3]"),
    -bdate
  )

# log transforms
offspring_long <- offspring_long %>%
  mutate(
    # cpd = log(cpd + 1),
    # dpd = log(dpd + 1),
    # pre_cpd = log(pre_cpd + 1),
    # pre_dpd = log(pre_dpd + 1)
    cpd = cpd + 1,
    dpd = dpd + 1,
    pre_cpd = pre_cpd + 1,
    pre_dpd = pre_dpd + 1
  )


# Fix all time variables to start of follow up ----------------------------

offspring_long <-
  offspring_long %>%
  mutate_at(date_vars, function(x) (x - offspring_long$date0))

offspring_long <-
  offspring_long %>%
  group_by(pid) %>%
  mutate(
    # approx_date = case_when(
    #   is.na(date) & !is.na(lag(date)) ~ lag(date) + 365.25 * 4,
    #   is.na(date) & !is.na(lag(date, 2)) ~ lag(date, 2) + 365.25 * 8,
    #   is.na(date) & !is.na(lag(date, 3)) ~ lag(date, 3) + 365.25 * 12,
    #   is.na(date) & !is.na(lag(date, 4)) ~ lag(date, 4) + 365.25 * 16,
    #   is.na(date) & !is.na(lag(date, 5)) ~ lag(date, 5) + 365.25 * 20,
    #   is.na(date) & !is.na(lag(date, 6)) ~ lag(date, 6) + 365.25 * 24,
    #   !is.na(date) ~ date
    # ),
    eventdate = pmin(chddate, strokedate, datedth),
    enddate = case_when(
      exam %in% 5:6 & !is.na(lead(date)) ~ lead(date) - 1,
      exam == 5 & is.na(lead(date)) ~ 1461,
      exam == 6 & is.na(lead(date)) ~ 2922,
      exam == 7 ~ 3652.5
    ),
    date = case_when(
      !is.na(date) ~ date,
      is.na(date) & fupat == "1_0_0" & exam == 6 ~ 1462,
      is.na(date) & fupat == "1_0_0" & exam == 7 ~ 2923,
      is.na(date) & fupat == "1_0_1" & exam == 6 ~ 1462,
      is.na(date) & fupat == "1_1_0" & exam == 7 ~ 2923,
    ),
    year = floor(date / 365.25),
    datedth = if_else(is.na(datedth), lastcon, datedth)
  ) %>%
  ungroup()


# Create variables for ASCVD ----------------------------------------------

offspring_long <- 
  offspring_long %>%
  mutate(
    ascvd = as.numeric(chd == 1 | stroke == 1),
    ascvddate = pmin(chddate, strokedate)
  )


# Censor individuals at the appropriate time ------------------------------

offspring_long <- 
  offspring_long %>%
  mutate(
    drop = case_when(
      # if you died or got CHD prior to interval then censor by dropping
      # subsequent exams
      datedth < date | chddate < date | strokedate < date ~ 1, 
      
      # if you have two or more consecutive misses drop after possibly carrying one forward
      str_detect(fupat, "^1_0_0") & exam > 6 ~ 1,
      # str_detect(fupat, "^1_1_0_0") & exam > 6 ~ 1,
      # str_detect(fupat, "^1_[01]_1_0_0") & exam > 7 ~ 1,
      # str_detect(fupat, "^1_[01]_[01]_1_0_0") & exam > 8 ~ 1,
      
      # if you died or got CHD during this interval then keep 
      (date <= datedth & datedth <= enddate) | (date <= chddate & chddate <= enddate) |  (date <= strokedate & strokedate <= enddate) ~ 0,
      
      # if you died or got CHD in a future interval then keep 
      (datedth > enddate) | (chddate > enddate) | (strokedate > enddate) ~ 0,
      
      # always keep exam 6
      exam == 6 ~ 0,
      
      # if you completed all exams then keep
      fupat == "1_1_1" ~ 0,
    )
  )

# censor by dropping censored observations 
offspring_long <- filter(offspring_long, drop == 0)


# Handle missing data -----------------------------------------------------

# last observation carry forward
offspring_long <-
  offspring_long %>%
  group_by(pid) %>%
  tidyr::fill(all_of(covs_tv[!covs_tv %in% "age"])) %>%
  ungroup()


# Create and clean outcomes -----------------------------------------------

offspring_long <- offspring_long %>%
  group_by(pid) %>%
  mutate(
    event_ascvd = case_when(
      # if never got ASCVD code as 0
      ascvd == 0 ~ 0,
      
      # if got ASCVD and event occurs within interval code as 1
      ascvd == 1 & (date <= ascvddate) & (ascvddate <= enddate) ~ 1,
      
      # if event occurs in future code as 0
      ascvd == 1 & ascvddate > enddate ~ 0
    ),
    event_chd = case_when(
      # if never got CHD code as 0
      chd == 0 ~ 0,
      
      # if got CHD and event occurs within interval code as 1
      chd == 1 & (date <= chddate) & (chddate <= enddate) ~ 1,
      
      # if event occurs in future code as 0
      chd == 1 & chddate > enddate ~ 0
    ),
    event_dth = case_when(
      # if never died code as 0
      is.na(dthrvwd) ~ 0,
      
      # if died and death occurred within interval code as 1
      dthrvwd == 1 & chd == 0 & (date <= datedth) & (datedth <= enddate) ~ 1,
      dthrvwd == 1 & chd == 1 & (date <= datedth) & (datedth <= enddate) & datedth < chddate ~ 1,
      dthrvwd == 1 & chd == 1 & (date <= datedth) & (datedth <= enddate) & datedth >= chddate ~ 0,
      
      # if event occurs in future code as 0
      dthrvwd == 1 & datedth > enddate ~ 0
    ),
    event_dth_ascvd = case_when(
      # if never died code as 0
      is.na(dthrvwd) ~ 0,
      
      # if died and death occurred within interval code as 1
      dthrvwd == 1 & ascvd == 0 & (date <= datedth) & (datedth <= enddate) ~ 1,
      dthrvwd == 1 & ascvd == 1 & (date <= datedth) & (datedth <= enddate) & datedth < ascvddate ~ 1,
      dthrvwd == 1 & ascvd == 1 & (date <= datedth) & (datedth <= enddate) & datedth >= ascvddate ~ 0,
      
      # if event occurs in future code as 0
      dthrvwd == 1 & datedth > enddate ~ 0
    ),
    event_cen = case_when(
      # if this is last time recorded and didn't die or get CHD define as censored
      time == last(time) & (event_ascvd == 0 | is.na(ascvd)) & event_dth_ascvd == 0 ~ 1,
      TRUE ~ 0
    ),
    enddate = pmin(enddate, ascvddate, datedth)
  ) %>%
  ungroup()

# defensive coding: stop if more than one event per person
stopifnot(
  offspring_long %>%
    group_by(pid) %>%
    summarise(events = sum(event_chd), .groups = "drop") %>%
    filter(events > 1) %>%
    nrow() == 0
)

stopifnot(
  offspring_long %>%
    group_by(pid) %>%
    summarise(events = sum(event_ascvd), .groups = "drop") %>%
    filter(events > 1) %>%
    nrow() == 0
)

stopifnot(
  offspring_long %>%
    group_by(pid) %>%
    summarise(events = sum(event_dth), .groups = "drop") %>%
    filter(events > 1) %>%
    nrow() == 0
)

stopifnot(
  offspring_long %>%
    group_by(pid) %>%
    summarise(events = sum(event_dth_ascvd), .groups = "drop") %>%
    filter(events > 1) %>%
    nrow() == 0
)
# verify distributions of follow up and events
offspring_long %>%
  group_by(pid) %>%
  summarise(
    start = min(exam),
    stop = max(exam),
    chd = sum(event_chd),
    died = sum(event_dth),
    .groups = "drop"
  ) %>%
  count(chd, died, start, stop)

# # A tibble: 18 x 5
#       chd  died start  stop     n
#     <dbl> <dbl> <dbl> <dbl> <int>
#   1     0     0     4     4     7
#   2     0     0     4     5    84
#   3     0     0     4     6    76
#   4     0     0     4     7    51
#   5     0     0     4     8   149
#   6     0     0     4     9  1540
#   7     0     1     4     4    28
#   8     0     1     4     5    40
#   9     0     1     4     6    32
#  10     0     1     4     7    61
#  11     0     1     4     8   152
#  12     0     1     4     9   191
#  13     1     0     4     4    49
#  14     1     0     4     5    52
#  15     1     0     4     6    60
#  16     1     0     4     7    89
#  17     1     0     4     8    79
#  18     1     0     4     9    32


# verify distributions of follow up and events
offspring_long %>%
  group_by(pid) %>%
  summarise(
    start = min(exam),
    stop = max(exam),
    ascvd = sum(event_ascvd),
    died = sum(event_dth_ascvd),
    .groups = "drop"
  ) %>%
  count(ascvd, died, start, stop)

# # A tibble: 17 x 5
#     ascvd  died start  stop     n
#     <dbl> <dbl> <dbl> <dbl> <int>
#   1     0     0     4     5    72
#   2     0     0     4     6    63
#   3     0     0     4     7    26
#   4     0     0     4     8   121
#   5     0     0     4     9  1523
#   6     0     1     4     4    28
#   7     0     1     4     5    37
#   8     0     1     4     6    31
#   9     0     1     4     7    57
#  10     0     1     4     8   135
#  11     0     1     4     9   186
#  12     1     0     4     4    56
#  13     1     0     4     5    67
#  14     1     0     4     6    74
#  15     1     0     4     7   118
#  16     1     0     4     8   124
#  17     1     0     4     9    54


# add lagged values of time-varying covariates
offspring_long <- offspring_long %>%
  group_by(pid) %>%
  mutate(across(all_of(covs_tv), lag, .names = "lag1_{col}", n = 1, default = 0)) %>%
  mutate(across(all_of(covs_tv), lag, .names = "lag2_{col}", n = 2, default = 0)) %>% 
  # mutate(
  #   lag1_age = replace_na(lag1_age, first(pre_age)),
  #   lag1_smk = replace_na(lag1_smk, first(eversmk)),
  #   #lag1_cpd = replace_na(lag1_cpd, first(pre_cpd)),
  #   #lag1_drk = replace_na(lag1_drk, as.numeric(first(pre_dpd) > 0)),
  #   #lag1_dpd = replace_na(lag1_dpd, first(pre_dpd)),
  #   lag1_bmi = replace_na(lag1_bmi, first(pre_bmi)),
  #   lag1_dm = replace_na(lag1_dm, first(pre_dm)),
  #   lag1_sbp = replace_na(lag1_sbp, first(pre_sbp)),
  #   #lag1_ldl = replace_na(lag1_ldl, first(pre_ldl)),
  #   lag1_tc = replace_na(lag1_tc, first(pre_tc)),
  #   lag1_hdl = replace_na(lag1_hdl, first(pre_hdl)),
  #   lag1_hrx = replace_na(lag1_hrx, first(pre_hrx)),
  #   lag1_liprx = replace_na(lag1_liprx, first(pre_liprx))
  # ) %>%
  ungroup()

# instead of carrying forward age (which doesn't make any sense) replace those
# with NA ages with average age difference between exam cycles
avg_age_diff <-
  offspring_long %>%
  filter(lag1_age != 0) %>%
  summarise(
    age = median(age - lag1_age, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pull(age) 

offspring_long <- offspring_long %>%
  mutate(
    age = if_else(is.na(age), lag1_age + avg_age_diff, age)
  ) %>% 
  group_by(pid) %>%
  mutate(
    # update lags
    lag1_age = lag(age, n = 1, default = 0),
    lag2_age = lag(age, n = 2, default = 0)
  ) %>%
  ungroup()

# offspring_long[, paste0("time_", 0:5)] <-
#   model.matrix(~ as.factor(time) - 1, data = offspring_long)

# add 10-year risk estimates from pooled cohort equations
offspring_long <- offspring_long %>%
  mutate(
    ascvd_10yr_accaha_risk = ascvd_10yr_accaha(
      race = "white",
      gender = if_else(sex == 1, "female", "male"), 
      age = age,
      totchol = tc, 
      hdl = hdl, 
      sbp = sbp,
      bp_med = hrx, 
      smoker = smk, 
      diabetes = dm
    ),
    ascvd_10yr_frs_risk = ascvd_10yr_frs(
      gender = if_else(sex == 1, "female", "male"), 
      age = age,
      totchol = tc, 
      hdl = hdl, 
      sbp = sbp,
      bp_med = hrx, 
      smoker = smk, 
      diabetes = dm
    ),
    ascvd_10yr_frs_simple_risk = ascvd_10yr_frs_simple(
      gender = if_else(sex == 1, "female", "male"), 
      age = age,
      bmi = bmi,
      sbp = sbp,
      bp_med = hrx, 
      smoker = smk, 
      diabetes = dm
    ),
    ascvd_10yr_accaha_risk_updated = ascvd_10yr_accaha(
      race = "white",
      gender = if_else(sex == 1, "female", "male"), 
      age = age,
      totchol = tc, 
      hdl = hdl, 
      sbp = sbp,
      bp_med = hrx, 
      smoker = smk, 
      diabetes = dm,      #    female, female-aa,      male,   male-aa       
      baseline_survival = c(0.9502099, 0.9502099, 0.8939252, 0.8939252)
    ),
    ascvd_10yr_frs_risk_updated = ascvd_10yr_frs(
      gender = if_else(sex == 1, "female", "male"), 
      age = age,
      totchol = tc, 
      hdl = hdl, 
      sbp = sbp,
      bp_med = hrx, 
      smoker = smk, 
      diabetes = dm,      #    female,      male      
      baseline_survival = c(0.8939252, 0.9502099)
    ),
    ascvd_10yr_frs_simple_risk_updated = ascvd_10yr_frs_simple(
      gender = if_else(sex == 1, "female", "male"), 
      age = age,
      bmi = bmi,
      sbp = sbp,
      bp_med = hrx, 
      smoker = smk, 
      diabetes = dm,      #    female,      male      
      baseline_survival = c(0.8939252, 0.9502099)
    )
  )

offspring_long <- select(
  offspring_long,
  all_of(covs_model),
  all_of(id_vars),
  all_of(dvs),
  starts_with("lag"),
  all_of(risk_scores),
  all_of(covs_time),
  "date",
  "enddate"
) %>% drop_na()

# update event_chd coding to conform with g-formula guide
offspring_long <- offspring_long %>%
  mutate(
    event_chd = if_else(event_dth == 1, NA_real_, event_chd),
    event_ascvd = if_else(event_dth_ascvd == 1, NA_real_, event_ascvd)
  )

