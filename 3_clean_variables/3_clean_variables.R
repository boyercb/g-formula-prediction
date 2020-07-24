
# Create and clean fixed and time-varying covariates ----------------------

analytic_long <- analytic_long %>%
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
    
    # exam indicators (0/1)
    exam = as.numeric(exam),
    time = exam - 4,
    time_f = factor(time)
  ) %>%
  # rename variables to final covariate names
  rename(
    pre_cpd = cpd3,
    pre_bmi = bmi3,
    pre_sbp = sbp3,
    pre_ldl = calc_ldl3,
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


# Fix all time variables to start of follow up ----------------------------

analytic_long <-
  analytic_long %>%
  mutate_at(date_vars, function(x) (x - analytic_long$date0))

analytic_long <-
  analytic_long %>%
  group_by(pid) %>%
  mutate(
    date = case_when(
      is.na(date) & !is.na(lag(date)) ~ lag(date) + 365.25 * 4,
      is.na(date) & !is.na(lag(date, 2)) ~ lag(date, 2) + 365.25 * 8,
      is.na(date) & !is.na(lag(date, 3)) ~ lag(date, 3) + 365.25 * 12,
      is.na(date) & !is.na(lag(date, 4)) ~ lag(date, 4) + 365.25 * 16,
      is.na(date) & !is.na(lag(date, 5)) ~ lag(date, 5) + 365.25 * 20,
      is.na(date) & !is.na(lag(date, 6)) ~ lag(date, 6) + 365.25 * 24,
      !is.na(date) ~ date
    ),
    edate = case_when(
      exam %in% 4:8 ~ lead(date) - 1,
      exam == 9 ~ date + 365.25 * 4
    ),
    datedth = if_else(is.na(datedth), lastcon, datedth)
  ) %>%
  ungroup()


# Censor individuals at the appropriate time ------------------------------

analytic_long <- 
  analytic_long %>%
  mutate(
    drop = case_when(
      # if you died or got CHD prior to interval then censor by dropping
      # subsequent exams
      datedth < date | chddate < date ~ 1, 
      
      # if you have two or more consecutive misses drop after possibly carrying one forward
      str_detect(fupat, "^1_0_0") & exam > 5 ~ 1,
      str_detect(fupat, "^1_1_0_0") & exam > 6 ~ 1,
      str_detect(fupat, "^1_1_1_0_0") & exam > 7 ~ 1,
      str_detect(fupat, "^1_1_1_1_0_0") & exam > 8 ~ 1,
      
      # if you died or got CHD during this interval then keep 
      (date <= datedth & datedth <= edate) | (date <= chddate & chddate <= edate) ~ 0,
      
      # if you died or got CHD in a future interval then keep 
      (datedth > edate) | (chddate > edate) ~ 0,
      
      # always keep exam 4
      exam == 4 ~ 0,
      
      # if you completed all exams then keep
      fupat == "1_1_1_1_1_1" ~ 0,
    )
  )

# censor by dropping censored observations 
analytic_long <- filter(analytic_long, drop == 0)


# Handle missing data -----------------------------------------------------

# last observation carry forward
analytic_long <- 
  analytic_long %>%
  group_by(pid) %>%
  fill(covs_tv) %>%
  ungroup()


# Create and clean outcomes -----------------------------------------------

analytic_long <- analytic_long %>%
  mutate(
    event_chd = case_when(
      # if never got CHD code as 0
      chd == 0 ~ 0,
      
      # if got CHD and event occurs within interval code as 1
      chd == 1 & (date <= chddate) & (chddate <= edate) ~ 1,
      
      # if event occurs in future code as 0
      chd == 1 & chddate > edate ~ 0
    ),
    event_dth = case_when(
      # if never died code as 0
      is.na(dthrvwd) ~ 0,
      
      # if died and death occured within interval code as 1
      dthrvwd == 1 & chd == 0 & (date <= datedth) & (datedth <= edate) ~ 1,
      dthrvwd == 1 & chd == 1 & (date <= datedth) & (datedth <= edate) & datedth < chddate ~ 1,
      dthrvwd == 1 & chd == 1 & (date <= datedth) & (datedth <= edate) & datedth >= chddate ~ 0,
      
      # if event occurs in future code as 0
      dthrvwd == 1 & datedth > edate ~ 0
    )
  )

# defensive coding: stop if more than one event per person
stopifnot(
  analytic_long %>%
    group_by(pid) %>%
    summarise(events = sum(event_chd)) %>%
    filter(events > 1) %>%
    nrow() == 0
)

stopifnot(
  analytic_long %>%
    group_by(pid) %>%
    summarise(events = sum(event_dth)) %>%
    filter(events > 1) %>%
    nrow() == 0
)

# verify distributions of follow up and events
analytic_long %>%
  group_by(pid) %>%
  summarise(
    start = min(exam),
    stop = max(exam),
    chd = sum(event_chd),
    died = sum(event_dth)
  ) %>%
  count(chd, died, start, stop)

# # A tibble: 11 x 5
#      chd  died start stop      n
#    <dbl> <dbl> <chr> <chr> <int>
#  1     0     0     4     5    73
#  2     0     0     4     6    61
#  3     0     0     4     7    24
#  4     0     0     4     8   138
#  5     0     0     4     9  1635
#  6     0     1     4     4    28
#  7     0     1     4     5    43
#  8     0     1     4     6    33
#  9     0     1     4     7    65
# 10     0     1     4     8   161
# 11     0     1     4     9   150
# 12     1     0     4     4    49
# 13     1     0     4     5    52
# 14     1     0     4     6    61
# 15     1     0     4     7    91
# 16     1     0     4     8    81
# 17     1     0     4     9    25

analytic_long <- select(
  analytic_long,
  covs_model,
  id_vars,
  dvs,
  time
)
