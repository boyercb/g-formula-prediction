# create offspring dataset
offspring <- select(csvs$vr_wkthru_ex09_1_1001d,
                   -starts_with("bg"),
                   -starts_with("creat"),
                   -starts_with("dbp"), 
                   -starts_with("dlvh"),
                   -starts_with("fasting"),
                   -starts_with("hgt"),
                   -starts_with("hip"),
                   -starts_with("trig"),
                   -starts_with("vent"),
                   -starts_with("waist"),
                   -starts_with("wgt"),
                   -starts_with("dmrx"))

# merge verified CVD survival data
offspring <- left_join(
  x = offspring, 
  y = csvs$vr_survcvd_2017_a_1194d, 
  by = id_vars
)

# merge verified stroke survival data
offspring <- left_join(
  x = offspring, 
  y = csvs$vr_svstk_2017_a_1196d, 
  by = id_vars
)

# merge verified all cause mortality data
offspring <- left_join(
  x = offspring, 
  y = csvs$vr_survdth_2017_a_1192d,
  by = id_vars
)

# merge verified diabetes diagnosis data
offspring <- left_join(
  x = offspring, 
  y = csvs$vr_diab_ex09_1_1002d,
  by = id_vars
)

# merge exam-specific data
for (i in 1:length(exam_list)) {
  offspring <- left_join(
    x = offspring, 
    y = select(exam_list[[i]], all_of(id_vars), exam_varlist[[i]]) %>%
      set_names(c(id_vars, names(exam_varlist[[i]]))),
    by = id_vars
  )
}

# define indicator for whether exam 4 (baseline) is complete
offspring$complete <- 
  offspring %>%
  select(
    age4,
    currsmk4,
    bmi4,
    curr_diab4,
    hrx4,
    liprx4,
    tc4,
    calc_ldl4,
    hdl4,
    sbp4,
    age5,
    currsmk5,
    bmi5,
    curr_diab5,
    hrx5,
    liprx5,
    tc5,
    calc_ldl5,
    hdl5,
    sbp5
  ) %>%
  select_if(function(x) any(is.na(x))) %>%
  is.na() %>%
  rowSums() == 0 

offspring$complete <- as.numeric(offspring$complete)

# define baseline date as date at exam 4 or mean if skipped exam 4
offspring$bdate <- 
  if_else(is.na(offspring$date4), mean(offspring$date4, na.rm = TRUE), offspring$date4)

# apply inclusion/exclusion criteria
offspring <-
  offspring %>% 
  filter(datedth > date5 | is.na(datedth)) %>%  # alive at 5th exam
  filter(pmin(cvddate, strokedate) > date5) %>% # no prior CVD diagnosis at 5th exam
  filter(age5 >= 40 & age5 <= 79) %>%           # age 40-79 at 5th exam
  filter(att5 == 1 & complete == 1)             # attended 5th exam and has complete baseline data

# create a variable representing follow up pattern
offspring <- unite(offspring, fupat, att5:att7, remove = FALSE)


rm(csvs, exam_list)
