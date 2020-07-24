# create analytic dataset
analytic <- select(framingham$vr_wkthru_ex09_1_1001d,
                   -starts_with("bg"),
                   -starts_with("creat"),
                   -starts_with("dbp"), 
                   -starts_with("dlvh"),
                   -starts_with("fasting"),
                   -starts_with("hdl"),
                   -starts_with("hgt"),
                   -starts_with("hip"),
                   -starts_with("tc"),
                   -starts_with("trig"),
                   -starts_with("vent"),
                   -starts_with("waist"),
                   -starts_with("wgt"),
                   -starts_with("dmrx"))

# merge verified CVD survival data
analytic <- left_join(
  x = analytic, 
  y = framingham$vr_survcvd_2014_a_1023d, 
  by = id_vars
)

# merge verified all cause mortality data
analytic <- left_join(
  x = analytic, 
  y = framingham$vr_survdth_2014_a_1025d,
  by = id_vars
)

# merge verified diabetes diagnosis data
analytic <- left_join(
  x = analytic, 
  y = framingham$vr_diab_ex09_1_1002d,
  by = id_vars
)

# merge exam-specific data
for (i in 1:length(exam_list)) {
 analytic <- left_join(
    x = analytic, 
    y = select(exam_list[[i]], id_vars, exam_varlist[[i]]) %>%
        set_names(c(id_vars, names(exam_varlist[[i]]))),
    by = id_vars
  )
}

# define indicator for whether exam 4 (baseline) is complete
analytic$complete4 <- 
  analytic %>%
  select(
    matches(".*4$"),
    cpd3,
    bmi3,
    sbp3,
    calc_ldl3,
    hrx3,
    liprx3,
    curr_diab3,
    beer_week3,
    wine_week3,
    liquor_week3,
    currsmk3
  ) %>%
  select_if(function(x) any(is.na(x))) %>%
  is.na() %>%
  rowSums() == 0 
analytic$complete4 <- as.numeric(analytic$complete4)
  
# define baseline date as date at exam 4 or mean if skipped exam 4
analytic$bdate <- 
  if_else(is.na(analytic$date4), mean(analytic$date4, na.rm = TRUE), analytic$date4)

# apply inclusion/exclusion criteria
analytic <-
  analytic %>% 
  filter(datedth > bdate | is.na(datedth)) %>%  # alive at 4th exam
  filter(cvddate > bdate) %>%                   # no prior CVD diagnosis at 4th exam
  filter(age4 <= 70) %>%                        # age 70 or younger at 4th exam
  filter(att4 == 1 & complete4 == 1)            # attended 4th exam and has complete baseline data
  
# create a variable representing follow up pattern
analytic <- unite(analytic, fupat, att4:att9, remove = FALSE)



