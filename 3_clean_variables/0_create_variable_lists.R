id_vars <- c("pid", "idtype")

exam2_vars <- c(
  "educ2" = "b43" # years of education
)

exam3_vars <- c(
  "beer_week3"    = "c83", # BEER-NUMBER PER WEEK
  "beer_days3"    = "c84", # BEER-DAYS PER WEEK DRINK
  "beer_limit3"   = "c85", # BEER-LIMIT AT ONE TIME
  "wine_week3"    = "c86", # WINE-NUMBER PER WEEK
  "wine_days3"    = "c87", # WINE-DAYS PER WEEK DRINK
  "wine_limit3"   = "c88", # WINE-LIMIT AT ONE TIME
  "liquor_week3"  = "c89", # COCKTAILS-NUMBER PER WEEK
  "liquor_days3"  = "c90", # COCKTAILS-DAYS PER WEEK DRINK
  "liquor_limit3" = "c91"  # COCKTAILS-LIMIT AT ONE TIME
)

exam4_vars <- c(
  "beer_week4"    = "d082", # BEER-NUMBER PER WEEK
  "beer_days4"    = "d083", # BEER-DAYS PER WEEK DRINK
  "beer_limit4"   = "d084", # BEER-LIMIT AT ONE TIME
  "wine_week4"    = "d085", # WINE-NUMBER PER WEEK
  "wine_days4"    = "d086", # WINE-DAYS PER WEEK DRINK
  "wine_limit4"   = "d087", # WINE-LIMIT AT ONE TIME
  "liquor_week4"  = "d088", # COCKTAILS-NUMBER PER WEEK
  "liquor_days4"  = "d089", # COCKTAILS-DAYS PER WEEK DRINK
  "liquor_limit4" = "d090", # COCKTAILS-LIMIT AT ONE TIME
  "marital"       = "d399"  # marital status
)

exam5_vars <- c(
  "beer_week5"    = "e310", # BEER-NUMBER PER WEEK
  "beer_days5"    = "e311", # BEER-DAYS PER WEEK DRINK
  "beer_limit5"   = "e312", # BEER-LIMIT AT ONE TIME
  "wine_week5"    = "e313", # WINE-NUMBER PER WEEK
  "wine_days5"    = "e314", # WINE-DAYS PER WEEK DRINK
  "wine_limit5"   = "e315", # WINE-LIMIT AT ONE TIME
  "liquor_week5"  = "e316", # COCKTAILS-NUMBER PER WEEK
  "liquor_days5"  = "e317", # COCKTAILS-DAYS PER WEEK DRINK
  "liquor_limit5" = "e318"  # COCKTAILS-LIMIT AT ONE TIME
)

exam6_vars <- c(
  "beer_week6"        = "f276", # BEER-NUMBER PER WEEK
  "beer_days6"        = "f277", # BEER-DAYS PER WEEK DRINK
  "beer_limit6"       = "f278", # BEER-LIMIT AT ONE TIME
  "white_wine_week6"  = "f279", # WHITE WINE-AVG. NUMBER PER WEEK
  "white_wine_days6"  = "f280", # WHITE WINE-NUMBER OF DAYS DRINK PER WEEK
  "white_wine_limit6" = "f281", # WHITE WINE-DRINKS AT ONE TIME
  "red_wine_week6"    = "f282", # RED WINE-AVG. NUMBER PER WEEK
  "red_wine_days6"    = "f283", # RED WINE-NUMBER OF DAYS DRINK PER WEEK
  "red_wine_limit6"   = "f284", # RED WINE-DRINKS AT ONE TIME
  "liquor_week6"      = "f285", # COCKTAILS-NUMBER PER WEEK
  "liquor_days6"      = "f286", # COCKTAILS-DAYS PER WEEK DRINK
  "liquor_limit6"     = "f287"  # COCKTAILS-LIMIT AT ONE TIME
)

exam7_vars <- c(
  "beer_week7"        = "g104", # BEER-NUMBER PER WEEK
  "beer_days7"        = "g105", # BEER-DAYS PER WEEK DRINK
  "beer_limit7"       = "g106", # BEER-LIMIT AT ONE TIME
  "white_wine_week7"  = "g107", # WHITE WINE-AVG. NUMBER PER WEEK
  "white_wine_days7"  = "g108", # WHITE WINE-NUMBER OF DAYS DRINK PER WEEK
  "white_wine_limit7" = "g109", # WHITE WINE-DRINKS AT ONE TIME
  "red_wine_week7"    = "g110", # RED WINE-AVG. NUMBER PER WEEK
  "red_wine_days7"    = "g111", # RED WINE-NUMBER OF DAYS DRINK PER WEEK
  "red_wine_limit7"   = "g112", # RED WINE-DRINKS AT ONE TIME
  "liquor_week7"      = "g113", # COCKTAILS-NUMBER PER WEEK
  "liquor_days7"      = "g114", # COCKTAILS-DAYS PER WEEK DRINK
  "liquor_limit7"     = "g115"  # COCKTAILS-LIMIT AT ONE TIME
)

exam8_vars <- c(
  "beer_week8"        = "h072", # BEER-NUMBER PER WEEK
  "wine_week8"        = "h075", # WINE-AVG. NUMBER PER WEEK
  "liquor_week8"      = "h078"  # COCKTAILS-NUMBER PER WEEK
)

exam9_vars <- c(
  "beer_week9"        = "j075", # BEER-NUMBER PER WEEK
  "wine_week9"        = "j077", # WINE-AVG. NUMBER PER WEEK
  "liquor_week9"      = "j079"  # COCKTAILS-NUMBER PER WEEK
)

psy_vars <- c(
  "edul3" = "py122", # years of schooling
  "educ3" = "py123"  # highest degree completed
)

exam_varlist <- list(
  exam2_vars,
  exam3_vars,
  exam4_vars,
  exam5_vars,
  exam6_vars,
  exam7_vars,
  exam8_vars,
  exam9_vars,
  psy_vars
)

exam_list <- list(
  csvs$ex1_2d_v3,
  csvs$ex1_3d_v1,
  csvs$ex1_4d_v1,
  csvs$ex1_5d_v1,
  csvs$ex1_6d_v1,
  csvs$ex1_7d_v2,
  csvs$e_exam_ex08_1_0005d,
  csvs$e_exam_ex09_1b_0844d,
  csvs$q_psych_ex03_1_0167d
)

date_vars <-
  c("chddate",
    "chfdate",
    "cvddate",
    "strokedate",
    "datedth",
    "date",
    "lastsoe",
    "lastatt",
    "lastcon")

dvs <- c(
  "event_chd",
  "event_ascvd",
  "event_dth",
  "event_dth_ascvd",
  "event_cen"
)

risk_scores <- c(
  "ascvd_10yr_accaha_risk",
  "ascvd_10yr_frs_risk",
  "ascvd_10yr_frs_simple_risk",
  "ascvd_10yr_accaha_risk_updated",
  "ascvd_10yr_frs_risk_updated",
  "ascvd_10yr_frs_simple_risk_updated"
)

covs_fixed <- c("sex")

covs_tv <- c(
  "age",
  "smk",
  "bmi",
  "dm",
  "hrx",
  "liprx",
  "tc",
  "ldl",
  "hdl",
  "sbp"
)

covs_lag1 <- paste0("lag1_", covs_tv[-1])
covs_lag2 <- paste0("lag2_", covs_tv[-1])

covs_binary <- c("smk", "dm", "hrx", "liprx")
vars_binary <- c("smk", "dm", "hrx", "liprx", dvs)
vars_cont <- c("age", "bmi", "tc", "ldl", "hdl", "sbp")

covs_time <- c(
  "time"
 # "time_0",
 # "time_1",
 # "time_2",
 # "time_3",
 # "time_4",
 # "time_5"
)

covs_refs <- c("")

covs_model <- c(covs_fixed, covs_tv)
covs_dvs <- covs_tv[!covs_tv %in% covs_refs & !grepl("dpd_", covs_tv)]

covs_ivs <- c(
  covs_model[!covs_model %in% c(covs_refs, "smk", "drk")], 
  "as.factor(time)", 
  "I(age0^2)"
)

