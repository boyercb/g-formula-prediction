# Define covariate models -------------------------------------------------

ascvd_interactions <- c(
  "age:tc",
  "age:hdl",
  "age:sbp",
  "age:hrx",
  "age:smk",
  "age:hrx:sbp",
  "hrx:sbp",
  "liprx:hdl",
  "liprx:tc",
  "age:liprx:hdl",
  "age:liprx:tc",
  "lag1_age:lag1_tc",
  "lag1_age:lag1_hdl",
  "lag1_age:lag1_sbp",
  "lag1_age:lag1_hrx",
  "lag1_age:lag1_smk",
  "lag1_age:lag1_hrx:lag1_sbp",
  "lag1_hrx:lag1_sbp",
  "lag1_liprx:lag1_hdl",
  "lag1_liprx:lag1_tc",
  "lag1_age:lag1_liprx:lag1_hdl",
  "lag1_age:lag1_liprx:lag1_tc",
  "lag2_age:lag2_tc",
  "lag2_age:lag2_hdl",
  "lag2_age:lag2_sbp",
  "lag2_age:lag2_hrx",
  "lag2_age:lag2_smk",
  "lag2_age:lag2_hrx:lag2_sbp",
  "lag2_hrx:lag2_sbp",
  "lag2_liprx:lag2_hdl",
  "lag2_liprx:lag2_tc",
  "lag2_age:lag2_liprx:lag2_hdl",
  "lag2_age:lag2_liprx:lag2_tc"
)

# initial experiment value of added lags ----------------------------------

covs_tv_ex1 <- covs_tv[!grepl("(bmi)|(ldl)|(liprx)", covs_tv)]
covs_lag1_ex1 <- covs_lag1[!grepl("(bmi)|(ldl)|(liprx)", covs_lag1)]
covs_lag2_ex1 <- covs_lag2[!grepl("(bmi)|(ldl)|(liprx)", covs_lag2)]
ascvd_interactions_ex1 <- ascvd_interactions[!grepl("(bmi)|(ldl)|(liprx)", ascvd_interactions)]

specifications_ex1 <- list(
  lag0 = list(
    "age"   = c("lag1_age"),
    "smk"   = c(covs_fixed, "poly(age, 2)"),
    "dm"    = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:2]),
    "hrx"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:3]),
    "tc"    = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:4]),
    "hdl"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:5]),
    "sbp"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:6]),
    "event_ascvd" = c(
      covs_fixed,
      "poly(age, 2)",
      covs_tv_ex1[2:7],
      ascvd_interactions_ex1[!grepl("lag", ascvd_interactions_ex1)]
    ),
    "event_dth_ascvd" = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:7])
  ),
  
  lag1 = list(
    "age"   = c("lag1_age"),
    "smk"   = c(covs_fixed, "poly(age, 2)", covs_lag1_ex1),
    "dm"    = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:2], covs_lag1_ex1[!grepl("dm", covs_lag1_ex1)]),
    "hrx"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:3], covs_lag1_ex1),
    "tc"    = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:4], covs_lag1_ex1),
    "hdl"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:5], covs_lag1_ex1),
    "sbp"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:6], covs_lag1_ex1),
    "event_ascvd" = c(
      covs_fixed,
      "poly(age, 2)",
      covs_tv_ex1[2:7],
      covs_lag1_ex1,
      "lag1_age",
      ascvd_interactions_ex1[!grepl("lag2", ascvd_interactions_ex1)]
    ),
    "event_dth_ascvd" = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:7], covs_lag1_ex1)
  ),
  
  lag2 = list(
    "age"   = c("lag1_age"),
    "smk"   = c(covs_fixed, "poly(age, 2)", covs_lag1_ex1, covs_lag2_ex1),
    "dm"    = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:2], covs_lag1_ex1[!grepl("dm", covs_lag1_ex1)], covs_lag2_ex1[!grepl("dm", covs_lag2_ex1)]),
    "hrx"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:3], covs_lag1_ex1, covs_lag2_ex1),
    "tc"    = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:4], covs_lag1_ex1, covs_lag2_ex1),
    "hdl"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:5], covs_lag1_ex1, covs_lag2_ex1),
    "sbp"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:6], covs_lag1_ex1, covs_lag2_ex1),
    "event_ascvd" = c(
      covs_fixed,
      "poly(age, 2)",
      covs_tv_ex1[2:7],
      covs_lag1_ex1,
      covs_lag2_ex1,
      "lag1_age",
      "lag2_age",
      ascvd_interactions_ex1
    ),
    "event_dth_ascvd" = c(covs_fixed, "poly(age, 2)", covs_tv_ex1[2:7], covs_lag1_ex1, covs_lag2_ex1)
  )
)

definitions_ex1 <-
  lapply(specifications_ex1,
         function(x)
           define_models(
             x,
             "event_ascvd",
             "event_dth_ascvd",
             time = TRUE,
             logvars = TRUE,
             sexints = TRUE
           ))

names(definitions_ex1) <- names(specifications_ex1)


# next experiment value of additional variables ---------------------------

covs_tv_ex2 <- list(
  covs_tv[!grepl("(ldl)|(liprx)", covs_tv)],
  covs_tv[!grepl("(tc)|(liprx)", covs_tv)],
  covs_tv[!grepl("(tc)", covs_tv)]
)
covs_lag1_ex2 <- list(
  covs_lag1[!grepl("(ldl)|(liprx)", covs_lag1)],
  covs_lag1[!grepl("(tc)|(liprx)", covs_lag1)],
  covs_lag1[!grepl("(tc)", covs_lag1)]
)
covs_lag2_ex2 <- list(
  covs_lag2[!grepl("(ldl)|(liprx)", covs_lag2)],
  covs_lag2[!grepl("(tc)|(liprx)", covs_lag2)],
  covs_lag2[!grepl("(tc)", covs_lag2)]
)
ascvd_interactions_ex2 <- list(
  ascvd_interactions[!grepl("(ldl)|(liprx)", ascvd_interactions)],
  ascvd_interactions[!grepl("(tc)|(liprx)", ascvd_interactions)],
  ascvd_interactions[!grepl("(tc)", ascvd_interactions)]
)

specifications_ex2 <- list(
  bmi = list(
    "age"   = c("lag1_age"),
    "smk"   = c(covs_fixed, "poly(age, 2)", covs_lag1_ex2[[1]], covs_lag2_ex2[[1]]),
    "bmi"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[1]][2:2], covs_lag1_ex2[[1]], covs_lag2_ex2[[1]]),
    "dm"    = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[1]][2:3], covs_lag1_ex2[[1]][!grepl("dm", covs_lag1_ex2[[1]])], covs_lag2_ex2[[1]][!grepl("dm", covs_lag2_ex2[[1]])]),
    "hrx"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[1]][2:4], covs_lag1_ex2[[1]], covs_lag2_ex2[[1]]),
    "tc"    = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[1]][2:5], covs_lag1_ex2[[1]], covs_lag2_ex2[[1]]),
    "hdl"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[1]][2:6], covs_lag1_ex2[[1]], covs_lag2_ex2[[1]]),
    "sbp"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[1]][2:7], covs_lag1_ex2[[1]], covs_lag2_ex2[[1]]),
    "event_ascvd" = c(
      covs_fixed,
      "poly(age, 2)",
      covs_tv_ex2[[1]][2:8],
      covs_lag1_ex2[[1]],
      covs_lag2_ex2[[1]],
      "lag1_age",
      "lag2_age",
      ascvd_interactions_ex2[[1]]
    ),
    "event_dth_ascvd" = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[1]][2:8], covs_lag1_ex2[[1]], covs_lag2_ex2[[1]])
  ),
  
  ldl = list(
    "age"   = c("lag1_age"),
    "smk"   = c(covs_fixed, "poly(age, 2)", covs_lag1_ex2[[2]], covs_lag2_ex2[[2]]),
    "bmi"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[2]][2:2], covs_lag1_ex2[[2]], covs_lag2_ex2[[2]]),
    "dm"    = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[2]][2:3], covs_lag1_ex2[[2]][!grepl("dm", covs_lag1_ex2[[2]])], covs_lag2_ex2[[2]][!grepl("dm", covs_lag2_ex2[[2]])]),
    "hrx"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[2]][2:4], covs_lag1_ex2[[2]], covs_lag2_ex2[[2]]),
    "ldl"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[2]][2:5], covs_lag1_ex2[[2]], covs_lag2_ex2[[2]]),
    "hdl"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[2]][2:6], covs_lag1_ex2[[2]], covs_lag2_ex2[[2]]),
    "sbp"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[2]][2:7], covs_lag1_ex2[[2]], covs_lag2_ex2[[2]]),
    "event_ascvd" = c(
      covs_fixed,
      "poly(age, 2)",
      covs_tv_ex2[[2]][2:8],
      covs_lag1_ex2[[2]],
      covs_lag2_ex2[[2]],
      "lag1_age",
      "lag2_age",
      ascvd_interactions_ex2[[2]]
    ),
    "event_dth_ascvd" = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[2]][2:8], covs_lag1_ex2[[2]], covs_lag2_ex2[[2]])
  ),
  
  liprx = list(
    "age"   = c("lag1_age"),
    "smk"   = c(covs_fixed, "poly(age, 2)", covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
    "bmi"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:2], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
    "dm"    = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:3], covs_lag1_ex2[[3]][!grepl("dm", covs_lag1_ex2[[3]])], covs_lag2_ex2[[3]][!grepl("dm", covs_lag2_ex2[[3]])]),
    "hrx"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:4], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
    "liprx" = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:5], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
    "ldl"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:6], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
    "hdl"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:7], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
    "sbp"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:8], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
    "event_ascvd" = c(
      covs_fixed,
      "poly(age, 2)",
      covs_tv_ex2[[3]][2:9],
      covs_lag1_ex2[[3]],
      covs_lag2_ex2[[3]],
      "lag1_age",
      "lag2_age",
      ascvd_interactions_ex2[[3]]
    ),
    "event_dth_ascvd" = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:9], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]])
  )
)

definitions_ex2 <-
  lapply(specifications_ex2,
         function(x)
           define_models(
             x,
             "event_ascvd",
             "event_dth_ascvd",
             time = TRUE,
             logvars = TRUE,
             sexints = TRUE
           ))

names(definitions_ex2) <- names(specifications_ex2)


# finally value of relaxing/changing functional form assumptions ----------

specifications_ex3 <- list(
  gam = list(
    "age"   = c("lag1_age"),
    "smk"   = c(covs_fixed, "poly(age, 2)", covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
    "bmi"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:2], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
    "dm"    = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:3], covs_lag1_ex2[[3]][!grepl("dm", covs_lag1_ex2[[3]])], covs_lag2_ex2[[3]][!grepl("dm", covs_lag2_ex2[[3]])]),
    "hrx"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:4], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
    "liprx" = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:5], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
    "ldl"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:6], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
    "hdl"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:7], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
    "sbp"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:8], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
    "event_ascvd" = c(
      covs_fixed,
      "poly(age, 2)",
      covs_tv_ex2[[3]][2:9],
      covs_lag1_ex2[[3]],
      covs_lag2_ex2[[3]],
      "lag1_age",
      "lag2_age",
      ascvd_interactions_ex2[[3]]
    ),
    "event_dth_ascvd" = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:9], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]])
  )
)

definitions_ex3 <-
  lapply(specifications_ex3,
         function(x)
           define_models(
             x,
             "event_ascvd",
             "event_dth_ascvd",
             time = TRUE,
             logvars = TRUE,
             sexints = TRUE
           ))

definitions_ex3$gam$X$ldl$link <- "identity"
definitions_ex3$gam$X$hdl$link <- "identity"
definitions_ex3$gam$X$sbp$link <- "identity"

names(definitions_ex3) <- names(specifications_ex3)


