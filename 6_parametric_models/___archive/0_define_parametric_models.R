# Define covariate models -------------------------------------------------

no_lags_specification <- list(
  "age"   = c("lag1_age"),
  "smk"   = c(covs_fixed, covs_tv[1]),
  "bmi"   = c(covs_fixed, covs_tv[1:2]),
  "dm"    = c(covs_fixed, covs_tv[1:3]),
  "hrx"   = c(covs_fixed, covs_tv[1:4]),
  "liprx" = c(covs_fixed, covs_tv[1:5]),
  "tc"    = c(covs_fixed, covs_tv[1:6]),
  "hdl"   = c(covs_fixed, covs_tv[1:7]),
  "sbp"   = c(covs_fixed, covs_tv[1:8]),
  "event_chd" = c(covs_fixed, covs_tv),
  "event_dth" = c(covs_fixed, covs_tv)
)

lag_1_specification <- list(
  "age"   = c("lag1_age"),
  "smk"   = c(covs_fixed, covs_tv[1], covs_lag1),
  "bmi"   = c(covs_fixed, covs_tv[1:2], covs_lag1),
  "dm"    = c(covs_fixed, covs_tv[1:3], covs_lag1[-3]),
  "hrx"   = c(covs_fixed, covs_tv[1:4], covs_lag1),
  "liprx" = c(covs_fixed, covs_tv[1:5], covs_lag1),
  "tc"    = c(covs_fixed, covs_tv[1:6], covs_lag1),
  "hdl"   = c(covs_fixed, covs_tv[1:7], covs_lag1),
  "sbp"   = c(covs_fixed, covs_tv[1:8], covs_lag1),
  "event_chd" = c(covs_fixed, covs_tv, covs_lag1),
  "event_dth" = c(covs_fixed, covs_tv, covs_lag1)
)

lag_2_specification <- list(
  "age"   = c("lag1_age"),
  "smk"   = c(covs_fixed, covs_tv[1], covs_lag1, covs_lag2),
  "bmi"   = c(covs_fixed, covs_tv[1:2], covs_lag1, covs_lag2),
  "dm"    = c(covs_fixed, covs_tv[1:3], covs_lag1[-3], covs_lag2[-3]),
  "hrx"   = c(covs_fixed, covs_tv[1:4], covs_lag1, covs_lag2),
  "liprx" = c(covs_fixed, covs_tv[1:5], covs_lag1, covs_lag2),
  "tc"    = c(covs_fixed, covs_tv[1:6], covs_lag1, covs_lag2),
  "hdl"   = c(covs_fixed, covs_tv[1:7], covs_lag1, covs_lag2),
  "sbp"   = c(covs_fixed, covs_tv[1:8], covs_lag1, covs_lag2),
  "event_chd" = c(covs_fixed, covs_tv, covs_lag1, covs_lag2),
  "event_dth" = c(covs_fixed, covs_tv, covs_lag1, covs_lag2)
)

X.models_lag0 <- define_models(no_lags_specification, type = "X")
Y.model_lag0 <- define_models(no_lags_specification, type = "Y")
D.model_lag0 <- define_models(no_lags_specification, type = "D")

X.models_lag1 <- define_models(lag_1_specification, type = "X")
Y.model_lag1 <- define_models(lag_1_specification, type = "Y")
D.model_lag1 <- define_models(lag_1_specification, type = "D")

X.models_lag2 <- define_models(lag_2_specification, type = "X")
Y.model_lag2 <- define_models(lag_2_specification, type = "Y")
D.model_lag2 <- define_models(lag_2_specification, type = "D")


analytic_long_std <- analytic_long
analytic_long_std[, c(vars_cont)] <- scale(analytic_long_std[, c(vars_cont)])

lasso_lag1_fits <- 
  map2(
    names(lag_1_specification),
    lag_1_specification,
    function(x, y) {
      print(x)
      print(which(grepl(paste0("(", paste(vars_cont, collapse = ")|("), ")"), y)))
      if (x %in% covs_tv) {
        cutpt <- 1
      } else {
        cutpt <- 0
      }
      select_covariates(
        outcome = x,
        covariates = if (x != "age") c(y, "time") else y,
        #contvars = if (x != "age") which(grepl(paste0("(", paste(vars_cont, collapse = ")|("), ")"), y)) else NULL,
        #degree = if (x != "age") 2 else NULL,
        data = subset(analytic_long_std, time >= cutpt & time <= 2 & !pid %in% censored_ids),
        logit = x %in% vars_binary #[-1], 
        #X.dependent = TRUE
      )
    }
  )  

lasso_lag1_specification <- 
  lapply(lasso_lag1_fits, get, x = "selected_covariates")

lasso_lag1_specification <- c(c("lag1_age"), c("lag1_smk"), lasso_lag1_specification[3:5], c("lag1_liprx"), lasso_lag1_specification[7:11])

names(lasso_lag1_specification) <- names(lag_1_specification)

# lasso_lag1_specification[["event_chd"]] <- 
#   unique(unlist(lasso_lag1_specification[-1]))

X.lasso_lag1 <- define_models(lasso_lag1_specification, type = "X", time = FALSE)
Y.lasso_lag1 <- define_models(lasso_lag1_specification, type = "Y", time = FALSE)
D.lasso_lag1 <- define_models(lasso_lag1_specification, type = "D", time = FALSE)

lasso_lag2_fits <- 
  map2(
    names(lag_2_specification)[-1],
    lag_2_specification[-1],
    function(x, y) {
      print(x)
      print(which(grepl(paste0("(", paste(vars_cont, collapse = ")|("), ")"), y)))
      if (x %in% covs_tv) {
        cutpt <- 1
      } else {
        cutpt <- 0
      }
      select_covariates(
        outcome = x,
        covariates = if (x != "age") c(y, "time") else y,
        #contvars = if (x != "age") which(grepl(paste0("(", paste(vars_cont, collapse = ")|("), ")"), y)) else NULL,
        #degree = if (x != "age") 2 else NULL,
        glmnet = TRUE,
        data = subset(analytic_long, time >= cutpt & time <= 2),
        logit = x %in% vars_binary, 
        #X.dependent = TRUE
      )
    }
  )  

lasso_lag2_specification <- 
  lapply(lasso_lag2_fits, function(y) {
    get(y, x = "selected_covariates")[-1]
  })

lasso_lag2_specification <- c(c("lag1_age"), lasso_lag2_specification)
names(lasso_lag2_specification) <- names(lag_2_specification)

# lasso_lag1_specification[["event_chd"]] <- 
#   unique(unlist(lasso_lag1_specification[-1]))

X.lasso_lag2 <- define_models(lasso_lag2_specification, type = "X", time = FALSE)
Y.lasso_lag2 <- define_models(lasso_lag2_specification, type = "Y", time = FALSE)
D.lasso_lag2 <- define_models(lasso_lag2_specification, type = "D", time = FALSE)

lasso_lag1_fits <- 
  map2(
    names(lag_1_specification)[-1],
    lag_1_specification[-1],
    function(x, y) {
      print(x)
      print(which(grepl(paste0("(", paste(vars_cont, collapse = ")|("), ")"), y)))
      if (x %in% covs_tv) {
        cutpt <- 1
      } else {
        cutpt <- 0
      }
      select_covariates(
        outcome = x,
        covariates = if (x != "age") c(y, "time") else y,
        #contvars = if (x != "age") which(grepl(paste0("(", paste(vars_cont, collapse = ")|("), ")"), y)) else NULL,
        #degree = if (x != "age") 2 else NULL,
        glmnet = TRUE,
        data = subset(analytic_long, time >= cutpt & time <= 2),
        logit = x %in% vars_binary
        #X.dependent = TRUE
      )
    }
  )  

lasso_lag1_specification <- 
  lapply(lasso_lag1_fits, function(y) {
    get(y, x = "selected_covariates")[-1]
  })

lasso_lag1_specification <- c(c("lag1_age"), lasso_lag1_specification)
names(lasso_lag1_specification) <- names(lag_1_specification)

X.lasso_lag1 <- define_models(lasso_lag1_specification, type = "X", time = FALSE)
Y.lasso_lag1 <- define_models(lasso_lag1_specification, type = "Y", time = FALSE)
D.lasso_lag1 <- define_models(lasso_lag1_specification, type = "D", time = FALSE)
# lasso_lag1_s



fit <- rlassologit(
  x = model.matrix(
    ~ (sex + log(age) + smk + log(bmi) + dm + hrx + liprx + log(tc) + log(hdl) + log(sbp))^2 - 1,
    data = na.omit(analytic_long)
    ),
  y = na.omit(analytic_long)[['event_chd']]
)

fit <- glinternet.cv(
  X = model.matrix(
    ~ sex + log(age) + smk + log(bmi) + dm + hrx + liprx + log(tc) + log(hdl) +
      log(sbp) + lag1_smk + log(lag1_bmi + 1) + lag1_dm + lag1_hrx + lag1_liprx + 
      log(lag1_tc + 1) + log(lag1_hdl + 1) + log(lag1_sbp + 1) + time - 1,
    data = na.omit(analytic_long)
  ),
  Y = na.omit(analytic_long)[['event_chd']],
  numLevels = c(2, 1, 2, 1, 2, 2, 2, 1, 1, 1, 2, 1, 2, 2, 2, 1, 1, 1, 1, 1),
  nFolds = 10,
  family = "binomial"
)
val.prob(predict(glm(
  event_chd ~ sex + log(age) + smk + log(bmi) + dm + hrx + liprx + log(tc) + log(hdl) +
    log(sbp) + sex:smk + sex:dm + sex:hrx + smk:dm + smk:liprx + log(age):log(bmi) + log(age):log(hdl) +
    log(tc):log(sbp) + log(hdl):log(sbp) + hrx:log(age) + liprx:log(age) + sex:log(tc) + dm:log(tc) + liprx:log(tc) + 
    liprx:log(hdl) + sex:log(sbp) + smk:log(sbp) + hrx:log(sbp),
  data = analytic_long,
  family = binomial(link = "logit")
), type = "response"),
na.omit(analytic_long)[['event_chd']]
)


lag_1_specification_w_ints <- 
  list(
    "age"   = c("lag1_age"),
    "smk"   = c(covs_fixed, covs_tv[1], covs_lag1),
    "bmi"   = c(covs_fixed, covs_tv[1:2], covs_lag1, "age:lag1_bmi", "age:time", "lag1_bmi:lag1_dm", "lag1_bmi:lag1_hdl", "lag1_hrx:lag1_liprx"),
    "dm"    = c(covs_fixed, covs_tv[1:3], covs_lag1[-3], "age:time", "lag1_bmi:lag1_hrx", "lag1_bmi:lag1_sbp", "lag1_hrx:lag1_liprx"),
    "hrx"   = c(covs_fixed, covs_tv[1:4], covs_lag1, "bmi:lag1_dm", "lag1_bmi:lag1_sbp", "lag1_dm:time", "lag1_hrx:lag1_tc", "lag1_hrx:time", "lag1_sbp:time"),
    "liprx" = c(covs_fixed, covs_tv[1:5], covs_lag1),
    "tc"    = c(covs_fixed, covs_tv[1:6], covs_lag1, "sex:lag1_tc", "sex:time", "liprx:lag1_bmi", "liprx:lag1_tc", "liprx:time", "lag1_tc:time"),
    "hdl"   = c(covs_fixed, covs_tv[1:7], covs_lag1, "sex:lag1_bmi", "sex:lag1_hdl", "sex:time", "bmi:lag1_hdl", "tc:lag1_hdl", "tc:time", "lag1_hdl:time"),
    "sbp"   = c(covs_fixed, covs_tv[1:8], covs_lag1, "sex:age", "age:bmi", "age:lag1_hdl", "age:time", "tc:lag1_hdl", "lag1_hrx:time"),
    "event_chd" = c(covs_fixed, covs_tv, covs_lag1, "dm:hrx"),
    "event_dth" = c(covs_fixed, covs_tv, covs_lag1)
  )

X.ints <- define_models(lag_1_specification_w_ints, type = "X")
Y.ints <- define_models(lag_1_specification_w_ints, type = "Y")
D.ints <- define_models(lag_1_specification_w_ints, type = "D")


# R2 = 0.45
ols(
  log(sbp) ~ rcs(age, 3) * sex * bmi + smk + hrx * (lag1_hrx + lag1_sbp * lag2_sbp) + tc * lag1_hdl + dm + catg(time),
  data = analytic_long
)

# C = 0.913
lrm(
  hrx ~ rcs(age, 3) * bmi + sex + lag1_sbp * lag2_sbp + lag1_hrx * lag2_hrx + dm + liprx, 
  data = analytic_long
)

# C = 0.932
lrm(
  liprx ~ rcs(age, 3) + bmi + sex + lag1_tc * lag2_tc + lag1_hdl * lag2_hdl + lag1_liprx * lag2_liprx + dm + hrx + lag1_hrx * lag1_sbp + lag2_hrx * lag2_sbp + catg(time), 
  data = analytic_long
)

# R2 = 0.664
ols(
  log(hdl) ~ rcs(age, 3) * sex * bmi + liprx * (lag1_liprx + lag1_tc * lag2_tc + lag1_hdl * lag2_hdl) + dm + catg(time),
  data = analytic_long
)

# R2 = 0.985
ols(
  log(tc) ~ rcs(age, 3) * sex * bmi + liprx * (lag1_liprx + tc + lag1_tc * lag2_tc + lag1_hdl * lag2_hdl) + dm + catg(time),
  data = analytic_long_std
)

# C = 0.807
lrm(
  dm ~ rcs(age, 3) * (rcs(bmi, 3) + rcs(lag1_bmi, 3) + sex) + smk + lag1_sbp * lag1_hrx + lag1_liprx * (lag1_tc + lag1_hdl) + catg(time), 
  data = analytic_long_std
)

lrm(
  event_chd ~ rcs(age, 3) * sex + smk + bmi + lag1_bmi + lag2_bmi + liprx * (tc + lag1_tc * lag2_tc + hdl + lag1_hdl * lag2_hdl) + lag1_liprx * lag2_liprx + dm + hrx * sbp + lag1_hrx * lag1_sbp + lag2_hrx * lag2_sbp + hrx * lag1_hrx * lag2_hrx + catg(time),
  data = subset(analytic_long_std, !pid %in% censored_ids)
)

lrm(
  event_dth ~ rcs(age, 3) * sex + smk + bmi + lag1_bmi + lag2_bmi + liprx * (tc + lag1_tc * lag2_tc + hdl + lag1_hdl * lag2_hdl) + lag1_liprx * lag2_liprx + dm + hrx * sbp + lag1_hrx * lag1_sbp + lag2_hrx * lag2_sbp + hrx * lag1_hrx * lag2_hrx + catg(time),
  data = subset(offspring_long, !pid %in% censored_ids)
)

lrm(
  smk ~ rcs(age, 3) * sex + lag1_smk + catg(time),
  data = analytic_long
)

ols(
  bmi ~ rcs(age, 3) * sex + catg(time) + lag1_bmi * lag2_bmi + lag1_hrx * lag1_dm + lag1_hdl + lag1_sbp + lag1_tc + lag1_liprx,
  data = analytic_long
)

X.models <-
  list (
    # logit model for probability of smoking
    "smk" = list(
      formula = smk ~ rcs(age, 3) * sex + lag1_smk + as.factor(time), 
      link = "logit",
      family = "binomial"
    ),
    
    # linear model of BMI
    "bmi" = list(
      formula = bmi ~ rcs(age, 3) * sex + lag1_bmi * lag2_bmi + 
        lag1_hrx * lag1_dm + lag1_hdl + lag1_sbp + lag1_tc + 
        lag1_liprx + as.factor(time),
      link = "identity",
      family = "normal"
    ),
    
    # logit model for diabetes (failure)
    "dm" = list(
      formula = dm ~ rcs(age, 3) * (rcs(bmi, 3) + rcs(lag1_bmi, 3) + sex) + smk +
        lag1_sbp * lag1_hrx + lag1_liprx * (lag1_tc + lag1_hdl) + as.factor(time),
      link = "logit",
      family = "binomial",
      restrict = list(
        subset = 'lag1_dm == 0',
        otherwise = 1
      )
    ),
    
    # logit model for hypertension meds
    "hrx" = list(
      formula = hrx ~ rcs(age, 3) * bmi + sex + lag1_sbp * lag2_sbp + 
        lag1_hrx * lag2_hrx + dm + liprx + as.factor(time),
      link = "logit",
      family = "binomial"
    ),
    
    # logit model for lipids meds
    "liprx" = list(
      formula = liprx ~ rcs(age, 3) + bmi + sex + lag1_tc * lag2_tc + 
        lag1_hdl * lag2_hdl + lag1_liprx * lag2_liprx + dm + hrx + 
        lag1_hrx * lag1_sbp + lag2_hrx * lag2_sbp + as.factor(time), 
      link = "logit",
      family = "binomial"
    ),
    
    # linear model for LDL cholesterol
    "tc" = list(
      formula = tc ~ rcs(age, 3) * sex * bmi + liprx * 
        (lag1_liprx + lag1_tc * lag2_tc + lag1_hdl * lag2_hdl) + 
        dm + as.factor(time),
      link = "log",
      family = "normal"
    ),
    
    # linear model for LDL cholesterol
    "hdl" = list(
      formula = hdl ~rcs(age, 3) * sex * bmi + liprx *
        (lag1_liprx + tc + lag1_tc * lag2_tc + lag1_hdl * lag2_hdl) + 
        dm + as.factor(time),
      link = "log",
      family = "normal"
    ),
    
    # linear model for systolic blood pressure
    "sbp" = list(
      formula = sbp ~ rcs(age, 3) * sex * bmi + smk + hrx * (lag1_hrx + lag1_sbp * lag2_sbp) +
        tc * lag1_hdl + dm + as.factor(time),
      link = "log",
      family = "normal"
    )
    
    
  )


# Define outcome model ----------------------------------------------------

# pooled logistic model of CVD events
Y.model <- 
  list(
    formula = event_chd ~ rcs(age, 3) * sex + smk + bmi + lag1_bmi + lag2_bmi +
      liprx * (tc + lag1_tc * lag2_tc + hdl + lag1_hdl * lag2_hdl) + 
      lag1_liprx * lag2_liprx + dm + hrx * sbp + lag1_hrx * lag1_sbp + 
      lag2_hrx * lag2_sbp + hrx * lag1_hrx * lag2_hrx + as.factor(time),
    link = "logit",
    family = "binomial"
  )



# Define competing event model --------------------------------------------

# pooled logistic model of death due to other causes
D.model <- 
  list(
    formula = event_dth ~ rcs(age, 3) * sex + smk + bmi + lag1_bmi + lag2_bmi +
      liprx * (tc + lag1_tc * lag2_tc + hdl + lag1_hdl * lag2_hdl) + 
      lag1_liprx * lag2_liprx + dm + hrx * sbp + lag1_hrx * lag1_sbp + 
      lag2_hrx * lag2_sbp + hrx * lag1_hrx * lag2_hrx + as.factor(time),
    link = "logit",
    family = "binomial"
  )
