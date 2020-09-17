# Define covariate models -------------------------------------------------

X.models <-
  list (
    # logit model for probability of smoking
    "smk" = list(
      formula = smk ~ sex + age0 + educ_1 + educ_2 + educ_3 + marital_1 +
        marital_2 + eversmk + pre_dpd + pre_bmi + pre_dm + pre_sbp + pre_cpd +
        pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
        lag1_sbp + lag1_ldl + lag1_hrx + lag1_liprx + as.factor(time), 
      link = "logit",
      family = "binomial"
    ),
    
    # log-linear model for number of cigarettes per day among smokers
    "cpd" = list(
      formula = cpd ~ sex + age0 + educ_1 + educ_2 + educ_3 + marital_1 + 
        marital_2 + eversmk + pre_dpd + pre_bmi + pre_dm + pre_sbp + pre_cpd +
        pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
        lag1_sbp + lag1_ldl + lag1_hrx + lag1_liprx + as.factor(time),
      link = "log",
      family = "normal",
      restrict = list(
        subset = 'smk == 1 & cpd > 0',
        otherwise = 0
      )
    ),
    
    # logit model for probability of smoking
    "drk" = list(
      formula = drk ~ sex + age0 + educ_1 + educ_2 + educ_3 + marital_1 + 
        marital_2 + eversmk + pre_dpd + pre_bmi + pre_dm + pre_sbp + pre_cpd +
        pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
        lag1_sbp + lag1_ldl + lag1_hrx + lag1_liprx + cpd + as.factor(time),
      link = "logit",
      family = "binomial"
    ),
    
    # log-linear model for number of drinks per day
    "dpd" = list(
      formula = dpd ~ sex + age0 + educ_1 + educ_2 + educ_3 + marital_1 + 
        marital_2 + eversmk + pre_dpd + pre_bmi + pre_dm + pre_sbp + pre_cpd +
        pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
        lag1_sbp + lag1_ldl + lag1_hrx + lag1_liprx + cpd + as.factor(time),
      link = "log",
      family = "normal",
      restrict = list(
        subset = 'drk == 1 & dpd > 0',
        otherwise = 0
      )
    ),
    
    # linear model of BMI
    "bmi" = list(
      formula = bmi ~ sex + age0 + educ_1 + educ_2 + educ_3 + marital_1 + 
        marital_2 + eversmk + pre_dpd + pre_bmi + pre_dm + pre_sbp + pre_cpd +
        pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
        lag1_sbp + lag1_ldl + lag1_hrx + lag1_liprx + cpd + dpd + as.factor(time),
      link = "identity",
      family = "normal"
    ),
    
    # logit model for diabetes (failure)
    "dm" = list(
      formula = dm ~ sex + age0 + educ_1 + educ_2 + educ_3 + marital_1 + 
        marital_2 + eversmk + pre_dpd + pre_bmi + pre_dm + pre_sbp + pre_cpd +
        pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi +
        lag1_sbp + lag1_ldl + lag1_hrx + lag1_liprx + cpd + dpd + bmi + as.factor(time),
      link = "logit",
      family = "binomial",
      restrict = list(
        subset = 'lag1_dm == 0',
        otherwise = 1
      )
    ),
    
    # linear model for systolic blood pressure
    "sbp" = list(
      formula = sbp ~ sex + age0 + educ_1 + educ_2 + educ_3 + marital_1 +
        marital_2 + eversmk + pre_dpd + pre_bmi + pre_dm + pre_sbp + pre_cpd +
        pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
        lag1_sbp + lag1_ldl + lag1_hrx + lag1_liprx + cpd + dpd + bmi + dm +
        as.factor(time),
      link = "identity",
      family = "normal"
    ),
    
    # linear model for LDL cholesterol
    "ldl" = list(
      formula = ldl ~ sex + age0 + educ_1 + educ_2 + educ_3 + marital_1 + 
        marital_2 + eversmk + pre_dpd + pre_bmi + pre_dm + pre_sbp + pre_cpd +
        pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
        lag1_sbp + lag1_ldl + lag1_hrx + lag1_liprx + cpd + dpd + bmi + dm + sbp + 
        as.factor(time),
      link = "identity",
      family = "normal"
    ),
    
    # logit model for hypertension meds
    "hrx" = list(
      formula = hrx ~ sex + age0 + educ_1 + educ_2 + educ_3 + marital_1 + 
        marital_2 + eversmk + pre_dpd + pre_bmi + pre_dm + pre_sbp + pre_cpd +
        pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
        lag1_sbp + lag1_ldl + lag1_hrx + lag1_liprx + cpd + dpd + bmi + dm + sbp + ldl + 
        as.factor(time),
      link = "logit",
      family = "binomial"
    ),
    
    # logit model for lipids meds
    "liprx" = list(
      formula = liprx ~ sex + age0 + educ_1 + educ_2 + educ_3 + marital_1 + marital_2 +
        eversmk + pre_dpd + pre_bmi + pre_dm + pre_sbp + pre_cpd +
        pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
        lag1_sbp + lag1_ldl + lag1_hrx + lag1_liprx + cpd + dpd + bmi + dm + sbp + ldl + 
        as.factor(time), 
      link = "logit",
      family = "binomial"
    )
  )


# Define outcome model ----------------------------------------------------

# pooled logistic model of CVD events
Y.model <- 
  list(
    formula = event_chd ~ sex + age0 + educ_1 + educ_2 + educ_3 +
      marital_1 + marital_2 + eversmk + pre_dpd + pre_bmi + pre_dm + pre_sbp +
      pre_cpd + pre_ldl + pre_hrx + pre_liprx + cpd + dpd + bmi +
      dm + sbp + ldl + hrx + liprx + as.factor(time),
    link = "logit",
    family = "binomial"
  )



# Define competing event model --------------------------------------------

# pooled logistic model of death due to other causes
D.model <- 
  list(
    formula = event_dth ~ sex + age0 + educ_1 + educ_2 + educ_3 +
      marital_1 + marital_2 + eversmk + pre_dpd + pre_bmi + pre_dm + pre_sbp +
      pre_cpd + pre_ldl + pre_hrx + pre_liprx + cpd + dpd + bmi +
      dm + sbp + ldl + hrx + liprx + as.factor(time),
    link = "logit",
    family = "binomial"
  )
