# Define covariate models -------------------------------------------------

covtypes <- c(
  "binary",
  "normal",
  "binary",
  "normal",
  "normal",
  "binary",
  "normal",
  "normal",
  "binary",
  "binary"
)

restrictions <- list(
  c('cpd', 'smk == 1 & cpd > 0', simple_restriction, 0),
  c('dpd', 'drk == 1 & dpd > 0', simple_restriction, 0),
  c('dm', 'lag1_dm == 0', simple_restriction, 1)
)

covparams <- list(
  covmodels = c(
    # logit model for probability of smoking 
    smk ~ sex + age0 + I(age0 ^ 2) + educ_1 + educ_2 + educ_3 + marital_1 + marital_2 +
      eversmk + pre_dpd_1 + pre_dpd_2 + pre_dpd_3 + pre_bmi + pre_dm +
      pre_sbp + pre_cpd_1 + pre_cpd_2 + pre_cpd_3 + pre_cpd_4 +
      pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
      lag1_sbp + lag1_ldl + lag1_hrx + lag1_liprx + as.factor(time), 
    
    # log-linear model for number of cigarettes per day among smokers
    cpd ~ sex + age0 + I(age0 ^ 2) + educ_1 + educ_2 + educ_3 + marital_1 + marital_2 +
      eversmk + pre_dpd_1 + pre_dpd_2 + pre_dpd_3 + pre_bmi + pre_dm +
      pre_sbp + pre_cpd_1 + pre_cpd_2 + pre_cpd_3 + pre_cpd_4 +
      pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
      lag1_sbp + lag1_ldl + lag1_hrx + lag1_liprx + as.factor(time), 
    
    # logit model for probability of smoking
    drk ~ sex + age0 + I(age0 ^ 2) + educ_1 + educ_2 + educ_3 + marital_1 + marital_2 +
      eversmk + pre_dpd_1 + pre_dpd_2 + pre_dpd_3 + pre_bmi + pre_dm +
      pre_sbp + pre_cpd_1 + pre_cpd_2 + pre_cpd_3 + pre_cpd_4 +
      pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
      lag1_sbp + lag1_ldl + lag1_hrx + cpd + as.factor(time), 
    
    # log-linear model for number of drinks per day
    dpd ~ sex + age0 + I(age0 ^ 2) + educ_1 + educ_2 + educ_3 + marital_1 + marital_2 +
      eversmk + pre_dpd_1 + pre_dpd_2 + pre_dpd_3 + pre_bmi + pre_dm +
      pre_sbp + pre_cpd_1 + pre_cpd_2 + pre_cpd_3 + pre_cpd_4 +
      pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
      lag1_sbp + lag1_ldl + lag1_hrx + cpd + as.factor(time), 
    
    # linear model of BMI 
    bmi ~ sex + age0 + I(age0 ^ 2) + educ_1 + educ_2 + educ_3 + marital_1 + marital_2 +
      eversmk + pre_dpd_1 + pre_dpd_2 + pre_dpd_3 + pre_bmi + pre_dm +
      pre_sbp + pre_cpd_1 + pre_cpd_2 + pre_cpd_3 + pre_cpd_4 +
      pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
      lag1_sbp + lag1_ldl + lag1_hrx + cpd + dpd + as.factor(time),
    
    # logit model for diabetes (failure) 
    dm ~ sex + age0 + I(age0 ^ 2) + educ_1 + educ_2 + educ_3 + marital_1 + marital_2 +
      eversmk + pre_dpd_1 + pre_dpd_2 + pre_dpd_3 + pre_bmi + pre_dm +
      pre_sbp + pre_cpd_1 + pre_cpd_2 + pre_cpd_3 + pre_cpd_4 +
      pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
      lag1_sbp + lag1_ldl + lag1_hrx + cpd + dpd + bmi + as.factor(time),
    
    # linear model for systolic blood pressure 
    sbp ~ sex + age0 + I(age0 ^ 2) + educ_1 + educ_2 + educ_3 + marital_1 + marital_2 +
      eversmk + pre_dpd_1 + pre_dpd_2 + pre_dpd_3 + pre_bmi + pre_dm +
      pre_sbp + pre_cpd_1 + pre_cpd_2 + pre_cpd_3 + pre_cpd_4 +
      pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
      lag1_sbp + lag1_ldl + lag1_hrx + cpd + dpd + bmi + dm + as.factor(time),
    
    # linear model for LDL cholesterol
    ldl ~ sex + age0 + I(age0 ^ 2) + educ_1 + educ_2 + educ_3 + marital_1 + marital_2 +
      eversmk + pre_dpd_1 + pre_dpd_2 + pre_dpd_3 + pre_bmi + pre_dm +
      pre_sbp + pre_cpd_1 + pre_cpd_2 + pre_cpd_3 + pre_cpd_4 +
      pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
      lag1_sbp + lag1_ldl + lag1_hrx + cpd + dpd + bmi + dm + sbp + as.factor(time),
    
    # logit model for hypertension meds
    hrx ~ sex + age0 + I(age0 ^ 2) + educ_1 + educ_2 + educ_3 + marital_1 + marital_2 +
      eversmk + pre_dpd_1 + pre_dpd_2 + pre_dpd_3 + pre_bmi + pre_dm +
      pre_sbp + pre_cpd_1 + pre_cpd_2 + pre_cpd_3 + pre_cpd_4 +
      pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
      lag1_sbp + lag1_ldl + lag1_hrx + cpd + dpd + bmi + dm + sbp + ldl + as.factor(time),
    
    # logit model for lipids meds
    liprx ~ sex + age0 + I(age0 ^ 2) + educ_1 + educ_2 + educ_3 + marital_1 + marital_2 +
      eversmk + pre_dpd_1 + pre_dpd_2 + pre_dpd_3 + pre_bmi + pre_dm +
      pre_sbp + pre_cpd_1 + pre_cpd_2 + pre_cpd_3 + pre_cpd_4 +
      pre_ldl + pre_hrx + pre_liprx + lag1_cpd + lag1_dpd + lag1_bmi + lag1_dm +
      lag1_sbp + lag1_ldl + lag1_liprx + cpd + dpd + bmi + dm + sbp + ldl + as.factor(time)
  ),
  covlink = c(
    "logit",
    "log",
    "logit",
    "log",
    "identity",
    "logit",
    "identity",
    "identity",
    "logit",
    "logit"
  )
)


# Define outcome model ----------------------------------------------------

ymodel <- 
  event_chd ~ sex + age0 + I(age0 ^ 2) + educ_1 + educ_2 + educ_3 +
  marital_1 + marital_2 + eversmk + pre_dpd_1 + pre_dpd_2 + pre_dpd_3 +
  pre_bmi + pre_dm + pre_sbp + pre_cpd_1 + pre_cpd_2 + pre_cpd_3 +
  pre_cpd_4 + pre_ldl + pre_hrx + pre_liprx + cpd + dpd + bmi +
  dm + sbp + ldl + hrx + liprx + as.factor(time) 


# Define competing event model --------------------------------------------

compevent_model <- 
  event_dth ~ sex + age0 + I(age0 ^ 2) + educ_1 + educ_2 + educ_3 +
  marital_1 + marital_2 + eversmk + pre_dpd_1 + pre_dpd_2 + pre_dpd_3 +
  pre_bmi + pre_dm + pre_sbp + pre_cpd_1 + pre_cpd_2 + pre_cpd_3 +
  pre_cpd_4 + pre_ldl + pre_hrx + pre_liprx + cpd + dpd + bmi +
  dm + sbp + ldl + hrx + liprx + as.factor(time) 


# Run gcomputation algorithm ----------------------------------------------

fit_pred <- 
  gformula_survival(
    obs_data = drop_na(analytic_long),
    id = "pid",
    time_name = "time",
    time_points = 4,
    covnames = covs_dvs, 
    covtypes = covtypes,
    covparams = covparams,
    outcome_name = "event_chd",
    ymodel = ymodel,
    compevent_name = "event_dth",
    compevent_model = compevent_model,
    restrictions = restrictions,
    basecovs = covs_fixed[!covs_fixed %in% covs_refs],
    histvars = list(covs_dvs),
    histories = c(lagged),
    nsamples = 0,
    sim_data_b = TRUE,
    nsimul = GFORM_SIM,
    seed = 4548076,
    show_progress = TRUE,
    model_fits = TRUE
  )


# predict risk on hold-out sample -----------------------------------------

newdata <- drop_na(analytic_long)

predrisk <- predict.gformula(
  object = fit_pred,
  obs_data = data.table::as.data.table(drop_na(analytic_long)),
  newdata = data.table::as.data.table(newdata),
  id = "pid",
  t0 = 0,
  covnames = covs_dvs, 
  covtypes = covtypes,
  covparams = covparams,
  outcome_name = "event_chd",
  ymodel = ymodel,
  compevent_name = "event_dth",
  compevent_model = compevent_model,
  restrictions = restrictions,
  basecovs = covs_fixed[!covs_fixed %in% covs_refs],
  histvars = list(covs_dvs),
  histories = c(lagged),
  nsamples = 0,
  nsimul = 500,
  seed = 4548076,
  show_progress = TRUE,
  model_fits = TRUE
)

Y.hat <- predrisk$Y.hat %>%
  pivot_longer(`0`:`3`, names_to = "time", values_to = "pred") %>%
  mutate(time = as.numeric(time))

analytic_hat <-
  left_join(analytic_long, Y.hat, by = c("pid", "time"))

stats <- map(0:3, function(t) {
  df <- filter(analytic_hat, time == t)
  rms::val.prob(df$pred, df$event_chd)
})

aucs <- map(0:3, function(t) {
  df <- filter(analytic_hat, time == t)
  pROC::roc(df$event_chd, df$pred)
})

map(aucs, pROC::ggroc)


# plot some sample trajectories -------------------------------------------

samples <- newdata[1:4,]

trajectories <- predict.gformula(
  object = fit_pred,
  obs_data = data.table::as.data.table(drop_na(analytic_long)),
  newdata = data.table::as.data.table(samples),
  id = "pid",
  t0 = 0,
  covnames = covs_dvs, 
  covtypes = covtypes,
  covparams = covparams,
  outcome_name = "event_chd",
  ymodel = ymodel,
  compevent_name = "event_dth",
  compevent_model = compevent_model,
  restrictions = restrictions,
  basecovs = covs_fixed[!covs_fixed %in% covs_refs],
  histvars = list(covs_dvs),
  histories = c(lagged),
  nsamples = 0,
  nsimul = 1000,
  seed = 4548076,
  model_fits = TRUE,
  return_sims = TRUE
)

ggplot(trajectories$sims, aes(x = time, y = survival)) +
  #facet_grid(~factor(id)) +
  geom_line(aes(group = factor(id)), alpha = 0.04) +
  stat_summary(geom = "line", fun = mean, color = "blue", size = 1.2) +
  stat_summary(geom = "point", fun = mean, color = "blue") +
  g_theme() +
  labs(
    x = "Follow up",
    y = "Predicted survival"
  ) +
  coord_cartesian(expand = F)

