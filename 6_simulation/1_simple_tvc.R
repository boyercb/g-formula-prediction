# simulate data -----------------------------------------------------------

sim <- datagen(5000)

sim_long <-
  sim %>% 
  pivot_longer(
    -id,
    names_to = c("variable", "time"),
    names_pattern = "([LAY])([0-9])"
  ) %>%
  pivot_wider(names_from = "variable") %>%
  mutate(time = as.numeric(time))

sim_w_lags <- mutate(
  sim_long, 
  across(c("L", "A"), lag, default = 0, .names = "lag1_{col}")
)


# define models -----------------------------------------------------------

X.models <-
  list(
    "L" = list(
      formula = L ~ lag1_L + lag1_A, 
      link = "logit",
      family = "binomial"
    ),
    "A" = list(
      formula = A ~ L + lag1_L + lag1_A, 
      link = "logit",
      family = "binomial"
    )
  )

D.model <- 
  list(
    formula = D ~ L + A + as.factor(time),
    link = "logit",
    family = "binomial"
  )

Y.model <- 
  list(
    formula = Y ~ L + A + as.factor(time),
    link = "logit",
    family = "binomial"
  )


# fit models --------------------------------------------------------------

Y.fit <- 
  fit_Y_model(
    Y.model,
    data = sim_w_lags
  )

X.fit <- 
  fit_X_models(
    X.models,
    data = sim_w_lags,
    time = "time"
  )

# D.fit <- 
#   fit_D_model(
#     D.model,
#     data = sim_w_lags
#   )


# run g-formula -----------------------------------------------------------

Y.hat_gformula <-
  gformula_mc(
    Y.fit,
    X.fit,
    #D.fit,
    data = sim_w_lags,
    id = "id",
    time = "time",
    hist.vars = c("lag1_A", "lag1_L"),
    hist.fun = "lag",
    mc.sims = 1000
  )


# run conventional model --------------------------------------------------

Y.fit_conventional <- 
  glm(Y1 ~ A0 + L0, data = sim, family = binomial(link = "logit"))

Y.hat_conventional <-
  predict(Y.fit_conventional, type = 'response')

sim$pred_c <- Y.hat_conventional


# run R package for g-formula ---------------------------------------------

r_pkg_fit <- 
  gformula_survival(
    obs_data = sim_long,
    id = "id",
    time_name = "time",
    time_points = 2,
    covnames = c("L", "A"), 
    covtypes = c("binary", "binary"),
    covparams = list(
      covmodels = c(
        L ~ lag1_L + lag1_A,
        A ~ L + lag1_A + lag1_L
      ),
      covlink = c(
        "logit",
        "logit"
      )
    ),
    outcome_name = "Y",
    ymodel = Y ~ L + A + as.factor(time),
    #compevent_name = "D",
    #compevent_model = D ~ L + A + as.factor(time),
    histvars = list(c("L", "A")),
    histories = c(lagged),
    nsamples = 0,
    sim_data_b = TRUE,
    seed = 4548076,
    show_progress = TRUE,
    model_fits = TRUE
  )

r_pkg_pred <- predict.gformula(
  object = r_pkg_fit,
  obs_data = data.table::as.data.table(sim_long),
  newdata = data.table::as.data.table(sim_long),
  id = "id",
  t0 = 0,
  covnames = c("L", "A"), 
  covtypes = c("binary", "binary"),
  covparams = list(
    covmodels = c(
      L ~ lag1_L + lag1_A,
      A ~ L + lag1_A + lag1_L
    ),
    covlink = c(
      "logit",
      "logit"
    )
  ),
  outcome_name = "Y",
  ymodel = Y ~ L + A + as.factor(time),
  #compevent_name = "D",
  #compevent_model = D ~ L + A + as.factor(time),
  histvars = list(c("L", "A")),
  histories = c(lagged),
  nsamples = 0,
  nsimul = 1000,
  seed = 4548076,
  show_progress = TRUE,
  return_sims = TRUE,
  model_fits = TRUE
)

Y.hat_r_pkg <- r_pkg_pred$sims %>%
  group_by(oldid, time) %>% 
  summarise(pred_r = mean(poprisk)) %>%
  rename(id = oldid)


# combine and validate ----------------------------------------------------

sim_w_pred <-
  cbind(sim_long, 'pred_g' = Y.hat_gformula)

sim_w_pred <-
  left_join(sim_w_pred, Y.hat_r_pkg, by = c("id", "time"))

sim_w_pred <- 
  left_join(sim, filter(sim_w_pred, time == 1))

rms::val.prob(sim_w_pred$pred_r, sim_w_pred$Y1)
rms::val.prob(sim_w_pred$pred_g, sim_w_pred$Y1)
rms::val.prob(sim_w_pred$pred_c, sim_w_pred$Y1)

