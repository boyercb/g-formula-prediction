
# globals -----------------------------------------------------------------

NSIMS <- 500

# scenario 1 --------------------------------------------------------------
# L1 and L2 are time-varying risk factors that determine whether treatment A is
# administered and there is treatment-confounder feedback. L5 is a time-varying
# risk factor that is unaffected by treatment. L3 and L4 are fixed baseline risk
# factors with L3 affecting whether future treatment is administered. Lagged
# effects, where present, are exponentially decreasing.

scenario1 <- list(
  alpha = matrix(
    # Int,               L1,               L2,            L3,            L4,              L5,                 A
    c(  0, c(0, 2^(-(1:3))), c(0, 2^(-(1:3))),     rep(0, 4),     rep(0, 4),       rep(0, 4), c(0, -1.5^(-(1:3))),  # L1
        0,       2^(-(1:4)), c(0, 2^(-(1:3))),     rep(0, 4),     rep(0, 4),       rep(0, 4), c(0, -1.5^(-(1:3))),  # L2
        0,        rep(0, 4),        rep(0, 4), c(0, 1, 0, 0),     rep(0, 4),       rep(0, 4),          rep(0, 4),  # L3
        0,        rep(0, 4),        rep(0, 4),     rep(0, 4), c(0, 1, 0, 0),       rep(0, 4),          rep(0, 4),  # L4
        0,        rep(0, 4),        rep(0, 4),     rep(0, 4),     rep(0, 4), c(0, 1.1, 0, 0),         rep(0, 4)), # L5
    nrow = 5,
    ncol = 25,
    byrow = TRUE
  ),
  gamma = NULL,
  delta = NULL,
  #                Int,          L1,         L2,               L3,              L4,         L5,                A
  beta =    c(log(0.1),  2^(-(1:4)), 3^(-(1:4)), -c(0, 0.5, 0, 0),       rep(0, 4),  rep(0, 4),  c(0, 2^(-(1:3)))),
  epsilon = c(log(0.05), 2^(-(1:4)), 3^(-(1:4)),  c(0.5, 0, 0, 0), c(0.5, 0, 0, 0), 1^(-(1:4)),       -3^(-(1:4))),
  
  #           L1   L2   L3   L4   L5
  error =  c(0.2, 0.2,   0,   0, 0.2),
  
  rho = 0.1,
  sigma = diag(rep(1, 5)),

  mc.sims = 100,
  sizes = c(500, 1500, 3000),
  p = 5,
  k = 4,
  
  treatment = NULL,
  shift = NULL,
  specification = c(
    "correct",
    "incorrect covariate",
    "incorrect outcome"
  )
)

factual_sims_df <- run_simulation(scenario1, NSIMS)

factual_sims_df %>%
  select(-train, -test, -truth, -p_function, -model) %>%
  readr::write_rds("../../3_data_encrypted/rds/factual_sims_df.rds")


# scenario 2 --------------------------------------------------------------
# same data generating process as simulation 1; however, now we add competing
# risk (D), which has same determinants as outcome, and we compare estimators
# under risk with and without elimination of competing events.

# change target for risk with elimination

scenario2 <- scenario1

scenario2$delta <- 
  #       Int,         L1,         L2,               L3,        L4,         L5,                A
  c(log(0.05), 2^(-(1:4)), 3^(-(1:4)), -c(0, 0.5, 0, 0), rep(0, 4),  rep(0, 4),  c(0, 2^(-(1:3))))

scenario2$specification <- c(
    "elimination risk", 
    "outcome-specific risk"
    )

compevent_sims_df <- run_simulation(scenario2, NSIMS)

compevent_sims_df %>%
  select(-train, -test, -truth, -p_function, -model) %>%
  readr::write_rds("../../3_data_encrypted/rds/compevent_sims_df.rds")


# scenario 3 --------------------------------------------------------------
    # counterfactual risk under initiation of sustained treatment (A)

scenario3 <- scenario1

scenario3$treatment <- "A"
scenario3$shift <- function(data, treatment, time) { 1 }

scenario3$specification <- c(
  "correct"
)

counterfactual_sims_df <- run_simulation(scenario3, NSIMS)


# counterfactual_sims_df %>%
#   select(-train, -test, -truth, -p_function, -model) %>%
#   readr::write_rds("../../3_data_encrypted/rds/counterfactual_sims_df.rds")
# 



# TODO:
# figure out what's going on with non-pooled version of g-formula (x)
# change target for risk with elimination  (x)
# consider adding gam-version of g-formula (x)
# use gaussian process to generate data

