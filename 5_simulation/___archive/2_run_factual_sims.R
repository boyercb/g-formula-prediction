plan(multisession, workers = 12)

# global simulation parameters --------------------------------------------

N <- 1000
K <- 4
N.SIMS <- 1000
MC.SIMS <- 50

BASE.MODEL <- list(
  formula = Y3 ~ X0 + W0,
  link = "logit",
  family = "binomial"
)

UPDATED.MODEL <- list(
  "base only" = list(
    formula = Y3 ~ X0 + W0,
    link = "logit",
    family = "binomial"
  ),
  "time 0 to 1" = list(
    formula = Y3 ~ X0 + W0 + X1 + W1,
    link = "logit",
    family = "binomial"
  ),
  "time 0 to 2" = list(
    formula = Y3 ~ X0 + W0 + X1 + W1 + X2 + W2,
    link = "logit",
    family = "binomial"
  ),
  "time 0 to 3" = list(
    formula = Y3 ~ X0 + W0 + X1 + W1 + X2 + W2 + X3 + W3,
    link = "logit",
    family = "binomial"
  )
)


# simulation 1: -----------------------------------------------------------

# simple markov model with no memory
with_progress({
  sim1 <- run_simulation(
    n = N,
    k = K,
    #          Int,  X, lag1_X, lag2_X, lag3_X, W, lag1_W, lag2_W, lag3_W
    a = c(       0,  0,      0,      0,      0, 0,      0,      0,      0),
    b = c(       0,  1,      0,      0,      0, 0,      0,      0,      0),
    e = c(log(0.1),  1,      0,      0,      0,-1,      0,      0,      0),
    sigma = 0.75,
    s0 = 1,
    n.sims = N.SIMS,
    mc.sims = MC.SIMS,
    X.models = list(
      "X" = list(formula = X ~ 1,
                 family = "normal"),
      "W" = list(
        formula = W ~ X,
        link = "logit",
        family = "binomial"
      )
    ),
    Y.model = list(
      formula = Y ~ X + W,
      link = "logit",
      family = "binomial"
    ),
    base.model = BASE.MODEL,
    updated.model = list(
      "base only" = list(
        formula = Y3 ~ X0 + W0,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 1" = list(
        formula = Y3 ~ X1 + W1,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 2" = list(
        formula = Y3 ~ X2 + W2,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 3" = list(
        formula = Y3 ~ X3 + W3,
        link = "logit",
        family = "binomial"
      )
    )
  )
})


# simulation 2: -----------------------------------------------------------

# covariates are determined by previous values, but outcome is still determined
# by most recent values only
with_progress({
  sim2 <- run_simulation(
    n = N,
    k = K,
    #          Int,  X, lag1_X, lag2_X, lag3_X, W, lag1_W, lag2_W, lag3_W
    a = c(       0,  0,      1,      0,      0, 0,   -0.5,      0,      0),
    b = c(       0,  1,    0.5,      0,      0, 0,      1,      0,      0),
    e = c(log(0.1),  1,      0,      0,      0,-1,      0,      0,      0),
    sigma = 0.75,
    s0 = 1,
    n.sims = N.SIMS,
    mc.sims = MC.SIMS,
    X.models = list(
      "X" = list(formula = X ~ lag1_X + lag1_W,
                 family = "normal"),
      "W" = list(
        formula = W ~ X + lag1_X + lag1_W,
        link = "logit",
        family = "binomial"
      )
    ),
    Y.model = list(
      formula = Y ~ X + W,
      link = "logit",
      family = "binomial"
    ),
    base.model = BASE.MODEL,
    updated.model = list(
      "base only" = list(
        formula = Y3 ~ X0 + W0,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 1" = list(
        formula = Y3 ~ X1 + W1,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 2" = list(
        formula = Y3 ~ X2 + W2,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 3" = list(
        formula = Y3 ~ X3 + W3,
        link = "logit",
        family = "binomial"
      )
    )
  )
})


# simulation 3: -----------------------------------------------------------

# outcome is now a function of last two values, equally
with_progress({
  sim3 <- run_simulation(
    n = N,
    k = K,
    #          Int,  X, lag1_X, lag2_X, lag3_X,   W, lag1_W, lag2_W, lag3_W
    a = c(       0,  0,      1,      0,      0,   0,   -0.5,      0,      0),
    b = c(       0,  1,    0.5,      0,      0,   0,      1,      0,      0),
    e = c(log(0.1),0.5,    0.5,      0,      0,-0.5,   -0.5,      0,      0),
    sigma = 0.75,
    s0 = 1,
    n.sims = N.SIMS,
    mc.sims = MC.SIMS,
    X.models = list(
      "X" = list(formula = X ~ lag1_X + lag1_W,
                 family = "normal"),
      "W" = list(
        formula = W ~ X + lag1_X + lag1_W,
        link = "logit",
        family = "binomial"
      )
    ),
    Y.model = list(
      formula = Y ~ X + W + lag1_X + lag1_W,
      link = "logit",
      family = "binomial"
    ),
    base.model = BASE.MODEL,
    updated.model = list(
      "base only" = list(
        formula = Y3 ~ X0 + W0,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 1" = list(
        formula = Y3 ~ X0 + W0 + X1 + W1,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 2" = list(
        formula = Y3 ~ X1 + W1 + X2 + W2,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 3" = list(
        formula = Y3 ~ X2 + W2 + X3 + W3,
        link = "logit",
        family = "binomial"
      )
    )
  )
})


# simulation 4: -----------------------------------------------------------

# outcome is now a function of last three values, equally
with_progress({
  sim4 <- run_simulation(
    n = N,
    k = K,
    #          Int,  X, lag1_X, lag2_X, lag3_X,    W, lag1_W, lag2_W, lag3_W
    a = c(       0,  0,      1,      0,      0,    0,   -0.5,      0,      0),
    b = c(       0,  1,    0.5,      0,      0,    0,      1,      0,      0),
    e = c(log(0.1),0.25,  0.25,   0.25,      0,-0.25,  -0.25,  -0.25,      0),
    sigma = 0.75,
    s0 = 1,
    n.sims = N.SIMS,
    mc.sims = MC.SIMS,
    X.models = list(
      "X" = list(formula = X ~ lag1_X + lag1_W,
                 family = "normal"),
      "W" = list(
        formula = W ~ X + lag1_X + lag1_W,
        link = "logit",
        family = "binomial"
      )
    ),
    Y.model = list(
      formula = Y ~ X + W + lag1_X + lag1_W + lag2_X + lag2_W,
      link = "logit",
      family = "binomial"
    ),
    base.model = BASE.MODEL,
    updated.model = list(
      "base only" = list(
        formula = Y3 ~ X0 + W0,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 1" = list(
        formula = Y3 ~ X0 + W0 + X1 + W1,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 2" = list(
        formula = Y3 ~ X0 + W0 + X1 + W1 + X2 + W2,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 3" = list(
        formula = Y3 ~ X1 + W1 + X2 + W2 + X3 + W3,
        link = "logit",
        family = "binomial"
      )
    )
  )
})


# simulation 5: -----------------------------------------------------------

# outcome is now a function of last four values, equally
with_progress({
  sim5 <- run_simulation(
    n = N,
    k = K,
    #          Int,  X, lag1_X, lag2_X, lag3_X,    W, lag1_W, lag2_W, lag3_W
    a = c(       0,  0,      1,      0,      0,    0,   -0.5,      0,      0),
    b = c(       0,  1,    0.5,      0,      0,    0,      1,      0,      0),
    e = c(log(0.1),0.25,  0.25,   0.25,   0.25,-0.25,  -0.25,  -0.25,  -0.25),
    sigma = 0.75,
    s0 = 1,
    n.sims = N.SIMS,
    mc.sims = MC.SIMS,
    X.models = list(
      "X" = list(formula = X ~ lag1_X + lag1_W,
                 family = "normal"),
      "W" = list(
        formula = W ~ X + lag1_X + lag1_W,
        link = "logit",
        family = "binomial"
      )
    ),
    Y.model = list(
      formula = Y ~ X + W + lag1_X + lag1_W + lag2_X + lag2_W + lag3_X + lag3_W,
      link = "logit",
      family = "binomial"
    ),
    base.model = BASE.MODEL,
    updated.model = list(
      "base only" = list(
        formula = Y3 ~ X0 + W0,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 1" = list(
        formula = Y3 ~ X0 + W0 + X1 + W1,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 2" = list(
        formula = Y3 ~ X0 + W0 + X1 + W1 + X2 + W2,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 3" = list(
        formula = Y3 ~ X0 + W0 + X1 + W1 + X2 + W2 + X3 + W3,
        link = "logit",
        family = "binomial"
      )
    )
  )
})


# simulation 6: -----------------------------------------------------------

# outcome is now a decaying function of last four values
with_progress({
  sim6 <- run_simulation(
    n = N,
    k = K,
    #          Int,  X, lag1_X, lag2_X, lag3_X,    W, lag1_W, lag2_W, lag3_W
    a = c(       0,  0,      1,      0,      0,    0,   -0.5,      0,      0),
    b = c(       0,  1,    0.5,      0,      0,    0,      1,      0,      0),
    e = c(log(0.1),  1,    0.5,   0.25,   0.125,  -1,  -0.75,  -0.50,  -0.25),
    sigma = 0.75,
    s0 = 1,
    n.sims = N.SIMS,
    mc.sims = MC.SIMS,
    X.models = list(
      "X" = list(formula = X ~ lag1_X + lag1_W,
                 family = "normal"),
      "W" = list(
        formula = W ~ X + lag1_X + lag1_W,
        link = "logit",
        family = "binomial"
      )
    ),
    Y.model = list(
      formula = Y ~ X + W + lag1_X + lag1_W + lag2_X + lag2_W + lag3_X + lag3_W,
      link = "logit",
      family = "binomial"
    ),
    base.model = BASE.MODEL,
    updated.model = list(
      "base only" = list(
        formula = Y3 ~ X0 + W0,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 1" = list(
        formula = Y3 ~ X0 + W0 + X1 + W1,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 2" = list(
        formula = Y3 ~ X0 + W0 + X1 + W1 + X2 + W2,
        link = "logit",
        family = "binomial"
      ),
      "time 0 to 3" = list(
        formula = Y3 ~ X0 + W0 + X1 + W1 + X2 + W2 + X3 + W3,
        link = "logit",
        family = "binomial"
      )
    )
  )
})

plan(sequential)


# save results ------------------------------------------------------------

factual_sim_results <- bind_rows(sim1,
                                 sim2,
                                 sim3,
                                 sim4,
                                 sim5,
                                 sim6,
                                 .id = "scenario") 

factual_sim_results <- factual_sim_results %>%
  mutate(
    cond = replace_na(cond, 1),
    sim = as.numeric(sim)
  )

write_rds(factual_sim_results, "../../3_data_encrypted/rds/factual_sim_results.rds")
