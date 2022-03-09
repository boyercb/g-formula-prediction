plan(multisession, workers = 6)

# global simulation parameters --------------------------------------------

N <- 1000
K <- 4
N.SIMS <- 1000
MC.SIMS <- 100

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
  cf_sim1 <- run_simulation(
    n = N,
    k = K,
    #          Int,  X, lag1_X, lag2_X, lag3_X, W, lag1_W, lag2_W, lag3_W
    a = c(       0,  0,      0,      0,      0, 0,      0,      0,      0),
    b = c(       0,  1,      0,      0,      0, 0,      0,      0,      0),
    e = c(log(0.1),  1,      0,      0,      0,-1,      0,      0,      0),
    sigma = 1,
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
    updated.model = UPDATED.MODEL
  )
})



plan(sequential)

