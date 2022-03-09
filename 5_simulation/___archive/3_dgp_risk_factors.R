datagen <- function(N, k, alpha, beta, gamma, delta, eta, sigma, output = "both") {
  X <- matrix(nrow = N, ncol = k)
  W <- matrix(nrow = N, ncol = k)
  C <- matrix(nrow = N, ncol = k)
  D <- matrix(nrow = N, ncol = k)
  Y <- matrix(nrow = N, ncol = k)
  
  colnames(X) <- paste0("X", 0:(k-1))
  colnames(W) <- paste0("W", 0:(k-1))
  colnames(C) <- paste0("C", 0:(k-1))
  colnames(D) <- paste0("D", 0:(k-1))
  colnames(Y) <- paste0("Y", 0:(k-1))
  
  for (i in 0:(k-1)) {
    if (i == 0) {
      lag1_X <- 0
      lag1_W <- 0
    } else {
      lag1_X <- ifelse(is.na(X[, i]), 0, X[, i])
      lag1_W <- ifelse(is.na(W[, i]), 0, W[, i])
    }
    
    X[, i + 1] <-
      rnorm(N, alpha[1] + alpha[2] * lag1_X + alpha[3] * lag1_W, sigma)
    W[, i + 1] <-
      rbinom(N, 1, plogis(beta[1] + beta[2] * X[, i + 1] + beta[3] * lag1_X + beta[4] * lag1_W))
    Y[, i + 1] <-
      rbinom(N, 1, plogis(eta[1] + eta[2] * X[, i + 1] + eta[3] * W[, i + 1] + eta[4] * lag1_X + eta[5] * lag1_W))
    
    if (!is.null(gamma) & i > 0) {
      C[, i + 1] <-
        rbinom(N, 1, plogis(gamma[1] + gamma[2] * X[, i + 1] + gamma[3] * W[, i + 1]))
    } else {
      C[, i + 1] <- 0
    }
    if (!is.null(delta)) {
      D[, i + 1] <- 
        rbinom(N, 1, plogis(delta[1] + delta[2] * X[, i + 1] + delta[3] * W[, i + 1]))
    } else {
      D[, i + 1] <- 0
    }
    
    if (i > 0) {
      X[, i + 1] <- ifelse(Y[, i] == 1 | D[, i] == 1 | C[, i + 1] == 1, NA, X[, i + 1])
      W[, i + 1] <- ifelse(Y[, i] == 1 | D[, i] == 1 | C[, i + 1] == 1, NA, W[, i + 1])
      C[, i + 1] <- ifelse(C[, i] == 1, 1, ifelse(Y[, i] == 1 | D[, i] == 1, NA, C[, i + 1]))
      D[, i + 1] <- ifelse(D[, i] == 1, 1, ifelse(Y[, i] == 1, 0, ifelse(C[, i + 1] == 1, NA, D[, i + 1])))
      Y[, i + 1] <- ifelse(Y[, i] == 1, 1, ifelse(D[, i + 1] == 1, 0, ifelse(C[, i + 1] == 1, NA, Y[, i + 1])))
    }
  }
  
  build_df <- function(N, X, W, Y, output) {
    switch(
      output,
      "long" = data.frame(
        id = rep(1:N, k),
        time = rep(0:(k-1), each = N),
        X = as.vector(X),
        W = as.vector(W),
        C = as.vector(C),
        D = as.vector(D),
        Y = as.vector(Y),
        lag1_X = c(rep(0, N), as.vector(X[, -k])),
        lag1_W = c(rep(0, N), as.vector(W[, -k]))
      ),
      "wide" = data.frame(
        id = 1:N,
        X, 
        W,
        C,
        D,
        Y
      )
    )
  }
  
  df <- 
    switch(
      output,
      "long" = build_df(N, X, W, Y, output = "long"),
      "wide" = build_df(N, X, W, Y, output = "wide"),
      "both" = list(
        "long" = build_df(N, X, W, Y, output = "long"),
        "wide" = build_df(N, X, W, Y, output = "wide")
      )
    )
  
  return(df)
}


df <- datagen(10000,      
              k = 5,
              alpha = c(0, 1.2, -1.5),
              beta = c(log(0.1), log(1.5), log(1.2), log(1.2)),
              gamma = NULL, #c(log(0.1), log(1.2), log(1.2)),
              delta = c(log(0.03), log(1.2), log(0.85)),
              eta = c(log(0.05), log(1.2), log(0.5), log(1.1), log(0.85)),
              sigma = 1)

datagen(5000,      
        k = 4,
        alpha = c(0, 1.2, -1.5),
        beta = c(log(0.25), log(3), log(3), log(2)),
        gamma = NULL,
        delta = NULL,
        eta = c(log(0.25), log(2), log(2)),
        sigma = 1)

# fit models for g-formula
Y.fit_g <- 
  fit_Y_model(
    model = list(
      formula = Y ~ X + W + lag1_X + lag1_W,
      link = "logit",
      family = "binomial"
    ),
    data = filter(df$long, D != 1)
  )

X.fit_g <-
  fit_X_models(
    models = list(
      "X" = list(
        formula = X ~ lag1_X + lag1_W, 
        family = "normal"
      ),
      "W" = list(
        formula = W ~ X + lag1_X + lag1_W, 
        link = "logit",
        family = "binomial"
      )
    ),
    data = df$long,
    time = "time"
  )

D.fit_g <-
  fit_Y_model(
    model = list(
      formula = D ~ X + W,
      link = "logit",
      family = "binomial"
    ),
    data = df$long
  )

Y.hat_g <-
  gformula_mc(
    Y.fit_g,
    X.fit_g,
    D.fit_g,
    data = filter(df$long, !is.na(X)),
    id = "id",
    time = "time",
    mc.start = 0,
    mc.stop = 4,
    hist.vars = c("lag1_X", "lag1_W"),
    hist.fun = "lag",
    mc.sims = 100,
    merge = FALSE
  )

t <- left_join(select(df$wide, Y4, id), filter(Y.hat_g, time == 4))

rms::val.prob(t$pred, t$Y4)


Y.fit_c_base <-
  glm(
    formula = reformulate(paste0(c("X", "W"), 0), paste0("Y", 4)),
    data = df$wide,
    family = binomial(link = "logit")
  )

Y.fit_c_base <-
  predict(Y.fit_c_base, type = 'response', newdata = df$wide)

grf::generate_causal_data()

