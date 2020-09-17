datagen <- function(N, k, alpha, beta, gamma, sigma, output = "both") {
  X <- matrix(nrow = N, ncol = k)
  W <- matrix(nrow = N, ncol = k)
  Y <- matrix(nrow = N, ncol = k)
  
  colnames(X) <- paste0("X", 1:k)
  colnames(W) <- paste0("W", 1:k)
  colnames(Y) <- paste0("Y", 1:k)
  
  for (i in 1:k) {
    if (i == 1) {
      lag1_X <- 0
      lag1_W <- 0
    } else {
      lag1_X <- ifelse(is.na(X[, i - 1]), 0, X[, i - 1])
      lag1_W <- ifelse(is.na(W[, i - 1]), 0, W[, i - 1])
    }
    
    X[, i] <-
      rnorm(N, alpha[1] + alpha[2] * lag1_X + alpha[3] * lag1_W, sigma)
    W[, i] <-
      rbinom(N, 1, plogis(beta[1] + beta[2] * X[, i] + beta[3] * lag1_X + beta[4] * lag1_W))
    Y[, i] <-
      rbinom(N, 1, plogis(gamma[1] + gamma[2] * X[, i] + gamma[3] * W[, i]))
    
    if (i > 1) {
      X[, i] <- ifelse(Y[, i - 1] == 1, NA, X[, i])
      W[, i] <- ifelse(Y[, i - 1] == 1, NA, W[, i])
      Y[, i] <- ifelse(Y[, i - 1] == 1, 1, Y[, i])
    }
  }
  
  build_df <- function(N, X, W, Y, output) {
    switch(
      output,
      "long" = data.frame(
        id = rep(1:N, k),
        time = rep(1:k, each = N),
        X = as.vector(X),
        W = as.vector(W),
        Y = as.vector(Y),
        lag1_X = c(rep(0, N), as.vector(X[, -k])),
        lag1_W = c(rep(0, N), as.vector(W[, -k]))
      ),
      "wide" = data.frame(
        id = 1:N,
        X, 
        W, 
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



