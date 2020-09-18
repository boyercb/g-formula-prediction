datagen <- function(N, k, alpha, beta, gamma, sigma, output = "both") {
  X <- matrix(nrow = N, ncol = k)
  W <- matrix(nrow = N, ncol = k)
  Y <- matrix(nrow = N, ncol = k)
  
  colnames(X) <- paste0("X", 0:(k-1))
  colnames(W) <- paste0("W", 0:(k-1))
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
      rbinom(N, 1, plogis(gamma[1] + gamma[2] * X[, i + 1] + gamma[3] * W[, i + 1]))
    
    if (i > 0) {
      X[, i + 1] <- ifelse(Y[, i] == 1, NA, X[, i + 1])
      W[, i + 1] <- ifelse(Y[, i] == 1, NA, W[, i + 1])
      Y[, i + 1] <- ifelse(Y[, i] == 1, 1, Y[, i + 1])
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



