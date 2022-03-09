# define data generation process ------------------------------------------

datagen <- function(N, 
                    a0 = 0,
                    a1 = 1,
                    a2 = -3, 
                    a3 = 0,
                    a4 = 0,
                    b0 = log(0.25),
                    b1 = log(2),
                    b2 = 0,
                    b3 = log(2),
                    c0 = log(0.25),
                    c1 = log(1.5),
                    c2 = log(0.5),
                    c3 = log(1.5),
                    c4 = log(0.5), 
                    c5 = 0,
                    sigma = 1, 
                    output = "long") {
  
  # time point 1
  L0 <- rnorm(N, a0, sigma)
  A0 <- rbinom(N, 1, plogis(b0 + b1 * L0))
  Y0 <- rbinom(N, 1, plogis(c0 + c1 * L0 + c2 * A0)) 
  
  # time point 2
  L1 <- rnorm(N, a0 + a1 * L0 + a2 * A0 + a3 * L0 * A0, sigma)
  if (b3 == Inf) {
    A1 <- A0
  } else {
    A1 <- rbinom(N, 1, plogis(b0 + b1 * L1 + b2 * L0 + b3 * A0))
  }
  Y1 <- rbinom(N, 1, plogis(c0 + c1 * L1 + c2 * A1 + c3 * L0 + c4 * A0 + c5 * A0 * A1))
  
  # adjust for survival
  L1 <- ifelse(Y0 == 1 , NA, L1)
  A1 <- ifelse(Y0 == 1 , NA, A1)
  Y1 <- ifelse(Y0 == 1, 1, Y1)

  build_df <- function(N, L0, A0, Y0, L1, A1, Y1, output) {
    switch(
      output,
      "long" = data.frame(
        id = rep(1:N, 2),
        time = rep(c(0, 1), each = N),
        L = c(L0, L1),
        A = c(A0, A1),
        Y = c(Y0, Y1),
        lag1_L = c(rep(0, N), L0),
        lag1_A = c(rep(0, N), A0)
      ),
      "wide" = data.frame(
        id = 1:N,
        L0 = L0,
        A0 = A0,
        Y0 = Y0,
        L1 = L1,
        A1 = A1,
        Y1 = Y1
      )
    )
  }
  
  # export simulated data in either wide or long format
  df <- 
    switch(
      output,
      "long" = build_df(N, L0, A0, Y0, L1, A1, Y1, output = "long"),
      "wide" = build_df(N, L0, A0, Y0, L1, A1, Y1, output = "wide"),
      "both" = list(
        "long" = build_df(N, L0, A0, Y0, L1, A1, Y1, output = "long"),
        "wide" = build_df(N, L0, A0, Y0, L1, A1, Y1, output = "wide")
      )
    )
  
  return(df)
}



