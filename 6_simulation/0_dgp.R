# define data generation process ------------------------------------------

datagen <- function(N, 
                    a0 = -1.5,
                    a1 = 1,
                    a2 = 1, 
                    a3 = 2,
                    a4 = 0.5,
                    # b0 = -3,
                    # b1 = -1,
                    # b2 = 0.75,
                    # b3 = 0.2,
                    c0 = -2,
                    c1 = -2,
                    c2 = 3,
                    c3 = 0.2,
                    d0 = 0,
                    d1 = 1,
                    d2 = 1) {
  
  # time point 1
  L0 <- rbinom(N, 1, plogis(d0))
  A0 <- rbinom(N, 1, plogis(a0 + a1 * L0))
  #D0 <- rbinom(N, 1, plogis(b0 + b1 * A1 + b2 * L1))
  Y0 <- rbinom(N, 1, plogis(c0 + c1 * A0 + c2 * L0)) 
  
  #Y0 <- ifelse(D0 == 1, NA, Y0)
  
  # time point 2
  L1 <- rbinom(N, 1, plogis(d0 + d1 * L0 + d2 * A0))
  A1 <- rbinom(N, 1, plogis(a0 + a1 * L1 + a2 * A0 + a3 * L0 + a4))
  #D1 <- rbinom(N, 1, plogis(b0 + b1 * A2 + b2 * L2 + b3))
  Y1 <- rbinom(N, 1, plogis(c0 + c1 * A1 + c2 * L1 + c3))
  
  # adjust for survival
  L1 <- ifelse(Y0 == 1 , NA, L1)
  A1 <- ifelse(Y0 == 1 , NA, A1)
  #D1 <- ifelse(Y0 == 1 | D0 == 1, NA, D2)
  Y1 <- ifelse(Y0 == 1, 1, Y1)
  #Y2 <- ifelse(D2 == 1, NA, Y2)
  
  # export simulated data
  data.frame(
      id = 1:N,
      L0 = L0, 
      A0 = A0,
      #D0 = D0,
      Y0 = Y0,
      L1 = L1,
      A1 = A1,
      #D1 = D1,
      Y1 = Y1
    )
}


