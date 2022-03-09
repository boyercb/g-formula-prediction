generate_sem_data <-
  function(n = 100,
           k = 2,
           p = 1,
           a = NULL,
           b = NULL,
           c = NULL,
           d = NULL,
           e = NULL,
           treatment = NULL,
           shift = NULL,
           sigma = diag(rep(1, p)),
           error = 0.1,
           long = FALSE) {
    
    C <- matrix(rep(0, n * k), nrow = n, ncol = k)
    L <- matrix(rep(0, n * k * p), nrow = n, ncol = k * p)
    A <- matrix(rep(0, n * k), nrow = n, ncol = k)
    D <- matrix(rep(0, n * k), nrow = n, ncol = k)
    Y <- matrix(rep(0, n * k), nrow = n, ncol = k)
    
    X <- cbind(1, L, A)
    
    if (is.null(dim(a))) {
      a <- t(a)
    }
    
    if (length(error) == 1) {
      error <- rep(error, p)
    }

    for (i in seq(1, k)) {
      if ((i + 1) <= k) {
        h <- c(i:1, (i + 1):k)
      } else {
        h <- c(i:1)
      }
      
      hL <- rep(seq(1, p), each = k) + (h - 1) * p + 1
      hA <- h + k * p + 1
      
      h <- c(1, hL, hA)
      
      iL <- seq(1, p) + (i - 1) * p + 1
      iA <- k * p + i + 1
      
      if (i > 1) {
        for (j in seq(1, p)) {
          if ("L" %in% treatment) {
            X[, iL[j]] <- shift(X[, h], "L", i)
          } else {
            X[, iL[j]] <- rnorm(n, X[, h] %*% a[j, ], error[j])
          }
        }
      } else {
        X[, 1:p + 1] <- rmvnorm(n, rep(0, p), sigma)
      }
      
      if ("A" %in% treatment) {
        X[, iA] <- shift(X[, h], "A", i)
      } else {
        X[, iA] <- rbinom(n, 1, plogis(X[, h] %*% b))
      }

      if (is.numeric(d)) {
        if ("D" %in% treatment) {
          D[, i] <- shift(X[, h], "D", i)
        } else {
          D[, i] <- rbinom(n, 1, plogis(X[, h] %*% d))
        }
      } else {
        D[, i] <- 0
      }
      
      if (is.numeric(c) & i > 1) {
        if ("C" %in% treatment) {
          C[, i] <- shift(X[, h], "C", i)
        } else {
          C[, i] <- rbinom(n, 1, plogis(X[, h] %*% c))
        }
      } else {
        C[, i] <- 0
      }
      
      Y[, i] <- rbinom(n, 1, plogis(X[, h] %*% e))
      
      
    }
    for (i in seq(1, k)) {
      iL <- seq(1, p) + (i - 1) * p + 1
      iA <- k * p + i + 1
      if (i > 1) {
        C[, i] <- ifelse(C[, i - 1], 1, ifelse(Y[, i - 1] | D[, i - 1] == 1, NA, C[, i]))
        for (j in seq(1, p)) {
          X[, iL[j]] <- ifelse(Y[, i - 1] | D[, i - 1] | C[, i], NA, X[, iL[j]])
        }
        X[, iA] <- ifelse(Y[, i - 1] | D[, i - 1] | C[, i], NA, X[, iA])
        D[, i] <- ifelse(D[, i - 1], 1, ifelse(Y[, i - 1] | C[, i], NA, D[, i]))
        Y[, i] <- ifelse(Y[, i - 1], 1, ifelse(D[, i] | C[, i], NA, Y[, i]))
      } else {
        Y[, i] <- ifelse(D[, i], NA, Y)
      }
      
    }
    
    mat <- cbind(C, X[, -1], D, Y)
    
    basenames <- list("C", paste0("L", 1:p, "_"), "A", "D", "Y")
    matnames <- lapply(basenames, function(x) paste0(x, rep(0:(k - 1), each = length(x))))
    colnames(mat) <- unlist(matnames)
    
    ordernames <- lapply(0:(k-1), function (x) paste0(unlist(basenames), as.character(x)))
    mat <- mat[, unlist(ordernames)]
    
    mat <- cbind("id" = 1:n, mat)
    
    if (long == TRUE) {
      mat_long <- cbind(
        rep(1:n, k),
        rep(0:(k-1), each = n),
        as.vector(C),
        sapply(1:p, function(x) as.vector(X[, (1:k - 1) * p + x + 1])),
        as.vector(X[, (ncol(X) - k + 1):ncol(X)]),
        as.vector(D),
        as.vector(Y),
        lapply(1:(k-1), function(l) {
          sapply(1:p, function(j) {
            as.vector(
              cbind(
                matrix(0, nrow = nrow(mat), ncol = l),
                X[, (1:(k-l) - 1) * p + j + 1]
              )
            )
          })
        }) %>% do.call("cbind", .),
        sapply(1:(k-1), function(l) {
          as.vector(
            cbind(
              matrix(0, nrow = nrow(mat), ncol = l), 
              X[, 1:(k-l) + k * p + 1]
            )
          )
        })
      )
      
      colnames(mat_long) <- c(
        "id",
        "time",
        "C",
        paste0("L", 1:p),
        "A",
        "D",
        "Y",
        paste0("lag", rep(1:(k - 1), each = p), "_L", rep(1:p, k-1)),
        paste0("lag", 1:(k - 1), "_A")
      )
      
      mat_long <- mat_long[order(mat_long[,1], mat_long[,2]),]
      
      return(list(
        "wide" = mat,
        "long" = mat_long
      ))
    } else {
      return(mat)
    }
    
  }


# tests
if (FALSE) {
  
  # generate factual data
  generate_sem_data(
    n = 100,
    k = 2,
    p = 1,
    #          Int, L1, lag1_L1,  A, lag1_A
    a = c(       0,  0,       1,  0,      1),
    b = c(       0,  1,     0.5,  1,    0.5),
    c = c(log(0.1),  1,     0.5,  1,    0.5),
    d = c(log(0.1),  1,     0.5,  1,    0.5),
    e = c(log(0.1),  1,     0.5,  1,    0.5)
  )
  
  # generate factual data with 4 time points
  generate_sem_data(
    n = 100,
    k = 4,
    p = 1,
    #          Int, L1, lag1_L1, lag2_L1, lag3_L1, A, lag1_A, lag2_A, lag3_A
    a = c(       0,  0,       1,       0,       0, 0,    0.5,      0,      0),
    b = c(       0,  1,     0.5,       0,       0, 0,      1,      0,      0),
    e = c(log(0.1),  1,     0.5,    0.25,   0.125, 1,    0.5,   0.25,  0.125),
  )
  
  # generate factual data with multiple risk factors
  generate_sem_data(
    n = 100,
    k = 2,
    p = 4,
    a = matrix(
      #          Int, L1, lag1_L1,  L2, lag1_L2,  L3, lag1_L3,  L4, lag1_L4,  A, lag1_A
      c(       0,  0,    0.25,   0,    0.25,   0,    0.25,   0,    0.25,  0,      1,  # L1
               0,  0,    0.25,   0,    0.25,   0,    0.25,   0,    0.25,  0,      1,  # L2
               0,  0,    0.25,   0,    0.25,   0,    0.25,   0,    0.25,  0,      1,  # L3
               0,  0,    0.25,   0,    0.25,   0,    0.25,   0,    0.25,  0,      1), # L4
      nrow = 4,
      ncol = 11,
      byrow = TRUE
    ),
    b = c(       0,  1,     0.5,   1,     0.5,   1,     0.5,   1,     0.5,  1,    0.5),
    e = c(log(0.1),  1,     0.5,   1,     0.5,   1,     0.5,   1,     0.5,  1,    0.5)
  )
  
  # generate factual data with long option
  generate_sem_data(
    n = 100,
    k = 4,
    p = 1,
    #          Int, L1, lag1_L1, lag2_L1, lag3_L1, A, lag1_A, lag2_A, lag3_A
    a = c(       0,  0,       1,       0,       0, 0,    0.5,      0,      0),
    b = c(       0,  1,     0.5,       0,       0, 0,      1,      0,      0),
    e = c(log(0.1),  1,     0.5,    0.25,   0.125, 1,    0.5,   0.25,  0.125),
    long = TRUE
  )
  
  generate_sem_data(
    n = 100,
    k = 2,
    p = 4,
    a = matrix(
      #      Int, L1, lag1_L1,  L2, lag1_L2,  L3, lag1_L3,  L4, lag1_L4,  A, lag1_A
      c(       0,  0,    0.25,   0,    0.25,   0,    0.25,   0,    0.25,  0,      1,  # L1
               0,  0,    0.25,   0,    0.25,   0,    0.25,   0,    0.25,  0,      1,  # L2
               0,  0,    0.25,   0,    0.25,   0,    0.25,   0,    0.25,  0,      1,  # L3
               0,  0,    0.25,   0,    0.25,   0,    0.25,   0,    0.25,  0,      1), # L4
      nrow = 4,
      ncol = 11,
      byrow = TRUE
    ),
    b = c(       0,  1,     0.5,   1,     0.5,   1,     0.5,   1,     0.5,  1,    0.5),
    e = c(log(0.1),  1,     0.5,   1,     0.5,   1,     0.5,   1,     0.5,  1,    0.5),
    long = TRUE
  )
  
  # generate counterfactual data
  generate_sem_data(
    k = 2,
    a = c(0, 0, 1, 0, 1),
    b = c(0, 1, 0.5, 1, 0.5),
    c = c(log(0.1), 1, 0.5, 1, 0.5),
    d = c(log(0.1), 1, 0.5, 1, 0.5),
    e = c(log(0.1), 1, 0.5, 1, 0.5),
    treatment = "A",
    shift = function(data, treatment, time) { 1 }
  )
  
  # generate counterfactual data (random treatment strategy)
  generate_sem_data(
    k = 2,
    a = c(0, 0, 1, 0, 1),
    b = c(0, 1, 0.5, 1, 0.5),
    c = c(log(0.1), 1, 0.5, 1, 0.5),
    d = c(log(0.1), 1, 0.5, 1, 0.5),
    e = c(log(0.1), 1, 0.5, 1, 0.5),
    treatment = "A",
    shift = function(data, treatment, time) { rbinom(nrow(data), 1, 0.5) }
  )
}

calculate_sem_truth <- function(n = 100,
                                k = 2,
                                p = 1,
                                a = NULL,
                                b = NULL,
                                c = NULL,
                                d = NULL,
                                e = NULL,
                                L0 = NULL,
                                A0 = NULL,
                                treatment = NULL,
                                shift = NULL,
                                time = 0,
                                sigma = diag(rep(1, p)),
                                error = 0.1) {
  
  C <- matrix(rep(0, n * k), nrow = n, ncol = k)
  L <- matrix(rep(0, n * k * p), nrow = n, ncol = k * p)
  A <- matrix(rep(0, n * k), nrow = n, ncol = k)
  D <- matrix(rep(0, n * k), nrow = n, ncol = k)
  Y <- matrix(rep(0, n * k), nrow = n, ncol = k)
  
  D_cum <- matrix(rep(0, n * k), nrow = n, ncol = k)
  Y_cum <- matrix(rep(0, n * k), nrow = n, ncol = k)
  
  X <- cbind(1, L, A)
  
  if (is.null(dim(a))) {
    a <- t(a)
  }
  
  if (is.null(dim(L0))) {
    L0 <- matrix(L0, ncol = 1)
  }
  
  if (is.null(dim(A0))) {
    A0 <- matrix(A0, ncol = 1)
  }
  
  for (i in seq(1, k)) {
    if ((i + 1) <= k) {
      h <- c(i:1, (i + 1):k)
    } else {
      h <- c(i:1)
    }
    
    hL <- rep(seq(1, p), each = k) + (h - 1) * p + 1
    hA <- h + k * p + 1
    
    h <- c(1, hL, hA)
    
    iL <- seq(1, p) + (i - 1) * p + 1
    iA <- k * p + i + 1

    if (i > time + 1) {
      for (j in seq(1, p)) {
        if ("L" %in% treatment) {
          X[, iL[j]] <- shift(X[, h], "L", i)
        } else {
          X[, iL[j]] <- X[, h] %*% a[j, ]
        }
      }
    } else {
      for (j in seq(1, p)) {
        X[, iL[j]] <- L0[, iL[j] - 1]
      }
    }
    
    if (i > time + 1) {
      if ("A" %in% treatment) {
        X[, iA] <- shift(X[, h], "A", i)
      } else {
        X[, iA] <- plogis(X[, h] %*% b)
      }
    } else {
      X[, iA] <- A0[, i]
    }
    
    if (i >= time + 1) {
      if (is.numeric(d)) {
        if ("D" %in% treatment) {
          D[, i] <- shift(X[, h], "D", i)
        } else {
          D[, i] <- plogis(X[, h] %*% d)
        }
      } else {
        D[, i] <- 0
      }
      
      if (is.numeric(c) & i > 1) {
        if ("C" %in% treatment) {
          C[, i] <- shift(X[, h], "C", i)
        } else {
          C[, i] <- plogis(X[, h] %*% c)
        }
      } else {
        C[, i] <- 0
      }
      
      Y[, i] <- plogis(X[, h] %*% e)
      if (i > 1) {
        if (is.numeric(d)) {
          D_cum[, i] <- D_cum[, i-1] + D[, i] * (1 - D_cum[, i-1])
          Y_cum[, i] <- Y_cum[, i-1] + Y[, i] * (1 - Y_cum[, i-1]) * (1 - D_cum[, i])
        } else{
          Y_cum[, i] <- Y_cum[, i-1] + Y[, i] * (1 - Y_cum[, i-1])
          D_cum[, i] <- 0
        }
      } else {
        if (is.numeric(d)) {
          D_cum[, i] <- D[, i]
          Y_cum[, i] <- Y[, i] * (1 - D[, i])
        } else{
          Y_cum[, i] <- Y[, i]
          D_cum[, i] <- 0
        }
      }
    } 
  }
  
  colnames(Y_cum) <- paste0("Y_cum", 0:(k-1))
  
  return(Y_cum)
  
}

if (FALSE) {
  # generate factual data
  df <- generate_sem_data(
    n = 5,
    k = 2,
    p = 1,
    #          Int, L1, lag1_L1,  A, lag1_A
    a = c(       0,  0,       1,  0,      1),
    b = c(       0,  1,     0.5,  1,    0.5),
    c = c(log(0.1),  1,     0.5,  1,    0.5),
    d = c(log(0.1),  1,     0.5,  1,    0.5),
    e = c(log(0.1),  1,     0.5,  1,    0.5)
  )
  
  # calculate truth
  calculate_sem_truth(
    n = 5,
    k = 2,
    p = 1,
    #          Int, L1, lag1_L1,  A, lag1_A
    a = c(       0,  0,       1,  0,      1),
    b = c(       0,  1,     0.5,  1,    0.5),
    c = c(log(0.1),  1,     0.5,  1,    0.5),
    d = c(log(0.1),  1,     0.5,  1,    0.5),
    e = c(log(0.1),  1,     0.5,  1,    0.5),
    L0 = df[, 3],
    A0 = df[, 4]
  )

  
  df <- generate_sem_data(
    n = 100,
    k = 2,
    p = 4,
    a = matrix(
      #      Int, L1, lag1_L1,  L2, lag1_L2,  L3, lag1_L3,  L4, lag1_L4,  A, lag1_A
      c(       0,  0,    0.25,   0,    0.25,   0,    0.25,   0,    0.25,  0,      1,  # L1
               0,  0,    0.25,   0,    0.25,   0,    0.25,   0,    0.25,  0,      1,  # L2
               0,  0,    0.25,   0,    0.25,   0,    0.25,   0,    0.25,  0,      1,  # L3
               0,  0,    0.25,   0,    0.25,   0,    0.25,   0,    0.25,  0,      1), # L4
      nrow = 4,
      ncol = 11,
      byrow = TRUE
    ),
    b = c(       0,  1,     0.5,   1,     0.5,   1,     0.5,   1,     0.5,  1,    0.5),
    e = c(log(0.1),  1,     0.5,   1,     0.5,   1,     0.5,   1,     0.5,  1,    0.5)
  )
  
  
  # calculate truth
  calculate_sem_truth(
    n = 100,
    k = 2,
    p = 4,
    a = matrix(
      #      Int, L1, lag1_L1,  L2, lag1_L2,  L3, lag1_L3,  L4, lag1_L4,  A, lag1_A
      c(       0,  0,    0.25,   0,    0.25,   0,    0.25,   0,    0.25,  0,      1,  # L1
               0,  0,    0.25,   0,    0.25,   0,    0.25,   0,    0.25,  0,      1,  # L2
               0,  0,    0.25,   0,    0.25,   0,    0.25,   0,    0.25,  0,      1,  # L3
               0,  0,    0.25,   0,    0.25,   0,    0.25,   0,    0.25,  0,      1), # L4
      nrow = 4,
      ncol = 11,
      byrow = TRUE
    ),
    b = c(       0,  1,     0.5,   1,     0.5,   1,     0.5,   1,     0.5,  1,    0.5),
    e = c(log(0.1),  1,     0.5,   1,     0.5,   1,     0.5,   1,     0.5,  1,    0.5),
    L0 = df[, 3:6],
    A0 = df[, 7]
  )
}


calc_position <- function(variable, matrix) {
  cnames <- colnames(matrix)
  which(cnames %in% variable)
}
