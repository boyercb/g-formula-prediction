library(cmdstanr)
library(posterior)
library(bayesplot)


df <-
  datagen(
    N = 100,
    k = 4,
    alpha = c(0, 2,-2),
    beta = c(log(0.25), log(3), log(3), log(2)),
    gamma = c(log(0.25), log(2), log(2)),
    sigma = 1,
    output = "long"
  )

dforig <- mutate(df, across(everything(), replace_na, 0))
df <- filter(df, !is.na(df$X))

stan_data <- list(
  N = 100,
  K = 4,
  NK = nrow(dforig),
  M = nrow(df),
  y = df$Y,
  w = df$W,
  x = df$X,
  lag1_w = df$lag1_W,
  lag1_x = df$lag1_X,
  worig = dforig$W,
  xorig = dforig$X
)

mod <- cmdstan_model("5_parametric_models/bayes.stan")
fit <- mod$sample(data = stan_data)

post <- as_draws_df(fit$draws())
