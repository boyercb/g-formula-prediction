offspring_long_std <- offspring_long

# standardize the continuous predictors
offspring_long_std[, c(vars_cont)] <-
  (1 / 2) * scale(offspring_long_std[, c(vars_cont)])

vars_cont_lag <- c(covs_lag1, covs_lag2)
vars_cont_lag <-  vars_cont_lag[str_remove(vars_cont_lag, "lag[12]_") %in% vars_cont]

for (i in seq_along(vars_cont_lag)) {
  offspring_long_std[offspring_long_std[[vars_cont_lag[i]]] > 0, vars_cont_lag[i]] <- 
    (1 / 2) * scale(offspring_long_std[offspring_long_std[[vars_cont_lag[i]]] > 0, vars_cont_lag[i]])
}

