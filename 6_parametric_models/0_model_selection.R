# smoking
fit <- glm(
  formula = smk ~ .,
  data = analytic_long_std[, c(covs_fixed, covs_tv[1], covs_lag1, covs_lag2, "time", "smk")],
  family = binomial(link = "logit")
)

test <- MASS::stepAIC(
  fit,
  scope = list(
    lower = ~ sex + lag1_smk + lag2_smk
  ),
  direction = "backward"
)

lrm(smk ~ sex + rcs(age, 3) + lag1_smk * lag2_smk + lag1_hrx + lag1_liprx + lag1_tc + 
      lag1_sbp + lag2_hdl + catg(time), data = analytic_long)



ga