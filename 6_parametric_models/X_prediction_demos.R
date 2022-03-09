# model specifications for demonstrations ---------------------------------

ascvd_interactions <- c(
  "age:tc",
  "age:hdl",
  "age:sbp",
  "age:hrx",
  "age:smk",
  "age:hrx:sbp",
  "hrx:sbp",
  "liprx:hdl",
  "liprx:tc",
  "age:liprx:hdl",
  "age:liprx:tc",
  "lag1_age:lag1_tc",
  "lag1_age:lag1_hdl",
  "lag1_age:lag1_sbp",
  "lag1_age:lag1_hrx",
  "lag1_age:lag1_smk",
  "lag1_age:lag1_hrx:lag1_sbp",
  "lag1_hrx:lag1_sbp",
  "lag1_liprx:lag1_hdl",
  "lag1_liprx:lag1_tc",
  "lag1_age:lag1_liprx:lag1_hdl",
  "lag1_age:lag1_liprx:lag1_tc",
  "lag2_age:lag2_tc",
  "lag2_age:lag2_hdl",
  "lag2_age:lag2_sbp",
  "lag2_age:lag2_hrx",
  "lag2_age:lag2_smk",
  "lag2_age:lag2_hrx:lag2_sbp",
  "lag2_hrx:lag2_sbp",
  "lag2_liprx:lag2_hdl",
  "lag2_liprx:lag2_tc",
  "lag2_age:lag2_liprx:lag2_hdl",
  "lag2_age:lag2_liprx:lag2_tc"
)

# initial experiment value of added lags ----------------------------------

covs_tv_ex1 <- covs_tv[!grepl("(bmi)|(ldl)|(liprx)", covs_tv)]
covs_lag1_ex1 <- covs_lag1[!grepl("(bmi)|(ldl)|(liprx)", covs_lag1)]
covs_lag2_ex1 <- covs_lag2[!grepl("(bmi)|(ldl)|(liprx)", covs_lag2)]
ascvd_interactions_ex1 <- ascvd_interactions[!grepl("(bmi)|(ldl)|(liprx)", ascvd_interactions)]

covs_tv_ex2 <- list(
  covs_tv[!grepl("(ldl)|(liprx)", covs_tv)],
  covs_tv[!grepl("(tc)|(liprx)", covs_tv)],
  covs_tv[!grepl("(tc)", covs_tv)]
)
covs_lag1_ex2 <- list(
  covs_lag1[!grepl("(ldl)|(liprx)", covs_lag1)],
  covs_lag1[!grepl("(tc)|(liprx)", covs_lag1)],
  covs_lag1[!grepl("(tc)", covs_lag1)]
)
covs_lag2_ex2 <- list(
  covs_lag2[!grepl("(ldl)|(liprx)", covs_lag2)],
  covs_lag2[!grepl("(tc)|(liprx)", covs_lag2)],
  covs_lag2[!grepl("(tc)", covs_lag2)]
)
ascvd_interactions_ex2 <- list(
  ascvd_interactions[!grepl("(ldl)|(liprx)", ascvd_interactions)],
  ascvd_interactions[!grepl("(tc)|(liprx)", ascvd_interactions)],
  ascvd_interactions[!grepl("(tc)", ascvd_interactions)]
)


# use single lag specification with main effect terms for now
lag_2_specification <- list(
  "age"   = c("lag1_age"),
  "smk"   = c(covs_fixed, "poly(age, 2)", covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
  "bmi"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:2], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
  "dm"    = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:3], covs_lag1_ex2[[3]][!grepl("dm", covs_lag1_ex2[[3]])], covs_lag2_ex2[[3]][!grepl("dm", covs_lag2_ex2[[3]])]),
  "hrx"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:4], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
  "liprx" = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:5], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
  "ldl"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:6], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
  "hdl"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:7], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
  "sbp"   = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:8], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]]),
  "event_ascvd" = c(
    covs_fixed,
    "poly(age, 2)",
    covs_tv_ex2[[3]][2:9],
    covs_lag1_ex2[[3]],
    covs_lag2_ex2[[3]],
    "lag1_age",
    "lag2_age",
    ascvd_interactions_ex2[[3]]
  ),
  "event_dth_ascvd" = c(covs_fixed, "poly(age, 2)", covs_tv_ex2[[3]][2:9], covs_lag1_ex2[[3]], covs_lag2_ex2[[3]])
)

# define models
models_lag2 <- define_models(
  lag_2_specification,
  "event_ascvd",
  "event_dth_ascvd",
  time = TRUE,
  logvars = FALSE,
  sexints = TRUE
)

# fit models
Y.fit_lag2 <- 
  fit_Y_model(
    models_lag2$Y,
    data = filter(offspring_long, !pid %in% censored_ids)
  )

X.fit_lag2 <- 
  fit_X_models(
    models_lag2$X,
    data = filter(offspring_long, !pid %in% censored_ids),
    time = "time"
  )

D.fit_lag2 <- 
  fit_D_model(
    models_lag2$D,
    data = filter(offspring_long, !pid %in% censored_ids)
  )



# helper function ---------------------------------------------------------


sim_intervention <- function(data, start, stop, intervention, treatment) {
  sim <- gformula_mc(
    Y.fit = Y.fit_lag2,
    X.fit = X.fit_lag2,
    D.fit = D.fit_lag2,
    data = data,
    id = "pid",
    time = "time",
    treatment = treatment,
    intervention = intervention,
    base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
    hist.vars = c(paste0("lag1_", covs_tv[!covs_tv %in% c("tc")]), paste0("lag2_", covs_tv[!covs_tv %in% c("tc")])),
    hist.fun = "lag",
    mc.sims = 500,
    mc.start = start,
    mc.stop = stop,
    bound.sims = FALSE,
    merge = FALSE,
    pool = TRUE
  )
  
}

sim_risk_factors <- function(data, start, stop) {
  sim <- gformula_mc(
    Y.fit = Y.fit_lag2,
    X.fit = X.fit_lag2,
    D.fit = D.fit_lag2,
    data = data,
    id = "pid",
    time = "time",
    base.covs = covs_fixed[!covs_fixed %in% c(covs_refs)],
    hist.vars = c(paste0("lag1_", covs_tv[!covs_tv %in% c("tc")]), paste0("lag2_", covs_tv[!covs_tv %in% c("tc")])),
    hist.fun = "lag",
    mc.sims = 500,
    mc.start = start,
    mc.stop = stop,
    bound.sims = FALSE,
    merge = FALSE,
    pool = TRUE
  )
  
}

plot_risk_factors <- function(data, add_points) {
  data %>%
    as_tibble() %>%
    arrange(pid, time, sim) %>%
    group_by(pid, time) %>%
    select(pid, time, covs_tv) %>%
    pivot_longer(covs_tv) %>%
    group_by(pid, name, time) %>%
    point_interval(
      value,
      .width = c(0.5, 0.8, 0.95),
      .point = mean,
      .interval = qi
    ) %>%
    mutate(
      .lower = ifelse(
        name %in% vars_binary,
        value + qnorm((1 - .width) / 2) * value * (1 - value),
        .lower
      ),
      .upper = ifelse(
        name %in% vars_binary,
        value + qnorm(1 - (1 - .width) / 2) * value * (1 - value),
        .upper
      ),
    ) %>%
    ungroup() %>%
    mutate(name = factor(
      name,
      levels = covs_tv,
      labels = c(
        expression("Age"),
        expression("Current~smoker~(0/1)"),
        expression("Body~mass~index~(kg/m^2)"),
        expression("Diabetes~(0/1)"),
        expression("Blood~pressure~medication~(0/1)"),
        expression("Lipid~lowering~medication~(0/1)"),
        expression("Total~cholesterol~(mg/dL)"),
        expression("HDL~cholesterol~(mg/dL)"),
        expression("Systolic~blood~pressure~(mmHg)")
      )
    )) %>% 
    ggplot(., aes(x = (time) * 4, y = value)) +
    facet_wrap(~name, scales = "free_y", labeller = label_parsed) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
    geom_point(
      data = add_points %>% 
        as_tibble() %>%
        arrange(pid, time) %>%
        group_by(pid, time) %>%
        select(pid, time, covs_tv) %>%
        pivot_longer(covs_tv) %>%
        ungroup() %>%
        mutate(name = factor(
          name,
          levels = covs_tv,
          labels = c(
            expression("Age"),
            expression("Current~smoker~(0/1)"),
            expression("Body~mass~index~(kg/m^2)"),
            expression("Diabetes~(0/1)"),
            expression("Blood~pressure~medication~(0/1)"),
            expression("Lipid~lowering~medication~(0/1)"),
            expression("Total~cholesterol~(mg/dL)"),
            expression("HDL~cholesterol~(mg/dL)"),
            expression("Systolic~blood~pressure~(mmHg)")
          )
        )
      )
    ) +
    scale_x_continuous(breaks = seq(0, 24, 4)) +
    scale_fill_brewer(name = "Simulation interval", type = "seq", palette = "Blues") +
    labs(
      x = "Years",
      y = "Predicted value"
    ) +
    slides_theme()
}

plot_ascvd <- function(data, add_line = NULL, color = "Reds") {
  df <- data %>%
    as_tibble() %>%
    arrange(pid, time, sim) %>%
    select(pid, time, poprisk) %>%
    pivot_longer(poprisk) %>%
    group_by(pid, name, time) %>%
    point_interval(
      value,
      .width = c(0.5, 0.8, 0.95),
      .point = mean,
      .interval = qi
    ) 
  
  p <- ggplot(df, aes(x = (time + 1) * 4, y = value, ymin = .lower, ymax = .upper)) +
    facet_wrap(~factor(name, label = "ASCVD"), scales = "free_y") +
    geom_lineribbon()
  
  if (!is.null(add_line)) {
    p <- p + geom_line(
      data = add_line %>%
        as_tibble() %>%
        arrange(pid, time, sim) %>%
        select(pid, time, poprisk) %>%
        pivot_longer(poprisk) %>%
        group_by(pid, name, time) %>%
        point_interval(
          value,
          .width = c(0.5, 0.8, 0.95),
          .point = mean,
          .interval = qi
        ),
      linetype = "longdash",
      color = "grey75",
      size = 1.25
    )
  }
  
  p + scale_x_continuous(breaks = seq(0, 24, 4)) +
    scale_y_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5)) +
    scale_fill_brewer(name = "Simulation interval",
                      type = "seq",
                      palette = color) +
    labs(x = "Years",
         y = "Predicted risk") +
    slides_theme()
}

# example 1: --------------------------------------------------------------

# a 60 year-old male smoker at baseline with BMI 30, no history of diabetes, no
# treatment history, and moderate risk factor levels

example1_df <- tibble(
  pid = 1,
  time = 0,
  sex = 1,
  age = 60,
  smk = 1,
  bmi = 30,
  dm = 0,
  hrx = 0,
  liprx = 0,
  ldl = 160,
  hdl = 40,
  tc = 200,
  sbp = 145,
  lag1_age = 0,
  lag1_smk = 0,
  lag1_bmi = 0,
  lag1_dm = 0,
  lag1_hrx = 0,
  lag1_liprx = 0,
  lag1_ldl = 0,
  lag1_tc = 0,
  lag1_hdl = 0,
  lag1_sbp = 0,
  lag2_age = 0,
  lag2_smk = 0,
  lag2_bmi = 0,
  lag2_dm = 0,
  lag2_hrx = 0,
  lag2_liprx = 0,
  lag2_ldl = 0,
  lag2_hdl = 0,
  lag2_tc = 0,
  lag2_sbp = 0,
  event_ascvd = 0,
  event_dth_ascvd = 0
)

example1_sims <- sim_risk_factors(example1_df, 0, 5)

pdf("9_results/figures/demo_example1_risk_factors.pdf", width = 8, height = 7)
plot_risk_factors(example1_sims, example1_df)
dev.off()

pdf("9_results/figures/demo_example1_ascvd.pdf", width = 8, height = 7)
plot_ascvd(example1_sims)
dev.off()


# example 2: --------------------------------------------------------------

# man in example 1 is diagnosed with diabetes at time 2
example2_df <- tibble(
  pid = 1,
  time = 1,
  sex = 1,
  age = 65,
  smk = 1,
  bmi = 30,
  dm = 1,
  hrx = 0,
  liprx = 0,
  ldl = 160,
  hdl = 40,
  sbp = 145,
  lag1_age = 60,
  lag1_smk = 1,
  lag1_bmi = 30,
  lag1_dm = 0,
  lag1_hrx = 0,
  lag1_liprx = 0,
  lag1_ldl = 160,
  lag1_hdl = 40,
  lag1_sbp = 145,
  lag2_age = 0,
  lag2_smk = 0,
  lag2_bmi = 0,
  lag2_dm = 0,
  lag2_hrx = 0,
  lag2_liprx = 0,
  lag2_ldl = 0,
  lag2_hdl = 0,
  lag2_sbp = 0,
  event_ascvd = 0,
  event_dth_ascvd = 0
)

example2_sims <- sim_risk_factors(example2_df, 1, 5)

pdf("9_results/figures/demo_example2_risk_factors.pdf", width = 8, height = 7)
plot_risk_factors(example2_sims, add_points = bind_rows(example1_df, example2_df))
dev.off()

pdf("9_results/figures/demo_example2_ascvd.pdf", width = 8, height = 7)
plot_ascvd(example2_sims, add_line = example1_sims)
dev.off()


# example 3: --------------------------------------------------------------

# if they also started statins and bp meds

example3_sims <-
  sim_intervention(example2_df, 1, 5, function(data, x) {
    x[, 1] <- 0
    x[, 2] <- 1
    x[, 3] <- 1
    x[, 4] <- 120
    x[, 5] <- 120
    x[, 6] <- 23
    x
  }, c("smk", "hrx", "liprx", "sbp", "ldl", "bmi"))

pdf("9_results/figures/demo_example3_risk_factors.pdf", width = 8, height = 7)
plot_risk_factors(example3_sims, add_points = bind_rows(example1_df, example2_df))
dev.off()

pdf("9_results/figures/demo_example3_ascvd.pdf", width = 8, height = 7)
plot_ascvd(example3_sims, add_line = example1_sims, color = "Greens")
dev.off()

