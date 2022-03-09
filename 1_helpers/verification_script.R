# verification that my gformula code matches the results of the gfoRmula package

# results using gfoRmula package from CRAN --------------------------------

shift1 <- function (newdf, pool, intvar, intvals, time_name, t) {
  if (t > 0) {
    data.table::set(newdf, j = intvar, value = intvals[[t + 1]])
  }
}


time_points <- 7

pkg_results <- 
  gformula_survival(
    obs_data = basicdata_nocomp,
    id = 'id',
    time_points = time_points,
    time_name = 't0',
    covnames = c('L1', 'L2', 'A'),
    covtypes = c('binary', 'normal', 'binary'),
    outcome_name = 'Y',
    histories = c(lagged),
    histvars = list(c('A', 'L1', 'L2')),
    covparams = list(
      covmodels = c(
        L1 ~ lag1_A + lag1_L1 + lag1_L2 + L3 + t0,
        L2 ~ lag1_A + L1 + lag1_L1 + lag1_L2 + L3 + t0,
        A ~ lag1_A + L1 + L2 + lag1_L1 + lag1_L2 + L3 + t0
      )
    ),
    ymodel = Y ~ A + L1 + L2 + L3 + lag1_A + lag1_L1 + lag1_L2 + t0,
    intvars = list('A', 'A'),
    interventions = list(list(c(shift1, c(0, rep(0, time_points - 1)))),
                         list(c(shift1, c(0, rep(1, time_points - 1))))),
    int_descript = c('Never treat', 'Always treat'),
    nsimul = length(unique(basicdata_nocomp$id)),
    basecovs = c('L3'),
    seed = 1234,
    model_fits = TRUE,
    sim_data_b = TRUE
  )


# results using my gformula_mc function  ----------------------------------

# create lagged variables
basicdata_nocomp_lags <- as.data.frame(basicdata_nocomp) %>%
  group_by(id) %>%
  mutate(
    lag1_A = lag(A, default = 0),
    lag1_L1 = lag(L1, default = 0),
    lag1_L2 = lag(L2, default = 0)
  ) %>%
  ungroup()

# fit covariate models
X.models <- fit_X_models(
  model = list(
    "L1" = list(
      formula = L1 ~ lag1_A + lag1_L1 + lag1_L2 + L3 + t0,
      link = "logit",
      family = "binomial"
    ),
    "L2" = list(
      formula = L2 ~ lag1_A + L1 + lag1_L1 + lag1_L2 + L3 + t0,
      link = "identity",
      family = "normal"
    ),
    "A" = list(
      formula = A ~ lag1_A + L1 + L2 + lag1_L1 + lag1_L2 + L3 + t0,
      link = "logit",
      family = "binomial"
    )
  ),
  data = basicdata_nocomp_lags,
  time = 't0'
)

# fit outcome models
Y.model <- fit_Y_model(
  model = list(
    formula = Y ~ A + L1 + L2 + L3 + lag1_A + lag1_L1 + lag1_L2 + t0,
    link = "logit",
    family = "binomial"
  ),
  data = basicdata_nocomp_lags
)

# run the monte carlo simulation
mc_results <- gformula_mc(
  Y.fit = Y.model,
  X.fit = X.models,
  data = basicdata_nocomp_lags,
  id = "id",
  time = "t0",
  base.covs = "L3",
  hist.vars = c(paste0("lag1_", c('A', 'L1', 'L2'))),
  hist.fun = "lag",
  mc.sims = 1,
  mc.start = 0,
  mc.stop = 6,
  merge = FALSE,
  pool = TRUE,
  seed = 56485
)

# run the monte carlo simulation
mc_never_treat <- gformula_mc(
  Y.fit = Y.model,
  X.fit = X.models,
  data = as.data.frame(basicdata_nocomp_lags),
  id = "id",
  time = "t0",
  base.covs = "L3",
  hist.vars = c(paste0("lag1_", c('A', 'L1', 'L2'))),
  hist.fun = "lag",
  treatment = "A",
  intervention = function(data, x) { 0 },
  mc.sims = 1,
  mc.start = 0,
  mc.stop = 6,
  merge = FALSE,
  seed = 1235
)

# run the monte carlo simulation
mc_always_treat <- gformula_mc(
  Y.fit = Y.model,
  X.fit = X.models,
  data = as.data.frame(basicdata_nocomp_lags),
  id = "id",
  time = "t0",
  base.covs = "L3",
  hist.vars = c(paste0("lag1_", c('A', 'L1', 'L2'))),
  hist.fun = "lag",
  treatment = "A",
  intervention = function(data, x) { 1 },
  mc.sims = 1,
  mc.start = 0,
  mc.stop = 6,
  merge = FALSE,
  seed = 1235
)


# verify that model fits are equivalent -----------------------------------

# Y
assertthat::are_equal(
  coef(Y.model), coef(pkg_results$fits$Y)
)

# L1
assertthat::are_equal(
  coef(X.models$L1), coef(pkg_results$fits$L1)
)

# L2
assertthat::are_equal(
  coef(X.models$L2), coef(pkg_results$fits$L2)
)

# A
assertthat::are_equal(
  coef(X.models$A), coef(pkg_results$fits$A)
)

# check simulations converge ----------------------------------------------

library(furrr)
library(future)

SIMS <- 10000

plan(multisession, workers = 4)

convergence_check <- 
  future_map_dfr(1:SIMS, function(x) {
    mc_result <- gformula_mc(
      Y.fit = Y.model,
      X.fit = X.models,
      data = basicdata_nocomp_lags,
      id = "id",
      time = "t0",
      base.covs = "L3",
      hist.vars = c(paste0("lag1_", c('A', 'L1', 'L2'))),
      hist.fun = "lag",
      mc.sims = 1,
      mc.start = 0,
      mc.stop = 6,
      merge = FALSE
    )
    
    mc_never_treat <- gformula_mc(
      Y.fit = Y.model,
      X.fit = X.models,
      data = as.data.frame(basicdata_nocomp_lags),
      id = "id",
      time = "t0",
      base.covs = "L3",
      hist.vars = c(paste0("lag1_", c('A', 'L1', 'L2'))),
      hist.fun = "lag",
      treatment = "A",
      intervention = function(data, x) { 0 },
      mc.sims = 1,
      mc.start = 0,
      mc.stop = 6,
      merge = FALSE
    )
    
    mc_always_treat <- gformula_mc(
      Y.fit = Y.model,
      X.fit = X.models,
      data = as.data.frame(basicdata_nocomp_lags),
      id = "id",
      time = "t0",
      base.covs = "L3",
      hist.vars = c(paste0("lag1_", c('A', 'L1', 'L2'))),
      hist.fun = "lag",
      treatment = "A",
      intervention = function(data, x) { 1 },
      mc.sims = 1,
      mc.start = 0,
      mc.stop = 6,
      merge = FALSE
    )
    
    time_points <- 7
    
    pkg_results <- 
      gformula_survival(
        obs_data = basicdata_nocomp,
        id = 'id',
        time_points = time_points,
        time_name = 't0',
        covnames = c('L1', 'L2', 'A'),
        covtypes = c('binary', 'normal', 'binary'),
        outcome_name = 'Y',
        histories = c(lagged),
        histvars = list(c('A', 'L1', 'L2')),
        covparams = list(
          covmodels = c(
            L1 ~ lag1_A + lag1_L1 + lag1_L2 + L3 + t0,
            L2 ~ lag1_A + L1 + lag1_L1 + lag1_L2 + L3 + t0,
            A ~ lag1_A + L1 + L2 + lag1_L1 + lag1_L2 + L3 + t0
          )
        ),
        ymodel = Y ~ A + L1 + L2 + L3 + lag1_A + lag1_L1 + lag1_L2 + t0,
        intvars = list('A', 'A'),
        interventions = list(list(c(shift1, rep(0, time_points))),
                             list(c(shift1, rep(1, time_points)))),
        int_descript = c('Never treat', 'Always treat'),
        nsimul = length(unique(basicdata_nocomp$id)),
        basecovs = c('L3'),
        seed = x + 153
      )
    
    c(
      "sim" = x,
      "mc_nc" = mean(filter(mc_result, t0 == 6) %>% pull(pred)),
      "mc_nt" = mean(filter(mc_never_treat, t0 == 6) %>% pull(pred)),
      "mc_at" = mean(filter(mc_always_treat, t0 == 6) %>% pull(pred)),
      "pkg_nc" = filter(pkg_results$result, `Interv.` == 0 & k == 6) %>% pull(`g-form risk`),
      "pkg_nt" = filter(pkg_results$result, `Interv.` == 1 & k == 6) %>% pull(`g-form risk`),
      "pkg_at" = filter(pkg_results$result, `Interv.` == 2 & k == 6) %>% pull(`g-form risk`)
    )
    
    },
    .progress = TRUE
  )

convergence_plot <- 
  pivot_longer(
    data = convergence_check,
    cols = -sim,
    names_to = "variable",
    values_to = "value"
  ) %>%
  separate(variable, c("source", "target")) %>%
  mutate(
    source = case_when(
      source == "mc" ~ "new conditional estimator",
      source == "pkg" ~ "gfoRmula package"
    ),
    target = case_when(
      target == "at" ~ "Always treat",
      target == "nt" ~ "Never treat",
      target == "nc" ~ "Natural course"
    )
  )

g <- ggplot(convergence_plot, aes(x = value, fill = source)) +
  facet_grid(~target, scales = "free_x") +
  geom_histogram(bins = 100, alpha = 0.5, position = "identity") +
  geom_vline(aes(xintercept = value),
             data = convergence_plot %>% 
               group_by(target, source) %>% 
               summarise(value = mean(value)),
             linetype = "dotted"
              ) +
  geom_text(
    aes(label = label, x = value_q90, y = y, color = source),
    data = convergence_plot %>% 
      group_by(target, source) %>% 
      summarise(mean_value = mean(value),
                value_q10 = quantile(value, 0.05),
                value_q90 = quantile(value, 0.95)) %>%
      mutate(
        y = case_when(
          source == "new conditional estimator" ~ 342.5,
          source == "gfoRmula package" ~ 327.5
        ),
        label = paste0("tilde(theta) ==", " ", specd(mean_value, 2))
      ) %>%
      group_by(target) %>%
      mutate(
        value_q10 = max(value_q10),
        value_q90 = max(value_q90)
      ),
    hjust = "left",
    parse = TRUE,
    size = 5
  ) + 
  coord_cartesian(expand = FALSE) +
  scale_fill_brewer(name = "", type = "qual", palette = "Set1") +
  scale_color_brewer(name = "", type = "qual", palette = "Set1", guide = FALSE) +
  labs(
    x = "\nestimate",
    y = "count\n"
  ) +
  slides_theme() +
  theme(
    text = element_text(size = 16)
  )

ggsave("9_results/figures/convergence_check.pdf", g, device = pdf, width = 11, height = 8, dpi = 600)
ggsave("9_results/figures/convergence_check.png", g, device = "png", dpi = 600, bg = "transparent")
