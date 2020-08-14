# simulate data -----------------------------------------------------------

sim <- datagen(2500)

# create long dataset
sim_long <-
  sim %>% 
  pivot_longer(
    -id,
    names_to = c("variable", "time"),
    names_pattern = "([LAY])([0-9])"
  ) %>%
  pivot_wider(names_from = "variable") %>%
  mutate(time = as.numeric(time))

# add lagged variables for L and A
sim_long <- mutate(
  sim_long, 
  across(c("L", "A"), lag, default = 0, .names = "lag1_{col}")
)

# 70/30 split into training and validation sets
train_rows <- sample(1:nrow(sim), size = floor(0.7 * nrow(sim)))

train <- sim[train_rows, ]
test <- sim[-train_rows, ] 

train_long <- left_join(train, sim_long, by = "id")
test_long <- left_join(test, sim_long, by = "id")


# define models -----------------------------------------------------------

X.models <-
  list(
    "L" = list(
      formula = L ~ lag1_L + lag1_A, 
      link = "logit",
      family = "binomial"
    ),
    "A" = list(
      formula = A ~ L + lag1_L + lag1_A, 
      link = "logit",
      family = "binomial"
    )
  )

D.model <- 
  list(
    formula = D ~ L + A + as.factor(time),
    link = "logit",
    family = "binomial"
  )

Y.model <- 
  list(
    formula = Y ~ L + A + as.factor(time),
    link = "logit",
    family = "binomial"
  )


# fit models --------------------------------------------------------------

Y.fit <- 
  fit_Y_model(
    Y.model,
    data = train_long
  )

X.fit <- 
  fit_X_models(
    X.models,
    data = train_long,
    time = "time"
  )

# D.fit <- 
#   fit_D_model(
#     D.model,
#     data = sim_w_lags
#   )


# run g-formula -----------------------------------------------------------

Y.hat_gformula <-
  gformula_mc(
    Y.fit,
    X.fit,
    #D.fit,
    data = test_long,
    id = "id",
    time = "time",
    hist.vars = c("lag1_A", "lag1_L"),
    hist.fun = "lag",
    mc.sims = 5000
  )


# run conventional model --------------------------------------------------

Y.fit_conventional <- 
  glm(Y1 ~ A0 + L0, data = train, family = binomial(link = "logit"))

Y.hat_conventional <-
  predict(Y.fit_conventional, type = 'response', newdata = test)

test$pred_c <- Y.hat_conventional


# combine and validate ----------------------------------------------------

test_w_pred <-
  cbind(test_long, 'pred_g' = Y.hat_gformula)

test_w_pred <- 
  left_join(test, filter(test_w_pred, time == 1))


test_plot_data <- 
  test_w_pred %>%
  mutate(bin_g = ntile(pred_g, 101),
         bin_c = ntile(pred_c, 101)) %>%
  pivot_longer(
    matches("_[gc]$"),
    names_to = c("metric", "calibration"),
    names_sep = "_"
  ) %>%
  pivot_wider(names_from = metric) %>%
  group_by(bin, calibration) %>%
  summarise(
    n = n(),
    # Get ests and CIs
    bin_pred = mean(pred),
    bin_prob = mean(Y1),
    se = sqrt((bin_prob * (1 - bin_prob)) / n),
    ul = bin_prob + 1.96 * se,
    ll = bin_prob - 1.96 * se
  ) %>%
  ungroup()

fit_lowess_g <- lowess(test_w_pred$pred_g, test_w_pred$Y1, iter = 0)
fit_lowess_c <- lowess(test_w_pred$pred_c, test_w_pred$Y1, iter = 0)
fit_lowess <- 
  bind_rows(fit_lowess_g, fit_lowess_c, .id = "id") %>%
  mutate(calibration = ifelse(id == 1, "g", "c"))

histogram_data <- 
  bind_rows(list("pred" = test_w_pred$pred_g),
            list("pred" = test_w_pred$pred_c),
            .id = "id") %>%
  mutate(calibration = ifelse(id == 1, "g", "c"))

stats_g <- rms::val.prob(test_w_pred$pred_g, test_w_pred$Y1)
stats_c <- rms::val.prob(test_w_pred$pred_c, test_w_pred$Y1)

stats <- bind_rows(stats_g, stats_c,.id = "id") %>% 
  mutate(calibration = ifelse(id == 1, "g", "c"))

pdf("9_results/figures/simulation_plot.pdf", width = 11.65, height = 6.59)
ggplot(test_plot_data, aes(x = bin_pred)) +
  facet_wrap(~ factor(calibration, labels = c("Conventional", "G-formula")), scales = "free_y") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  geom_abline(aes(
    linetype = "ideal",
    color = "ideal",
    alpha = "ideal",
    slope = 1,
    intercept = 0
  ),
  show.legend = FALSE) +
  geom_smooth(
    aes(y = bin_prob,
        color = "linear calibration",
        linetype = "linear calibration", 
        alpha = "linear calibration"),
    method = "lm",
    se = FALSE,
    formula = y ~ -1 + x,
    size = 1.3
  ) +
  geom_line(
    aes(
      x = x,
      y = y,
      color = "nonparametric calibration",
      linetype = "nonparametric calibration",
      alpha = "nonparametric calibration"
    ),
    data = fit_lowess,
    size = 1.3,
  ) +
  geom_text(aes(
    x = 0.15,
    y = 0.85,
    label = paste0(
      "C = ", specd(`C (ROC)`, 3), 
      "\nslope = ", specd(Slope, 3),
      "\nintercept = ", specd(Intercept, 3), 
      "\nR2 = ", specd(R2, 3),
      "\nE90 = ", specd(E90, 3)
      )
  ),
  data = stats,
  hjust = "left",
  vjust = "top",
  family = "Helvetica",
  size = 3.5) +
  scale_color_manual(
    name = "",
    labels = c("ideal", "linear calibration", "nonparametric calibration"),
    values = c(scales::alpha("black", 0.5), scales::alpha("red", 0.5), scales::alpha("blue", 0.5))
  ) +
  scale_linetype_manual(
    name = "",
    labels = c("ideal", "linear calibration", "nonparametric calibration"),
    values = c("dotted", "solid", "solid")
  ) +
  scale_alpha_manual(
    name = "",
    labels = c("ideal", "linear calibration", "nonparametric calibration"),
    values = c(0.5, 0.5, 0.5)
  ) +
  xlab("\nPredicted Probability") +
  ylab("Observed Probability\n") +
  coord_cartesian(expand = FALSE) +
  slides_theme() +
  theme(panel.spacing = unit(1.5, "lines"), 
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = NA)# get rid of legend panel bg  
        )
dev.off()

