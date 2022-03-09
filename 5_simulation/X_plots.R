factual_sims_df %>%
  unnest(rmspe) %>%
  ggplot(., aes(x = factor(n), y = rmspe, color = estimator, shape = estimatop)) + 
  facet_grid(~ time) +
  geom_point() +
  ggpubr::theme_cleveland()


factual_sims_df %>%
  unnest(rmspe) %>%
  group_by(estimator, time, n) %>%
  summarise(
    mean = mean(rmspe),
    lwr = quantile(rmspe, 0.10),
    upr = quantile(rmspe, 0.90)
  )  %>%
  ggplot(.,
         aes(
           x = factor(n),
           y = mean,
           color = estimator,
           shape = estimator
         )) +
  facet_grid( ~ time) +
  geom_pointrange(aes(ymin = lwr, ymax = upr),
                  size = 1,
                  position = position_dodge(width = 0.75)) +
  ggpubr::theme_pubr()


factual_sims_df %>%
  unnest(bias) %>%
  group_by(estimator, time, n) %>%
  summarise(
    mean = mean(bias),
    lwr = quantile(bias, 0.10),
    upr = quantile(bias, 0.90)
  )  %>%
  ggplot(.,
         aes(
           x = factor(n),
           y = mean,
           color = estimator,
           shape = estimator
         )) +
  facet_grid( ~ time) +
  geom_pointrange(aes(ymin = lwr, ymax = upr),
                  size = 1,
                  position = position_dodge(width = 0.75)) +
  ggpubr::theme_pubr()

factual_sims_df %>%
  unnest(rmspe) %>%
  group_by(estimator, time, n) %>%
  summarise(
    mean = mean(rmspe),
    lwr = quantile(rmspe, 0.10),
    upr = quantile(rmspe, 0.90)
  )  %>%
  ggplot(.,
         aes(
           x = time,
           y = mean,
           color = estimator,
           shape = estimator,
           group = estimator
         )) +
  facet_grid( ~ factor(n)) +
  scale_color_brewer(palette = "Set2") +
  geom_pointrange(aes(ymin = lwr, ymax = upr),
                  size = 1,
                  position = position_dodge(width = 0.75)) +
  ggpubr::theme_transparent()


factual_sims_df %>% 
  unnest(rmspe) %>% 
  select(sim, n, estimator, specification, time, rmspe) %>%
  ggplot(.,
         aes(
           x = factor(n),
           y = rmspe,
           color = estimator,
           shape = estimator,
           group = estimator
         )) +
  facet_grid(specification ~ time, scales = "free_y") +
  scale_color_brewer(palette = "Set2") +
  geom_point(position = position_dodge(width = 0.75), size = 2)

factual_sims_df %>% 
  unnest(bias) %>% 
  select(sim, n, estimator, specification, time, bias) %>%
  ggplot(.,
         aes(
           x = factor(n),
           y = bias,
           color = estimator,
           shape = estimator,
           group = estimator
         )) +
  facet_grid(specification ~ time, scales = "free_y") +
  scale_color_brewer(palette = "Set2") +
  geom_point(position = position_dodge(width = 0.75), size = 2)
