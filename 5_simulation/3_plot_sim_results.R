
# plot function -----------------------------------------------------------

plot_rmspe_comparisons <- function(data) {
  plot_df <- data %>%
    unnest(rmspe) %>%
    select(sim, n, estimator, specification, time, rmspe) %>%
    group_by(estimator, specification, time, n) %>%
    summarise(
      mean = mean(rmspe),
      lwr = quantile(rmspe, 0.10),
      upr = quantile(rmspe, 0.90),
      .groups = 'drop'
    )  %>%
    mutate(time = factor(time, labels = c("1 exam", "2 exams", "3 exams", "4 exams")))
  
    ggplot(plot_df,
           aes(
             x = factor(n),
             y = mean,
             color = estimator,
             shape = estimator,
             fill = estimator,
             group = estimator
           )) +
    facet_grid(specification ~ time, scales = "free_y") +
    scale_shape_manual(name = "", values = 21:24) +
    scale_color_brewer(name = "", palette = "RdPu") +
    scale_fill_brewer(name = "", palette = "RdPu") +
    geom_pointrange(aes(ymin = lwr, ymax = upr),
                    size = 0.5,
                    position = position_dodge(width = 0.75)) +
    ggpubr::theme_pubclean() +
    labs(
      x = "\nTraining sample size (N)",
      y = "RMSPE"
    ) +
    theme(legend.position = "bottom")
}


# generate plots ----------------------------------------------------------

plot_rmspe_comparisons(factual_sims_df)
plot_rmspe_comparisons(compevent_sims_df)
plot_rmspe_comparisons(counterfactual_sims_df)
