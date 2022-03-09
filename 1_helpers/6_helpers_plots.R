latex_theme <-
  function() {
    theme_bw() +
      theme(
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = '#eeeeee'),
        strip.background = element_blank(),
        legend.position = "bottom",
        text = element_text(family = "Palatino")
      )
  }

slides_theme <-
  function() {
    theme_bw() +
      theme(
        #strip.background = element_blank(),
        legend.position = "bottom",
        text = element_text(family = "Helvetica"),
        panel.spacing = unit(1.5, "lines"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", color = NA)
      )
  }
