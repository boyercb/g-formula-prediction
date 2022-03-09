baseline_mean <- 
  offspring_long %>%
  select(all_of(covs_fixed)) %>%
  summarise_all(mean, na.rm = T) %>%
  pivot_longer(everything()) %>%
  rename(mean = value)

baseline_sd <- 
  offspring_long %>%
  select(all_of(covs_fixed)) %>%
  summarise_all(sd, na.rm = T) %>%
  pivot_longer(everything()) %>%
  rename(sd = value)

baseline_vars <-
  tibble(
    name = baseline_mean$name,
    label = c(
      "Female sex",
      "Age at baseline",
      "None of below",
      "High school/junior college",
      "Bachelor's degree",
      "Postgraduate degree",
      "Single",
      "Married",
      "Divorced, widowed, or separated",
      "Ever smoked",
      "Drinks per day",
      "Body mass index",
      "Diabetes mellitus",
      "Systolic blood pressure",
      "Cigarettes per day",
      "LDL-cholesterol",
      "Blood pressure medication",
      "Anti-cholesterol medication"
    ),
    exam = c(
      rep("4", 2),
      rep("3", 4),
      rep("4", 3),
      rep("3", 9)
    ),
    yrs = c(
      rep("1987–1991", 2),
      rep("1984-1987", 4),
      rep("1987–1991", 3),
      rep("1984-1987", 9)
    ),
    type = c("bin",
             "cont",
             rep("bin", 8),
             "cont",
             "cont",
             "bin",
             "cont",
             "cont",
             "cont",
             "bin",
             "bin")
  )

baseline_summary <-
  left_join(baseline_vars, baseline_mean, by = "name") %>%
  left_join(baseline_sd, by = "name") %>%
  mutate(
    mean = if_else(type == "bin", specd(mean * 100, 1), specd(mean, 1)),
    sd = specd(sd, 1),
    label = if_else(type == "cont", paste0(label, ", mean (SD)"), label),
    value = if_else(type == "bin", as.character(mean), paste0(mean, " (", sd, ")"))
  ) %>%
  select(label, name, exam, yrs, value)
  
sink("9_results/tables/baseline.tex")
options(knitr.kable.NA = '')
kable(
  x = baseline_summary,
  format = "latex",
  digits = 1,
  align = "lcccc",
  col.names = c("Characteristic", "Variable", "Examination cycle", "Years", "%"),
  booktabs = TRUE,
  linesep = ""
)
sink()