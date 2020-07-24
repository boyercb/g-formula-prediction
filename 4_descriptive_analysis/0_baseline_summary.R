baseline_mean <- 
  analytic_long %>%
  select(covs_fixed) %>%
  summarise_all(mean, na.rm = T) %>%
  pivot_longer(everything()) %>%
  rename(mean = value)

baseline_sd <- 
  analytic_long %>%
  select(covs_fixed) %>%
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
      "No drinks per day",
      "1 to <2 drinks per day",
      "2 to <4 drinks per day",
      "4 or more drinks per day",
      "Body mass index",
      "Diabetes mellitus",
      "Systolic blood pressure",
      "No cigarretes per day",
      "1 cigarretes per day",
      "2 to 5 cigarretes per day",
      "5 to 24 or cigarretes per day",
      "25 or more cigarretes per day",
      "LDL-cholesterol",
      "Blood pressure medication",
      "Anti-cholesterol medication"
    ),
    exam = c(
      rep("4", 2),
      rep("3", 4),
      rep("4", 3),
      rep("3", 16)
    ),
    yrs = c(
      rep("1987–1991", 2),
      rep("1984-1987", 4),
      rep("1987–1991", 3),
      rep("1984-1987", 16)
    ),
    type = c("bin",
             "cont",
             rep("bin", 12),
             "cont",
             "bin",
             "cont",
             rep("bin", 5),
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