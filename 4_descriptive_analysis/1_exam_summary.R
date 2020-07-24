exam_mean <- 
  analytic_long %>%
  mutate(
    cpd = replace(cpd, smk == 0, NA),
    dpd = replace(dpd, drk == 0, NA)
  ) %>%
  select(covs_tv, time) %>%
  group_by(time) %>%
  summarise_all(mean, na.rm = T) %>%
  pivot_longer(-time) %>%
  rename(mean = value)

exam_sd <- 
  analytic_long %>%
  mutate(
    cpd = replace(cpd, smk == 0, NA),
    dpd = replace(cpd, drk == 0, NA)
  ) %>%
  select(covs_tv, time) %>%
  group_by(time) %>%
  summarise_all(sd, na.rm = T) %>%
  pivot_longer(-time) %>%
  rename(sd = value)


exam_vars <-
  tibble(
    name = unique(exam_mean$name),
    label = c(
      "Current smoker",
      "Cigarettes per day among smokers",
      "Current drinker",
      "Drinks per day among drinkers",
      "Body mass index (kg/m$^2$)",
      "Diabetes mellitus",
      "Systolic blood pressure (mmHg)",
      "LDL-cholesterol",
      "Blood pressure medication",
      "Lipid lowering medication"
    ),
    type = c("bin",
             "cont",
             "bin",
             "cont",
             "cont",
             "bin",
             "cont",
             "cont",
             "bin",
             "bin")
  )

exam_summary <-
  left_join(exam_vars, exam_mean, by = "name") %>%
  left_join(exam_sd, by = c("name", "time")) %>%
  mutate(
    mean = if_else(type == "bin", specd(mean * 100, 1), specd(mean, 1)),
    sd = specd(sd, 1),
    label = if_else(type == "cont", paste0(label, ", mean (SD)"), paste0(label, ", (\\%)")),
    value = if_else(type == "bin", as.character(mean), paste0(mean, " (", sd, ")"))
  ) %>%
  select(label, name, time, value) %>%
  pivot_wider(names_from = "time", values_from = "value", names_prefix = "exam_") 

sink("9_results/tables/exam_summary.tex")
options(knitr.kable.NA = '')
kable(
  x = exam_summary,
  format = "latex",
  digits = 1,
  align = "lccccc",
  col.names = c(
    "Characteristic",
    "Variable",
    "\\shortstack{4th exam \\\\ (1987–1991)}",
    "\\shortstack{5th exam \\\\ (1991–1994)}",
    "\\shortstack{6th exam \\\\ (1994–1998)}",
    "\\shortstack{7th exam \\\\ (1998–2001)}"
  ),
  escape = FALSE,
  booktabs = TRUE,
  linesep = ""
)
sink()