exam_mean <- 
  offspring_long %>%
  # mutate(
  #   cpd = replace(cpd, smk == 0, NA),
  #   dpd = replace(dpd, drk == 0, NA)
  # ) %>%
  select(all_of(covs_tv), all_of(dvs), time) %>%
  group_by(time) %>%
  summarise(
    across(all_of(covs_tv), mean, na.rm = T), 
    across(all_of(dvs), sum, na.rm = T),
    .groups = "drop"
    ) %>%
  pivot_longer(-time) %>%
  rename(mean = value)

exam_sd <- 
  offspring_long %>%
  # mutate(
  #   cpd = replace(cpd, smk == 0, NA),
  #   dpd = replace(cpd, drk == 0, NA)
  # ) %>%
  select(all_of(covs_tv), all_of(dvs), time) %>%
  group_by(time) %>%
  summarise(
    across(everything(), sd, na.rm = T),
    .groups = "drop"
    ) %>%
  pivot_longer(-time) %>%
  rename(sd = value)


exam_vars <-
  tibble(
    name = unique(exam_mean$name),
    label = c(
      "Age",
      "Current smoker",
      # "Cigarettes per day among smokers",
      # "Current drinker",
      # "Drinks per day among drinkers",
      "Body mass index (kg/m$^2$)",
      "Diabetes mellitus",
      "Blood pressure medication",
      "Lipid lowering medication",
      "Total cholesterol (mg/dL)",
      "LDL cholesterol (mg/dL)",
      "HDL cholesterol (mg/dL)",
      "Systolic blood pressure (mmHg)",
      "Coronary heart disease events (Y)",
      "Atherosclerotic cardiovascular disease events (Y)",
      "non-CHD deaths (D)",
      "non-ASCVD deaths (D)",
      "Lost to follow up (C)"
    ),
    type = c("cont",
             "bin",
             # "cont",
             # "bin",
             # "cont",
             "cont",
             "bin",
             "bin",
             "bin",
             "cont",
             "cont",
             "cont",
             "cont",
             "count",
             "count",
             "count",
             "count",
             "count")
  )

exam_summary <-
  left_join(exam_vars, exam_mean, by = "name") %>%
  left_join(exam_sd, by = c("name", "time")) %>%
  mutate(
    mean = case_when(
      type == "bin" ~ specd(mean * 100, 1), 
      type == "cont"~ specd(mean, 1),
      type == "count"~ specd(mean, 0)),
    sd = specd(sd, 1),
    label = case_when(
      type == "cont" ~ paste0(label, ", mean (SD)"), 
      type == "bin" ~ paste0(label, ", (\\%)"),
      type == "count" ~ label
    ),
    value = case_when(
      type == "bin" ~ as.character(mean), 
      type == "cont" ~ paste0(mean, " (", sd, ")"),
      type == "count" ~ as.character(mean)
      )
  ) %>%
  select(label, name, time, value) %>%
  pivot_wider(names_from = "time", values_from = "value", names_prefix = "exam_") 

sink("9_results/tables/exam_summary.tex")
options(knitr.kable.NA = '')
kable(
  x = exam_summary,
  format = "latex",
  digits = 1,
  align = "lcccc",
  col.names = c(
    "Characteristic (Z)",
    "Variable",
    # "\\shortstack{4th exam \\\\ (1987–1991)}",
    "\\shortstack{5th exam \\\\ (1991–1994)}",
    "\\shortstack{6th exam \\\\ (1994–1998)}",
    "\\shortstack{7th exam \\\\ (1998–2001)}"
    # "\\shortstack{8th exam \\\\ (2001–2004)}",
    # "\\shortstack{9th exam \\\\ (2004–2008)}"
  ),
  escape = FALSE,
  booktabs = TRUE,
  linesep = ""
) %>% kable_styling() %>%
  row_spec(length(covs_tv), hline_after = TRUE) %>%
  print()
sink()
