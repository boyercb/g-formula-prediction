gform_table_ldl <- 
  fit_ldl$result %>% 
  filter(k == 3) %>% 
  mutate(`Interv.` = case_when(
    `Interv.` == 0 ~ "natural course",
    `Interv.` == 1 ~ "low LDL (<130 mg/dL) for 16 years",
    `Interv.` == 2 ~ "moderate LDL (130 mg/dL to <160 mg/dL) for 16 years",
    `Interv.` == 3 ~ "high LDL (160 mg/dL to <190 mg/dL) for 16 years",
    `Interv.` == 4 ~ "very high LDL (>190 mg/dL) for 16 years"
  ),
  r_95 = paste0("(", specd(`Risk lower 95% CI`, 3), ", ", specd(`Risk upper 95% CI`, 3), ")"),
  rr_95 = paste0("(", specd(`RR lower 95% CI`, 3), ", ", specd(`RR upper 95% CI`, 3), ")"),
  rd_95 = paste0("(", specd(`RD lower 95% CI`, 3), ", ", specd(`RD upper 95% CI`, 3), ")")
  ) %>%
  select(
    `Interv.`,
    `NP risk`,
    `g-form risk`,
    r_95,
    `Risk ratio`,
    rr_95, 
    `Risk difference`,
    rd_95
  )

gform_table_sbp <- 
  fit_sbp$result %>% 
  filter(k == 3) %>% 
  mutate(`Interv.` = case_when(
    `Interv.` == 0 ~ "natural course",
    `Interv.` == 1 ~ "low SBP (<120 mmHg) for 16 years",
    `Interv.` == 2 ~ "prehypertension (120 mmHg to <140 mg/dL) for 16 years",
    `Interv.` == 3 ~ "stage 1 hypertension (140 mmHg to <160 mmHg) for 16 years",
    `Interv.` == 4 ~ "stage 2 hypertension (>160 mmHg) for 16 years"
  ),
  r_95 = paste0("(", specd(`Risk lower 95% CI`, 3), ", ", specd(`Risk upper 95% CI`, 3), ")"),
  rr_95 = paste0("(", specd(`RR lower 95% CI`, 3), ", ", specd(`RR upper 95% CI`, 3), ")"),
  rd_95 = paste0("(", specd(`RD lower 95% CI`, 3), ", ", specd(`RD upper 95% CI`, 3), ")")
  ) %>%
  select(
    `Interv.`,
    `NP risk`,
    `g-form risk`,
    r_95,
    `Risk ratio`,
    rr_95, 
    `Risk difference`,
    rd_95
  )

sink("9_results/tables/gformula.tex")
options(knitr.kable.NA = '')
kable(
  bind_rows(gform_table_ldl, gform_table_sbp), 
  digits = 3,
  align = "lccccccc",
  col.names = c(
    "Intervention",
    "Nonparametric risk",
    "G-formula risk",
    "95% CI",
    "Risk ratio",
    "95% CI",
    "Risk difference",
    "95% CI"
  ),
  format = "latex",
  booktabs = TRUE
) %>% 
  kable_styling() %>%
  group_rows(
    start_row = 1,
    5,
    group_label = "Interventions on LDL-cholesterol",
    bold = FALSE,
    italic = TRUE
  ) %>%
  group_rows(
    start_row = 6,
    10,
    group_label = "Interventions on systolic blood pressure",
    bold = FALSE,
    italic = TRUE
  ) %>% 
print()
sink()