# create variables representing age and date at time 0 (start of follow up)
analytic$age0 <- analytic$age4
analytic$date0 <- analytic$date4

# pivot to appropriate person-time format for g-formula
analytic_long <- 
  analytic %>%
  pivot_longer(
    cols = matches(".*[4-9]"),
    names_to = c("variable", "exam"),
    names_pattern = "(.*)([4-9])"
  ) %>%
  pivot_wider(
    names_from = "variable", 
    values_from = "value"
  )

# limit to exams 4 through 7
analytic_long <- filter(analytic_long, exam >= 4 & exam <= 7)
