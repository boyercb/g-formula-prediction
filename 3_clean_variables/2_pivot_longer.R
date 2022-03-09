# create variables representing age and date at time 0 (start of follow up)
offspring$age0 <- offspring$age5
offspring$date0 <- offspring$date5

# pivot to appropriate person-time format for g-formula
offspring_long <- 
  offspring %>%
  pivot_longer(
    cols = matches(".*[5-9]"),
    names_to = c("variable", "exam"),
    names_pattern = "(.*)([5-9])"
  ) %>%
  pivot_wider(
    names_from = "variable", 
    values_from = "value"
  )

# limit to exams 4 through 9
#offspring_long <- filter(offspring_long, exam >= 5 & exam <= 7)
offspring_long <- filter(offspring_long, exam >= 4 & exam <= 9)