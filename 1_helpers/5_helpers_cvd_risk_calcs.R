#' ACC/AHA 2013 ASCVD risk score
#'
#' Computes 10-year risk for hard ASCVD event (defined as first occurrence of
#' non-fatal myocardial infarction (MI), congestive heart disease (CHD) death,
#' or fatal or nonfatal stroke).
#'
#' @param race patient race (white, aa)
#' @param gender patient gender (male, female)
#' @param age patient age (years)
#' @param totchol Total cholesterol (mg/dL)
#' @param hdl HDL cholesterol (mg/dL)
#' @param sbp Systolic blood pressure (mm Hg)
#' @param bp_med Patient is on a blood pressure medication (1=Yes, 0=No)
#' @param smoker Current smoker (1=Yes, 0=No)
#' @param diabetes Diabetes (1=Yes, 0=No)
#' @param ... Additional predictors can be passed and will be ignored
#'
#'
#' @return Estimated 10-Y Risk for hard ASCVD (percent)
#'
#' @export
#'
#' @examples
#' library(CVrisk)
#' ascvd_10yr_accaha(
#'   race = "aa", gender = "male", age = 55,
#'   totchol = 213, hdl = 50, sbp = 140,
#'   bp_med = 0, smoker = 0, diabetes = 0
#' )
#' @references
#' Goff, David C., et al. "2013 ACC/AHA guideline on the assessment of
#' cardiovascular risk: a report of the American College of
#' Cardiology/American Heart Association Task Force on Practice
#' Guidelines." Journal of the American College of Cardiology 63.25
#' Part B (2014): 2935-2959.
ascvd_10yr_accaha <- function(race = "white", gender = c("male", "female"),
                             age, totchol, hdl, sbp,
                             bp_med, smoker, diabetes, 
                             baseline_survival = c(0.9665, 0.9533, 0.9144, 0.8954),  ...) {
  if (any(!race %in% c("aa", "white") | missing(race))) {
    stop("race must be either 'aa' or 'white'")
  }
  
  if (any(!gender %in% c("male", "female") | missing(gender))) {
    stop("gender must be either 'male' or 'female'")
  }
  
  if (any(!is.numeric(age) | age < 1 | age > 120 | missing(age))) {
    stop("age must be a valid numeric value'")
  }
  
  if (any(!is.numeric(totchol) | totchol < 1 | totchol > 999 | missing(totchol))) {
    stop("totchol must be a valid numeric value'")
  }
  
  
  d <- tribble(
    ~race, ~gender, ~ln_age, ~ln_age_squared, ~ln_totchol, ~ln_age_totchol, ~ln_hdl, ~ln_age_hdl, ~ln_treated_sbp, ~ln_age_treated_sbp, ~ln_untreated_sbp, ~ln_age_untreated_sbp, ~smoker, ~ln_age_smoker, ~diabetes, ~group_mean, ~baseline_survival,
    "white", "female", -29.799, 4.884, 13.54, -3.114, -13.578, 3.149, 2.019, 0, 1.957, 0, 7.574, -1.665, 0.661, -29.18, baseline_survival[1],
    "aa", "female", 17.114, 0, 0.94, 0, -18.92, 4.475, 29.291, -6.432, 27.82, -6.087, 0.691, 0, 0.874, 86.61, baseline_survival[2],
    "white", "male", 12.344, 0, 11.853, -2.664, -7.99, 1.769, 1.797, 0, 1.764, 0, 7.837, -1.795, 0.658, 61.18, baseline_survival[3],
    "aa", "male", 2.469, 0, 0.302, 0, -0.307, 0, 1.916, 0, 1.809, 0, 0.549, 0, 0.645, 19.54, baseline_survival[4]
  )
  
  r <- mapply(
    FUN = function(race, gender)
      which(d$race == race & d$gender == gender),
    race = race,
    gender = gender
  )
  
  pooled_coef <- d[r, ]

  sbp_treated <- ifelse(bp_med == 1, sbp, 1)
  sbp_untreated <- ifelse(bp_med == 0, sbp, 1)

  indv_sum <- log(age) * pooled_coef$ln_age +
    log(age)^2 * pooled_coef$ln_age_squared +
    log(totchol) * pooled_coef$ln_totchol +
    log(age) * log(totchol) * pooled_coef$ln_age_totchol +
    log(hdl) * pooled_coef$ln_hdl +
    log(age) * log(hdl) * pooled_coef$ln_age_hdl +
    log(sbp_treated) * pooled_coef$ln_treated_sbp +
    log(sbp_treated) * log(age) * pooled_coef$ln_age_treated_sbp +
    log(sbp_untreated) * pooled_coef$ln_untreated_sbp +
    log(sbp_untreated) * log(age) * pooled_coef$ln_age_untreated_sbp +
    smoker * pooled_coef$smoker +
    smoker * log(age) * pooled_coef$ln_age_smoker +
    diabetes * pooled_coef$diabetes

  risk_score <- (1 - (pooled_coef$baseline_survival ^ exp(indv_sum - pooled_coef$group_mean)))

  return(risk_score)
}



#' Framingham 2008 ASCVD risk score (with lab measurement)
#'
#' Computes 10-year risk for ASCVD event (coronary death, myocardial
#' infarction (MI), coronary insufficiency, angina, ischemic stroke,
#' hemorrhagic stroke, transient ischemic attack, peripheral artery disease,
#' or heart failure).
#'
#' @param gender patient gender (male, female)
#' @param age patient age (years), between 30 and 74
#' @param hdl HDL cholesterol (mg/dL)
#' @param totchol Total cholesterol (mg/dL)
#' @param sbp Systolic blood pressure (mm Hg)
#' @param bp_med Patient is on a blood pressure medication (1=Yes, 0=No)
#' @param smoker Current smoker (1=Yes, 0=No)
#' @param diabetes Diabetes (1=Yes, 0=No)
#' @param ... Additional predictors can be passed and will be ignored
#'
#' @return Estimated 10-Y Risk for hard ASCVD event (percent)
#'
#' @export
#'
#' @examples
#' library(CVrisk)
#' ascvd_10y_frs(
#'   gender = "male", age = 55,
#'   hdl = 50, totchol = 213, sbp = 140,
#'   bp_med = 0, smoker = 0, diabetes = 0
#' )
#'
#' # 16.7
#' @references
#' D’agostino, R.B., Vasan, R.S., Pencina, M.J., Wolf, P.A., Cobain, M.,
#' Massaro, J.M. and Kannel, W.B., 2008. General cardiovascular risk
#' profile for use in primary care: the Framingham Heart Study.
#' Circulation, 117(6), pp.743-753.
ascvd_10yr_frs <- function(gender = c("male", "female"),
                          age, hdl, totchol, sbp,
                          bp_med, smoker, diabetes,
                          baseline_survival = c(0.88936, 0.95012), ...) {
  if (any(!gender %in% c("male", "female") | missing(gender))) {
    stop("gender must be either 'male' or 'female'")
  }
  
  if (any(!is.numeric(age) | missing(age))) {
    stop("age must be a valid numeric value'")
  }
  
  if (any(age < 1 | age > 120)) {
    return(NA)
  }
  
  
  # retrieve model coefficients
  d <- tribble(
    ~gender,
    ~ln_age,
    ~ln_totchol,
    ~ln_hdl,
    ~ln_untreated_sbp,
    ~ln_treated_sbp,
    ~smoker,
    ~diabetes,
    ~group_mean,
    ~baseline_survival,
    "male",
    3.06117,
    1.12370,
    -0.93263,
    1.93303,
    1.99881,
    0.65451,
    0.57367,
    23.9802,
    baseline_survival[1],
    "female",
    2.32888,
    1.20904,
    -0.70833,
    2.76157,
    2.82263,
    0.52873,
    0.69154,
    26.1931,
    baseline_survival[2]
  )
  
  r <- mapply(
    FUN = function(race, gender)
      which(d$gender == gender),
    gender = gender
  )
  
  model_coef <- d[r, ]
  
  sbp_treated <- ifelse(bp_med == 1, sbp, 1)
  sbp_untreated <- ifelse(bp_med == 0, sbp, 1)
  
  indv_sum <- log(age) * model_coef$ln_age +
    log(hdl) * model_coef$ln_hdl +
    log(totchol) * model_coef$ln_totchol +
    log(sbp_treated) * model_coef$ln_treated_sbp +
    log(sbp_untreated) * model_coef$ln_untreated_sbp +
    smoker * model_coef$smoker +
    diabetes * model_coef$diabetes
  
  risk_score <- (1 - (model_coef$baseline_survival ^ exp(indv_sum - model_coef$group_mean)))
  
  return(risk_score)  
}

#' Framingham 2008 ASCVD risk score (no lab measurement)
#'
#' Computes 10-year risk for ASCVD event (coronary death, myocardial
#' infarction (MI),coronary insufficiency, angina, ischemic stroke,
#' hemorrhagic stroke, transient ischemic attack, peripheral artery
#' disease, or heart failure).
#'
#' @param gender patient gender (male, female)
#' @param age patient age (years), between 30 and 74
#' @param bmi Body mass index (kg/m2)
#' @param sbp Systolic blood pressure (mm Hg)
#' @param bp_med Patient is on a blood pressure medication (1=Yes, 0=No)
#' @param smoker Current smoker (1=Yes, 0=No)
#' @param diabetes Diabetes (1=Yes, 0=No)
#' @param ... Additional predictors can be passed and will be ignored
#'
#' @return Estimated 10-Y Risk for hard ASCVD (percent)
#'
#' @export
#'
#' @examples
#' library(CVrisk)
#' ascvd_10y_frs_simple(
#'   gender = "male", age = 55,
#'   bmi = 30, sbp = 140,
#'   bp_med = 0, smoker = 0, diabetes = 0
#' )
#'
#' # 16.7
#' @references
#' D’agostino, R.B., Vasan, R.S., Pencina, M.J., Wolf, P.A., Cobain, M.,
#' Massaro, J.M. and Kannel, W.B., 2008. General cardiovascular risk
#' profile for use in primary care: the Framingham Heart Study.
#' Circulation, 117(6), pp.743-753.
ascvd_10yr_frs_simple <- function(gender = c("male", "female"),
                                 age, bmi, sbp,
                                 bp_med, smoker, diabetes,
                                 baseline_survival = c(0.88431, 0.94833), ...) {
  
  if (any(!gender %in% c("male", "female") | missing(gender))) {
    stop("gender must be either 'male' or 'female'")
  }
  
  if (any(!is.numeric(age) | missing(age))) {
    stop("age must be a valid numeric value'")
  }
  
  if (any(age < 1 | age > 120)) {
    return(NA)
  }
  
  # retrieve model coefficients
  d <- tribble(
    ~gender,
    ~ln_age,
    ~ln_bmi,
    ~ln_untreated_sbp,
    ~ln_treated_sbp,
    ~smoker,
    ~diabetes,
    ~group_mean,
    ~baseline_survival,
    "male",
    3.11296,
    0.79277,
    1.85508,
    1.92672,
    0.70953,
    0.53160,
    23.9388,
    baseline_survival[1],
    "female",
    2.72107,
    0.51125,
    2.81291,
    2.88267,
    0.61868,
    0.77763,
    26.0145,
    baseline_survival[2]
  )
  
  r <- mapply(
    FUN = function(race, gender)
      which(d$gender == gender),
    gender = gender
  )
  
  model_coef <- d[r, ]
  
  sbp_treated <- ifelse(bp_med == 1, sbp, 1)
  sbp_untreated <- ifelse(bp_med == 0, sbp, 1)
  
  indv_sum <- log(age) * model_coef$ln_age +
    log(bmi) * model_coef$ln_bmi +
    log(sbp_treated) * model_coef$ln_treated_sbp +
    log(sbp_untreated) * model_coef$ln_untreated_sbp +
    smoker * model_coef$smoker +
    diabetes * model_coef$diabetes
  
  risk_score <- (1 - (model_coef$baseline_survival ^ exp(indv_sum - model_coef$group_mean)))
  
  return(risk_score)
}
