# Define global analysis parameters ---------------------------------------

CSV_ROOT_DIR <- "../../3_data_encrypted/Framingham_Offspring_2020c/Datasets/CSV/"
GFORM_SIM <- 100
GFORM_BSAMP <- 0

# Packages and helper functions -------------------------------------------

source("0_packages/load_packages.R")

source("1_helpers/1_helpers_generic.R")

source("1_helpers/2_helpers_parameteric_models.R")

source("1_helpers/3_helpers_gformula.R")

source("1_helpers/4_helpers_lasso.R")

source("1_helpers/5_helpers_cvd_risk_calcs.R")

source("1_helpers/6_helpers_plots.R")


# Load data ---------------------------------------------------------------

source("2_load_data/load_data.R")


# Data cleaning and preparation -------------------------------------------

source("3_clean_variables/0_create_variable_lists.R")

source("3_clean_variables/1_merge_datasets.R")

source("3_clean_variables/2_pivot_longer.R")

source("3_clean_variables/3_clean_variables.R")

source("3_clean_variables/4_pivot_wider.R")

source("3_clean_variables/5_scale_and_center.R")


# Descriptive analysis ----------------------------------------------------

#source("4_descriptive_analysis/0_baseline_summary.R")

source("4_descriptive_analysis/1_exam_summary.R")

source("4_descriptive_analysis/2_descriptive_plots.R")


# Simulations -------------------------------------------------------------

source("5_simulation/0_sem.R")

source("5_simulation/1_sim_functions.R")

source("5_simulation/2_run_sims.R")

source("5_simulation/3_plot_sim_results.R")


# Run parametric models ---------------------------------------------------

source("6_parametric_models/0_define_parametric_models.R")

source("6_parametric_models/1_fit_parametric_models.R")

source("6_parametric_models/2_model_diagnostics.R")

source("6_parametric_models/3_prediction.R")

source("6_parametric_models/4_performance.R")

source("6_parametric_models/5_validation.R")




