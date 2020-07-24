# Define global analysis parameters ---------------------------------------

CSV_ROOT_DIR <- "../../3_data_encrypted/Framingham_Offspring_2018b/Datasets/CSV/"
GFORM_SIM <- 10000
GFORM_BSAMP <- 10

# Packages and helper functions -------------------------------------------

source("0_packages/load_packages.R")

source("1_helpers/1_helpers_generic.R")

source("1_helpers/2_helpers_gformula.R")

source("1_helpers/3_helpers_plots.R")

#source("1_helpers/02_helpers_estimators.R")

#source("1_helpers/03_helpers_bootstrap.R")

#source("1_helpers/04_helpers_tables.R")

#source("1_helpers/05_helpers_plots.R")


# Load data ---------------------------------------------------------------

source("2_load_data/load_data.R")


# Data cleaning and preparation -------------------------------------------

source("3_clean_variables/0_create_variable_lists.R")

source("3_clean_variables/1_merge_datasets.R")

source("3_clean_variables/2_pivot_longer.R")

source("3_clean_variables/3_clean_variables.R")

source("3_clean_variables/4_pivot_wider.R")


# Descriptive analysis ----------------------------------------------------

source("4_descriptive_analysis/0_baseline_summary.R")

source("4_descriptive_analysis/1_exam_summary.R")


# Run parametric models ---------------------------------------------------

source("5_parametric_models/1_gformula_ldl.R")

source("5_parametric_models/2_gformula_sbp.R")

source("5_parametric_models/3_gformula_table.R")

#source("5_parametric_models/gcomp_subs.R")


