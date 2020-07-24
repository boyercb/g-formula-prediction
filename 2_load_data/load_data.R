csv_list <- c(
  "ex1_1d_v3.csv",
  "ex1_2d_v3.csv",
  "ex1_3d_v1.csv",
  "ex1_4d_v1.csv",
  "ex1_5d_v1.csv",
  "ex1_6d_v1.csv",
  "ex1_7d_v2.csv",
  "q_psych_ex03_1_0167d.csv",
  "vr_diab_ex09_1_1002d.csv",
  "vr_survcvd_2014_a_1023d.csv",
  "vr_survdth_2014_a_1025d.csv",
  "vr_wkthru_ex09_1_1001d.csv"
)

framingham <- lapply(csv_list, function (x) {
  df <- read_csv(get_data(x))
  names(df) <- tolower(names(df))
  return(df)
}) 

names(framingham) <- gsub(".csv", "", csv_list)
