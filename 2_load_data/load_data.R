csv_list <- c(
  "ex1_1d_v3.csv",
  "ex1_2d_v3.csv",
  "ex1_3d_v1.csv",
  "ex1_4d_v1.csv",
  "ex1_5d_v1.csv",
  "ex1_6d_v1.csv",
  "ex1_7d_v2.csv",
  "e_exam_ex08_1_0005d.csv", 
  "e_exam_ex09_1b_0844d.csv",
  "q_psych_ex03_1_0167d.csv",
  "vr_diab_ex09_1_1002d.csv",
  "vr_survcvd_2017_a_1194d.csv",
  "vr_survdth_2017_a_1192d.csv",
  "vr_svstk_2017_a_1196d.csv",
  "vr_wkthru_ex09_1_1001d.csv"
)

csvs <- lapply(csv_list, function (x) {
  df <- read_csv(get_data(x), show_col_types = FALSE)
  names(df) <- tolower(names(df))
  return(df)
}) 

names(csvs) <- gsub(".csv", "", csv_list)
