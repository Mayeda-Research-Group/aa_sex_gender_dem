# Tables
# Created by: Yingyan Wu
# Sep.2.2024

#---- Package loading + options ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "plyr", "haven", "labelled", 
       "gtsummary")
# Do not load "summarytools" as it contains conflicted functions as dplyr
#No scientific notation
options(scipen = 999)

# Paths
source(here::here("01_R", "0_paths.R"))

#---- Load the data ----
load(paste0(path_to_analytic_data, "analysis_data/aa_adrd_tte.RData"))

#---- Table 1 ----
data_for_table_1 <- aa_adrd_tte_selected %>%
  mutate(current_work_status = case_when(current_work_status == 1 ~ 1,
                                         current_work_status %in% c(2, 3) ~ 0),
         
         sizeofhh_2 = case_when(sizeofhh_3 == 1 ~ 1,
                                sizeofhh_3 > 1 ~ 0)) %>%
  select(survey_age, female, ethnicity_rev, usaborn_rev, edu_4, 
         current_work_status, income_pp, 
         # employ_working, employ_retired, employ_not_workforce, 
         marital_2, sizeofhh_2,  
         generalhealth_3, sr_bmi, 
         sr_hyp, sr_diabetes, sr_stroke, sr_depress, 
         # sr_dem_kin, 
         main_dem_v1_fu_time, main_dem_v1_end_type) %>%
  set_variable_labels(survey_age = "Survey age, years",
                      female = "Female",
                      ethnicity_rev = "Race/ethnicity",
                      usaborn_rev = "US-born",
                      edu_4 = "Education attainment",
                      current_work_status = "Currently working", # "Work status"
                      marital_2 = "Married or living as married",
                      sizeofhh_2 = "Household size: living alone",
                      income_pp = "Household-adjusted income, USD",
                      generalhealth_3 = "Self-rated health",
                      sr_bmi = "Self-reported BMI",
                      sr_depress = "Self-reported depression",
                      sr_hyp = "Self-reported hypertension",
                      sr_stroke = "Self-reported stroke",
                      sr_diabetes = "Self-reported diabetes",
                      main_dem_v1_end_type = "End of follow-up event",
                      # sr_dem_kin = "Family member conditions: Dementia/Alzheimer's Disease",
                      main_dem_v1_fu_time = "Follow up time") %>%
  set_value_labels(edu_4 = c("<= Grade 11" = 1, 
                             "High school diploma" = 2,
                             "Some college" = 3,
                             "College graduate" = 4),
                   female = c("Women" = 1, "Men" = 0),
                   marital_2 = c("Yes" = 1, "No" = 0),
                   # sizeofhh_3 = c("Live alone" = 1, "Two" = 2, "Three or more" = 3),
                   sizeofhh_2 = c("Yes" = 1, "No" = 0),
                   # current_work_status = 
                   #   c("Current working" = 1, "Current retired" = 2, 
                   #     "Not in workforce" = 3),
                   current_work_status = c("Yes" = 1, "No" = 0),
                   generalhealth_3 = c("Excellent/Very Good" = 3,
                                       "Good" = 2,
                                       "Fair/Poor" = 1),
                   main_dem_v1_end_type = 
                     c("Dementia" = "DEMENTIA",
                       "Administratively Censored" = "ADMIN CENSORED",
                       "Death" = "DEATH",
                       "End of Membership" = "END OF MEMBERSHIP"))
# ethnicity_rev = c(
#   "South Asian" = 1,
#   "Chinese" = 2,
#   "Japanese" = 3,
#   "Korean" = 4,
#   "Filipino" = 5,
#   "Vietnamese" = 6,
#   "Other Southeast Asian" = 7,
#   "Any Pacific Islander" = 8,
#   "Non-Latino White" = 9,
#   "Multiple Asian ethn" = 10))
# sr_dem_kin = c("Yes" = 1, "No" = 0))

table_1_res <- data_for_table_1 %>%
  modify_if(is.labelled, to_factor) %>%
  tbl_strata(
    strata = ethnicity_rev,
    .tbl_fun =
      ~ .x %>%
      tbl_summary(type = all_continuous() ~ "continuous",
                  digits = list(all_continuous() ~ 1,
                                all_categorical() ~ c(0, 0),
                                c("income_pp") ~ 0), # change income to 0 digit
                  statistic = list(all_categorical() ~ "{n} ({p})",
                                   all_continuous() ~ "{mean} ({sd})"),
                  by = female,
                  missing = "ifany",
                  missing_text = "Missing") %>%
      # add_overall %>%
      modify_header(all_stat_cols() ~ 
                      "**{level}**, \nN = {n} ({style_percent(p, digits = 0)}%)")) %>%
  bold_labels()

table_1_res %>% as_flex_table()

table_1_res %>% 
  as_hux_xlsx(here::here("03_tables", "table1.xlsx"))
table_1_res %>% 
  as_flex_table() %>% 
  flextable::save_as_docx(path = here::here("03_tables", "table1.docx"))
