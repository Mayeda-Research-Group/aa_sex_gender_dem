# Data construction for sex/gender & dementia paper
# Created by: Yingyan Wu
# Date: 09.02.2024

#---- Package loading + options ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "plyr", "haven", "labelled", 
       "gtsummary", "janitor")

options(scipen = 999)

#---- Load the datasets ----
source(here::here("01_R", "0_paths.R"))

aa_adrd_tte <- read_sas(paste0(path_to_raw_data, 
                               "aa_adrd_cardiometabolic_tte.sas7bdat")) %>%
  clean_names()

aa_adrd_time_to_event <- read_sas(paste0(path_to_raw_data, 
                                         "aa_adrd_time_to_event.sas7bdat")) %>%
  clean_names()

#---- Education ----
aa_adrd_tte %<>%
  mutate(edu_4 = case_when(education_rev %in% c(1, 2) ~ 1,
                           education_rev == 3 ~ 2,
                           education_rev == 4 ~ 3,
                           education_rev %in% c(5, 6) ~ 4))
#---- Marital status ----
# "Never Married" = 1, "Married or living as married" = 2,
# "Separated/Divorced" = 3, "Widowed" = 4
aa_adrd_tte %<>%
  mutate(marital_2 = case_when(maritalstatus == 2 ~ 1,
                               maritalstatus %in% c(1, 3, 4) ~ 0))

#---- Household size ----
aa_adrd_tte %<>%
  mutate(sizeofhh_3 = case_when(sizeofhh == 1 ~ 1,
                                sizeofhh == 2 ~ 2,
                                sizeofhh >= 3 ~ 3))

#---- General health ----
# "Excellent" = 5, "Very Good" = 4, "Good" = 3, "Fair" = 2, "Poor" = 1
aa_adrd_tte %<>%
  mutate(generalhealth_3 = case_when(generalhealth %in% c(1, 2) ~ 1,
                                     generalhealth == 3 ~ 2,
                                     generalhealth %in% c(4, 5) ~ 3))

#---- Currently work ----
aa_adrd_tte %>% select(starts_with("employment")) %>%
  apply(., 2, table, useNA = "ifany")

aa_adrd_tte %<>%
  mutate(employ_working = case_when(
    employment_full_time_employed == 1 | employment_part_time_employed == 1 ~ 1,
    rowSums(across(starts_with("employment"))) == 0 ~ NA_real_,
    employment_full_time_employed == 0 & employment_part_time_employed == 0 ~ 0),
    employ_retired = case_when(
      employ_working == 1 ~ 0,
      rowSums(across(starts_with("employment"))) == 0 ~ NA_real_,
      TRUE ~ employment_retired),
    employ_not_workforce = case_when(
      employ_working == 1 | employ_retired == 1 ~ 0,
      employment_disabled == 1 | employment_full_time_student == 1 |
        employment_homemaker == 1 | employment_unemployed == 1 |
        employment_other == 1 ~ 1,
      rowSums(across(starts_with("employment"))) == 0 ~ NA_real_,
      TRUE ~ 0))

# Sanity check
aa_adrd_tte %>% select(starts_with("employment")) %>%
  apply(., 2, table, useNA = "ifany")
with(aa_adrd_tte, table(employ_working, employ_retired, useNA = "ifany"))
with(aa_adrd_tte, table(employ_working, employ_not_workforce, useNA = "ifany"))
with(aa_adrd_tte, table(employ_not_workforce, employ_retired, useNA = "ifany"))


aa_adrd_tte %<>%
  mutate(current_work_status = case_when(
    employ_working == 1 ~ 1,
    employ_retired == 1 ~ 2,
    employ_not_workforce == 1 ~ 3,
    rowSums(across(starts_with("employ_")), na.rm = T) == 0 ~ NA_real_))

with(aa_adrd_tte, table(current_work_status, useNA = "ifany"))

#---- Smoking status ----
# Ever/never smoker
aa_adrd_tte %<>%
  mutate(smoking_2 = case_when(smoking_status %in% c(2, 3) ~ 1,
                               smoking_status == 1 ~ 0))
#---- Censor 90 + ----
# We are using the imputed 90 + age and change the censored 90+ type to admin censored
aa_adrd_tte %<>%
  mutate(main_dem_v1_end_type =
           case_when(main_dem_v1_end_type == "CENSORED 90+" ~ "ADMIN CENSORED",
                     TRUE ~ main_dem_v1_end_type))

#---- Event variables ----
aa_adrd_tte %<>% 
  mutate(
    end_type = case_when(
      main_dem_v1_end_type == "DEMENTIA" ~ "Dementia",
      main_dem_v1_end_type == "ADMIN CENSORED" ~ "Administratively Censored",
      main_dem_v1_end_type == "DEATH" ~ "Death",
      main_dem_v1_end_type == "END OF MEMBERSHIP" ~ "End of Membership"),
    event_endmem = ifelse(end_type == "End Of Membership", 1, 0),
    event_death = ifelse(end_type == "Death", 1, 0),
    event_dem = ifelse(end_type == "Dementia", 1, 0))

#---- Filter participants ----
# n = 4535 not in Main dementia definition sample
aa_adrd_tte %>% filter(!(main_dem_v1_sample == 1)) %>% nrow()
aa_adrd_tte %<>% filter(main_dem_v1_sample == 1) # n = 180394

# everyone aged 60+ at survey
aa_adrd_tte %<>% mutate(survey_age_r = round(survey_age))
aa_adrd_tte %>% filter(!survey_age_r >= 60) %>% nrow()

# n = 20917 not in main Asian/white group: South Asian, Chinese, Japanese, Filipino, White
with(aa_adrd_tte, table(ethnicity_rev, main_dem_v1_end_dem_flag, useNA = "ifany"))
table(aa_adrd_tte$ethnicity_rev, useNA = "ifany")
aa_adrd_tte %>% filter(!ethnicity_rev %in% c(1, 2, 3, 5, 9)) %>% nrow()
aa_adrd_tte_selected <- aa_adrd_tte %>% 
  filter(ethnicity_rev %in% c(1, 2, 3, 5, 9)) %>% # n = 159477
  mutate(ethnicity_rev = factor(ethnicity_rev,
                                levels = c(2, 5, 3, 1, 9),
                                labels = c("Chinese", "Filipino", "Japanese", 
                                           "South Asian", "Non-Latino White")))

#---- Save the data ----
save(aa_adrd_tte, aa_adrd_tte_selected,
     file = paste0(path_to_analytic_data, "analysis_data/aa_adrd_tte.RData"))

#---- Check the 90+ imputation for women ----
# longest follow up for women is < 13 years
# check the fu year
with(aa_adrd_tte_selected %>% 
       filter(female == 1, 
              main_dem_v1_fu_time > 13, main_dem_v1_end_type == "DEMENTIA"),
     tapply(main_dem_v1_fu_time, ethnicity_rev, summary))

aa_adrd_tte_selected %>% 
  filter(female == 1, 
         main_dem_v1_fu_time > 13) %>% 
  ggplot(aes(x = main_dem_v1_fu_time)) +
  geom_density(aes(color = ethnicity_rev))


with(aa_adrd_tte_selected %>% 
       filter(female == 1, 
              main_dem_v1_fu_time > 13), 
     table(ethnicity_rev, main_dem_v1_end_type), useNA = "ifany")

aa_adrd_tte_selected %>% 
  filter(female == 1, 
         main_dem_v1_fu_time > 13) %>%
  ggplot(aes(x = main_dem_v1_end_age, color = ethnicity_rev)) +
  geom_density(aes(y = after_stat(count)))

aa_adrd_tte_selected %>% 
  filter(female == 1, 
         main_dem_v1_fu_time > 13, ethnicity_rev != "White") %>%
  ggplot(aes(x = main_dem_v1_end_age, color = ethnicity_rev)) +
  geom_density(aes(y = after_stat(count)))

with(aa_adrd_tte_selected %>% 
       filter(female == 1), table(ethnicity_rev))

with(aa_adrd_tte_selected %>% 
       filter(female == 1, 
              main_dem_v1_fu_time >= 13, main_dem_v1_end_type == "DEMENTIA"), 
     hist(main_dem_v1_fu_time))


with(aa_adrd_tte_selected %>% filter(female == 1, main_dem_v1_end_age > 90), 
     table(main_dem_v1_90flag))

aa_adrd_tte_selected %>% filter(main_dem_v1_90flag == 0) %>% 
  select(survey_age, contains("main_dem_v1")) %>%
  view()

with(aa_adrd_tte_selected %>% filter(main_dem_v1_90flag == 0, main_dem_v1_end_age > 90) %>% 
       select(survey_age, contains("main_dem_v1")), table(main_dem_v1_end_type))
