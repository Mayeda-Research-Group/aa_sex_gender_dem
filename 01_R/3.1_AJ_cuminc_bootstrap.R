# cumulative incidence for dementia and dementia free death -- cluster script
# Created by: Yingyan Wu
# Data: 09.17.2024

#---- Argument ----
# Read in the arguments listed in the:
# R CMD BATCH --no-save --no-restore "--args scenario_num=$SGE_TASK_ID"  
## expression:
args = (commandArgs(TRUE))

# Check to see if arguments are passed and set default values if not.
# If so, parse the arguments. (I only have one argument here.)
if (length(args) == 0) {
  print("No arguments supplied.")
  ##supply default values
  scenario_num <- 1
} else {
  for (i in 1:length(args)) {
    eval(parse(text = args[[i]])) # parse each argument (only have 1 here)
  }
}

#---- Package loading ----
library(tidyverse)
library(magrittr)
library(survival)
library(broom)
# No scientific notation
options(scipen = 999)

#---- Load the data ----
# Paths (local)
# source(here::here("01_R", "0_paths.R"))
# load(paste0(path_to_analytic_data, "analysis_data/aa_adrd_tte.RData"))

# Paths cluster
project_folder <- paste0("/u/home/y/yingyanw/aa_sex_gender_dem/")
load(paste0(project_folder, "data/aa_adrd_tte.RData"))

#---- Prepare the data ----
aa_adrd_tte_selected %<>% 
  mutate(
    # factorize the strata variables
    female = as.factor(female),
    usaborn_rev = as.factor(usaborn_rev),
    # Set missing as a level for education for now
    edu_4_miss = case_when(is.na(edu_4) ~ "Missing",
                           TRUE ~ as.character(edu_4)),
    # create multiple events for dementia and death for competing event analysis
    event_multi =  factor(case_when(event_dem == 1 & event_death == 0 ~ 1,
                                    event_death == 1 ~ 2, # dementia free death
                                    TRUE ~ 0),
                          labels = c("s(0)", 
                                     "dementia", "death")),
    starttime = 0)

#---- Cumulative incidence using Aalen-Johansen estimator (left truncation, competing risk) ----
AJ_boot <- function(data, AJ_form, timescale = "age", 
                    weight = FALSE, bootstrap = TRUE){
  # data = aa_adrd_tte_selected
  # AJ_form = Surv(survey_age, main_dem_v1_end_age, event_multi, type = "mstate") ~
  #   strata(ethnicity_rev, female)
  # # # , usaborn_rev)
  # weight = FALSE
  # timescale = "age"
  # i_boot = 1
  # bootstrap = i_boot != 0
  
  ##---- 0. Data preparation ----
  # Filter according to race/ethnicity groups
  # Obtain sampling with replacement ids within the specific race/ethnicity groups
  if (bootstrap == TRUE) {
    set.seed(i_boot)
    boot_ids <- sample(unique(data$subjid), 
                       length(unique(data$subjid)), 
                       replace = TRUE)
  } else {
    boot_ids <- unique(data$subjid)
  }
  
  data_boot <- tibble(subjid = boot_ids,
                      newid = 1:length(boot_ids)) %>% 
    left_join(data, by = "subjid", relationship = "many-to-many") %>%
    select(-subjid)
  
  ##---- 1. Weight ----
  # Weight for education
  if (weight == TRUE) {
    tryCatch({
      female_num_mod <<- glm(female ~ ethnicity_rev,
                             family = binomial(link = "logit"),
                             data = data_boot)
      w_wt_num_1 <<- c(w_wt_num_1, 0)
    }, warning = function(w) {
      w_wt_num_1 <<- c(w_wt_num_1, 1)
      female_num_mod <<- glm(female ~ ethnicity_rev,
                             family = binomial(link = "logit"),
                             data = data_boot)
      
    })
    
    tryCatch({
      female_denom_mod <<- glm(as.formula(weight_edu_formula),
                               family = binomial(link = "logit"),
                               data = data_boot)
      w_wt_num_1 <<- c(w_wt_num_1, 0)
    }, warning = function(w) {
      w_wt_denom_1 <<- c(w_wt_denom_1, 1)
      female_denom_mod <<- glm(as.formula(weight_edu_formula),
                               family = binomial(link = "logit"),
                               data = data_boot)
      
    })
    
    data_boot %<>%
      mutate(pred_num = predict(female_num_mod , newdata = data_boot, 
                                type = "response"),
             pred_denom = predict(female_denom_mod, newdata = data_boot, 
                                  type = "response"),
             pmnum = case_when(female == 1 ~ pred_num,
                               female == 0 ~ 1 - pred_num),
             pmdenom = case_when(female == 1 ~ pred_denom,
                                 female == 0 ~ 1 - pred_denom),
             ipmw_unstblz = 1/pmdenom,
             ipmw_stblz = pmnum/pmdenom)
  }
  
  ##---- 2. AJ fit ----
  if (weight == FALSE){
    tryCatch({
      AJ_fit <<- survfit(formula = AJ_form, data = data_boot,
                         cluster = newid, id = newid)
      w_unw_AJ_1 <<- c(w_unw_AJ_1, 0)
    }, warning = function(w) {
      w_unw_AJ_1 <<- c(w_unw_AJ_1, 1)
      AJ_fit <<- survfit(formula = AJ_form, data = data_boot,
                         cluster = newid, id = newid)
    })
  } else if (weight == TRUE){
    tryCatch({
      AJ_fit <<- survfit(formula = AJ_form, data = data_boot,
                         cluster = newid, id = newid, weights = ipmw_stblz)
      w_wtd_AJ_1 <<- c(w_wtd_AJ_1, 0)
    }, warning = function(w) {
      w_wtd_AJ_1 <<- c(w_wtd_AJ_1, 1)
      AJ_fit <<- survfit(formula = AJ_form, data = data_boot,
                         cluster = newid, id = newid, weights = ipmw_stblz)
    })
  }
  
  ##---- 3. Tidy the AJ fit ----
  strata_temp <- gsub(".*strata\\((.*)\\).*", "\\1", AJ_form)[3]
  strata <- strsplit(strata_temp, ", ") %>% unlist()
  AJ_1y_tib <-
    broom::tidy(AJ_fit) %>%
    dplyr::rename(!!sym(timescale) := time) %>%
    mutate(strata = gsub(".*strata\\(.*\\)=(.*)", "\\1", strata)) %>%
    separate(strata, sep = ", ", into = strata) %>%
    # YW 09.26.2024: add dementia free death to the table
    filter(state %in% c("dementia", "death")) %>%
    mutate(!!sym(paste0(timescale, "_cut")) :=
             # Create age/time cut as [t, t+1)
             cut(!!sym(timescale),
                 breaks = case_when(timescale == "age" ~ list(seq(60, 103, 1)),
                                    timescale == "time" ~ list(seq(0, 19, 1)))[[1]],
                 labels = case_when(timescale == "age" ~ list(seq(60, 102, 1)),
                                    timescale == "time" ~ list(seq(0, 18, 1)))[[1]],
                 right = FALSE)) %>%
    group_by_at(c(strata, "state", paste0(timescale, "_cut"))) %>%
    # Only include the first row of observation available at each [t, t+1)
    # As long as participants are not at t+1 age, they are considered as aged t
    slice_head() %>%
    ungroup() %>%
    select(paste0(timescale, "_cut"), cuminc = estimate, all_of(strata), state) %>%
    mutate(cuminc = cuminc*100) %>%
    pivot_wider(names_from = female,
                values_from = c("cuminc"),
                names_prefix = "cuminc_female_") %>%
    mutate(RR = cuminc_female_1/cuminc_female_0,
           RD = cuminc_female_1 - cuminc_female_0) %>%
    pivot_wider(values_from = c("cuminc_female_1", "cuminc_female_0", "RR", "RD"),
                names_from = paste0(timescale, "_cut"),
                names_glue = paste0("{.value}_", timescale, "_{", timescale, 
                                    "_cut}"))
  
  return(AJ_1y_tib)
}

#---- Scenario tibble ----
n_bootstrap_per_job = 100
n_bootstrap = 1100
scenario_tib <- expand_grid(
  boot_i = c(paste0(seq(1, n_bootstrap - n_bootstrap_per_job + 1, 
                        by = n_bootstrap_per_job), "-",
                    seq(n_bootstrap_per_job, n_bootstrap, 
                        by = n_bootstrap_per_job)))) %>%
  separate(boot_i, into = c("boot_i_start", "boot_i_end"), sep = "-") %>%
  mutate_at(vars(boot_i_start), ~ifelse(boot_i_start == 1, 0, .))

# Make sure that scenario_num is not greater than all the scenario index
if (scenario_num > nrow(scenario_tib)) {
  print("out of bound senario number")
  # supply default value
  scenario_num <- 1
}

# print values just to make sure:
print(scenario_num)
scenario <- scenario_tib[scenario_num, ]
boot_i_start <- scenario$boot_i_start
boot_i_end <- scenario$boot_i_end

print(boot_i_start)
print(boot_i_end)

w_wt_num_1 <- w_wt_denom_1 <- w_unw_AJ_1 <- w_wtd_AJ_1 <- numeric()

##---- Model formulas ----
AJ_formula <- Surv(survey_age, main_dem_v1_end_age, event_multi, type = "mstate") ~
  strata(ethnicity_rev, female)
# # Time as time scale
# AJ_form = Surv(starttime, main_dem_v1_fu_time, event_multi, type = "mstate") ~
#   strata(ethnicity_rev, female)

# # YW 09.26.2024: drop the weight
# weight_edu_formula <- "female ~ edu_4_miss*ethnicity_rev"

#---- Bootstrap ----
# for (i_boot in c(0, 902, 211, 212)){
for (i_boot in boot_i_start:boot_i_end) {
  print(i_boot)
  # boot_i_start <- 0
  # i_boot <- 2 # point estimate, using full dataset, no resampling
  # i_boot <- 902 # and so on: bootstrap
  # i_boot <- 212 # and so on: bootstrap
  start <- Sys.time() # record starting time
  
  crude_tib_i <- AJ_boot(aa_adrd_tte_selected, AJ_form = AJ_formula, 
                         timescale = "age",
                         weight = FALSE, bootstrap = i_boot != 0)
  
  res_boot_i <- crude_tib_i %>%
    mutate(i_boot = i_boot,
           scenario_num = scenario_num)
  
  # # YW 09.26.2024: drop the weighted estimate
  # wtd_tib_i <- AJ_boot(aa_adrd_tte_selected, AJ_form = AJ_formula, 
  #                      timescale = "age",
  #                      weight = TRUE, bootstrap = i_boot != 0)
  # res_boot_i <- rbind(crude_tib_i %>% mutate(weight = "Crude"),
  #                     wtd_tib_i %>% mutate(weight = "Weighted")) %>%
  #   relocate(ethnicity_rev, weight) %>%
  #   mutate(i_boot = i_boot,
  #          scenario_num = scenario_num)
  
  end <- Sys.time() # record ending time
  
  w_boot_i <- tibble(
    i_boot = i_boot,
    scenario = scenario_num,
    duration = difftime(end, start, units = 'mins'),
    w_wt_num_1 = sum(w_wt_num_1), 
    w_wt_denom_1= sum(w_wt_denom_1),
    w_unw_AJ_1 = sum(w_unw_AJ_1),
    w_wtd_AJ_1 = sum(w_wtd_AJ_1))
  
  # collect output for a specific ethnicity in a given scenario
  if (i_boot == boot_i_start) {
    out <- res_boot_i
    warnings <- w_boot_i
  } else {
    out <- rbind(out, res_boot_i)
    warnings <- rbind(warnings, w_boot_i)
  }
}

# res_col <- colnames(res_boot_i)
# out_col <- colnames(out)
# out_col[which(!out_col %in% res_col)]

saveRDS(out, paste0(project_folder, "AJ_bootstrap/",
                    paste0("bootstrap_AJ_scenario_", scenario_num, ".RDS")))

saveRDS(warnings, paste0(project_folder, "warnings/",
                         paste0("warning_summary_scenario_", scenario_num, ".RDS")))
