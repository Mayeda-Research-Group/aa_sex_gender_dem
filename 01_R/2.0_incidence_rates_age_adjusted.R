# Incidence rate age adjusted and age specified
# Created by: Yingyan Wu
# Creation date: 9.2.2024

#---- Package loading + options ----
rm(list = ls())

if (!require("pacman")) {
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}
p_load("here", "tidyverse", "dplyr","formattable", "writexl", "readxl", "gridExtra",
       "magrittr")

#No scientific notation
options(scipen = 999)

# ----load the dataset & functions ----
source(here::here("01_R", "0_paths.R"))
load(paste0(path_to_analytic_data, "analysis_data/aa_adrd_tte.RData"))

#---- Incidence rate function ----
# modified function from AA diabetes and AA APOE project
calculate_pys <- function(data, age_cat, fu_start, fu_end, event_flag) {
  
  age_cat_labels <- c()
  for (i in 1:(length(age_cat) - 1)) {
    # initialize the age interval
    cat_start <- age_cat[i]
    cat_end <- age_cat[i+1]
    # create the variable names for py and case contribution
    age_int <- paste0(cat_start, "_", cat_end)
    vars <- paste0(c("py_", "case_"), age_int)
    py_var <- vars[1]
    case_var <- vars[2]
    
    # collect the age intervals for setting up the IR output table
    age_cat_labels <- c(age_cat_labels, age_int)
    
    # calculate py and case contribution by age category
    data <- data %>% 
      mutate(
        !! py_var := case_when(
          # start of fu after age category or end of fu before age category: 
          # contribute 0 pys
          get(fu_start) > cat_end | get(fu_end) <= cat_start ~ 0, 
          
          # start of fu before age category and end of fu during age category: 
          get(fu_start) <= cat_start & get(fu_end) <= cat_end ~ 
            get(fu_end) - cat_start, 
          
          # start of fu before age category and end of fu after age category: 
          get(fu_start) <= cat_start & get(fu_end) > cat_end ~ 
            cat_end - cat_start, 
          
          # start and end of fu during age category:
          get(fu_start) <= cat_end & get(fu_end) <= cat_end ~ 
            get(fu_end) - get(fu_start), 
          
          # start of fu during age category and end of fu after age cateogory:
          get(fu_start) <= cat_end & get(fu_end) > cat_end ~ 
            cat_end - get(fu_start), 
          
          TRUE ~ NA_real_),  
        
        # if end of fu is during age category, case contribution is the dem flag. 
        # otherwise, case contribution is zero
        !! case_var := ifelse(
          get(fu_end) <= cat_end & get(fu_end) > cat_start, 
          get(event_flag), 0
        )
      )
  }
  
  return(list(data = data, age_cat_labels = age_cat_labels))
}

calculate_ir <- function(data, age_cat_labels, std_pop, exposure_treated, 
                         per_pys = 1000) {
  # # # testing
  # # data <- sex_gender_pys$data
  # data <- sex_gender_pys_death$data
  # age_cat_labels <- sex_gender_pys$age_cat_labels
  # std_pop <- c(16817924, 12435263, 9278166, 7317795, 5743327, 5493433)
  # # using 2010 US census population as the standard population
  # # set up age-specific IR table
  # per_pys = 1000
  # exposure_treated <- c("female", 1)
  
  data <- data %>%
    filter(get(exposure_treated[1]) == exposure_treated[2]) %>%
    split(., .$ethnicity_rev)
  
  out <- lapply(data, function(x) {
    # #Test
    # x <- data
    # using 2000 US census population as the standard population
    # set up age-specific IR table
    age_spec_ir <- tibble(
      age_range = age_cat_labels, 
      pop_count = std_pop, # std pop - count by age category
      pop_total = sum(pop_count), # std pop total 
      pop_prop = pop_count / pop_total, # std pop prop by age category
      unw_pys = NA, 
      unw_cases = NA)
    
    for (i in 1:length(age_cat_labels)) {
      # sum up pys and cases by age category
      cat <- age_cat_labels[i]
      age_spec_ir[i, "unw_pys"] <- sum(x[, paste0("py_", cat)]) 
      age_spec_ir[i, "unw_cases"] <- sum(x[, paste0("case_", cat)])
      
    }
    
    # age specific IRs
    age_spec_ir <- age_spec_ir %>% 
      mutate(
        unw_adj_ir = unw_cases / unw_pys,
        unw_adj_ir_se = sqrt(unw_cases / unw_pys^2),
        unw_adj_ir_CI_l = unw_adj_ir - 1.96 * unw_adj_ir_se,
        unw_adj_ir_CI_u = unw_adj_ir + 1.96 * unw_adj_ir_se
      ) 
    
    # age adjusted IR
    ir_out <- age_spec_ir %>% 
      mutate(
        # weight with standard pop prop
        unw_adj_ir_pop = unw_cases / unw_pys * pop_prop,
        unw_adj_ir_pop_var = unw_cases / unw_pys^2 * pop_prop^2,
        unw_pys_squared = unw_pys^2) %>% 
      dplyr::summarise(
        across(c(unw_cases, unw_pys, unw_pys_squared), sum), 
        crude_ir = unw_cases/unw_pys,
        crude_ir_var = unw_cases/unw_pys_squared,
        adj_ir = sum(unw_adj_ir_pop),
        adj_ir_var = sum(unw_adj_ir_pop_var)) %>% 
      mutate(
        adj_ir_se = sqrt(adj_ir_var),
        adj_ir_l = adj_ir - 1.96 * adj_ir_se,
        adj_ir_u = adj_ir + 1.96 * adj_ir_se,
        crude_ir_se = sqrt(crude_ir_var),
        crude_ir_l = crude_ir - 1.96 * crude_ir_se,
        crude_ir_u = crude_ir + 1.96 * crude_ir_se) 
    
    # format the numbers to 1 decimal point and per 1000 pys
    out_1 <- age_spec_ir %>% 
      select(age_range, unw_cases, unw_pys,
             unw_adj_ir, unw_adj_ir_CI_l, unw_adj_ir_CI_u) %>%
      mutate(across(contains("ir"), function(x) x*per_pys))
    
    out_2 <- ir_out %>% 
      select(unw_cases, unw_pys, 
             crude_ir, crude_ir_l, crude_ir_u, adj_ir, adj_ir_l, adj_ir_u) %>% 
      mutate(across(contains("ir"), function(x) x*per_pys))
    
    return(list("age_spec" = out_1, "age_adj" = out_2))
  }
  )
  
  return(out)
}

#---- IR tables ----
# Applying function to wide dataset
# age_cat <- c(59, 65, 70, 75, 80, 85, 120)
age_cat <- c(59, 75, 80, 85, 90, 95, 120)
sex_gender_pys <- 
  calculate_pys(aa_adrd_tte_selected, age_cat, "survey_age", "main_dem_v1_end_age",
                "main_dem_v1_end_dem_flag")

# Standard population from 2000 census
# Table 1 from: https://www2.census.gov/library/publications/decennial/2000/briefs/c2kbr01-12.pdf
# corresponding to 60-75, 70-75, 75-80, 80-85, 85-90, 90-95, 95+
std_pop <- c(29196433, 7415813, 4945367, 2789818, 1112531, 337238)
# # Standard population from 2000 census
# # corresponding to 60-65, 65-70, 70-75, 75-80, 80-85, 85+
# std_pop <- c(10805447, 9533545, 8857441, 7415813, 4945367, 4239587) 
# # Standard population from 2010 census
# std_pop <- c(16817924, 12435263, 9278166, 7317795, 5743327, 5493433)

IR_temp <- lapply(c(0, 1), function(x) {
  calculate_ir(sex_gender_pys$data,
               sex_gender_pys$age_cat_labels, std_pop, 
               exposure_treated = exposure_treated <- c("female", x),
               per_pys = 1000)
})

names(IR_temp) <- c("male", "female")

# Save
save(IR_temp, file = paste0(path_to_analytic_data, "analysis_data/IR_list.RData"))

##---- functions to formatting tables ----
collect_age_adj_ir <- function(IR_data, sex) {
  # # Test
  # IR_data <- IR_temp
  # ethns <- levels(aa_adrd_tte_selected$ethnicity_rev)
  age_adj <- lapply(
    ethns, function(x) {IR_data[[sex]][[x]]$age_adj %>%
        mutate(across(contains("ir"), 
                      function(x) sprintf(x, fmt = "%.1f")))}
  )
  age_adj <- do.call(rbind, age_adj) %>% 
    mutate(
      age_crude_ir = paste0(crude_ir, " (", crude_ir_l, ", ", crude_ir_u, ")"),
      age_adj_ir = paste0(adj_ir, " (", adj_ir_l, ", ", adj_ir_u, ")")) %>%
    rename_all(function(x) paste0(x, "_", sex)) %>%
    mutate(ethn = ethns) %>% 
    relocate(ethn)
  
  return(age_adj)
}

# function to collect age-specific rates and calculate IRR into tables
collect_age_spec_ir_irr <- function(IR_data) {
  # # Test
  # IR_data <- IR_temp
  # ethns <- levels(aa_adrd_tte_selected$ethnicity_rev)
  
  for (sex in names(IR_data)){
    temp <- lapply(ethns, function(x) {
      IR_data[[sex]][[x]]$age_spec %>% 
        mutate(across(contains("ir"), function(x) sprintf(x, fmt = "%.1f"),
                      .names = "{.col}_rounded"),
               age_spec_ir = ifelse(
                 unw_cases < 5, "< 5 events", 
                 paste0(unw_adj_ir_rounded, " (", unw_adj_ir_CI_l_rounded, 
                        ", ", unw_adj_ir_CI_u_rounded, ")")),
               ethn = NA) %>% 
        select(ethn, age_range, unw_cases, unw_adj_ir, unw_pys, age_spec_ir) %>%
        rename_at(vars(contains(c("unw_", "_ir"))), 
                  function(x) paste0(x, "_", sex)) %>%
        add_row(ethn = x, .before = 1)
    }
    )
    assign(paste0("age_spec_", sex), temp)
  }
  
  age_spec <-
    do.call(rbind, age_spec_female) %>%
    cbind(do.call(rbind, age_spec_male) %>% select(-ethn, -age_range)) %>%
    mutate(unw_adj_irr = 
             case_when(unw_cases_male >= 5 & unw_cases_female >= 5 ~ 
                         unw_adj_ir_female/unw_adj_ir_male,
                       TRUE ~ NA_real_),
           se_ln_irr = sqrt(1/unw_cases_male + 1/unw_cases_female),
           unw_adj_irr_l = exp(log(unw_adj_irr) - 1.96*se_ln_irr),
           unw_adj_irr_u = exp(log(unw_adj_irr) + 1.96*se_ln_irr),
           age_spec_irr = sprintf(unw_adj_irr, fmt = "%.1f"),
           age_spec_irr = case_when(
             !is.na(unw_adj_irr) ~ 
               paste0(sprintf(unw_adj_irr, fmt = "%.1f"), " (", 
                      sprintf(unw_adj_irr_l, fmt = "%.1f"), ", ", 
                      sprintf(unw_adj_irr_u, fmt = "%.1f"), ")"),
             is.na(unw_adj_irr) & is.na(age_spec_ir_female) ~ NA,
             TRUE ~ "< 5 events"),
           age_range = case_when(age_range == "59_75" ~ "<= 74 years old",
                                 age_range == "75_80" ~ "75-79 years old",
                                 age_range == "80_85" ~ "80-84 years old",
                                 age_range == "85_90" ~ "85-89 years old",
                                 age_range == "90_95" ~ "90-94 years old",
                                 age_range == "95_120" ~ "95+ years")) %>%
    # age_range = case_when(age_range == "59_65" ~ "60-64 years",
    #                       age_range == "65_70" ~ "65-69 years",
    #                       age_range == "70_75" ~ "70-74 years",
    #                       age_range == "75_80" ~ "75-79 years",
    #                       age_range == "80_85" ~ "80-84 years",
    #                       age_range == "85_120" ~ "85+ years")) %>%
    select(ethn, age_range, contains(c("unw_cases", "unw_pys", "age_spec_ir", 
                                       "age_spec_irr")))
  
  return(age_spec)
}

##---- Format Tables ----
ethns <- levels(aa_adrd_tte_selected$ethnicity_rev)

age_spec_ir_tib <- collect_age_spec_ir_irr(IR_temp) %>%
  set_colnames(c("Race/ethnicity", "Age range", 
                 paste0("Cases ", c("Women", "Men")), 
                 paste0("Person-Years ",  c("Women", "Men")),
                 paste0("IR (95% CI) ",  c("Women", "Men")),
                 "IRR (95% CI)"))

male_age_adj_ir <- collect_age_adj_ir(IR_temp, "male")
female_age_adj_ir <- collect_age_adj_ir(IR_temp, "female")

age_adj_ir_tib <- cbind(female_age_adj_ir, male_age_adj_ir %>% 
                          select(-"ethn")) %>%
  select(ethn, contains("unw_cases"), contains("pys"), 
         contains("age_crude_ir"), contains("age_adj_ir")) %>%
  set_colnames(c("Race/ethnicity", 
                 paste0("Cases ",  c("Women", "Men")), 
                 paste0("Person-Years ",  c("Women", "Men")),
                 paste0("Crude IR (95% CI) ",  c("Women", "Men")),
                 paste0("Age-standardized IR (95% CI) ",  c("Women", "Men"))))

##---- Save the table ----
write_xlsx(list("age_spec" = age_spec_ir_tib, "crude_adj" = age_adj_ir_tib), 
           here::here("03_tables", "raw_tables", 
                      "table_age_ir_femalexethn.xlsx"))
#---- IR Figure ----
male_age_adj_ir <- collect_age_adj_ir(IR_temp, "male") %>%
  mutate(sex = "Men") %>%
  rename_at(vars(contains("male")), ~str_replace(., "_male", "")) %>%
  select(ethn, sex, starts_with("adj_ir"))

female_age_adj_ir <- collect_age_adj_ir(IR_temp, "female") %>%
  mutate(sex = "Women") %>%
  rename_at(vars(contains("male")), ~str_replace(., "_female", "")) %>%
  select(ethn, sex, starts_with("adj_ir"))

age_adj_ir_forplot <- rbind(male_age_adj_ir, female_age_adj_ir) %>%
  mutate(across(c(adj_ir, adj_ir_l, adj_ir_u), as.numeric)) %>%
  mutate(sex = factor(sex, levels = c("Women", "Men")),
         ethn = factor(ethn, levels = c("Chinese", "Filipino", "Japanese", 
                                        "South Asian", "Non-Latino White")))

color_palette <- c("#7B5AA3", "#ED9DB2", "#54B663", "#B69833", "#4B4B4B")
age_adj_ir_forplot %>%
  ggplot(aes(x = ethn, y = adj_ir, 
             fill = ethn, alpha = sex)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(
    stat = "identity",
    aes(ymin = adj_ir_l, ymax = adj_ir_u), width = .2,
    position = position_dodge(.9)) + 
  scale_y_continuous(breaks = seq(0, 30, by = 5)) +
  scale_alpha_discrete(range = c(0.5, 1)) +
  #theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust = 1)) +
  theme_bw() + 
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 11),
    strip.text = element_text(size = 11)
  ) +
  scale_fill_manual(values = color_palette) + 
  labs(x = "Race/Ethnicity", y = "Age-adjusted incidence rate per 1,000 person-years")

ggsave(here::here("02_figs", "fig_age_adj_ir_ethnxsex_2000.png"),
       device = "png", width = 7, height = 5, units = "in", dpi = 300)

#---- Death table -----
sex_gender_pys_death <- 
  calculate_pys(aa_adrd_tte_selected, age_cat, "survey_age", "main_dem_v1_end_age",
                "death_flag")
# # Standard population from 2000 census
# # corresponding to 60-65, 65-70, 70-75, 75-80, 80-85, 85+
# std_pop <- c(10805447, 9533545, 8857441, 7415813, 4945367, 4239587) 
# # Standard population from 2010 census
# std_pop <- c(16817924, 12435263, 9278166, 7317795, 5743327, 5493433)

IR_temp_death <- lapply(c(0, 1), function(x) {
  calculate_ir(sex_gender_pys_death$data,
               sex_gender_pys_death$age_cat_labels, std_pop, 
               exposure_treated = exposure_treated <- c("female", x),
               per_pys = 1000)
})

names(IR_temp_death) <- c("male", "female")
IR_temp_death

age_spec_ir_tib_death <- collect_age_spec_ir_irr(IR_temp_death) 
# South Asian 95+ < 5 events
# South Asian 90-94, female 4 events

# #---- OLD ----
# # Age specific IR table
# male_age_spec_ir <- collect_age_spec_ir(IR_temp, "male")
# female_age_spec_ir <- collect_age_spec_ir(IR_temp, "female")
# 
# age_spec_ir_tib <- cbind(female_age_spec_ir, male_age_spec_ir %>% 
#                            select(-"ethn", -"age_range")) %>%
#   relocate(ethn, age_range, contains("unw_cases"), contains("pys"), contains("ir")) %>%
#   set_colnames(c("Race/ethnicity", "Age range", 
#                  paste0("Cases ", c("Women", "Men")), 
#                  paste0("Person-Years ",  c("Women", "Men")),
#                  paste0("IR (95% CI) ",  c("Women", "Men"))))
# collect_age_spec_ir <- function(IR_data, sex) {
#   # # Test
#   # IR_data <- IR_temp
#   # ethns <- levels(aa_adrd_tte_selected$ethnicity_rev)
#   age_spec <- lapply(
#     ethns, function(x) {
#       IR_data[[sex]][[x]]$age_spec %>% 
#         mutate(across(contains("ir"), 
#                       function(x) sprintf(x, fmt = "%.1f"))) %>%
#         mutate(age_spec_ir = ifelse(
#           unw_cases < 5, "< 5 events", 
#           paste0(unw_adj_ir, " (", unw_adj_ir_CI_l, ", ", unw_adj_ir_CI_u, ")")),
#           ethn = NA) %>% 
#         select(ethn, age_range, unw_cases, unw_pys, age_spec_ir) %>%
#         add_row(ethn = x, .before = 1)
#     }
#   )
#   
#   age_spec <- do.call(rbind, age_spec) %>% 
#     mutate(age_range = case_when(age_range == "59_65" ~ "60-64 years",
#                                  age_range == "65_70" ~ "65-69 years",
#                                  age_range == "70_75" ~ "70-74 years",
#                                  age_range == "75_80" ~ "75-79 years",
#                                  age_range == "80_85" ~ "80-84 years",
#                                  age_range == "85_120" ~ "85+ years")) %>%
#     rename_at(vars(contains(c("unw", "ir"))), function(x) paste0(x, "_", sex))
#   
#   return(age_spec)
# }