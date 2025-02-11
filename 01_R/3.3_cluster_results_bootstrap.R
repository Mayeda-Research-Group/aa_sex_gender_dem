# Bootstrap results (from cluster)
# Created by: Yingyan Wu
# Sep.18.2024
# Update log:
# Sep.26: Update with no IPCW weighted results and include dementia free death results

#---- Package loading + options ----
rm(list = ls())
if (!require("pacman")) {
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}
p_load("here", "tidyverse", "broom", "magrittr", "ggh4x")

# No scientific notation
options(scipen = 999)

# Paths
source(here::here("01_R", "0_paths.R"))

#---- Load the results ----
AJ_results <- tibble()
warnings <- tibble()
for (s in 1:11){
  AJ_results <- bind_rows(
    AJ_results, readRDS(paste0(path_to_analytic_data, 
                               "hoffman2/AJ_bootstrap/",
                               "bootstrap_AJ_scenario_", s, ".RDS")))
  warnings <- bind_rows(
    warnings, readRDS(paste0(path_to_analytic_data, 
                             "hoffman2/warnings/",
                             "warning_summary_scenario_", s, ".RDS")))
}

AJ_results %<>%
  mutate(ethnicity_rev = factor(ethnicity_rev, 
                                levels = c("Chinese", "Filipino", "Japanese", 
                                           "South Asian", "Non-Latino White")))

#---- Warnings ----
warnings %>%
  select(starts_with("w_")) %>%
  apply(., 2, table, useNA = "ifany")

# No warnings for crude AJ estimates
# Weighting were not run at this time but the placeholder warning numbers were collected as 0s

#----- RR RD ----
# There are NANs due to dividing by 0. -999 or 999 were randomly assigned 
# so that they lie on the extreme ends of the bootstrap distribution
# There are also Infs due to 0 divide by 0. -999 or 999 will be randomly assigned
# so that they lie on the extreme ends of the bootstrap distribution
# There are also NAs due to the bootstrapping, and cases were not found at that age
# (YW checked with the raw output by bootstrapping for several i_boots)
# na.rm = T for now. 

AJ_results_edit <- AJ_results %>% 
  rowwise() %>%
  mutate(across(starts_with(c("cuminc", "RD", "RR")), 
                ~ifelse(is.nan(.x) | is.infinite(.x), 
                        sample(c(-999, 999), size = 1, replace = TRUE), .x))) %>%
  group_by(ethnicity_rev, state) %>%
  slice_head(n = 1001) %>% # slice on first 1000 bootstraps and one point estimate
  ungroup()

# No Inf or NaNs for RR at age 75, 80, 85, 90 and 95
AJ_results %>%
  select(ethnicity_rev, paste0("RR_age_", seq(75, 95, by = 5))) %>%
  group_by(ethnicity_rev) %>%
  summarise_all(~sum(is.na(.x) | is.infinite(.x) | is.nan(.x)))

AJ_CI <- AJ_results_edit %>% filter(i_boot != 0) %>%
  select(-i_boot, -scenario_num) %>%
  group_by(ethnicity_rev, state) %>%
  reframe(across(everything(), ~quantile(.x, c(0.025, 0.975), na.rm = T))) %>%
  mutate(estimate =  rep(c("p2.5th", "p97.5th"), 
                         length(unique(AJ_results_edit$ethnicity_rev)) *
                           length(unique(AJ_results_edit$state)))) %>%
  ungroup() %>%
  pivot_longer(cols = -c(estimate, ethnicity_rev, state), 
               names_to = c(".value", "age"),
               names_pattern = "(.*)_age_(.*)") %>%
  mutate(age = as.numeric(age)) %>%
  arrange(age, estimate)

AJ_pe <- AJ_results_edit %>% filter(i_boot == 0) %>%
  mutate(estimate = "pe") %>%
  pivot_longer(cols = -c(estimate, ethnicity_rev, state), 
               names_to = c(".value", "age"),
               names_pattern = "(.*)_age_(.*)") %>%
  mutate(age = as.numeric(age)) %>%
  arrange(age, estimate)

AJ_pe_CI <- bind_rows(AJ_pe, AJ_CI)

#---- Format the results ----
AJ_75_95 <- AJ_pe_CI %>%
  filter(age %in% seq(75, 95, by = 5)) %>%
  # Change cuminc and RD to 1 decimal but RR to 2 decimal
  mutate(across(contains(c("cuminc", "RD")), \(x) sprintf(x, fmt = '%#.1f')),
         RR = sprintf(RR, fmt = '%#.2f')) %>%
  # mutate(across(`cuminc_female_1`:`RD`, \(x) sprintf(x, fmt = '%#.2f'))) %>%
  pivot_wider(names_from = estimate, values_from = `cuminc_female_1`:`RD`) %>%
  mutate(cuminc_Women_CI = paste0(cuminc_female_1_pe, " (", cuminc_female_1_p2.5th,
                                  ", ", cuminc_female_1_p97.5th, ")"),
         cuminc_Men_CI = paste0(cuminc_female_0_pe, " (", cuminc_female_0_p2.5th,
                                ", ", cuminc_female_0_p97.5th, ")"),
         RR_CI = paste0(RR_pe, " (", RR_p2.5th, ", ", RR_p97.5th, ")"),
         RD_CI = paste0(RD_pe, " (", RD_p2.5th, ", ", RD_p97.5th, ")")) %>%
  pivot_wider(names_from = `state`, values_from = `cuminc_female_1_pe`:`RD_CI`,
              names_glue = "{state} {.value}") %>%
  select(ethnicity_rev, age, ends_with("CI")) %>%
  relocate(ethnicity_rev, age, starts_with("dementia"), starts_with("death")) %>%
  arrange(ethnicity_rev, age)
  # mutate(across(ends_with("CI"), 
  #               # South Asian, age 95+, < 5 events for both dementia and deaths 
  #               ~case_when(ethnicity_rev == "South Asian" & age == 95 ~
  #                            "< 5 events",
  #                          TRUE ~ .x)))

# Check South Asian 95 estimates
AJ_75_95 %>%
  filter(age == 95, ethnicity_rev == "South Asian")

AJ_75_95_dem <- AJ_75_95 %>% select(-contains("death"))
AJ_75_95_death <- AJ_75_95 %>% select(-contains("dementia"))

writexl::write_xlsx(list("Dementia" = AJ_75_95_dem, 
                         "Death" = AJ_75_95_death), 
                    here::here("03_tables", "raw_tables", 
                               "table_RR_RD_AJ_dem_ethns.xlsx"))

#---- Figure of the results ----
hline_tib <- tibble(scale = c("Risk Ratio", "Risk Difference (%)"), ref = c(1, 0))

RR_upper_limit <- 2.5 # For dementia
RD_upper_limit <- 25 # For dementia

AJ_pe_CI %>%
  filter(age == 95, ethnicity_rev == "South Asian", 
         state == "dementia")

AJ_75_95_toplot <- AJ_pe_CI %>%
  filter(age %in% seq(75, 95, by = 5)) %>%
  select(-contains("cuminc")) %>%
  pivot_longer(RR:RD, names_to = "scale") %>%
  pivot_wider(names_from = estimate, values_from = value) %>%
  mutate(
    # across(c(pe, p2.5th, p97.5th), \(x) 
    #        # drop SA's 95 estimates as < 5 events for dementia and death
    #        case_when(age == 95 & ethnicity_rev == "South Asian" ~ NA_real_,
    #                  TRUE ~ x)),
    
    # markpos = case_when(
    #   state == "dementia" &
    #     scale == "RR" & age %in% c(90, 95) & ethnicity_rev == "South Asian" ~ 1.7,
    #   state == "dementia" &
    #     scale == "RD" & age %in% c(90, 95) & ethnicity_rev == "South Asian" ~ 10,
    #   state == "death" &
    #     scale == "RR" & age %in% c(90, 95) & ethnicity_rev == "South Asian" ~ 0.8,
    #   state == "death" &
    #     scale == "RD" & age %in% c(90, 95) & ethnicity_rev == "South Asian" ~ -5),
    upper = case_when(scale == "RR" & p97.5th > RR_upper_limit ~ NA_real_,
                      scale == "RD" & p97.5th > RD_upper_limit ~ NA_real_,
                      TRUE ~ p97.5th),
    upper_arrow_pos = case_when(
      scale == "RR" & p97.5th > RR_upper_limit ~ RR_upper_limit,
      scale == "RD" & p97.5th > RD_upper_limit ~ RD_upper_limit),
    scale = case_when(scale == "RD" ~ "Risk Difference (%)",
                      scale == "RR" ~ "Risk Ratio"))

color_palette <- c("#7B5AA3", "#ED9DB2", "#54B663", "#B69833", "#4B4B4B")
for (s in c("dementia", "death")){
  p <- AJ_75_95_toplot %>%
    filter(state == s) %>%
    ggplot(aes(x = age, y = pe, color = ethnicity_rev)) +
    # alpha = age)) +
    geom_point() +
    # # Do not show anything for SA's 90 and 95 estimates
    # geom_point(aes(y = markpos), shape = 4, show.legend = FALSE) + 
    geom_errorbar(aes(ymin = p2.5th, ymax = upper), width = 1) +
    geom_segment(aes(x = age, xend = age, y = p2.5th, yend = upper_arrow_pos),
                 arrow = arrow(length = unit(0.2, "cm"))) +
    geom_text(aes(x = age, y = upper_arrow_pos * 1.05, 
                  label = sprintf(p97.5th, fmt = '%#.2f')), size = 3) +
    facet_grid(rows = vars(scale),
               cols = vars(ethnicity_rev),
               scales = "free", switch = "y") +
    # scale_y_continuous(
    #   breaks = case_when(str_detect(scale , "Ratio") ~ seq(1, 2.6, by = 0.2), 
    #                     str_detect(scale , "Difference") ~ seq(0, 18, by = 2)),
    #   position = "right") +
    scale_x_continuous(breaks = seq(75, 95, by = 5), limits = c(74, 96.5)) +
    ggh4x::facetted_pos_scales(
      y = case_when(
        s == "dementia" ~ 
          list(scale_y_continuous(breaks = seq(-4, 25, by = 2), position = "right"),
               # log scale for RR
               scale_y_continuous(transform = "log", 
                                  limits = c(0.4, 2.6), 
                                  breaks = seq(0.4, 2.6, by = 0.2),
                                  minor_breaks = seq(1, 2.6, by = 0.2),
                                  position = "right")),
        s == "death" ~ 
          list(scale_y_continuous(breaks = seq(-25, 5, by = 2), position = "right"),
               # log scale for RR
               scale_y_continuous(transform = "log", 
                                  limits = c(0.35, 1.6),
                                  breaks = seq(0.4, 1.6, by = 0.2), 
                                  position = "right")))) +
    geom_hline(data = hline_tib, aes(yintercept = ref), linetype = "dashed",
               color = "grey30") +
    scale_color_manual(values = color_palette) + 
    theme_bw() +
    labs(x = "Age", y = element_blank()) +
    guides(color = "none", alpha = "none")
  assign(paste0("p_", s), p)
}

p_dementia
p_death

ggsave(p_dementia, file = here::here("02_figs", "RD_RR_dem_crude_all_ages.png"), 
       device = "png", width = 9, height = 7, units = "in", dpi = 300)
ggsave(p_death, file = here::here("02_figs", "RD_RR_death_crude_all_ages.png"),
       device = "png", width = 9, height = 7, units = "in", dpi = 300)

#----- OLD ----
# # Check NAN and NAs in the results for RD and RR
# AJ_unwtd_results <- AJ_results %>%
#   filter(weight == "Crude")
# 
# colSums(is.na(AJ_unwtd_results))
# missing_RD <- 
#   AJ_unwtd_results %>%
#   mutate(across(contains("RD"), ~case_when(.x == 0 ~ NA_real_,
#                                            TRUE ~ .x))) %>%
#   select(ethnicity_rev, contains(c("RD"))) %>%
#   group_by(ethnicity_rev) %>%
#   summarise_all(~sum(is.na(.x)))
# 
# missing_RR <- AJ_unwtd_results %>%
#   select(ethnicity_rev, contains(c("RR"))) %>%
#   group_by(ethnicity_rev) %>%
#   summarise_all(~sum(is.na(.x)))
# 
# diffdf::diffdf(missing_RD %>% set_colnames(colnames(missing_RR)), missing_RR)
# # No diff. RD 0s are from 0 cumincs and that will give RR NAs
# 
# AJ_crude_CI <- AJ_results_edit %>% filter(i_boot != 0, weight == "Crude") %>%
#   select(-i_boot, -scenario_num, -weight) %>%
#   group_by(ethnicity_rev, weight) %>%
#   reframe(across(everything(), ~quantile(.x, c(0.025, 0.975), na.rm = T))) %>%
#   mutate(estimate =  rep(c("p2.5th", "p97.5th"), 
#                          length(unique(AJ_results_crude_edit$ethnicity_rev)))) %>%
#   pivot_longer(cols = -c(estimate, ethnicity_rev), 
#                names_to = c(".value", "age"),
#                names_pattern = "(.*)_age_(.*)") %>%
#   mutate(age = as.numeric(age)) %>%
#   arrange(age, estimate)
# 
# AJ_crude_pe <- AJ_results_edit %>%
#   filter(i_boot == 0, weight == "Crude") %>%
#   mutate(estimate = "pe") %>%
#   pivot_longer(cols = -c(estimate, ethnicity_rev), 
#                names_to = c(".value", "age"),
#                names_pattern = "(.*)_age_(.*)") %>%
#   mutate(age = as.numeric(age)) %>%
#   arrange(age, estimate)
# 
# AJ_crude_results <- bind_rows(AJ_crude_pe, AJ_crude_CI)