# Descriptive statistics for baseline and cumulative incidence
# Created by: Yingyan Wu
# Data: 09.02.2024

#---- Package loading + options ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "readr", "tidyverse", "magrittr", "plyr", "haven", "labelled", 
       "gtsummary", "survival", "broom", "ggsurvfit", "survminer", "cowplot",
       "ggsci")
# Do not load "summarytools" as it contains conflicted functions as dplyr

#No scientific notation
options(scipen = 999)

# ----load the dataset ----
source(here::here("01_R", "0_paths.R"))
load(paste0(path_to_analytic_data, "analysis_data/aa_adrd_tte.RData"))

aa_adrd_tte_selected %<>% 
  mutate(female = factor(female, levels = c(1, 0), labels = c("Women", "Men")),
         edu_4_miss = case_when(is.na(edu_4) ~ "Missing",
                                TRUE ~ as.character(edu_4)))

#----Density plot of baseline age by race/ethnicity and sex/gender----
aa_adrd_tte_selected %>%
  ggplot(aes(x = survey_age, group = female, color = female)) +
  geom_density() +
  facet_grid(~ethnicity_rev, scales = "free_y") +
  scale_color_manual(name = NULL,
                     values = c("#3271AE", "#E5A84B"),
                     labels = c("Women", "Men")) +
  labs(
    # title = "Density plot of baseline age stratified by sex/gender and race/ethnicity",
    x = "Age",
    y = "Frequency",
    color = "Sex/gender") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

ggsave(here::here("02_figs", "fig_baseline_age_ethn.png"),
       device = "png", width = 7, height = 5, units = "in", dpi = 300)

# Density plot of end age
aa_adrd_tte_selected %>%
  ggplot(aes(x = main_dem_v1_end_age, group = female, color = female)) +
  geom_density() +
  facet_grid(~ethnicity_rev, scales = "free_y") +
  scale_color_manual(name = NULL,
                     values = c("#3271AE", "#E5A84B"),
                     labels = c("Women", "Men")) +
  labs(
    # title = "Density plot of baseline age stratified by sex/gender and race/ethnicity",
    x = "Age",
    y = "Frequency",
    color = "Sex/gender") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

# Counts of end ages over 90 by gender and race
aa_adrd_tte_selected %>%
  select(ethnicity_rev, female, main_dem_v1_end_age) %>%
  filter(main_dem_v1_end_age >= 90) %>%
  mutate(end_age = case_when(round(main_dem_v1_end_age) >= 95 ~ "95 +",
                             TRUE ~ paste0(round(main_dem_v1_end_age)))) %>%
  select(-main_dem_v1_end_age) %>%
  group_by(ethnicity_rev, female) %>%
  count()


#---- count of deaths after dementia ----
aa_adrd_tte_selected %<>%
  mutate(death_after_dem_flag = 
           case_when(main_dem_v1_end_dem_flag == 1 & death_flag == 1 ~ 1,
                     main_dem_v1_end_dem_flag == 1 & death_flag == 0 ~ 0)) 


death_after_dem_tib <- aa_adrd_tte_selected %>%
  group_by(ethnicity_rev, female) %>%
  dplyr::summarize(n = n(),
                   n_dem = sum(main_dem_v1_end_dem_flag),
                   prop_dem = n_dem/n*100,
                   n_death_after_dem = sum(death_after_dem_flag, na.rm = T),
                   prop_death_after_dem = n_death_after_dem/n*100) %>%
  mutate(n_dem_prop = paste0(n_dem, " (", round(prop_dem, 1), "%)"),
         n_death_after_dem_prop = paste0(n_death_after_dem, " (", 
                                         round(prop_death_after_dem, 1), "%)")) %>%
  select(ethnicity_rev, female, n, n_dem_prop, n_death_after_dem_prop) %>%
  set_colnames(c("Race/ethnicity", "Sex/gender", "N", "Dementia cases (%)", 
                 "Death cases after dementia (%)"))

write_xlsx(death_after_dem_tib, here::here("03_tables", "raw_tables", 
                                           "table_death_counts_after_dem.xlsx"))

#---- Cumulative incidence using Aalen-Johansen estimator (left truncation, competing risk) ----
aa_adrd_tte_selected %<>%
  mutate(event_multi =  factor(case_when(event_dem == 1 & event_death == 0 ~ 1,
                                         event_death == 1 ~ 2,
                                         TRUE ~ 0),
                               labels = c("s(0)", 
                                          "dementia", "death")),
         starttime = 0)

# # Check the age
# summary(aa_adrd_tte_selected$main_dem_v1_end_age)
# aa_adrd_tte_selected %>%
#   select(event_multi, main_dem_v1_end_age, survey_age, main_dem_v1_end_dem_age,
#          main_dem_v1_dem_age, death_age, membership_end_age, cens90_surv_age) %>% 
#   filter(event_multi == "dementia", main_dem_v1_end_dem_age < death_age,
#          main_dem_v1_end_dem_age != main_dem_v1_end_age) %>%
#   view()

#---- Age as time scale ----
##---- Crude model ----
crude_AJ_fit <- 
  survfit(Surv(survey_age, main_dem_v1_end_age, event_multi, type = "mstate") ~ 
            strata(ethnicity_rev, female), data = aa_adrd_tte_selected,
          cluster = subjid, id = subjid)

crude_AJ_tib <- broom::tidy(crude_AJ_fit) %>%
  dplyr::rename(age = time) %>%
  mutate(strata = stringr::str_remove(strata, "strata\\(ethnicity_rev, female\\)="),
         ethnicity = factor(str_split_i(strata, ",", 1), 
                            levels = levels(aa_adrd_tte_selected$ethnicity_rev)),
         female = factor(str_split_i(strata, ", ", 2),
                         levels = levels(aa_adrd_tte_selected$female)),
         state = case_when(state == "(s0)" ~ "Alive without dementia",
                           state == "death" ~ "Death without dementia",
                           state == "dementia" ~ "Dementia")) 

##---- Cumulative incidence plot Dementia only ----
crude_AJ_tib %>% 
  # We decided to keep the cumulative incidence for South Asians until 95
  # # truncate at age 90 for South Asians
  # filter(!(age >= 90 & ethnicity == "South Asian")) %>%
  mutate(across(c("estimate", "conf.low", "conf.high"), function(x) x*100)) %>%
  # cut the y axis at 50, top coded conf.high at 50
  mutate(conf.high_tc = case_when(conf.high >= 50 ~ 50,
                                  TRUE ~ conf.high)) %>%
  filter(state == "Dementia", age <= 95) %>%
  ggplot(aes(x = age, y = estimate, group = female)) +
  geom_line(aes(color = female)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high_tc, fill = female), alpha = 0.2) +
  facet_grid(cols = vars(ethnicity), scales = "free", space = "free") +
  scale_color_manual(name = NULL,
                     values = c("#3271AE", "#E5A84B"),
                     labels = c("Women", "Men")) +
  scale_fill_manual(name = NULL,
                    values = c("#3271AE", "#E5A84B"),
                    labels = c("Women", "Men")) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Age (years)",
       y = "Cumulative incidence (%)") 

crude_AJ_tib %>% 
  mutate(across(c("estimate", "conf.low", "conf.high"), function(x) x*100)) %>%
  # cut the y axis at 50, top coded conf.high at 50
  mutate(conf.high_tc = case_when(conf.high >= 50 ~ 50,
                                  TRUE ~ conf.high)) %>%
  filter(state == "Dementia", age <= 95) %>%
  filter(conf.high_tc != conf.high)

# OLD area version
# crude_AJ_tib %>%
#   ggplot(aes(x= age, y = estimate, fill = state)) + geom_area() +
#   facet_grid(rows = vars(female), cols = vars(ethnicity)) +
#   cowplot::theme_cowplot() +
#   ggsci::scale_fill_jco() +
#   theme(legend.position = "bottom") +
#   labs(x = "Age (years)",
#        y = "Cumulative incidence (%)",
#        fill = NULL)

ggsave(here::here("02_figs", "fig_cuminc_AJ_dem_ethns_withCI.png"),
       device = "png", width = 7, height = 5, units = "in", dpi = 300)

# color palette
color_palette <- c("#7B5AA3", "#ED9DB2", "#54B663", "#B69833", "#4B4B4B")
crude_AJ_tib %>%
  filter(state == "Dementia", age > 75, age <= 95) %>%
  ggplot(aes(age, estimate)) +
  geom_line(aes(color = ethnicity)) +
  facet_grid(rows = vars(female), scales = "free") +
  scale_color_manual(values = color_palette) + 
  # geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = ethnicity), alpha = 0.2) +
  labs(
    title = "Cumulative incidence curves",
    x = "Age (years)",
    y = "Cumulaitve incidence",
    colour = "Race/Ethnicity",
    fill = "Race/Ethnicity") +
  # scale_x_continuous(limits = c(0, 18), 
  #                    breaks = seq(0, 18, 2)) +
  # scale_y_continuous(limits = c(0, 1),  breaks=seq(0, 1, 0.2)) +
  theme_bw()

ggsave(here::here("02_figs", "fig_cuminc_curve_AJ_dem_ethns.png"),
       device = "png", width = 9, height = 7, units = "in", dpi = 300)

##---- Cumulative incidence plot death only ----
crude_AJ_tib %>% 
  # # truncate at age 90 for South Asians
  # filter(!(age >= 90 & ethnicity == "South Asian")) %>%
  mutate(across(c("estimate", "conf.low", "conf.high"), function(x) x*100)) %>%
  # # cut the y axis at 50, top coded conf.high at 50
  # mutate(conf.high_tc = case_when(conf.high >= 50 ~ 50,
  #                                 TRUE ~ conf.high)) %>%
  filter(state == "Death without dementia", age <= 95) %>%
  ggplot(aes(x = age, y = estimate, group = female)) +
  geom_line(aes(color = female), linetype = "dashed") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = female), alpha = 0.2) +
  facet_grid(cols = vars(ethnicity), scales = "free", space = "free") +
  scale_color_manual(name = NULL,
                     values = c("#3271AE", "#E5A84B"),
                     labels = c("Women", "Men")) +
  scale_fill_manual(name = NULL,
                    values = c("#3271AE", "#E5A84B"),
                    labels = c("Women", "Men")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Age (years)",
       y = "Cumulative incidence (%)") 

ggsave(here::here("02_figs", "fig_cuminc_AJ_death_ethns_withCI.png"),
       device = "png", width = 7, height = 5, units = "in", dpi = 300)

##---- Cumulative incidence plot dem and death ----
crude_AJ_tib %>% 
  filter(state %in% c("Dementia", "Death without dementia"), age <= 95) %>%
  ggplot(aes(x = age, y = estimate, group = interaction(state, female))) + 
  geom_line(aes(color = female, linetype = state)) +
  facet_grid(cols = vars(ethnicity)) +
  scale_color_manual(name = NULL,
                     values = c("#3271AE", "#E5A84B"),
                     labels = c("Women", "Men")) +
  scale_linetype_manual(name = NULL,
                        values = c(2, 1),
                        labels = c("Death without dementia", "Dementia")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Age (years)",
       y = "Cumulative incidence (%)") 

ggsave(here::here("02_figs", "fig_cuminc_AJ_dem_death_ethns.png"),
       device = "png", width = 7, height = 5, units = "in", dpi = 300)

# #---- OLD ----
# ###----- RR RD ----
# crude_AJ_1y_tib <- crude_AJ_tib %>%
#   filter(state == "Dementia") %>%
#   select(ethnicity, age, female, estimate, conf.high, conf.low) %>%
#   filter(!is.na(conf.high) & !is.na(conf.low)) %>%
#   mutate(age_cut = cut(age, breaks = seq(60, 103, 1),
#                        labels = seq(60, 102, 1),
#                        right = FALSE)) %>% 
#   group_by(ethnicity, female, age_cut) %>%
#   slice_head() %>%
#   ungroup() %>%
#   arrange(age_cut, female) %>% 
#   select(ethnicity, age_cut, cuminc = estimate, female)
# 
# # Sanity check
# # # Manually checking the ages with NA estimate for cuminc
# # with(crude_AJ_1y_tib, table(female, age_cut, ethnicity, useNA = "ifany"))
# # crude_AJ_tib %>% filter(state == "Dementia", ethnicity == "Chinese", age < 63) %>% view()
# # crude_AJ_tib %>% filter(state == "Dementia", ethnicity == "Chinese", age >100, age < 102) %>% view()
# # crude_AJ_tib %>% filter(state == "Dementia", ethnicity == "Filipino", age < 62) %>% view()
# # crude_AJ_tib %>% filter(state == "Dementia", ethnicity == "Filipino", age > 99, age < 103) %>% view()
# 
# crude_AJ_RR_RD_tib <- crude_AJ_1y_tib %>%
#   mutate(cuminc = cuminc*100) %>%
#   pivot_wider(id_cols = c("ethnicity", "age_cut"),
#               names_from = female,
#               values_from = c("cuminc"),
#               names_prefix = "cuminc_") %>%
#   na.omit() %>%
#   mutate(RR = cuminc_Women/cuminc_Men,
#          RD = cuminc_Women - cuminc_Men) %>%
#   filter(age_cut %in% seq(75, 95, 5)) %>% 
#   select(ethnicity, age_cut, cuminc_Women, cuminc_Men, RR, RD) %>%
#   arrange(ethnicity, age_cut)
# 
# writexl::write_xlsx(crude_AJ_RR_RD_tib , here::here("03_tables", "raw_tables", 
#                                                     "table_crude_RR_RD_AJ_dem_ethns.xlsx"))
# 
# #---- Weight by education ----
# female_num_mod <- glm(female ~ ethnicity_rev,
#                       family = binomial(link = "logit"),
#                       data = aa_adrd_tte_selected)
# 
# weight_edu_formula <- "female ~ edu_4_miss*ethnicity_rev"
# female_denom_mod <- glm(as.formula(weight_edu_formula),
#                         family = binomial(link = "logit"),
#                         data = aa_adrd_tte_selected)
# # summary(female_denom_mod)
# 
# # # test if the weight is the same with stratified by ethnicity 
# # model_Chinese <- glm("female ~ edu_4_miss",
# #                      family = binomial(),
# #                      data = aa_adrd_tte_selected %>% filter(ethnicity_rev == "Chinese"))
# # summary(model_Chinese)
# # 
# # model_Filipino <- glm("female ~ edu_4_miss",
# #                      family = binomial(),
# #                      data = aa_adrd_tte_selected %>% filter(ethnicity_rev == "Filipino"))
# # summary(model_Filipino)
# # # It's the same with the full model with ethnicity interaction terms
# 
# aa_adrd_tte_selected %<>%
#   mutate(
#     pred_num = predict(female_num_mod , newdata = aa_adrd_tte_selected, 
#                        type = "response"),
#     pred_denom = predict(female_denom_mod, newdata = aa_adrd_tte_selected, 
#                          type = "response"),
#     pmnum = case_when(female == "Women" ~ 1 - pred_num,
#                       female == "Men" ~ pred_num),
#     pmdenom = case_when(female == "Women" ~ 1 - pred_denom,
#                         female == "Men" ~ pred_denom),
#     ipmw_unstblz = 1/pmdenom,
#     ipmw_stblz = pmnum/pmdenom)
# 
# aa_adrd_tte_selected %>% 
#   select(ethnicity_rev, contains("ipmw")) %>% 
#   split(.$ethnicity_rev) %>% map(summary)
# 
# # weighted models
# # crude_cox_mod <- coxph(Surv(survey_age, main_dem_v1_end_age, event_dem) ~ 
# #                          ethnicity_rev*female, data = aa_adrd_tte_selected)
# # wtd_cox_mod <- coxph(Surv(survey_age, main_dem_v1_end_age, event_dem) ~ 
# #                        ethnicity_rev*female, data = aa_adrd_tte_selected,
# #                      weights = ipmw_stblz)
# ##---- Wtd age as time scale ----
# wtd_AJ_fit <- 
#   survfit(Surv(survey_age, main_dem_v1_end_age, event_multi, type = "mstate") ~ 
#             strata(ethnicity_rev, female), data = aa_adrd_tte_selected,
#           cluster = subjid, id = subjid, weights = ipmw_stblz)
# 
# wtd_AJ_tib <- broom::tidy(wtd_AJ_fit) %>%
#   dplyr::rename(age = time) %>%
#   mutate(strata = stringr::str_remove(strata, "strata\\(ethnicity_rev, female\\)="),
#          ethnicity = factor(str_split_i(strata, ",", 1), 
#                             levels = levels(aa_adrd_tte_selected$ethnicity_rev)),
#          female = factor(str_split_i(strata, ", ", 2),
#                          levels = levels(aa_adrd_tte_selected$female)),
#          state = case_when(state == "(s0)" ~ "Alive without dementia",
#                            state == "death" ~ "Death without dementia",
#                            state == "dementia" ~ "Dementia")) 
# 
# ###---- Cumulative incidence plot ----
# wtd_AJ_tib %>% 
#   filter(state != "Alive without dementia") %>%
#   ggplot(aes(x = age, y = estimate)) +
#   # group = interaction(state, female))) + 
#   geom_line(aes(color = female, 
#                 linetype = state)) +
#   facet_grid(cols = vars(ethnicity)) +
#   scale_color_manual(name = NULL,
#                      values = c("salmon", "cyan3"),
#                      labels = c("Women", "Men")) +
#   scale_linetype_manual(name = NULL,
#                         values = c(2, 1),
#                         labels = c("Death without dementia", "Dementia")) +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   labs(x = "Age (years)",
#        y = "Cumulative incidence (%)") 
# 
# # wtd_AJ_tib %>% 
# #   ggplot(aes(age, estimate, fill = state)) + geom_area() + 
# #   facet_grid(rows = vars(female), cols = vars(ethnicity)) +
# #   cowplot::theme_cowplot() +
# #   ggsci::scale_fill_jco() +
# #   theme(legend.position = "bottom") +
# #   labs(x = "Age (years)",
# #        y = "Cumulative incidence (%)", 
# #        fill = NULL) 
# 
# ggsave(here::here("02_figs", "fig_cuminc_AJ_dem_ethns_wtd_edu.png"),
#        device = "png", width = 9, height = 7, units = "in", dpi = 300)
# 
# ###----- RR RD ----
# wtd_AJ_1y_tib <- wtd_AJ_tib %>%
#   filter(state == "Dementia") %>%
#   select(ethnicity, age, female, estimate, conf.high, conf.low) %>%
#   filter(!is.na(conf.high) & !is.na(conf.low)) %>%
#   mutate(age_cut = cut(age, breaks = seq(60, 103, 1),
#                        labels = seq(60, 102, 1),
#                        right = FALSE)) %>% 
#   group_by(ethnicity, female, age_cut) %>%
#   slice_head() %>%
#   ungroup() %>%
#   arrange(age_cut, female) %>% 
#   select(ethnicity, age_cut, cuminc = estimate, female)
# 
# # Sanity check
# # # Manually checking the ages with NA estimate for cuminc
# # with(wtd_AJ_1y_tib, table(female, age_cut, ethnicity, useNA = "ifany"))
# # wtd_AJ_tib %>% filter(state == "Dementia", ethnicity == "Chinese", age < 63) %>% view()
# # wtd_AJ_tib %>% filter(state == "Dementia", ethnicity == "Chinese", age >100, age < 102) %>% view()
# # wtd_AJ_tib %>% filter(state == "Dementia", ethnicity == "Filipino", age < 62) %>% view()
# # wtd_AJ_tib %>% filter(state == "Dementia", ethnicity == "Filipino", age > 99, age < 103) %>% view()
# 
# wtd_AJ_RR_RD_tib <- wtd_AJ_1y_tib %>%
#   mutate(cuminc = cuminc*100) %>%
#   pivot_wider(id_cols = c("ethnicity", "age_cut"),
#               names_from = female,
#               values_from = c("cuminc"),
#               names_prefix = "cuminc_") %>%
#   na.omit() %>%
#   mutate(RR = cuminc_Women/cuminc_Men,
#          RD = cuminc_Women - cuminc_Men) %>%
#   filter(age_cut %in% seq(75, 95, 5)) %>% 
#   select(ethnicity, age_cut, cuminc_Women, cuminc_Men, RR, RD) %>%
#   arrange(ethnicity, age_cut)
# 
# writexl::write_xlsx(wtd_AJ_RR_RD_tib , here::here("03_tables", "raw_tables", 
#                                                   "table_wtd_RR_RD_AJ_dem_ethns.xlsx"))
# 
# 
# ##---- All RR RD ----
# AJ_RR_RD_tib <- crude_AJ_RR_RD_tib %>% rename_at(vars(contains(c("cuminc", "RR", "RD"))),
#                                                  function(x) paste0("Unweighted ", x)) %>%
#   left_join(wtd_AJ_RR_RD_tib %>% rename_at(vars(contains(c("cuminc", "RR", "RD"))),
#                                            function(x) paste0("Weighted ", x)),
#             by = c("ethnicity", "age_cut"))
# 
# writexl::write_xlsx(AJ_RR_RD_tib , here::here("03_tables", "raw_tables", 
#                                               "table_RR_RD_AJ_dem_ethns.xlsx"))

# #---- OLD ----
# #---- Time as time scale ----
# ##---- Crude model ----
# crude_AJ_fit_time_time <- 
#   survfit(Surv(time = starttime, time2 = main_dem_v1_fu_time, event_multi, type = "mstate") ~ 
#             strata(ethnicity_rev, female), data = aa_adrd_tte_selected,
#           cluster = subjid, id = subjid)
# 
# crude_AJ_tib_time <- broom::tidy(crude_AJ_fit_time) %>%
#   mutate(strata = stringr::str_remove(strata, "strata\\(ethnicity_rev, female\\)="),
#          ethnicity = factor(str_split_i(strata, ",", 1), 
#                             levels = levels(aa_adrd_tte_selected$ethnicity_rev)),
#          female = factor(str_split_i(strata, ", ", 2),
#                          levels = levels(aa_adrd_tte_selected$female)),
#          state = case_when(state == "(s0)" ~ "Alive without dementia",
#                            state == "death" ~ "Death without dementia",
#                            state == "dementia" ~ "Dementia")) 
# 
# ###---- Cumulative incidence plot ----
# crude_AJ_tib_time %>% 
#   ggplot(aes(time, estimate, fill = state)) + geom_area() + 
#   facet_grid(rows = vars(female), cols = vars(ethnicity)) +
#   cowplot::theme_cowplot() +
#   ggsci::scale_fill_jco() +
#   theme(legend.position = "bottom") +
#   labs(x = "Follow up time (years)",
#        y = "Cumulative incidence (%)", 
#        fill = NULL) 
# 
# ggsave(here::here("02_figs", "fig_cuminc_AJ_dem_ethns_time.png"),
#        device = "png", width = 9, height = 7, units = "in", dpi = 300)
# 
# crude_AJ_tib_time %>%
#   filter(state == "Dementia") %>%
#   ggplot(aes(time, estimate)) +
#   geom_line(aes(color = ethnicity)) +
#   facet_grid(rows = vars(female), scales = "free") +
#   # geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = ethnicity), alpha = 0.2) +
#   labs(
#     title = "Cumulative incidence curves",
#     x = "Follow up time (years)",
#     y = "Cumulaitve incidence",
#     colour = "Race/Ethnicity",
#     fill = "Race/Ethnicity") +
#   # scale_x_continuous(limits = c(0, 18), 
#   #                    breaks = seq(0, 18, 2)) +
#   # scale_y_continuous(limits = c(0, 1),  breaks=seq(0, 1, 0.2)) +
#   theme_bw()
# ggsave(here::here("02_figs", "fig_cuminc_curve_dem_ethns*gender_time.png"),
#        device = "png", width = 9, height = 7, units = "in", dpi = 300)
# 
# ###----- RR RD ----
# crude_AJ_1y_tib_time <- crude_AJ_tib_time %>%
#   filter(state == "Dementia") %>%
#   select(ethnicity, time, female, estimate, conf.high, conf.low) %>%
#   filter(!is.na(conf.high) & !is.na(conf.low)) %>%
#   mutate(time_cut = cut(time, breaks = seq(0, 19, 1),
#                         labels = seq(0, 18, 1),
#                         right = FALSE)) %>% 
#   group_by(ethnicity, female, time_cut) %>%
#   slice_head() %>%
#   ungroup() %>%
#   arrange(time_cut, female) %>% 
#   select(ethnicity, time_cut, cuminc = estimate, female)
# 
# crude_AJ_RR_RD_tib_time <- crude_AJ_1y_tib_time %>%
#   mutate(cuminc = cuminc*100) %>%
#   pivot_wider(id_cols = c("ethnicity", "time_cut"),
#               names_from = female,
#               values_from = c("cuminc"),
#               names_prefix = "cuminc_") %>%
#   na.omit() %>%
#   mutate(RR = cuminc_Women/cuminc_Men,
#          RD = cuminc_Women - cuminc_Men) %>%
#   filter(time_cut %in% seq(0, 18, 5)) %>%
#   select(ethnicity, time_cut, cuminc_Women, cuminc_Men, RR, RD) %>%
#   arrange(ethnicity, time_cut)
# 
# writexl::write_xlsx(crude_AJ_RR_RD_tib_time, 
#                     here::here("03_tables", "raw_tables", 
#                                "table_crude_RR_RD_AJ_dem_ethns_time.xlsx"))
# # cuminc tib
# crude_AJ_cuminc_tib <- crude_AJ_tib %>%
#   filter(state == "Dementia") %>%
#   mutate(age_grp = cut(time,
#                        breaks = c(60, 65, 70, 75, 80, 85, 90, 95, 120),
#                        labels = c(60, 65, 70, 75, 80, 85, 90, 95),
#                        right = FALSE)) %>%
#   group_by(age_grp, strata) %>%
#   slice_min(n = 1, order_by = time) %>%
#   ungroup() %>%
#   filter(age_grp %in% seq(70, 95, by = 5)) %>%
#   dplyr::select(ethnicity, age = time, female, estimate) %>%
#   mutate(age = floor(age),
#          estimate = estimate*100) %>%
#   pivot_wider(names_from = "female", values_from = "estimate", 
#               names_prefix = "cuminc_") %>%
#   arrange(ethnicity) %>%
#   relocate(ethnicity, age, cuminc_Women, cuminc_Men)
# 
# write_xlsx(crude_AJ_cuminc_tib, here::here("03_tables", "raw_tables", 
#                                     "table_cuminc_AJ_dem_ethns.xlsx"))

# estimate_AJ <- function(data) {
#   # Test
#   # data <- aa_adrd_tte_selected %>% filter(ethnicity_rev == "South Asian")
#   require(survival)
#   require(broom)
#   
#   # because of bootstrap resampling, the original ID's are no longer unique
#   # need to create new ID's for survfit
#   data <- data %>% mutate(new_id = 1:n())
#   
#   crude_AJ_fit <- 
#     survfit(Surv(survey_age, main_dem_v1_end_age, event_multi) ~ female, 
#             id = new_id, data = data)
#   # For ggcompetingrisk plot's facet labels
#   names(crude_AJ_fit$strata) <- c("Women", "Men")
#   
#   crude_AJ_tib <- broom::tidy(crude_AJ_fit) %>% 
#     filter(state == "dementia") %>% 
#     mutate(age_grp = cut(time, 
#                          breaks = c(60, 65, 70, 75, 80, 85, 90, 120), 
#                          labels = c(60, 65, 70, 75, 80, 85, 90),
#                          right = FALSE)) %>% 
#     group_by(age_grp, strata) %>% 
#     slice_min(n = 1, order_by = time) %>% 
#     ungroup() %>% 
#     # Order for sex/gender: women then men
#     # arrange(time, desc(strata)) %>%
#     filter(age_grp != 60) %>%
#     select(time, strata, estimate)
#   
#   # plot cumulative incidence
#   crude_AJ_plot <- ggcompetingrisks(crude_AJ_fit) +
#     facet_wrap(~rev(strata))
#   cowplot::theme_cowplot() +
#     ggsci::scale_fill_jco() +
#     labs(x = "Age")
#   
#   
#   
#   out <- list(crude_AJ_fit = crude_AJ_fit, crude_AJ_tib = crude_AJ_tib,
#               crude_AJ_plot = crude_AJ_plot)
#   
#   return(out)
# }
# 
# ethns <- levels(aa_adrd_tte_selected$ethnicity_rev)
# 
# crude_AJ_list <- lapply(ethns, function(x) 
#   estimate_AJ(aa_adrd_tte_selected %>% filter(ethnicity_rev == x)))
# names(crude_AJ_list) <- ethns
# 
# for (ethn in ethns){
#   p <- crude_AJ_list[[ethn]]$crude_AJ_plot +
#     theme(plot.background = element_rect(fill = "white", color = NA)) +
#     labs(title = paste0("Cumulative incidence of dementia:", ethn))
#   ggsave(here::here("02_figs", paste0("fig_cuminc_AJ_dem_", ethn, ".png")),
#          device = "png", width = 9, height = 7, units = "in", dpi = 300)
# }

# #----Time as time scale Weight by education ----
# female_num_mod <- glm(female ~ ethnicity_rev,
#                       family = binomial(link = "logit"),
#                       data = aa_adrd_tte_selected)
# 
# weight_edu_formula <- "female ~ edu_4_miss*ethnicity_rev  + survey_age*ethnicity_rev"
# female_denom_mod <- glm(as.formula(weight_edu_formula),
#                         family = binomial(link = "logit"),
#                         data = aa_adrd_tte_selected)
# # summary(female_denom_mod)
# 
# # # test if the weight is the same with stratified by ethnicity 
# # model_Chinese <- glm("female ~ edu_4_miss",
# #                      family = binomial(),
# #                      data = aa_adrd_tte_selected %>% filter(ethnicity_rev == "Chinese"))
# # summary(model_Chinese)
# # 
# # model_Filipino <- glm("female ~ edu_4_miss",
# #                      family = binomial(),
# #                      data = aa_adrd_tte_selected %>% filter(ethnicity_rev == "Filipino"))
# # summary(model_Filipino)
# # # It's the same with the full model with ethnicity interaction terms
# 
# aa_adrd_tte_selected %<>%
#   mutate(
#     pred_num = predict(female_num_mod , newdata = aa_adrd_tte_selected, 
#                        type = "response"),
#     pred_denom = predict(female_denom_mod, newdata = aa_adrd_tte_selected, 
#                          type = "response"),
#     pmnum = case_when(female == "Women" ~ 1 - pred_num,
#                       female == "Men" ~ pred_num),
#     pmdenom = case_when(female == "Women" ~ 1 - pred_denom,
#                         female == "Men" ~ pred_denom),
#     ipmw_unstblz = 1/pmdenom,
#     ipmw_stblz = pmnum/pmdenom)
# 
# aa_adrd_tte_selected %>% 
#   select(ethnicity_rev, contains("ipmw")) %>% 
#   split(.$ethnicity_rev) %>% map(summary)
# 
# # weighted models
# # crude_cox_mod <- coxph(Surv(survey_age, main_dem_v1_end_age, event_dem) ~ 
# #                          ethnicity_rev*female, data = aa_adrd_tte_selected)
# # wtd_cox_mod <- coxph(Surv(survey_age, main_dem_v1_end_age, event_dem) ~ 
# #                        ethnicity_rev*female, data = aa_adrd_tte_selected,
# #                      weights = ipmw_stblz)
# 
# ##---- Wtd model ----
# wtd_AJ_fit_time <- 
#   survfit(Surv(time = starttime, time2 = main_dem_v1_fu_time, event_multi, type = "mstate") ~ 
#             strata(ethnicity_rev, female), data = aa_adrd_tte_selected,
#           cluster = subjid, id = subjid, weights = ipmw_stblz)
# 
# wtd_AJ_tib_time <- broom::tidy(wtd_AJ_fit_time) %>%
#   mutate(strata = stringr::str_remove(strata, "strata\\(ethnicity_rev, female\\)="),
#          ethnicity = factor(str_split_i(strata, ",", 1), 
#                             levels = levels(aa_adrd_tte_selected$ethnicity_rev)),
#          female = factor(str_split_i(strata, ", ", 2),
#                          levels = levels(aa_adrd_tte_selected$female)),
#          state = case_when(state == "(s0)" ~ "Alive without dementia",
#                            state == "death" ~ "Death without dementia",
#                            state == "dementia" ~ "Dementia")) 
# 
# ###---- Cumulative incidence plot ----
# wtd_AJ_tib_time %>% 
#   ggplot(aes(time, estimate, fill = state)) + geom_area() + 
#   facet_grid(rows = vars(female), cols = vars(ethnicity)) +
#   cowplot::theme_cowplot() +
#   ggsci::scale_fill_jco() +
#   theme(legend.position = "bottom") +
#   labs(x = "Follow up time (years)",
#        y = "Cumulative incidence (%)", 
#        fill = NULL) 
# 
# # ggsave(here::here("02_figs", "fig_cuminc_AJ_dem_ethns_wtd_edu_time.png"),
# #        device = "png", width = 9, height = 7, units = "in", dpi = 300)
# 
# ###----- RR RD ----
# wtd_AJ_1y_tib_time <- wtd_AJ_tib_time %>%
#   filter(state == "Dementia") %>%
#   select(ethnicity, time, female, estimate, conf.high, conf.low) %>%
#   filter(!is.na(conf.high) & !is.na(conf.low)) %>%
#   mutate(time_cut = cut(time, breaks = seq(0, 19, 1),
#                         labels = seq(0, 18, 1),
#                         right = FALSE)) %>% 
#   group_by(ethnicity, female, time_cut) %>%
#   slice_head() %>%
#   ungroup() %>%
#   arrange(time_cut, female) %>% 
#   select(ethnicity, time_cut, cuminc = estimate, female)
# 
# wtd_AJ_RR_RD_tib_time <- wtd_AJ_1y_tib_time %>%
#   mutate(cuminc = cuminc*100) %>%
#   pivot_wider(id_cols = c("ethnicity", "time_cut"),
#               names_from = female,
#               values_from = c("cuminc"),
#               names_prefix = "cuminc_") %>%
#   na.omit() %>%
#   mutate(RR = cuminc_Women/cuminc_Men,
#          RD = cuminc_Women - cuminc_Men) %>%
#   filter(time_cut %in% seq(0, 18, 5)) %>%
#   select(ethnicity, time_cut, cuminc_Women, cuminc_Men, RR, RD) %>%
#   arrange(ethnicity, time_cut)
# 
# ##---- All RR RD ----
# AJ_RR_RD_tib_time <- crude_AJ_RR_RD_tib_time %>% 
#   rename_at(vars(contains(c("cuminc", "RR", "RD"))),
#             function(x) paste0("Unweighted ", x)) %>%
#   left_join(wtd_AJ_RR_RD_tib_time %>% 
#               rename_at(vars(contains(c("cuminc", "RR", "RD"))),
#                         function(x) paste0("Weighted ", x)),
#             by = c("ethnicity", "time_cut"))
# 
# writexl::write_xlsx(AJ_RR_RD_tib_time, here::here("03_tables", "raw_tables", 
#                                                   "table_RR_RD_AJ_dem_ethns_time.xlsx"))