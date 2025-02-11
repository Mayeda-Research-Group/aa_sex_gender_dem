# aa_sex_gender_dem
This repo contains all the code necessary (folder 01_R) to replicate data construction and analyses in the paper: Sex/gender differences in lifetime dementia risk among  Asian American ethnic groups and non-Latino White older adults in California

## Main scripts include: 

* Scripts 0: specifying paths to data folders. consolidate datasets from various sources.
* Script 1 constructs the time-to-event dataset in a wide format for analysis.
* Script 2.0: Calculates age-specific and age-adjusted IR
* Script 2.1: Generate descriptive stats and figures for baseline age, count of deaths, and cumulative incidence (Aalen–Johansen estimators).
* Script 3.1: Script to run Aalen–Johansen estimation with bootstrap on the Hoffman cluster. The script contains a function from bootstrapping data, fitting the weight model for education, estimating AJ estimators, and tidying the results. A tibble containing job scenarios was also set in the script.
* Script 3.2: Submission script for submitting script 3.1 to Hoffman.
* Script 3.3: Reads in the bootstrap results. Obtain confidence intervals for RD and RR and generate plots and tables.
* Script 9: Generates descriptive tables: Table 1stratified by sex/gender and race/ethnicity.
