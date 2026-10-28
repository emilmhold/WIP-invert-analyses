# ---
# title: clean 2023 WIP data
# author: Emily H
# created: "2026-10-21"
# last edited: 2026-10-27
# inputs: [input_data.csv]
# outputs: [results.csv]
# notes: 
#   This script prepares the 2023 invert data for the taxonomic resolution
#   algorithm.
# ---

# 1. Setup ----
rm(list=ls())
gc()

## 1.1 Load packages ----
# Include comments for what the packages are used for
# Specify the package version for reproducibility.
# Example:
library(readxl)
library(abmi.themes)
library(tidyverse)

## 1.3 Import data ----
# lookup
lookup <- read_csv("0_data/external/Look-up Table_Invert Taxonomy_2023-09-28.csv")
str(lookup)

# coarse data
coarse_data <- read_excel("0_data/external/WIP_2023_Coarse Level Community Analysis.xlsx") |>
  filter(!is.na(`Group/Sample Number`)) |>
  filter(row_number() < which(`Group/Sample Number` == "Primary Organisms")) |>
  select(!any_of(c("Mean", "Min", "Max", "Total"))) |>
  mutate(across(-`Group/Sample Number`, \(x) na_if(trimws(x), "no sample"))) |>
  mutate(`Group/Sample Number` = str_replace_all(`Group/Sample Number`, " ", ".")) |>
  column_to_rownames(var = "Group/Sample Number") |>
  t() |>
  as.data.frame()
str(coarse_data)

# fine data 
fine_2023 <- read_excel("0_data/external/WIP_2023_Data Matrix_DPAC-SG.xlsx", sheet = 1) |>
  select(!"Row Number") |>
  column_to_rownames(var = "Site") |>
  mutate(across(where(is.numeric) | where(is.logical), \(x) replace(x, is.na(x), 0))) |>
  rename_with(\(x) gsub(" ", ".", x, fixed = TRUE)) |>
  rename_with(\(x) gsub("..", ".", x, fixed = TRUE)) 
str(fine_2023)

# UMOS data
umos_2023 <- read_excel("0_data/external/WIP_2023_Data Matrix_DPAC-SG.xlsx", sheet = 2) |>
  select(!c('Row Number','Total UMOS')) |>
  column_to_rownames(var = "Site") |>
  mutate(across(where(is.numeric) | where(is.logical), \(x) replace(x, is.na(x), 0))) |>
  rename_with(\(x) gsub(" ", ".", x, fixed = TRUE)) |>
  rename_with(\(x) gsub("..", ".", x, fixed = TRUE)) 
str(umos_2023)

marchant_counts <- read_excel("0_data/external/WIP_2023_Coarse Level Community Analysis.xlsx") |>
  filter(`Group/Sample Number` == "Marchant Cells Counted") |>
  select(-`Group/Sample Number`) |>
  pivot_longer(
    everything(),
    names_to = "Sample",
    values_to = "Marchant.Cells.Counted",
    values_transform = as.numeric
  )

# 2. Check for species not in lookup table ----
nm <- lookup$Analysis_Name
need_lookup <- list(
  fine_not_in = setdiff(names(fine_2023), nm),
  umos_not_in = setdiff(names(umos_2023), nm)
)
need_lookup_df <- tibble::enframe(need_lookup, name = "source", value = "column_name") |>
  tidyr::unnest(column_name)

#export list for export review
write_csv(need_lookup_df, "2_pipeline/2023 taxa for expert review.csv")

# 3. Export cleaned data ----
# Read by "resolve 2023 WIP data.R"
write_rds(coarse_data, "0_data/processed/WIP_2023_coarse_clean.rds")
write_rds(marchant_counts, "0_data/processed/WIP_2023_marchant_counts.rds")
write_rds(fine_2023, "0_data/processed/WIP_2023_fine_clean.rds")
write_rds(umos_2023, "0_data/processed/WIP_2023_umos_clean.rds")

# End of script ----
