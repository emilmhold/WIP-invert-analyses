# ---
# title: clean 2025 WIP data
# author: Emily H
# created: "2026-10-29"
# last edited: 2026-10-29
# inputs: [input_data.csv]
# outputs: [results.csv]
# notes: 
#   This script prepares the 2025 invert data for the taxonomic resolution
#   algorithm. This is separate from the other data prep script because the
#   data came in a different format.
# ---

# 1. Setup ----
rm(list=ls())
gc()

## 1.1 Load packages ----
# Include comments for what the packages are used for
# Specify the package version for reproducibility.
# Example:
library(readxl)
library(tidyverse)

## 1.3 Import data ----
# lookup
lookup <- read_excel("0_data/external/Look-up Table_Invert Taxonomy_2024-08-27.xlsx")
str(lookup)

# fine data 
## spring data
spring_fine_2025 <- read_excel("0_data/external/Advanced ID_WGP_2025_AB_spring.xlsx") |>
  rename_with(\(x) gsub(" ", ".", x, fixed = TRUE)) |>
  mutate(Sample.Number = if_else(
    is.na(Sample.Number), NA_character_, str_c(Sample.Number, "_spring")
  )) |>
  select(Sample.Number, Coarse.Sorting.Group, `Operational.Taxonomic.Unit.(OTU)`, Advanced.ID.Count) |>
  rename(OTU = `Operational.Taxonomic.Unit.(OTU)`,
         Site = Sample.Number) |>
  filter(!Coarse.Sorting.Group %in% c("UMOS", "Homogenized Sample"))

spring_doubles <- spring_fine_2025 |>
  add_count(Site, OTU) |>
  filter(n > 1)
write.csv(spring_doubles, "2_pipeline/2025 double observations for expert review.csv")

spring_fine_2025 <- spring_fine_2025 |>
  pivot_wider(id_cols = Site, names_from = OTU, values_from = Advanced.ID.Count,
              values_fill = 0, values_fn = sum) |>
  rename_with(\(x) gsub(" ", ".", x, fixed = TRUE)) |>
  rename_with(\(x) gsub("..", ".", x, fixed = TRUE))
str(spring_fine_2025)

## summer data
summer_fine_2025 <- read_excel("0_data/external/Advanced ID_WGP_2025_AB_summer_CORRECTED-7Oct2026.xlsx") |>
  rename_with(\(x) gsub(" ", ".", x, fixed = TRUE)) |>
  mutate(Sample.Number = if_else(
    is.na(Sample.Number), NA_character_, str_c(Sample.Number, "_summer")
  )) |>
  select(Sample.Number, Coarse.Sorting.Group, `Operational.Taxonomic.Unit.(OTU)`, Advanced.ID.Count) |>
  rename(OTU = `Operational.Taxonomic.Unit.(OTU)`,
         Site = Sample.Number) |>
  filter(!Coarse.Sorting.Group %in% c("UMOS", "Homogenized Sample"))|>
  pivot_wider(id_cols = Site, names_from = OTU, values_from = Advanced.ID.Count,
              values_fill = 0, values_fn = sum) |>
  rename_with(\(x) gsub(" ", ".", x, fixed = TRUE)) |>
  rename_with(\(x) gsub("..", ".", x, fixed = TRUE))
str(summer_fine_2025)

## join datasets
fine_2025 <- bind_rows(spring_fine_2025, summer_fine_2025) |>
  mutate(across(-Site, \(x) replace_na(x, 0))) |>
  column_to_rownames("Site") |>
  select(!UID)

# UMOS data
## spring data
spring_umos_2025 <- read_excel("0_data/external/Advanced ID_WGP_2025_AB_spring.xlsx") |>
  rename_with(\(x) gsub(" ", ".", x, fixed = TRUE)) |>
  mutate(Sample.Number = if_else(
    is.na(Sample.Number), NA_character_, str_c(Sample.Number, "_spring")
  )) |>
  select(Sample.Number, Coarse.Sorting.Group, `Operational.Taxonomic.Unit.(OTU)`, Advanced.ID.Count) |>
  rename(OTU = `Operational.Taxonomic.Unit.(OTU)`,
         Site = Sample.Number) |>
  filter(Coarse.Sorting.Group == "UMOS") |>
  pivot_wider(id_cols = Site, names_from = OTU, values_from = Advanced.ID.Count,
              values_fill = 0, values_fn = sum) |>
  rename_with(\(x) gsub(" ", ".", x, fixed = TRUE)) |>
  rename_with(\(x) gsub("..", ".", x, fixed = TRUE))
str(spring_umos_2025)

## summer data
summer_umos_2025 <- read_excel("0_data/external/Advanced ID_WGP_2025_AB_summer_CORRECTED-7Oct2026.xlsx") |>
  rename_with(\(x) gsub(" ", ".", x, fixed = TRUE)) |>
  mutate(Sample.Number = if_else(
    is.na(Sample.Number), NA_character_, str_c(Sample.Number, "_summer")
  )) |>
  select(Sample.Number, Coarse.Sorting.Group, `Operational.Taxonomic.Unit.(OTU)`, Advanced.ID.Count) |>
  rename(OTU = `Operational.Taxonomic.Unit.(OTU)`,
         Site = Sample.Number) |>
  filter(Coarse.Sorting.Group == "UMOS") |>
  pivot_wider(id_cols = Site, names_from = OTU, values_from = Advanced.ID.Count,
              values_fill = 0, values_fn = sum) |>
  rename_with(\(x) gsub(" ", ".", x, fixed = TRUE)) |>
  rename_with(\(x) gsub("..", ".", x, fixed = TRUE))
str(summer_umos_2025)

## join datasets
umos_2025 <- bind_rows(spring_umos_2025, summer_umos_2025) |>
  mutate(across(-Site, \(x) replace_na(x, 0))) |>
  column_to_rownames("Site") 

# Marchant counts
## spring
spring_marchant_counts_2025 <- read_excel("0_data/external/Advanced ID_WGP_2025_AB_spring.xlsx") |>
  rename_with(\(x) gsub(" ", ".", x, fixed = TRUE)) |>
  mutate(Sample.Number = if_else(
    is.na(Sample.Number), NA_character_, str_c(Sample.Number, "_spring")
  )) |>
  filter(!Coarse.Sorting.Group %in% c("UMOS", "Homogenized Sample")) |>
  select(Sample.Number, Number.of.Cells.Counted.in.Marchant.Box) |>
  rename(Marchant.Cells.Counted = Number.of.Cells.Counted.in.Marchant.Box) |>
  distinct() |>
  rename(Site = Sample.Number)
str(spring_marchant_counts_2025)

## summer
summer_marchant_counts_2025 <- read_excel("0_data/external/Advanced ID_WGP_2025_AB_summer_CORRECTED-7Oct2026.xlsx") |>
  rename_with(\(x) gsub(" ", ".", x, fixed = TRUE)) |>
  mutate(Sample.Number = if_else(
    is.na(Sample.Number), NA_character_, str_c(Sample.Number, "_summer")
  )) |>
  filter(!Coarse.Sorting.Group %in% c("UMOS", "Homogenized Sample")) |>
  select(Sample.Number, Number.of.Cells.Counted.in.Marchant.Box) |>
  rename(Marchant.Cells.Counted = Number.of.Cells.Counted.in.Marchant.Box) |>
  distinct() |>
  rename(Site = Sample.Number)
str(summer_marchant_counts_2025)

## join datasets
marchant_counts_2025 <- rbind(spring_marchant_counts_2025, summer_marchant_counts_2025)

# 2. Check for species not in lookup table ----
nm <- lookup$Analysis_Name
need_lookup <- list(
  fine_not_in = setdiff(names(fine_2025), nm),
  umos_not_in = setdiff(names(umos_2025), nm)
)
need_lookup_df <- tibble::enframe(need_lookup, name = "source", value = "column_name") |>
  tidyr::unnest(column_name)

#export list for export review
write_csv(need_lookup_df, "2_pipeline/2025 taxa for expert review.csv")

# 3. Export cleaned data ----
# Read by "resolve WIP data.R"
write_rds(coarse_data, "0_data/processed/WIP_2025_coarse_clean.rds")
write_rds(marchant_counts, "0_data/processed/WIP_2025_marchant_counts.rds")
write_rds(fine_2025, "0_data/processed/WIP_2025_fine_clean.rds")
write_rds(umos_2025, "0_data/processed/WIP_2025_umos_clean.rds")

# End of script ----
