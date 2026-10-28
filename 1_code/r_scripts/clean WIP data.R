# ---
# title: Clean WIP data (any year)
# author: Emily H
# created: "2026-10-27"
# inputs: [0_data/external/WIP_<year>_Data Matrix_DPAC-SG.xlsx,
#          0_data/external/WIP_<year>_Coarse Level Community Analysis.xls(x),
#          0_data/external/Look-up Table_Invert Taxonomy_<date>.csv,
#          1_code/reference_tables/*.csv]
# outputs: [0_data/processed/WIP_<year>_*.rds,
#           2_pipeline/<year> taxa for expert review.csv]
# notes:
#   Year-generic version of "clean 2023 WIP data.R". To run a new year, put the
#   year's workbooks in 0_data/external, set wip_year (and lookup_file) in 1.1,
#   and run from the project root. Then run "resolve WIP data.R" with the same
#   settings.
#   The script stops (before writing anything) if sample IDs are duplicated.
#   Fix real ID errors in 1_code/reference_tables/sample_id_corrections.csv.
#   Taxa that are not in the lookup are listed in the expert review file; the
#   resolve script will not run until each one has a crosswalk row
#   (1_code/reference_tables/taxa_crosswalk.csv) or a lookup row.
# ---

# 1. Setup ----

## 1.1 Settings ----
# Keep wip_year and lookup_file the same in "resolve WIP data.R"
wip_year <- 2023
lookup_file <- "0_data/external/Look-up Table_Invert Taxonomy_2023-09-28.csv"

# Sheet names in the year's workbooks
sheet_fine <- "RAW Counts"        # DPAC-SG workbook, main (fine-level) counts
sheet_umos <- "UMOS Counts"       # DPAC-SG workbook, unique/mature search counts
sheet_coarse <- "Coarse Group Raw" # coarse workbook, counts by coarse group

## 1.2 Load packages ----
library(readxl) # read the Excel inputs
library(tidyverse)

## 1.3 Source helpers ----
source("1_code/r_scripts/functions/WIP pipeline helpers.R")

## 1.4 Locate inputs ----
find_input <- function(pattern, dir = "0_data/external") {
  hit <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(hit) != 1) {
    stop("Expected exactly one file in ", dir, " matching ", pattern,
         " but found ", length(hit))
  }
  hit
}
coarse_file <- find_input(sprintf("^WIP_%s_Coarse Level Community Analysis\\.xlsx?$", wip_year))
dpac_file <- find_input(sprintf("^WIP_%s_Data Matrix_DPAC-SG\\.xlsx$", wip_year))
stopifnot("lookup_file does not exist" = file.exists(lookup_file))

## 1.5 Import data ----
# lookup and crosswalk
lookup <- read_invert_lookup(lookup_file)
name_crosswalk <- read_taxa_crosswalk()

# coarse workbook (used for coarse_data and marchant_counts)
coarse_sheet <- read_excel(coarse_file, sheet = sheet_coarse)
stopifnot(
  "Coarse sheet has no 'Group/Sample Number' column" = "Group/Sample Number" %in% names(coarse_sheet),
  "Coarse sheet has no 'Primary Organisms' row" = "Primary Organisms" %in% coarse_sheet$`Group/Sample Number`,
  "Coarse sheet has no 'Marchant Cells Counted' row" = "Marchant Cells Counted" %in% coarse_sheet$`Group/Sample Number`
)

# coarse data
coarse_data <- coarse_sheet |>
  filter(!is.na(`Group/Sample Number`)) |>
  filter(row_number() < which(`Group/Sample Number` == "Primary Organisms")) |>
  select(!any_of(c("Mean", "Min", "Max", "Total"))) |>
  mutate(across(-`Group/Sample Number`, \(x) na_if(trimws(x), "no sample"))) |>
  mutate(`Group/Sample Number` = str_replace_all(`Group/Sample Number`, " ", ".")) |>
  column_to_rownames(var = "Group/Sample Number") |>
  t() |>
  as.data.frame()
str(coarse_data)

# Marchant cells counted per sample (includes the Mean/Min/Max/Total summary
# columns of the workbook; the resolve script removes them)
marchant_counts <- coarse_sheet |>
  filter(`Group/Sample Number` == "Marchant Cells Counted") |>
  select(-`Group/Sample Number`) |>
  pivot_longer(
    everything(),
    names_to = "Sample",
    values_to = "Marchant.Cells.Counted",
    values_transform = as.numeric
  )

# Sample ID corrections for this year (see notes above)
id_corrections <- read_csv(
  file.path(reference_dir, "sample_id_corrections.csv"),
  col_types = cols(year = "i", row_number = "i", .default = "c")
) |>
  filter(year == wip_year)

# fine and UMOS data: one row per sample, one column per taxon
read_dpac_sheet <- function(sheet, drop_cols = character()) {
  raw <- read_excel(dpac_file, sheet = sheet)

  # Corrections are matched on Row Number and Site so a stale entry cannot
  # be applied to the wrong row
  for (i in seq_len(nrow(id_corrections))) {
    hit <- which(raw$`Row Number` == id_corrections$row_number[i] &
                   raw$Site == id_corrections$site[i])
    if (length(hit) != 1) {
      stop("Sample ID correction for row ", id_corrections$row_number[i], " (",
           id_corrections$site[i], ") matches ", length(hit), " rows in sheet '", sheet, "'")
    }
    raw$Site[hit] <- id_corrections$corrected_site[i]
  }

  # CHECK: sample IDs become row names, so they must be present and unique
  if (anyNA(raw$Site)) {
    stop("Missing Site in sheet '", sheet, "', Row Number: ",
         paste(raw$`Row Number`[is.na(raw$Site)], collapse = ", "))
  }
  if (anyDuplicated(raw$Site)) {
    dup <- unique(raw$Site[duplicated(raw$Site)])
    stop("Duplicated Site in sheet '", sheet, "': ", paste(dup, collapse = ", "),
         " (Row Number ", paste(raw$`Row Number`[raw$Site %in% dup], collapse = ", "),
         "). Check against the coarse workbook; if it is a labelling error, add the fix to ",
         file.path(reference_dir, "sample_id_corrections.csv"))
  }

  out <- raw |>
    select(!any_of(c("Row Number", drop_cols))) |>
    column_to_rownames(var = "Site") |>
    mutate(across(where(is.numeric) | where(is.logical), \(x) replace(x, is.na(x), 0))) |>
    rename_with(\(x) gsub(" ", ".", x, fixed = TRUE)) |>
    rename_with(\(x) gsub("..", ".", x, fixed = TRUE))

  # CHECK: every column should be a count
  not_numeric <- names(out)[!vapply(out, is.numeric, logical(1))]
  if (length(not_numeric) > 0) {
    stop("Non-numeric column(s) in sheet '", sheet, "': ", paste(not_numeric, collapse = ", "))
  }
  out
}

fine_clean <- read_dpac_sheet(sheet_fine)
str(fine_clean)

umos_clean <- read_dpac_sheet(sheet_umos, drop_cols = "Total UMOS")
str(umos_clean)

# 2. Check for taxa not in lookup table ----
# Every fine / UMOS column that is not an Analysis_Name in the lookup is listed
# for expert review. handled_by says what the resolve script will do with it:
#   crosswalk      already has a row in taxa_crosswalk.csv (status shows whether
#                  the taxonomist has confirmed it)
#   analysis group is a group name (Oligochaeta, Amphipoda, ...) that matches as is
#   NEEDS REVIEW   new: send to the taxonomist, then add a crosswalk row
taxa_review <- imap(list(fine_not_in = fine_clean, umos_not_in = umos_clean), \(x, source) {
  flag_taxa_names(names(x), lookup, name_crosswalk) |>
    mutate(source = source, n_specimens = unname(colSums(x)))
}) |>
  list_rbind() |>
  filter(!in_lookup) |>
  left_join(name_crosswalk, by = c("column_name" = "raw_name")) |>
  mutate(n_specimens = if_else(column_name == row_total_col, NA, n_specimens)) |>
  select(source, column_name, handled_by, lookup_name, confidence, status, note, n_specimens) |>
  # NEEDS REVIEW first within each source
  arrange(source, handled_by != "NEEDS REVIEW")
count(taxa_review, source, handled_by)

if (any(taxa_review$handled_by == "NEEDS REVIEW")) {
  message(n_distinct(taxa_review$column_name[taxa_review$handled_by == "NEEDS REVIEW"]),
          " taxon name(s) in the ", wip_year, " data are not in the lookup or the crosswalk. ",
          "Send them to the taxonomist and add crosswalk rows before running the resolve script.")
}

# 3. Export cleaned data ----
# Read by "resolve WIP data.R"
dir.create("0_data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("2_pipeline", showWarnings = FALSE, recursive = TRUE)

write_csv(taxa_review, sprintf("2_pipeline/%s taxa for expert review.csv", wip_year), na = "")

processed_path <- \(name) sprintf("0_data/processed/WIP_%s_%s.rds", wip_year, name)
write_rds(coarse_data, processed_path("coarse_clean"))
write_rds(marchant_counts, processed_path("marchant_counts"))
write_rds(fine_clean, processed_path("fine_clean"))
write_rds(umos_clean, processed_path("umos_clean"))

# End of script ----
