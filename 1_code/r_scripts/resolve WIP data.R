# ---
# title: Resolve WIP data (any year)
# author: Emily H
# created: "2026-10-27"
# inputs: [0_data/processed/WIP_<year>_*.rds,
#          0_data/external/Look-up Table_Invert Taxonomy_<date>.csv,
#          1_code/reference_tables/*.csv]
# outputs: [3_output/data/WIP <year> resolved taxa_*.xlsx,
#           2_pipeline/WIP <year> *.csv / *.rds,
#           2_pipeline/<year> taxa crosswalk for expert review.csv]
# notes:
#   Year-generic version of "resolve 2023 WIP data.R". Run "clean WIP data.R"
#   first with the same wip_year and lookup_file, then run this from the
#   project root.
#   This script applies a taxonomic resolution algorithm to aquatic
#     invertebrate data identified to varying resolutions.
#   It stops if a taxon name is not in the lookup or the crosswalk
#     (1_code/reference_tables/taxa_crosswalk.csv), or if a sample has no
#     Marchant count and drop_samples_without_marchant is FALSE.
#   Open decisions are tagged DECISION:, data gaps DATA GAP: and checks CHECK:.
#   They are all listed in AGENTS.md.
# ---

# 1. Setup ----

## 1.1 Load packages ----
library(openxlsx)
library(tidyverse)
# Also required but not attached here: mefa4 (groupSums), data.table and
# vegan (attached by the functions sourced below)

## 1.2 Settings ----
# Keep wip_year and lookup_file the same in "clean WIP data.R"
wip_year <- 2023
lookup_file <- "0_data/external/Look-up Table_Invert Taxonomy_2023-09-28.csv"

# DECISION: Marchant standardisation. The CABIN pipeline and the 2023 Excel
# workflow ("RAW 100 + UMOS" sheet) both scale fine counts to 100 cells and do
# not scale UMOS counts.
mc_standard <- 100

# DECISION: samples with fine + UMOS data but no Marchant count in the coarse
# workbook. Counts for them can be supplied in
# 1_code/reference_tables/marchant_counts_supplement.csv (the 16 Hattie 5 C/D
# "invert"/"plant" samples in 2023 are there, back-calculated).
#   use_marchant_supplement = FALSE ignores that file.
#   drop_samples_without_marchant = TRUE drops samples still without a count;
#   FALSE stops with the list of samples.
use_marchant_supplement <- TRUE
drop_samples_without_marchant <- FALSE

# DECISION: which algorithm / regional rule. Jenet chose APTC-frequent for CABIN;
# the 2023 Excel workflow used DPAC-SG. Running all four until confirmed.
run_grid <- expand_grid(
  method = c("aptc_abmi", "dpac_abmi"),
  region = c("frequent", "abundant")
) |>
  mutate(run = paste(method, region, sep = "_"))

# DATA GAP: no climate / natural region data for the WIP sites in this repo, so
# regional fallbacks are study-wide (strata = NULL), like the Chao1 x Chao2
# step in the Excel workflow.
strata_primary <- NULL
strata_secondary <- NULL

## 1.3 Source functions ----
source("1_code/r_scripts/functions/Taxonomic ambiguous assignment abmi_clean.R")
source("1_code/r_scripts/functions/WIP pipeline helpers.R")

## 1.4 Import data ----
processed_path <- \(name) sprintf("0_data/processed/WIP_%s_%s.rds", wip_year, name)
stopifnot(
  "Cleaned data not found for wip_year; run 'clean WIP data.R' first" =
    all(file.exists(processed_path(c("coarse_clean", "marchant_counts", "fine_clean", "umos_clean")))),
  "lookup_file does not exist" = file.exists(lookup_file)
)

# coarse data (QA only; the algorithm does not use it)
coarse_data <- read_rds(processed_path("coarse_clean"))

# Marchant cells counted per sample
marchant_counts <- read_rds(processed_path("marchant_counts"))

# fine data
fine_clean <- read_rds(processed_path("fine_clean"))

# UMOS data
umos_clean <- read_rds(processed_path("umos_clean"))

# 2. Set up lookup table ----
# First step: USE "Analysis_Grouping" to aggregate "Analysis_Name";
# but keep taxon names with "Distribute to Children" in the Analysis_Grouping as is in raw data.
# read_invert_lookup() adds Taxon_Analysis for this (see the helpers file).

## EM'S NOTE: I need to update this to align with the primary and secondary taxonomic groups
# specified by the WIP project ----

# DECISION: lookup version is set in lookup_file (1.1). 2023-09-28 and
# 2025-07-28 give identical grouping and rank for every 2023 taxon; the 2025
# file adds taxa not in the 2023 data. Check new taxa in later years against
# the version used.
Invert_lookup <- read_invert_lookup(lookup_file)

## 2.1 Lookup checks ----
# CHECK: Analysis_Name is duplicated in the lookup. match() takes the first row,
# so for ambiguous parents the first row's Rank decides which children are used.
lookup_dupes <- Invert_lookup |>
  mutate(lookup_row = row_number()) |>
  filter(Analysis_Name %in% Analysis_Name[duplicated(Analysis_Name)]) |>
  filter(Analysis_Name %in% c(names(fine_clean), names(umos_clean))) |>
  select(lookup_row, Analysis_Name, Final_ID, Rank, Analysis_Grouping)
lookup_dupes

# 3. Match taxon names to the lookup ----

## 3.1 Name crosswalk ----
# DATA GAP: names in the data that are not Analysis_Name values in the lookup are
# mapped with 1_code/reference_tables/taxa_crosswalk.csv. Rows with status
# "provisional" still need taxonomist sign-off. lookup_name empty = leave out of
# the analysis.
#   confidence: high   = spelling variant, or a lookup row for the same ID exists
#               medium = follows a decision already made in the CABIN pipeline
#               low    = judgement call, needs expert review
name_crosswalk <- read_taxa_crosswalk()

# CHECK: every name is in the lookup, is an analysis group, or has a crosswalk row
taxa_flags <- flag_taxa_names(unique(c(names(fine_clean), names(umos_clean))),
                              Invert_lookup, name_crosswalk)
needs_review <- filter(taxa_flags, handled_by == "NEEDS REVIEW")
if (nrow(needs_review) > 0) {
  stop(nrow(needs_review), " taxon name(s) are not in the lookup or the crosswalk: ",
       paste(needs_review$column_name, collapse = ", "),
       ".\nThey are listed in '2_pipeline/", wip_year, " taxa for expert review.csv' ",
       "(written by 'clean WIP data.R'). Once the taxonomist has decided, add a row for each to ",
       file.path(reference_dir, "taxa_crosswalk.csv"), " (or add the taxon to the lookup).")
}

# Crosswalk rows used by this year's data, with specimen counts, for review
crosswalk_used <- name_crosswalk |>
  filter(raw_name %in% c(names(fine_clean), names(umos_clean))) |>
  mutate(
    n_fine = replace_na(unname(colSums(fine_clean)[raw_name]), 0),
    n_umos = replace_na(unname(colSums(umos_clean)[raw_name]), 0),
    # the row total is not a taxon
    across(c(n_fine, n_umos), \(x) if_else(raw_name == row_total_col, NA, x))
  )
count(crosswalk_used, status, confidence)

## 3.2 Apply crosswalk ----
# CHECK: the row total should equal the sum of every other fine column
stopifnot(
  "Fine data has no row total column" = row_total_col %in% names(fine_clean),
  "Row total column does not equal the sum of the other fine columns" =
    all(rowSums(select(fine_clean, -all_of(row_total_col))) == fine_clean[[row_total_col]])
)

# Rename columns using the crosswalk; columns that end up with the same name are summed
apply_crosswalk <- function(x, crosswalk) {
  hit <- match(names(x), crosswalk$raw_name)
  new_name <- ifelse(is.na(hit), names(x), crosswalk$lookup_name[hit])
  keep <- !is.na(new_name)
  out <- rowsum(t(as.matrix(x[, keep, drop = FALSE])), group = new_name[keep], reorder = FALSE)
  list(
    data = as.data.frame(t(out)),
    left_out = colSums(x[, !keep, drop = FALSE])
  )
}
fine_xw <- apply_crosswalk(fine_clean, name_crosswalk)
umos_xw <- apply_crosswalk(umos_clean, name_crosswalk)

# Specimens left out because the crosswalk is NA (excluding the row total)
bind_rows(fine = fine_xw$left_out, umos = umos_xw$left_out, .id = "source") |>
  select(-any_of(row_total_col))

## 3.3 Aggregate to Taxon_Analysis ----
# Names that are already analysis groups (Oligochaeta, Hirudinea, Hydrachnidia,
# Amphipoda) match on Taxon_Analysis because they are not Analysis_Name values
to_taxon_analysis <- function(x, lookup) {
  i <- match(x, lookup$Analysis_Name)
  out <- lookup$Taxon_Analysis[i]
  is_group <- is.na(i) & x %in% lookup$Taxon_Analysis
  out[is_group] <- x[is_group]
  out
}

group_to_analysis <- function(x, lookup) {
  taxon_gp <- to_taxon_analysis(names(x), lookup)
  if (anyNA(taxon_gp)) {
    stop("Not in lookup: ", paste(names(x)[is.na(taxon_gp)], collapse = ", "))
  }
  as.data.frame(mefa4::groupSums(as.matrix(x), 2, by = as.character(taxon_gp), na.rm = TRUE))
}
fine_in <- group_to_analysis(fine_xw$data, Invert_lookup)
umos_in <- group_to_analysis(umos_xw$data, Invert_lookup)

# DECISION: "Do.not.analyze" is removed from main and UMOS data. The CABIN script
# only removed it from UMOS.
do_not_analyze <- list(fine = fine_in$Do.not.analyze, umos = umos_in$Do.not.analyze) |>
  map_dbl(\(x) sum(x, na.rm = TRUE))
do_not_analyze
fine_in <- select(fine_in, -any_of("Do.not.analyze"))
umos_in <- select(umos_in, -any_of("Do.not.analyze"))

# 4. Marchant cell standardisation ----

## 4.1 Marchant counts ----
# DATA GAP (2023): no Marchant counts in the repo for the 16 Hattie 5 C/D samples.
# Values in marchant_counts_supplement.csv were back-calculated from the 2023
# Excel workflow ("RAW 100 + UMOS" sheet in WIP_2023_Data Matrix_DPAC-SG.xlsx):
# raw count / (scaled count - UMOS) * 100. The same calculation reproduces the
# coarse-workbook counts exactly for the other 79 samples. Confirm against the
# source counts.
mc_supplement <- read_csv(
  file.path(reference_dir, "marchant_counts_supplement.csv"),
  col_types = cols(year = "i", Sample = "c", Marchant.Cells.Counted = "d", .default = "c")
) |>
  filter(year == wip_year, use_marchant_supplement)

# marchant_counts also carries the Mean / Min / Max / Total summary columns from the workbook
mc_workbook <- marchant_counts |>
  filter(!Sample %in% c("Mean", "Min", "Max", "Total"))

# CHECK: the supplement fills gaps; it must not disagree with a workbook count
mc_conflicts <- mc_workbook |>
  filter(!is.na(Marchant.Cells.Counted), Sample %in% mc_supplement$Sample)
if (nrow(mc_conflicts) > 0) {
  stop("marchant_counts_supplement.csv has samples that already have a count in the coarse workbook: ",
       paste(mc_conflicts$Sample, collapse = ", "))
}

mc_table <- mc_workbook |>
  filter(!Sample %in% mc_supplement$Sample) |>
  bind_rows(select(mc_supplement, Sample, Marchant.Cells.Counted))

# Samples without a Marchant count
no_marchant_count <- setdiff(
  rownames(fine_in),
  mc_table$Sample[!is.na(mc_table$Marchant.Cells.Counted)]
)
if (length(no_marchant_count) > 0) {
  if (!drop_samples_without_marchant) {
    stop(length(no_marchant_count), " sample(s) have fine data but no Marchant count: ",
         paste(no_marchant_count, collapse = ", "),
         ".\nAdd counts to ", file.path(reference_dir, "marchant_counts_supplement.csv"),
         " or set drop_samples_without_marchant <- TRUE.")
  }
  message("Dropping ", length(no_marchant_count), " sample(s) without a Marchant count: ",
          paste(no_marchant_count, collapse = ", "))
}

# Samples used in the run
run_samples <- setdiff(rownames(fine_in), no_marchant_count)

# CHECK: sample lists across the cleaned datasets
sample_checks <- list(
  fine_not_umos     = setdiff(rownames(fine_in), rownames(umos_in)),
  umos_not_fine     = setdiff(rownames(umos_in), rownames(fine_in)),
  fine_not_coarse   = setdiff(rownames(fine_in), rownames(coarse_data)),
  coarse_not_fine   = setdiff(rownames(coarse_data), rownames(fine_in)),
  no_marchant_count = no_marchant_count
)
sample_checks

# The UMOS table is indexed by the fine samples below, so a missing sample would become an NA row
missing_umos <- setdiff(run_samples, rownames(umos_in))
if (length(missing_umos) > 0) {
  stop("Samples in the fine data but not in the UMOS data: ", paste(missing_umos, collapse = ", "))
}

# CHECK: fine totals vs coarse totals
coarse_vs_fine <- tibble(
  sample = intersect(rownames(coarse_data), rownames(fine_clean)),
  coarse_total = rowSums(coarse_data[sample, ] |> mutate(across(everything(), as.numeric)), na.rm = TRUE),
  fine_total = fine_clean[sample, row_total_col]
) |>
  filter(coarse_total != fine_total)
coarse_vs_fine

## 4.2 Scale fine data ----
mc <- mc_table$Marchant.Cells.Counted[match(run_samples, mc_table$Sample)]
stopifnot(!anyNA(mc), all(mc > 0))

fine_bal <- fine_in[run_samples, ]
fine_bal <- round(fine_bal / mc * mc_standard)
umos_bal <- umos_in[run_samples, ]

# 5. Run the algorithm ----
# Match main and auxiliary (UMOS) data so both have the same sites and taxa
mmdata <- fn_matchdata(data = fine_bal, aux.data = umos_bal)
fine_bal <- mmdata$data
umos_bal <- mmdata$aux.data

# CHECK: NAs prevent the algorithm from running
c(fine = anyNA(fine_bal), umos = anyNA(umos_bal))

# CHECK: every taxon is in the lookup as Taxon_Analysis (the functions stop otherwise)
setdiff(colnames(fine_bal), Invert_lookup$Taxon_Analysis)

results <- pmap(run_grid, \(method, region, run) {
  match.fun(method)(
    data = fine_bal, aux.data = umos_bal, taxalookup = Invert_lookup,
    strata.primary = strata_primary, strata.secondary = strata_secondary,
    method.region = region, reg.prop = FALSE
  )
}) |>
  set_names(run_grid$run)

# 6. Check results ----

## 6.1 Ambiguous parents ----
# Number of children found for each parent (the same for every run)
ap_info <- results[[1]]$AP.info |>
  mutate(n_children = as.integer(`Number of children`)) |>
  select(Parent, n_children) |>
  arrange(n_children)
ap_info

# CHECK: parents with no children keep their own (scaled) counts
colSums(results[[1]]$Nochild.data) |>
  keep(\(x) x > 0)

## 6.2 Mass balance ----
# CHECK: total counts in = total counts out. Main data should balance (DPAC can
# be off by rounding). Out.combined is short by the UMOS counts of parents with
# 0 or 2+ children, because the functions never add the redistributed UMOS
# counts back (X.Abu.spp is calculated but not used). Out.aux also contains the
# main data's no-child parents, so use Out.main and Out.combined only.
# DECISION: leave the functions as they are, or fix them (see AGENTS.md).
aux_parents_dropped <- results[[1]]$Raw.aux[, ap_info$Parent[ap_info$n_children != 1], drop = FALSE] |>
  colSums() |>
  keep(\(x) x > 0)
aux_parents_dropped

mass_balance <- imap(results, \(x, run) {
  tibble(
    run = run,
    main_in = sum(x$Raw.main), main_out = sum(x$Out.main),
    max_site_diff_main = max(abs(rowSums(x$Out.main) - rowSums(x$Raw.main))),
    combined_in = sum(x$Raw.main) + sum(x$Raw.aux), combined_out = sum(x$Out.combined),
    combined_gap = combined_in - combined_out
  )
}) |>
  list_rbind()
mass_balance

# 7. Export results ----
dir.create("2_pipeline", showWarnings = FALSE, recursive = TRUE)
dir.create("3_output/data", showWarnings = FALSE, recursive = TRUE)

write_csv(crosswalk_used, sprintf("2_pipeline/%s taxa crosswalk for expert review.csv", wip_year), na = "")

# Each workbook has one sheet per output (Out.main, Out.aux, Out.combined, ...)
walk2(results, names(results), \(x, run) {
  write.xlsx(x, rowNames = TRUE,
             file = sprintf("3_output/data/WIP %s resolved taxa_%s.xlsx", wip_year, run))
})
write_rds(results, sprintf("2_pipeline/WIP %s resolved taxa results.rds", wip_year))
write_csv(mass_balance, sprintf("2_pipeline/WIP %s resolution mass balance.csv", wip_year))

# End of script ----
