# ---
# title: Resolve 2023 WIP data
# author: Emily H
# created: "2026-10-21"
# inputs: [0_data/processed/WIP_2023_*.rds, 0_data/external/Look-up Table_Invert Taxonomy_2023-09-28.csv]
# outputs: [3_output/data/WIP 2023 resolved taxa_*.xlsx, 2_pipeline/WIP 2023 *.csv]
# notes: 
#   This script applies a taxonomic resolution algorithm to aquatic 
#     invertebrate data identified to varying resolutions.
#   Run "clean 2023 WIP data.R" first (it writes the cleaned inputs).
#   Open decisions are tagged DECISION:, data gaps DATA GAP: and checks CHECK:.
#   They are all listed in AGENTS.md.
# ---

# 1. Setup ----

## 1.1 Load packages ----
library(openxlsx)
library(abmi.themes)
library(tidyverse)
library(backports)
# Also required but not attached here: mefa4 (groupSums), data.table and
# vegan (attached by the functions sourced below)

## 1.2 Source functions ----
source("1_code/r_scripts/functions/Taxonomic ambiguous assignment abmi_clean.R")

## 1.3 Import data ----
# coarse data (QA only; the algorithm does not use it)
coarse_data <- read_rds("0_data/processed/WIP_2023_coarse_clean.rds")

# Marchant cells counted per sample
marchant_counts <- read_rds("0_data/processed/WIP_2023_marchant_counts.rds")

# fine data 
fine_2023 <- read_rds("0_data/processed/WIP_2023_fine_clean.rds")

# UMOS data
umos_2023 <- read_rds("0_data/processed/WIP_2023_umos_clean.rds")

## 1.4 Analysis settings ----
# DECISION: Marchant standardisation. The CABIN pipeline and the 2023 Excel
# workflow ("RAW 100 + UMOS" sheet) both scale fine counts to 100 cells and do
# not scale UMOS counts.
mc_standard <- 100

# DECISION: the 16 Hattie 5 C/D "invert"/"plant" samples have fine + UMOS data but
# no coarse data and no Marchant count in the 2023 coarse workbook. TRUE keeps
# them using the back-calculated counts in section 4.1; FALSE drops them.
include_hattie5_CD <- TRUE

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

# 2. Set up lookup table ----
#First step: USE "Analysis_Grouping" to aggregate "Analysis_Name";
#but keep taxon names with "Distribute to Children" in the Analysis_Grouping as is in raw data
#handle this in LOOK up table, add Taxon_Analysis

## EM'S NOTE: I need to update this to align with the primary and secondary taxonomic groups
# specified by the WIP project ----

# DECISION: using the 2023-09-28 lookup that is in this repo (the clean script
# uses it too). The 2025-07-28 version only adds 17 taxa that are not in the
# 2023 data and does not change grouping or rank for any 2023 taxon.
Invert_lookup <-  read.csv("0_data/external/Look-up Table_Invert Taxonomy_2023-09-28.csv",
                           fileEncoding = "latin1") |>
  # The file has trailing non-breaking spaces (0xA0) in some Phylum / Order /
  # Family / Subfamily cells, which would split taxa when the algorithm matches on rank
  mutate(across(where(is.character), \(x) trimws(gsub("\u00a0", " ", x, fixed = TRUE)))) |>
#Add a column with analysis name for ambiguous parents designated as "distribute to children" in analysis grouping
  mutate(Taxon_Analysis = Analysis_Grouping)
parent.id <- which(Invert_lookup$Taxon_Analysis=="Distribute.to.Children")
Invert_lookup$Taxon_Analysis[parent.id ] <- Invert_lookup$Analysis_Name[parent.id]

## 2.1 Lookup checks ----
# CHECK: Analysis_Name is duplicated in the lookup. match() takes the first row,
# so for ambiguous parents the first row's Rank decides which children are used.
lookup_dupes <- Invert_lookup |>
  mutate(lookup_row = row_number()) |>
  filter(Analysis_Name %in% Analysis_Name[duplicated(Analysis_Name)]) |>
  filter(Analysis_Name %in% c(names(fine_2023), names(umos_2023))) |>
  select(lookup_row, Analysis_Name, Final_ID, Rank, Analysis_Grouping)
lookup_dupes

# 3. Match taxon names to the lookup ----

## 3.1 Name crosswalk ----
# DATA GAP: 28 names in the fine data (26 in UMOS) are not Analysis_Name values
# in the lookup (see "2023 taxa for expert review.csv"). The mappings below are
# PROVISIONAL and need taxonomist sign-off. NA = leave out of the analysis.
#   high   = spelling variant, or a lookup row for the same ID exists
#   medium = follows a decision already made in the CABIN pipeline
#   low    = my judgement call, needs expert review
name_crosswalk <- tribble(
  ~raw_name,                  ~lookup_name,          ~confidence, ~note,
  "Lepidoptrea.UID",          "Lepidoptera.UID",     "high",   "typo",
  "Ceratopogonide",           "Ceratopogonidae.UID", "high",   "typo; coarse-level count",
  "Ylodes.ruteri",            "Ylodes.reuteri",      "high",   "typo",
  "Liodessus.affinus",        "Liodessus.affinis",   "high",   "typo",
  "Neoporus.superiorus",      "Neoporus.superioris", "high",   "typo",
  "Peltodytes",               "Peltodytes.sp.",      "high",   "genus-level name",
  "Chironomidae",             "Chironomidae.UID",    "high",   "coarse-level count; lookup parent, expect no children so kept as is",
  "Tipulidae",                "Tipulidae.UID",       "high",   "coarse-level count",
  "Culicidae",                "Culicidae.UID",       "high",   "coarse-level count",
  "Chaoboridae",              "Chaoboridae.UID",     "high",   "coarse-level count",
  "Other.Diptera",            "Other.Diptera.UID",   "high",   "coarse-level count",
  "Other:.Conchostraca",      "Conchostraca.UID",    "high",   "coarse-level count",
  "Total.Primary.Organisms",  NA,                    "high",   "row total, not a taxon (checked in 3.2)",
  "Coenagrion/Enallagma.sp.", "Coenagrionidae.UID",  "medium", "CABIN precedent; lookup Final_ID 'Coenagrion/Enallagma' maps here",
  "Hygrotus.sp.",             "Hygrotus.sp.other",   "medium", "CABIN precedent",
  "Terrestrial",              "Do.not.analyze",      "medium", "terrestrial specimens; dropped in 3.3",
  "Hygrotus.marklini",        "Hygrotus.sp.other",   "low",    "species not in lookup; extends the CABIN Hygrotus.sp. rule",
  "Hygrotus.sellatus",        "Hygrotus.sp.other",   "low",    "species not in lookup; extends the CABIN Hygrotus.sp. rule",
  "Berosus.hatchi",           "Berosus.sp.",         "low",    "species not in lookup",
  "Liodessus.adjectus.?",     "Liodessus.sp.",       "low",    "uncertain ID, species not in lookup",
  "Stictotarsus.?",           "Stictotarsus.sp.",    "low",    "uncertain ID",
  "Planorbella.morph",        "Planorbella.sp.",     "low",    "morphotype; assumed genus-level",
  "Sigara/Callicorixa.sp.",   "Corixidae.UID",       "low",    "genus pair; family-level parent is the closest lookup row",
  "Saldidae.UID",             NA,                    "low",    "no lookup entry; left out until one is added"
)
write_csv(name_crosswalk, "2_pipeline/2023 taxa crosswalk for expert review.csv")

## 3.2 Apply crosswalk ----
# CHECK: Total.Primary.Organisms should equal the sum of every other fine column
stopifnot(all(rowSums(select(fine_2023, -Total.Primary.Organisms)) == fine_2023$Total.Primary.Organisms))

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
fine_xw <- apply_crosswalk(fine_2023, name_crosswalk)
umos_xw <- apply_crosswalk(umos_2023, name_crosswalk)

# Specimens left out because the crosswalk is NA (excluding the row total)
bind_rows(fine = fine_xw$left_out, umos = umos_xw$left_out, .id = "source") |>
  select(-any_of("Total.Primary.Organisms"))

## 3.3 Aggregate to Taxon_Analysis ----
# Names that are already analysis groups (Oligochaeta, Hirudinea, Hydrachnidia,
# Amphipoda) match on Taxon_Analysis because they are not Analysis_Name values
to_taxon_analysis <- function(x) {
  i <- match(x, Invert_lookup$Analysis_Name)
  out <- Invert_lookup$Taxon_Analysis[i]
  is_group <- is.na(i) & x %in% Invert_lookup$Taxon_Analysis
  out[is_group] <- x[is_group]
  out
}

group_to_analysis <- function(x) {
  taxon_gp <- to_taxon_analysis(names(x))
  if (anyNA(taxon_gp)) {
    stop("Not in lookup: ", paste(names(x)[is.na(taxon_gp)], collapse = ", "))
  }
  as.data.frame(mefa4::groupSums(as.matrix(x), 2, by = as.character(taxon_gp), na.rm = TRUE))
}
fine_in <- group_to_analysis(fine_xw$data)
umos_in <- group_to_analysis(umos_xw$data)

# DECISION: "Do.not.analyze" is removed from main and UMOS data. The CABIN script
# only removed it from UMOS.
do_not_analyze <- list(fine = fine_in$Do.not.analyze, umos = umos_in$Do.not.analyze) |>
  map_dbl(\(x) sum(x, na.rm = TRUE))
do_not_analyze
fine_in <- select(fine_in, -any_of("Do.not.analyze"))
umos_in <- select(umos_in, -any_of("Do.not.analyze"))

# 4. Marchant cell standardisation ----

## 4.1 Marchant counts ----
# DATA GAP: no Marchant counts in the repo for the 16 Hattie 5 C/D samples.
# Values below were back-calculated from the 2023 Excel workflow ("RAW 100 + UMOS"
# sheet in WIP_2023_Data Matrix_DPAC-SG.xlsx): raw count / (scaled count - UMOS) * 100.
# The same calculation reproduces the coarse-workbook counts exactly for the
# other 79 samples. Confirm against the source counts.
mc_hattie5_CD <- tribble(
  ~Sample,                  ~Marchant.Cells.Counted,
  "Hattie 5-C1 spr invert", 6,
  "Hattie 5-C1 spr plant",  17,
  "Hattie 5-C2 spr invert", 6,
  "Hattie 5-C2 spr plant",  31,
  "Hattie 5-C3 spr invert", 5,
  "Hattie 5-C3 spr plant",  18,
  "Hattie 5-C4 spr invert", 5,
  "Hattie 5-C4 spr plant",  100,
  "Hattie 5-D1 spr invert", 22,
  "Hattie 5-D1 spr plant",  31,
  "Hattie 5-D2 spr invert", 21,
  "Hattie 5-D2 spr plant",  100,
  "Hattie 5-D3 spr invert", 10,
  "Hattie 5-D3 spr plant",  35,
  "Hattie 5-D4 spr invert", 5,
  "Hattie 5-D4 spr plant",  100
)

mc_table <- marchant_counts |>
  # marchant_counts also carries the Mean / Min / Max / Total summary columns from the workbook
  filter(!Sample %in% c("Mean", "Min", "Max", "Total")) |>
  bind_rows(mc_hattie5_CD)

# Samples used in the run
run_samples <- rownames(fine_in)
if (!include_hattie5_CD) {
  run_samples <- setdiff(run_samples, mc_hattie5_CD$Sample)
}

# CHECK: sample lists across the cleaned datasets
list(
  fine_not_umos     = setdiff(rownames(fine_in), rownames(umos_in)),
  umos_not_fine     = setdiff(rownames(umos_in), rownames(fine_in)),
  fine_not_coarse   = setdiff(rownames(fine_in), rownames(coarse_data)),
  coarse_not_fine   = setdiff(rownames(coarse_data), rownames(fine_in)),
  no_marchant_count = setdiff(run_samples, mc_table$Sample[!is.na(mc_table$Marchant.Cells.Counted)])
)

# CHECK: fine totals vs coarse totals
coarse_vs_fine <- tibble(
  sample = intersect(rownames(coarse_data), rownames(fine_2023)),
  coarse_total = rowSums(coarse_data[sample, ] |> mutate(across(everything(), as.numeric)), na.rm = TRUE),
  fine_total = fine_2023[sample, "Total.Primary.Organisms"]
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
Taxon_Lookup <- Invert_lookup
mmdata <- fn_matchdata(data = fine_bal, aux.data = umos_bal)
fine_bal <- mmdata$data
umos_bal <- mmdata$aux.data

# CHECK: NAs prevent the algorithm from running
c(fine = anyNA(fine_bal), umos = anyNA(umos_bal))

# CHECK: every taxon is in the lookup as Taxon_Analysis (the functions stop otherwise)
setdiff(colnames(fine_bal), Taxon_Lookup$Taxon_Analysis)

results <- pmap(run_grid, \(method, region, run) {
  match.fun(method)(
    data = fine_bal, aux.data = umos_bal, taxalookup = Taxon_Lookup,
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
# Each workbook has one sheet per output (Out.main, Out.aux, Out.combined, ...)
walk2(results, names(results), \(x, run) {
  write.xlsx(x, rowNames = TRUE,
             file = paste0("3_output/data/WIP 2023 resolved taxa_", run, ".xlsx"))
})
write_rds(results, "2_pipeline/WIP 2023 resolved taxa results.rds")
write_csv(mass_balance, "2_pipeline/WIP 2023 resolution mass balance.csv")

# End of script ----
