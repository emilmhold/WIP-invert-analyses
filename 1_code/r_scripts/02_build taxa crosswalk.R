# ---
# title: Build taxa crosswalk
# author: Emily H
# created: "2026-10-29"
# last edited: 2026-11-10
# inputs: [0_data/processed/WIP_<year>_fine_clean[_]<date>.rds (latest dated, else undated),
#          0_data/processed/WIP_<year>_umos_clean[_]<date>.rds (latest dated, else undated),
#          0_data/external/Look-up Table_Invert Taxonomy_<date>.xlsx (lookup_file),
#          1_code/reference_tables/taxa_crosswalk_11.05.2026.csv (annotated)]
# outputs: [1_code/reference_tables/taxa_crosswalk_11.05.2026.csv (rows appended),
#           1_code/reference_tables/taxa_unresolved<date>.csv]
# notes:
#   Proposes crosswalk rows for the taxon names in the fine and UMOS data of the
#   years in wip_years, and appends them to taxa_crosswalk_11.05.2026.csv. Run
#   "clean WIP data.R" first (it writes the cleaned data this script reads), then
#   this script, then "resolve WIP data.R".
#   Only names that are not already handled are considered: not an Analysis_Name
#   or analysis group in the lookup, and not already in the crosswalk. Existing
#   rows are never changed, so decisions made by hand or by the taxonomist stay.
#   Proposals come from the rules in section 2 (rank suffix, close spelling,
#   genus fallback). They are appended as status "provisional" with a note that
#   starts "auto-proposed", and need taxonomist sign-off like the other rows.
#   Names no rule can map (name pairs such as A/B.sp., names with no genus in the
#   lookup) are listed but NOT added: an empty lookup_name means "leave out", so
#   resolve keeps stopping on them until someone adds a row by hand.
#   The crosswalk is the annotated file (taxonomist notes in extra columns).
#   New rows are appended at the bottom of it with the same column order; the
#   annotation columns stay empty for them. Hand-added rows may have an empty
#   confidence or status.
#   Re-run "clean WIP data.R" afterwards to refresh the expert review file.
#   Run from the project root.
# ---

# 1. Setup ----
rm(list=ls())
gc()

## 1.1 Settings ----
# Years whose cleaned fine and UMOS data are scanned for names
wip_years <- c(2023, 2024, 2025)

# Same lookup as in "clean WIP data.R" and "resolve WIP data.R"
lookup_file <- "0_data/external/Look-up Table_Invert Taxonomy_2026-11-09.xlsx"

# TRUE: print the proposals but do not write the crosswalk
dry_run <- FALSE

## 1.2 Load packages and helpers ----
library(tidyverse)
source("1_code/r_scripts/functions/WIP pipeline helpers.R")

# Annotated crosswalk (taxonomist notes in extra columns). Read as is; missing
# names are appended at the bottom, existing rows and columns are not touched.
crosswalk_file <- file.path(reference_dir, "taxa_crosswalk_11.05.2026.csv")

# 2. Proposal rules ----
# Applied in this order to a name that is not handled yet; the first rule that
# gives a lookup name wins.
#   0. Row total column -> left out (the only automatic NA).
#   1. Name pairs ("A/B") -> no proposal.
#   2. Rank suffix: add ".sp." or ".UID" if that is an Analysis_Name (high).
#      "Other:." prefixes are dropped first ("Other:.Conchostraca").
#   3. Close spelling: one lookup name within 2 edits (and 20% of its length),
#      same genus for species names (medium).
#   4. Genus fallback: species-level name with no lookup row -> Genus.sp. (low).
# A trailing "?" (uncertain ID) caps the confidence at low.
base_name <- \(x) sub("\\.(sp\\.|UID)$", "", x)
genus_of <- \(x) sub("\\..*", "", x)

propose_mapping <- function(raw_name, lookup_names) {
  uncertain <- grepl("\\?$", raw_name)
  x <- sub("^Other:\\.", "", sub("\\.?\\?$", "", raw_name))
  pick <- function(target, confidence, why) {
    tibble(
      raw_name = raw_name,
      lookup_name = target,
      confidence = if (uncertain) "low" else confidence,
      note = paste0("auto-proposed: ", why, if (uncertain) "; uncertain ID in raw name" else "")
    )
  }

  if (raw_name == row_total_col) return(pick(NA_character_, "high", "row total, not a taxon"))
  if (grepl("/", raw_name, fixed = TRUE)) return(NULL)

  # Rank suffix
  hit <- intersect(paste0(x, c(".sp.", ".UID")), lookup_names)
  if (length(hit) > 0) return(pick(hit[1], "high", "rank suffix added"))

  # Close spelling
  b <- base_name(x)
  ub <- unique(base_name(lookup_names))
  d <- adist(b, ub)[1, ]
  dotted <- grepl(".", b, fixed = TRUE)
  ok <- d <= 2 & d / nchar(b) <= 0.2 &
    (grepl(".", ub, fixed = TRUE) == dotted) &
    (!dotted | genus_of(ub) == genus_of(b))
  best <- which(ok & d == suppressWarnings(min(d[ok])))
  if (length(best) == 1) {
    # Prefer the same rank suffix as the raw name, then the bare name, .sp., .UID
    raw_suffix <- substring(x, nchar(b) + 1)
    target <- intersect(paste0(ub[best], c(raw_suffix, "", ".sp.", ".UID")), lookup_names)[1]
    return(pick(target, "medium", sprintf("close spelling of %s (%d edit%s)",
                                          target, d[best], if (d[best] == 1) "" else "s")))
  }

  # Genus fallback
  if (dotted && paste0(genus_of(b), ".sp.") %in% lookup_names) {
    return(pick(paste0(genus_of(b), ".sp."), "low", "name not in lookup; assigned to the genus"))
  }
  NULL
}

# 3. Find names that need a decision ----
stopifnot("lookup_file does not exist" = file.exists(lookup_file))
lookup <- read_invert_lookup(lookup_file)

# Most recent dated file for the year (WIP_<year>_<name>[_]<yyyy-mm-dd>.rds);
# the 2025 files have no underscore before the date
processed_path <- function(year, name) {
  files <- list.files("0_data/processed",
                      pattern = sprintf("^WIP_%s_%s_?\\d{4}-\\d{2}-\\d{2}\\.rds$", year, name),
                      full.names = TRUE)
  if (length(files) == 0) return(NA_character_)
  sort(files, decreasing = TRUE)[1]
}
needed_files <- c(outer(wip_years, c("fine_clean", "umos_clean"), Vectorize(processed_path)))
stopifnot("Cleaned data not found; run 'clean WIP data.R' for each year in wip_years first" =
            !anyNA(needed_files))
message("Reading: ", paste(basename(needed_files), collapse = ", "))

# One row per name, year and source, with specimen counts
data_names <- expand_grid(year = wip_years, source = c("fine_clean", "umos_clean")) |>
  pmap(\(year, source) {
    x <- read_rds(processed_path(year, source))
    tibble(year = year, source = source, column_name = names(x), n_specimens = unname(colSums(x)))
  }) |>
  list_rbind()

crosswalk_empty <- tibble(raw_name = character(), lookup_name = character(),
                          confidence = character(), status = character(), note = character())
existing <- if (file.exists(crosswalk_file)) read_taxa_crosswalk(crosswalk_file) else crosswalk_empty

to_decide <- flag_taxa_names(unique(data_names$column_name), lookup, existing) |>
  filter(handled_by == "NEEDS REVIEW") |>
  pull(column_name)
message(length(to_decide), " name(s) in the ", paste(wip_years, collapse = ", "),
        " data are not in the lookup or the crosswalk.")

# 4. Propose rows ----
proposals <- map(to_decide, propose_mapping, lookup_names = lookup$Analysis_Name) |>
  list_rbind() |>
  mutate(status = "provisional") |>
  # keeps the columns when there are no proposals
  (\(x) bind_rows(crosswalk_empty, x))()

name_totals <- data_names |>
  summarise(n_specimens = sum(n_specimens), years = paste(sort(unique(year)), collapse = ", "),
            .by = column_name)

# Proposed rows, with specimen counts
proposals |>
  left_join(name_totals, by = c("raw_name" = "column_name")) |>
  arrange(confidence, desc(n_specimens))

# DATA GAP: names no rule could map. They are not added to the crosswalk, so
# resolve stops on them until a row is added by hand (ask the taxonomist).
unresolved <- name_totals |>
  filter(column_name %in% setdiff(to_decide, proposals$raw_name)) |>
  arrange(desc(n_specimens))
unresolved
if (nrow(unresolved) > 0) {
  message(nrow(unresolved), " name(s) have no proposal: add a row for each by hand.")
}

# 5. Checks and write ----
crosswalk_new <- bind_rows(existing, proposals)

# CHECK: the same format rules read_taxa_crosswalk() enforces, plus allowed
# values. Rows added by hand to the annotated file may have confidence or status
# still empty, so only filled-in values are checked.
stopifnot(
  "Duplicated raw_name" = !anyDuplicated(crosswalk_new$raw_name),
  "confidence must be high, medium or low" =
    all(crosswalk_new$confidence %in% c("high", "medium", "low", NA)),
  "status must be provisional or confirmed" =
    all(crosswalk_new$status %in% c("provisional", "confirmed", NA))
)
n_incomplete <- sum(is.na(crosswalk_new$confidence) | is.na(crosswalk_new$status))
if (n_incomplete > 0) {
  message(n_incomplete, " crosswalk row(s) have an empty confidence or status (hand-added).")
}

# CHECK: every lookup_name must exist in the lookup (as an Analysis_Name or a
# Taxon_Analysis group); resolve stops on any that do not
bad_targets <- crosswalk_new |>
  filter(!is.na(lookup_name), !lookup_name %in% c(lookup$Analysis_Name, lookup$Taxon_Analysis))
if (nrow(bad_targets) > 0) {
  stop("lookup_name not in the lookup (", basename(lookup_file), "): ",
       paste(bad_targets$lookup_name, collapse = ", "))
}

if (nrow(proposals) == 0) {
  message("No new crosswalk rows; ", crosswalk_file, " is unchanged.")
} else if (dry_run) {
  message("dry_run: ", nrow(proposals), " proposed row(s) not written.")
} else {
  # Append only the new rows, with the columns in the file's order (extra
  # annotation columns stay empty), so the existing lines are left byte for byte
  file_cols <- names(read_csv(crosswalk_file, n_max = 0, show_col_types = FALSE))
  new_rows <- proposals |>
    mutate(!!!set_names(rep(list(NA_character_), length(setdiff(file_cols, names(proposals)))),
                        setdiff(file_cols, names(proposals)))) |>
    select(all_of(file_cols))
  write_csv(new_rows, crosswalk_file, na = "", append = TRUE, eol = "\r\n")
  message("Appended ", nrow(proposals), " row(s) to ", crosswalk_file,
          " (", nrow(crosswalk_new), " in total). Review them, then re-run 'clean WIP data.R'.")
}

# 6. Save unresolved names ----
# One file shared across years: an open to-do list of names that still need a
# hand-made crosswalk row. Names seen in earlier runs are kept as they are, new
# names are appended, and names that now have a crosswalk row are dropped.
unresolved_file <- file.path(reference_dir, paste0("taxa_unresolved", Sys.Date(), ".csv"))

unresolved_empty <- tibble(raw_name = character(), years = character(), n_specimens = numeric())
unresolved_existing <- if (file.exists(unresolved_file)) {
  read_csv(unresolved_file, col_types = "ccd", na = "")
} else {
  unresolved_empty
}

unresolved_new <- unresolved |>
  transmute(raw_name = column_name, years, n_specimens) |>
  filter(!raw_name %in% unresolved_existing$raw_name)

unresolved_all <- bind_rows(unresolved_existing, unresolved_new) |>
  filter(!raw_name %in% crosswalk_new$raw_name)

if (nrow(unresolved_new) == 0 && nrow(unresolved_all) == nrow(unresolved_existing)) {
  message("No change to ", unresolved_file)
} else if (dry_run) {
  message("dry_run: ", nrow(unresolved_new), " unresolved name(s) not written.")
} else {
  dir.create(reference_dir, showWarnings = FALSE, recursive = TRUE)
  write_csv(unresolved_all, unresolved_file, na = "")
  message("Added ", nrow(unresolved_new), " name(s) to ", unresolved_file,
          " (", nrow(unresolved_all), " in total).")
}

# End of script ----
