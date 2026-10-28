# ---
# title: WIP pipeline helpers
# notes:
#   Shared by "clean WIP data.R" and "resolve WIP data.R" so that the names
#   flagged for expert review by the first script are exactly the names the
#   second script refuses to run on.
#   Needs readr and tibble (both attached by the tidyverse in the calling scripts).
# ---

# Column in the fine matrix that is a row total, not a taxon
row_total_col <- "Total.Primary.Organisms"

# Reference tables (tracked in git; edit these rather than the scripts)
reference_dir <- "1_code/reference_tables"

# Read the invertebrate taxonomy lookup and add Taxon_Analysis -------------
# Taxon_Analysis = Analysis_Grouping, except that parents marked
# "Distribute.to.Children" keep their own Analysis_Name.
read_invert_lookup <- function(file) {
  lookup <- read.csv(file, fileEncoding = "latin1")
  # Trailing non-breaking spaces (0xA0) in some Phylum / Order / Family /
  # Subfamily cells would split taxa when the algorithm matches on rank
  is_chr <- vapply(lookup, is.character, logical(1))
  lookup[is_chr] <- lapply(lookup[is_chr], \(x) trimws(gsub("\u00a0", " ", x, fixed = TRUE)))

  lookup$Taxon_Analysis <- lookup$Analysis_Grouping
  parent <- which(lookup$Taxon_Analysis == "Distribute.to.Children")
  lookup$Taxon_Analysis[parent] <- lookup$Analysis_Name[parent]
  lookup
}

# Read the name crosswalk --------------------------------------------------
# One row per raw name that is not an Analysis_Name in the lookup.
# lookup_name = NA means "leave out of the analysis".
read_taxa_crosswalk <- function(file = file.path(reference_dir, "taxa_crosswalk.csv")) {
  xw <- readr::read_csv(file, col_types = readr::cols(.default = "c"), na = c("", "NA"))
  needed <- c("raw_name", "lookup_name", "confidence", "status", "note")
  if (!all(needed %in% names(xw))) {
    stop("Crosswalk is missing column(s): ", paste(setdiff(needed, names(xw)), collapse = ", "))
  }
  if (anyDuplicated(xw$raw_name)) {
    stop("Duplicated raw_name in crosswalk: ",
         paste(unique(xw$raw_name[duplicated(xw$raw_name)]), collapse = ", "))
  }
  xw
}

# Flag taxon names against the lookup and crosswalk ------------------------
# handled_by (same precedence the resolve script uses):
#   crosswalk      raw name has a crosswalk row (applied first)
#   lookup         name is an Analysis_Name in the lookup
#   analysis group name is a Taxon_Analysis group (Oligochaeta, Amphipoda, ...)
#   NEEDS REVIEW   none of the above: send to the taxonomist, then add a
#                  crosswalk row (or a lookup row)
flag_taxa_names <- function(x, lookup, crosswalk) {
  handled_by <- dplyr::case_when(
    x %in% crosswalk$raw_name        ~ "crosswalk",
    x %in% lookup$Analysis_Name      ~ "lookup",
    x %in% lookup$Taxon_Analysis     ~ "analysis group",
    .default = "NEEDS REVIEW"
  )
  tibble::tibble(
    column_name = x,
    in_lookup   = x %in% lookup$Analysis_Name,
    handled_by  = handled_by
  )
}
