# Metadata: WIP aquatic invertebrate resolved-taxa outputs

## 1. What these files are

Aquatic invertebrate counts from the ABMI / DUC Wetland Invertebrate (WIP) samples, with
**ambiguous taxa resolved to a common taxonomic level** so that samples can be compared.
Many specimens can only be identified to a coarse level (for example `Chironomidae.UID`,
`Corixinae.UID`). The resolution algorithms (Ermias Azeria's `aptc_abmi` and `dpac_abmi`)
reassign each such "parent" count to its identified "children" using the taxonomy lookup table.

One workbook = one year x one algorithm run. File names:

```
WIP <year> resolved taxa_<method>_<region>.xlsx
```

| Part | Values | Meaning |
|---|---|---|
| `year` | 2023, 2024, 2025 | Sampling year |
| `method` | `aptc_abmi` | Assign Parents To Children: at a site that has children, the whole parent count goes to the single most numerous child at that site |
| | `dpac_abmi` | Distribute Parents Among Children: the parent count is split among the children at that site in proportion to their counts |
| `region` | `frequent` | If a site has no children for a parent, use the child that is most widespread (found in the most sites) across the data (Chao2-like) |
| | `abundant` | Same fallback, but use the child with the highest total abundance (Chao1-like) |

All four combinations are exported because the final choice is still TBD.
The 2023 WIP Excel workflow was DPAC-SG; the CABIN pipeline used APTC-frequent.

Fallback regions: every WIP site is treated as Parkland natural region (one stratum). No
natural subregions are available, so the fallback uses all sites.

## 2. Sheets in each workbook

Every workbook has the same 12 sheets. **Rows are samples, columns are taxa**, except `AP.info`.
The first column has no header; it holds the sample ID (it was written from R row names).

| Sheet | Use it for | Contents |
|---|---|---|
| `Out.combined` | **Main analysis table** | Resolved fine data plus resolved UMOS. Equals `Out.main + Out.aux` |
| `Out.main` | | Resolved fine (main) data only, no UMOS |
| `Out.aux` | | Resolved UMOS data only |
| `Raw.main` | Audit | Input to the algorithm: fine counts after the name crosswalk and scaling, before resolution. Contains both children and ambiguous parents. **Not unscaled raw counts** |
| `Raw.aux` | Audit | Input UMOS counts (not scaled), same columns as `Raw.main` |
| `Raw.child` | Audit | Fine counts of taxa that were already resolved (non-parent columns of `Raw.main`) |
| `New.child` | Audit | Fine counts added to children by redistributing parents. `Raw.child + New.child` plus `Nochild.data` = `Out.main` |
| `Raw.child.aux` | Audit | Same as `Raw.child`, for UMOS |
| `New.child.aux` | Audit | Same as `New.child`, for UMOS |
| `AP.info` | QA | One row per ambiguous parent found in the data: `Parent` and the number of children found for it (`Number of children`; R's `read.xlsx` shows it as `Number.of.children`). Has a blank first column (row names). It is not samples x taxa |
| `Nochild.data` | QA | One column per ambiguous parent, filled only for parents with **no children**; fine counts that could not be redistributed and stay under the parent's name |
| `Nochild.data.aux` | QA | Same, for UMOS |

Why the column counts differ between sheets: `Raw.main` / `Raw.aux` have parents and children;
`Out.*` have no redistributed parents (only children, plus the no-child parents kept as they are).

## 3. Column (taxon) names

Names are the `Taxon_Analysis` values from the taxonomy lookup (`Look-up Table_Invert Taxonomy_<date>.xlsx`).

- `Genus.species`: species level (for example `Lymnaea.stagnalis`).
- `Genus.sp.`: genus level, not identified to species.
- `<Family or higher>.UID`: unidentified at that level (for example `Chironomidae.UID`, `Hemiptera.UID`).
- Names such as `Oligochaeta`, `Amphipoda`, `Hydrachnidia`, `Hirudinida` are analysis groups with no finer resolution.
- Raw names that were not in the lookup were mapped to a lookup name with the taxon crosswalk
  (section 5), so a column can hold specimens recorded under several raw names.

Ambiguous parents (lookup `Analysis_Grouping` = `Distribute.to.Children`) disappear from `Out.*`
unless they have no children; the numbers of columns are in section 4.

## 4. Quick facts per year

Counts are totals over all samples. "UMOS" = unique mature organisms search (a separate
search for mature specimens; its counts are not scaled).

| | 2023 | 2024 | 2025 |
|---|---|---|---|
| Samples (rows) | 95 | 100 | 100 |
| Taxa in `Out.*` | 77 | 92 | 86 |
| Ambiguous parents in data (`AP.info` rows) | 42 | 45 | 40 |
| of which no children (kept as is) | 17 | 12 | 11 |
| Fine (main) total, scaled | 671,078 | 891,607 | 280,815 |
| UMOS total | 398 | 520 | 554 |
| `Out.combined` total (APTC) | 671,476 | 892,127 | 281,369 |
| Most abundant taxon | `Chironomidae.UID` (186,349) | `Chironomidae.UID` (428,699) | `Chironomidae.UID` (105,224) |

Sample IDs:

- 2023: `Hattie 1-A spr`, wetland - replicate (A/B) - season (`spr`, `sum`). The Hattie 5 C/D samples are "invert/plant" samples.
- 2024: `Armitage-01-spr`, site - number - season.
- 2025: `218613-21-1_spring` / `_summer`. **TO CONFIRM**: meaning of the numeric parts.

## 5. How the numbers were produced

Scripts (run from the project root, in this order):
`1_code/r_scripts/01_generic clean WIP data.R` (2025: `01b_clean 2025 WIP data.R`), then
`02_build taxa crosswalk.R` if there are new names, then `03_generic resolve WIP data.R`.

1. **Clean**: read the fine, UMOS and Marchant-count tables; stop on missing or duplicated sample IDs.
2. **Name crosswalk**: names not in the lookup are mapped using
   `1_code/reference_tables/taxa_crosswalk_11.05.2026.csv` (columns: `raw_name`, `lookup_name`,
   `confidence` high / medium / low, `status` provisional / confirmed). A blank `lookup_name` means
   the specimens are left out (for example `Terrestrial`, row totals).
3. **Aggregate** raw names to `Taxon_Analysis` and drop `Do.not.analyze`.
4. **Scale fine counts to 100 Marchant cells**: `round(count / cells counted x 100)`. UMOS counts are not scaled.
5. **Resolve** parents with the algorithm in `Taxonomic ambiguous assignment abmi_UMOSfix.R`
   (a draft fix of the original, see section 7).
6. Export one workbook per run.

## 6. How to read and use the data

- **Values are standardised counts, not raw counts and not densities.** Fine counts were scaled
  to a 100-cell Marchant subsample and rounded. A sample where only 5 cells were counted is
  multiplied by 20, so rare taxa in such samples are inflated.
- Values should be whole numbers. Totals match the input exactly for APTC. For DPAC, rounding
  makes totals differ slightly from the input (2023: +9 specimens over 95 samples, at most 1 per sample;
  2024: +2; 2025: 0). See `2_pipeline/WIP <year> resolution mass balance.csv`.
- Use `Out.combined` for community analyses. It includes UMOS specimens.
- Taxa listed in `Nochild.data` (for example `Chironomidae.UID`, `Enallagma.sp.` in 2025) are parents with no
  children in the lookup, so their counts stay at the parent level. `Chironomidae.UID` is the largest part of the
  total in every year (about 28% in 2023, 48% in 2024, 37% in 2025), so be careful with totals and richness that include it.
- **Choice of run changes the numbers.** In 2025, compared with `aptc_abmi_frequent`, `aptc_abmi_abundant`
  differs in 70 cells of `Out.combined`, `dpac_abmi_frequent` in 98 and `dpac_abmi_abundant` in 160
  (of 8,600). Parents with a single child send 100% to that child under every run. Use the same run for all years.
- Samples are not pooled: A/B replicates (2023) are separate rows.

Reading in R:

```r
library(openxlsx)
x <- read.xlsx("3_output/data/WIP 2025 resolved taxa_aptc_abmi_frequent.xlsx",
               sheet = "Out.combined", rowNames = TRUE)
# readxl: the first column is the sample ID; rename it, for example
# readxl::read_xlsx(path, sheet = "Out.combined") |> dplyr::rename(sample = 1)
```