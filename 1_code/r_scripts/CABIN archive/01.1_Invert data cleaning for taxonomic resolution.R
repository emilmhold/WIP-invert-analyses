## Clean and prep aquatic invert data for Ermias' algorithm
## Author: Emily H
## Created: May 22, 2025
## Last edited: July 3, 2025

#install.packages("readxl")
#install.packages("tidyverse")
#install.packages("abmi.themes")

library(readxl)
library(abmi.themes)
library(tidyverse)

setwd("~/GitHub/New aquatic invertebrate analyses")

#### import and clean 2018 CABIN coarse data ####
coarse.2018 <- read_csv("data/Invert Sorting_Coarse Summary_2018-CABIN.csv")
coarse.2018 <- coarse.2018[,1:29]
coarse.2018 <- coarse.2018 %>%
    filter(!Group %in% c("Other =") & !is.na(Group)) %>% #get rid of non-taxa row values
    column_to_rownames("Group") %>%
    t() %>% #transpose
    as.data.frame() %>%
    dplyr::select(-c("Damaged Organisms", "Primary Organisms", "Total Organisms", "Unique/Mature Search", "Other:", "Secondary Organisms", "Groups")) %>% #remove extraneous variables
    rownames_to_column(var = "Site") %>% #to match Ermias' data formatting
    mutate(Year = "2018") %>%
    mutate(SiteYear = paste0(Site, "_", Year)) %>%
    rename(MC_counted = "Marchant Cells Counted") %>% #create and rename columns to match
    dplyr::select(Site, Year, SiteYear, MC_counted, everything()) %>% #reorder columns to match
    mutate(across(4:37, as.numeric)) %>% #correct variables
    mutate(Year = as.integer(Year))
rownames(coarse.2018) <- coarse.2018$SiteYear #set rownames
str(coarse.2018)

##import and clean CABIN fine data
fine.2018 <- read.csv("data/Invert Sorting_Advanced ID_2018-CABIN-Complete +Chironomidae (2).csv") %>%
    dplyr::select(c(Site.Number, Confirmed.Genus...Species, Number.of.Specimens)) %>% #select columns I want
    rename(Site = Site.Number,
           Confirmed.Genus.Species = Confirmed.Genus...Species,
           specimens = Number.of.Specimens) %>% #rename for ease of use
    filter(Confirmed.Genus.Species != "") %>% #filter for advanced IDs to species
    mutate(specimens = as.numeric(specimens)) %>%
    pivot_wider(names_from = Confirmed.Genus.Species, values_from = specimens, values_fn = sum) %>% #reshape into wide format
    filter(Site != "") %>% #remove an empty row
    mutate(Year = "2018", .after = 1) %>% #add year data
    mutate(SiteYear = paste0(Site, "_", Year), .after = 2) %>% #create identifying column
    rename_with(~ str_replace_all(., " ", ".")) %>% #replace spaces with .
    mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>% # replace blanks with zeroes
    dplyr::select(
        -all_of(names(.)[4:102]),                        # remove target columns
        all_of(sort(names(.)[4:102])),                   # re-add them alphabetically
        everything())  %>%                                  # keep remaining columns
    mutate(Year = as.integer(Year))
rownames(fine.2018) <- fine.2018$SiteYear #set rownames
str(fine.2018)

##import and clean CABIN UMOS data
umos.2018 <- read.csv("data/Invert Sorting_Advanced ID_2018-CABIN-Complete +Chironomidae (2).csv") %>%
    filter(Coarse.Sorting.Group == "UMOS") %>%
    dplyr::select(Site.Number, Family...........Sub.family.for.Chironomidae., Genus...Species, Number.of.Specimens) %>%
    rename(Site = Site.Number,
           Family = Family...........Sub.family.for.Chironomidae.,
           Genus.Species = Genus...Species,
           specimens = Number.of.Specimens) %>%
    filter(Family != "") %>% #remove an empty row
    mutate(specimens = as.numeric(specimens)) %>%
    mutate(New.genus.species = if_else(Genus.Species == "UID", paste0(Family, " UID"), Genus.Species)) %>%
    filter(New.genus.species != "UID UID") %>% #remove because this is not helpful
    dplyr::select(-c(Family, Genus.Species)) %>%
    pivot_wider(names_from = New.genus.species, values_from = specimens, values_fn = sum) %>% #reshape into wide format
    mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>% # replace blanks with zeroes
    mutate(Year = "2018", .after = 1) %>% #add year data
    mutate(SiteYear = paste0(Site, "_", Year), .after = 2) %>% #create identifying column
    rename_with(~ str_replace_all(., " ", ".")) %>% #replace spaces with .
    mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>% # replace blanks with zeroes
    dplyr::select(
        -all_of(names(.)[4:80]),                        # remove target columns
        all_of(sort(names(.)[4:80])),                   # re-add them alphabetically
        everything()) %>%
    mutate(Year = as.integer(Year))
rownames(umos.2018) <- umos.2018$SiteYear #set rownames
str(umos.2018)

## pull out CABIN site names for 2018
site.names <- coarse.2018 %>%
    dplyr::select(Site, Year) %>%
    rename(CABIN.Site = Site) %>%
    mutate(Site = sub("^CABIN-", "", CABIN.Site)) %>% #remove CABIN prefix
    dplyr::select(!CABIN.Site)
str(site.names)

#### import and clean 2018 ABMI data ####
##import ABMI 2018 coarse data
species.2018 <- read.csv("data/2018 umos aquatic invert data.csv")
MCN.2018 <- read.csv("data/Marchant_cell_counts.csv") %>%
    filter(Year == 2018) %>%
    rename(Site = Site.Name,
           MC_counted = COUNT...)

##Coarse data
ABMI.coarse.2018 <- species.2018 %>%
    dplyr::select(ABMI.Site, Coarse.Group, Specimen.Count) %>%
    rename(Site = ABMI.Site) %>%
    filter(Coarse.Group != "UMOS") %>%
    pivot_wider(names_from = Coarse.Group, values_from = Specimen.Count, values_fn = sum) %>%
    mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>% #replace blanks with zeroes
    left_join(MCN.2018, by = "Site") %>%
    mutate(Year = ifelse(is.na(Year), "2018", Year)) %>% #make sure all the year values are filled in
    mutate(SiteYear = paste0(Site, "_", Year)) %>% #create identifying column
    mutate(Year = as.integer(Year)) %>%
    dplyr::select(Site, Year, SiteYear, MC_counted, everything()) %>%
    filter(Site %in% site.names$Site) #include only sites sampled with the CABIN protocol
rownames(ABMI.coarse.2018) <- ABMI.coarse.2018$SiteYear
str(ABMI.coarse.2018)
#
# ABMI.coarse.2018 <- read.csv("data/W07A_ all coarse data.csv") %>%
#     dplyr::select(ABMI.Site, Year, Coarse.Group, Coarse.Group.Count) %>%
#     filter(Year == 2018) %>%
#     filter(ABMI.Site %in% site.names$ABMI.Site) %>%
#     mutate(Site = paste0("W",ABMI.Site), .before = 1)%>%
#     dplyr::select(!ABMI.Site) %>%
#     pivot_wider(names_from = Coarse.Group, values_from = Coarse.Group.Count, values_fn = sum) %>%
#     mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>% #replace blanks with zeroes
#     mutate(SiteYear = paste0(Site, "_", Year), .after = 2) %>%
#     left_join(MCN.2018, by = c("Site", "Year"))
# str(ABMI.coarse.2018)

##import ABMI 2018 fine data
# ABMI.fine.2018 <- read.csv("data/Invert Sorting_Advanced ID_2018 ABMI Non-targets (1).csv") %>%
#     dplyr::select(c(Site.Number, Genus...Species.1, Number.of.Specimens)) %>% #select columns I want
#     rename(Site = Site.Number,
#            Confirmed.Genus.Species = Genus...Species.1,
#            specimens = Number.of.Specimens) %>% #rename for ease of use
#     filter(Confirmed.Genus.Species != "") %>% #filter for advanced IDs to species
#     mutate(specimens = as.numeric(specimens)) %>%
#     pivot_wider(names_from = Confirmed.Genus.Species, values_from = specimens, values_fn = sum) %>% #reshape into wide format
#     filter(Site != "") %>% #remove an empty row
#     mutate(Year = "2018", .after = 1) %>% #add year data
#     mutate(SiteYear = paste0(Site, "_", Year), .after = 2) %>% #create identifying column
#     rename_with(~ str_replace_all(., " ", ".")) %>% #replace spaces with .
#     mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>% # replace blanks with zeroes
#     dplyr::select(!UID) %>%
#     mutate(Year = as.integer(Year))
# str(ABMI.fine.2018) ##note: only 22 species? Check in with Cheryl

##ABMI 2018 Fine data
ABMI.fine.2018 <- species.2018 %>%
    dplyr::select(ABMI.Site, Genus, Species, Specimen.Count) %>%
    rename(Site = ABMI.Site) %>%
    mutate(Species.ID = paste0(Genus,".",Species)) %>%
    dplyr::select(-c(Genus, Species)) %>%
    pivot_wider(names_from = Species.ID, values_from = Specimen.Count, values_fn = sum) %>%
    mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>% #replace blanks with zeroes
    mutate(Year = 2018, .after = Site) %>%
    mutate(SiteYear = paste0(Site, "_", Year), .after = Year) %>% #create identifying column
    dplyr::select(!.UID) %>%
    mutate(Year = as.integer(Year)) %>%
    filter(Site %in% site.names$Site) #include only sites sampled with the CABIN protocol
rownames(ABMI.fine.2018) <- ABMI.fine.2018$SiteYear
str(ABMI.fine.2018)

# ##import ABMI 2018 UMOS data
# ABMI.umos.2018 <- read.csv("data/2018 umos aquatic invert data.csv") %>%
#     dplyr::select(ABMI.Site, Coarse.Group, Genus, Species, Specimen.Count) %>%
#     rename(Site = ABMI.Site) %>%
#     filter(Coarse.Group == "UMOS") %>%
#     mutate(Species.ID = paste0(Genus,".",Species)) %>%
#     dplyr::select(-c(Genus, Species, Coarse.Group)) %>%
#     pivot_wider(names_from = Species.ID, values_from = Specimen.Count, values_fn = sum) %>%
#     mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>% #replace blanks with zeroes
#     mutate(Year = 2022, .after = 1) %>%
#     mutate(SiteYear = paste0(Site, "_", Year), .after = 2) %>% #create identifying column
#     filter(Site %in% site.names$Site) #filter for only CABIN sites
# rownames(ABMI.umos.2018) <- ABMI.umos.2018$SiteYear
# str(ABMI.umos.2018)

##UMOS data
ABMI.umos.2018 <- species.2018 %>%
    dplyr::select(ABMI.Site, Coarse.Group, Genus, Species, Specimen.Count) %>%
    rename(Site = ABMI.Site) %>%
    filter(Coarse.Group == "UMOS") %>%
    mutate(Species.ID = paste0(Genus,".",Species)) %>%
    dplyr::select(-c(Genus, Species, Coarse.Group)) %>%
    filter(Species.ID != c(".UID")) %>% #remove unidentified specimens
    pivot_wider(names_from = Species.ID, values_from = Specimen.Count, values_fn = sum) %>%
    mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>% #replace blanks with zeroes
    mutate(Year = 2018, .after = Site) %>%
    mutate(SiteYear = paste0(Site, "_", Year), .after = Year) %>% #create identifying column
    mutate(Year = as.integer(Year)) %>%
    filter(Site %in% site.names$Site) #include only sites sampled with the CABIN protocol
rownames(ABMI.umos.2018) <- ABMI.umos.2018$SiteYear
str(ABMI.umos.2018)

#### import 2022 data ####
species.2022 <- read.csv("data/aquatic invert data_2022.csv")
MCN.2022 <- read.csv("data/Marchant_cell_counts.csv") %>%
    filter(Year == 2022) %>%
    rename(Site = Site.Name,
           MC_counted = COUNT...)

##Coarse data
coarse.2022 <- species.2022 %>%
    dplyr::select(ABMI.Site, Coarse.Group, Specimen.Count) %>%
    rename(Site = ABMI.Site) %>%
    filter(Coarse.Group != "UMOS") %>%
    pivot_wider(names_from = Coarse.Group, values_from = Specimen.Count, values_fn = sum) %>%
    mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>% #replace blanks with zeroes
    left_join(MCN.2022, by = "Site") %>%
    mutate(SiteYear = paste0(Site, "_", Year), .after = 2) %>% #create identifying column
    dplyr::select(Site, Year, SiteYear, MC_counted, everything())
rownames(coarse.2022) <- coarse.2022$SiteYear
str(coarse.2022)

##Fine data
fine.2022 <- species.2022 %>%
    dplyr::select(ABMI.Site, Genus, Species, Specimen.Count) %>%
    rename(Site = ABMI.Site) %>%
    mutate(Species.ID = paste0(Genus,".",Species)) %>%
    dplyr::select(-c(Genus, Species)) %>%
    pivot_wider(names_from = Species.ID, values_from = Specimen.Count, values_fn = sum) %>%
    mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>% #replace blanks with zeroes
    mutate(Year = 2022) %>%
    mutate(SiteYear = paste0(Site, "_", Year), .after = 2) #create identifying column
rownames(fine.2022) <- fine.2022$SiteYear
str(fine.2022)

##UMOS data
umos.2022 <- species.2022 %>%
    dplyr::select(ABMI.Site, Coarse.Group, Genus, Species, Specimen.Count) %>%
    rename(Site = ABMI.Site) %>%
    filter(Coarse.Group == "UMOS") %>%
    mutate(Species.ID = paste0(Genus,".",Species)) %>%
    dplyr::select(-c(Genus, Species, Coarse.Group)) %>%
    filter(Species.ID != c(".UID")) %>% #remove unidentified specimens
    pivot_wider(names_from = Species.ID, values_from = Specimen.Count, values_fn = sum) %>%
    mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>% #replace blanks with zeroes
    mutate(Year = 2022, .after = 1) %>%
    mutate(SiteYear = paste0(Site, "_", Year), .after = 2) #create identifying column
rownames(umos.2022) <- umos.2022$SiteYear
str(umos.2022)

#### import 2023 data ####
species.2023 <- read.csv("data/aquatic invert data_2023.csv")
str(species.2023)
MCN.2023 <- read.csv("data/Marchant_cell_counts.csv") %>%
    filter(Year == 2023) %>%
    rename(Site = Site.Name,
           MC_counted = COUNT...)

##Coarse data
coarse.2023 <- species.2023 %>%
    dplyr::select(ABMI.Site, Coarse.Group, Specimen.Count) %>%
    rename(Site = ABMI.Site) %>%
    filter(Coarse.Group != "UMOS") %>%
    pivot_wider(names_from = Coarse.Group, values_from = Specimen.Count, values_fn = sum) %>% #reshape into wide format
    mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>% #replace blanks with zeroes
    left_join(MCN.2023, by = "Site") %>%
    mutate(SiteYear = paste0(Site, "_", Year), .after = 2) %>% #create identifying column
    dplyr::select(Site, Year, SiteYear, MC_counted, everything())
rownames(coarse.2023) <- coarse.2023$SiteYear
str(coarse.2023)

##Fine data
fine.2023 <- species.2023 %>%
    dplyr::select(ABMI.Site, Genus, Species, Specimen.Count) %>%
    rename(Site = ABMI.Site) %>%
    mutate(Species.ID = paste0(Genus,".",Species)) %>%
    dplyr::select(-c(Genus, Species)) %>%
    pivot_wider(names_from = Species.ID, values_from = Specimen.Count, values_fn = sum) %>% #reshape into wide format
    mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>% #replace blanks with zeroes
    mutate(Year = 2023, .after = 1) %>%
    mutate(SiteYear = paste0(Site, "_", Year), .after = 2) #create identifying column
rownames(fine.2023) <- fine.2023$SiteYear
str(fine.2023)

##UMOS data
umos.2023 <- species.2023 %>%
    dplyr::select(ABMI.Site, Coarse.Group, Genus, Species, Specimen.Count) %>%
    rename(Site = ABMI.Site) %>%
    filter(Coarse.Group == "UMOS") %>%
    mutate(Species.ID = paste0(Genus,".",Species)) %>%
    dplyr::select(-c(Genus, Species, Coarse.Group)) %>%
    pivot_wider(names_from = Species.ID, values_from = Specimen.Count, values_fn = sum) %>% #reshape into wide format
    mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>% #replace blanks with zeroes
    mutate(Year = 2023, .after = 1) %>% #mark year
    mutate(SiteYear = paste0(Site, "_", Year), .after = 2) #create identifying column
rownames(umos.2023) <- umos.2023$SiteYear
str(umos.2023)

#### join dataframes ####
##coarse data
coarse.data.1 <- bind_rows(coarse.2018, ABMI.coarse.2018, coarse.2022, coarse.2023) %>%
    mutate(across(-c(Site, Year, SiteYear, MC_counted), ~replace_na(., 0))) %>%
    mutate(Site.ID = str_extract(Site, "\\d+")) # Extract the numeric part from all Site names (e.g., "W1598" or "CABIN-W1598" -> "1598")

# Identify Site.IDs that have both a CABIN and non-CABIN version
coarse.paired_ids <- coarse.data.1 %>%
    mutate(IsCABIN = str_starts(Site, "CABIN-")) %>%
    distinct(Site.ID, IsCABIN) %>%
    group_by(Site.ID) %>%
    filter(n() == 2) %>%
    pull(Site.ID)

# Filter the original data to keep only rows with these paired Site.IDs
coarse.data <- coarse.data.1 %>%
    filter(Site.ID %in% coarse.paired_ids) %>%
    dplyr::select(!Site.ID)
rownames(coarse.data) <- coarse.data$SiteYear
str(coarse.data)

##fine data
fine.data.1 <- bind_rows(fine.2018, ABMI.fine.2018, fine.2022, fine.2023) %>%
    mutate(across(-c(Site, Year, SiteYear), ~replace_na(., 0))) %>%
    mutate(Site.ID = str_extract(Site, "\\d+"), .after = SiteYear) # Extract the numeric part from all Site names (e.g., "W1598" or "CABIN-W1598" -> "1598")
str(fine.data.1)

# Identify Site.IDs that have both a CABIN and non-CABIN version
fine.paired_ids <- fine.data.1 %>%
    mutate(IsCABIN = str_starts(Site, "CABIN-")) %>%
    distinct(Site.ID, IsCABIN) %>%
    group_by(Site.ID) %>%
    filter(n() == 2) %>%
    pull(Site.ID)

# Filter the original data to keep only rows with these paired Site.IDs
fine.data2 <- fine.data.1 %>%
    filter(Site.ID %in% fine.paired_ids) %>%
    dplyr::select(!Site.ID)
rownames(fine.data2) <- fine.data2$SiteYear
str(fine.data2)

##UMOS
UMOS.data.1<- bind_rows(umos.2018, ABMI.umos.2018, umos.2022, umos.2023) %>%
    mutate(across(-c(Site, Year, SiteYear), ~replace_na(., 0))) %>%
    #dplyr::select(!.UID) %>%
    mutate(Site.ID = str_extract(Site, "\\d+"), .after = SiteYear) # Extract the numeric part from all Site names (e.g., "W1598" or "CABIN-W1598" -> "1598")
str(UMOS.data.1)

# Identify Site.IDs that have both a CABIN and non-CABIN version
UMOS.paired_ids <- UMOS.data.1 %>%
    mutate(IsCABIN = str_starts(Site, "CABIN-")) %>%
    distinct(Site.ID, IsCABIN) %>%
    group_by(Site.ID) %>%
    filter(n() == 2) %>%
    pull(Site.ID)

# Filter the original data to keep only rows with these paired Site.IDs
UMOS.data.2 <- UMOS.data.1 %>%
    filter(Site.ID %in% UMOS.paired_ids) %>%
    dplyr::select(!Site.ID)
rownames(UMOS.data.2) <- UMOS.data.2$SiteYear
str(UMOS.data.2)

#### check sites ####
# Helper function to process each dataset
process_dataset <- function(df, dataset_name) {
    df %>%
        mutate(base_site = sub("^CABIN-", "", Site),
               source = ifelse(grepl("^CABIN-", Site), "CABIN", "ABMI")) %>%
        select(base_site, Year, source) %>%
        distinct() %>%
        pivot_wider(names_from = source,
                    values_from = source,
                    values_fn = length,
                    values_fill = 0) %>%
        mutate(across(c(ABMI, CABIN), ~ .x > 0)) %>%
        rename_with(~ paste0("has_", dataset_name, "_", tolower(.x)), c("ABMI", "CABIN"))
}

# Process each dataset
coarse_summary <- process_dataset(coarse.data, "coarse")
fine_summary   <- process_dataset(fine.data, "fine")
UMOS_summary   <- process_dataset(UMOS.data, "UMOS")

# Merge all together
summary_df <- reduce(list(coarse_summary, fine_summary, UMOS_summary),
                     full_join, by = c("base_site", "Year"))

# Replace NAs with FALSE
summary_df[is.na(summary_df)] <- FALSE
str(summary_df)

# remove all rows where all values are false
summary_df_filtered <- summary_df %>%
    filter(!if_all(ends_with("_cabin"), ~ .x == FALSE))
str(summary_df_filtered) #105 sites
# export
write.csv(summary_df_filtered, "output/CABIN sites and data availability_updated June 17 2025.csv", row.names = FALSE)

#### clean up fine data to match the lookup table ####
fine.data <- fine.data2 %>%
    rename_with(~ gsub("/", ".", .x)) %>% #replace '/' with '.' to match names in the lookup table
    dplyr::select(!UID) %>%
    mutate(Helobdella.cf.stagnalis = rowSums(across(c("Helobdella.cf..stagnalis", "Helobdella.cf. stagnalis")), na.rm = TRUE)) %>% #collapse duplicate columns
    mutate(Limnephilus.cf.externus = rowSums(across(c("Limnephilus.cf..externus", "Limnephilus.cf. externus")), na.rm = TRUE)) %>%
    dplyr::select(!c("Helobdella.cf..stagnalis", "Helobdella.cf. stagnalis", "Limnephilus.cf..externus", "Limnephilus.cf. externus")) %>% #drop columns
    rename_with(~ gsub(" ", ".", .x)) %>% #replace spaces with '.' to match names in lookup table
    mutate(Cricotopus.sp.C.RPH = rowSums(across(c("Cricotopus.sp.C.RPH", "Cricotopus.sp..C.RPH")), na.rm = TRUE)) %>% #collapse duplicate columns
    mutate(Psectrocladius.sp.M.RPH = rowSums(across(c("Psectrocladius..sp..M.RPH", "Psectrocladius.sp..M.RPH")), na.rm = TRUE)) %>% #collapse duplicate columns
    dplyr::select(!c(Cricotopus.sp..C.RPH, Psectrocladius..sp..M.RPH, Psectrocladius.sp..M.RPH)) %>% #remove duplicate columns
    rename_with(~ gsub("\\.+", ".", .x)) %>% #replace '..' with '.' to match names in lookup table
    mutate(Bezzia.Palpomyia = rowSums(across(c("Bezzia.Palpomyia.sp.", ".Bezzia.Palpomyia.sp.", ".Bezzia.Palpomyia.Probezzia", ".Bezzia.Palpomyia", ".Bezzia.Palpomyia.Probezzia")), na.rm = TRUE),
           Bezzia.Probezzia.Palpomyia = rowSums(across(c("Bezzia.Probezzia.Palpomyia.sp.", ".Bezzia.Probezzia.Palpomyia.sp.", ".Bezzia.Probezzia.Palpomyia", "Bezzia.Palpomyia")), na.rm = TRUE),
           Ceratopogon.Culicoides.Stilobezzia = rowSums(across(c("Ceratopogon.Culicoides.Stilobezzia.sp.", ".Ceratopogon.Culicoides.Stilobezzia.sp.", ".Ceratopogon.Culicoides.Stilobezzia")), na.rm = TRUE),
           Helius.sp. = rowSums(across(c("Helius.sp.", ".Helius.sp.")), na.rm = TRUE),
           Lasioseius.sp. = rowSums(across(c("Lasioseius.sp.", ".Lasioseius.sp.", "Lasioseius.sp.")), na.rm = TRUE),
           Pericoma.Telmatoscopus = rowSums(across(c("Pericoma.Telmatoscopus.sp.", ".Pericoma.Telmatoscopus.sp.")), na.rm = TRUE),
    )  %>% #collapse duplicate columns
    dplyr::select(!c(Bezzia.Palpomyia.sp., .Bezzia.Palpomyia.sp., .Bezzia.Palpomyia, Bezzia.Probezzia.Palpomyia.sp., .Bezzia.Probezzia.Palpomyia.sp., .Bezzia.Palpomyia.Probezzia,
                     .Bezzia.Probezzia.Palpomyia, Ceratopogon.Culicoides.Stilobezzia.sp., .Ceratopogon.Culicoides.Stilobezzia.sp.,
                     .Ceratopogon.Culicoides.Stilobezzia, .Helius.sp., Lasioseius.sp., Lasioseius.sp., Pericoma.Telmatoscopus.sp., .Pericoma.Telmatoscopus.sp., Bezzia.Palpomyia,
                     )) %>%
    rename_with(~ gsub("^\\.", "", .x)) %>% #remove '.' at the start of column names
    dplyr::select(!c(SNI, SNI.Terrestrial, UID)) %>%
    mutate(Haliplus.cribrarius = rowSums(across(c("Haliplus.cribrarius", "Haliplus.cribarius")), na.rm = TRUE),
           Promenetus.exacuous = rowSums(across(c("Promenetus.exacuous", "Promenetus.exacuous.megas", "Promenetus.exacuous.exacuous")), na.rm = TRUE),
           Lymnaea.stagnalis = rowSums(across(c("Lymnaea.stagnalis", "Lymnaea.stagnalis.jugularis")), na.rm = TRUE),
           Endochironomus.nigricans = rowSums(across(c("Endochironomus.nigricans", "Endochrinomus.nigricans")), na.rm = TRUE),
           Psectrocladius.sp.L.RPH = rowSums(across(c("Psectrocladius.sp.L.RPH", "Psectrocaldius.sp.L.RPH")), na.rm = TRUE),
           Thienemannimyia.group = rowSums(across(c("Thienemannimyia.group", "Thienemannimyia.Group")), na.rm = TRUE),
           Cenocorixa.bifida = rowSums(across(c("Cenocorixa.bifida", "Cenocorixa.bifida.bifida")), na.rm = TRUE),
           Cordulia.shurtleffii = rowSums(across(c("Cordulia.shurtleffii", "Cordulia.shurtlefii")), na.rm = TRUE),
           Porolohmannella.violacea = rowSums(across(c("Porolohmannella.violacea", "Porolohmanella.violacea")), na.rm = TRUE),
           Hydroporinae.UID = rowSums(across(c("Hydroporinae.UID", "Hydroporinae")), na.rm = TRUE)
    ) %>%
    dplyr::select(!c(Haliplus.cribarius, Promenetus.exacuous.megas, Lymnaea.stagnalis.jugularis, Endochrinomus.nigricans,
                     Psectrocaldius.sp.L.RPH, Thienemannimyia.Group, Promenetus.exacuous.exacuous,
                     Cenocorixa.bifida.bifida, Cordulia.shurtlefii, Porolohmanella.violacea, Hydroporinae)) %>%
    rename(Ferrissia.californica = Ferrissia.fragilis, #fix misspellings/subspecies names
           Caenis.cf.diminuta = Caenis.cf.diminuta,
           Hygrotus.sp.other = Hygrotus.sp.,
           Stagnicola.catascopium = Stagnicola.catascopium.catascopium,
           Mystacides.sp.other = Mystacides.sp.,
           Donacia.sp. = Donaciinae,
           "Psectrotanypus.(Derotanypus).sp." = "Psectrotanypus.(Derotanypus)",
           Coenagrionidae.UID = Coenagrion.Enallagma,
           Polypedilum.scaelum.group = Polypedilum.scalaenum.group,
           Eukiefferiella.sp. = Eukiefferiella,
           Hydrobaenus.sp. = Hydrobaenus
           ) %>%
    mutate(Caenis.cf.diminuta = rowSums(across(c("Caenis.cf.diminuta", "Caenis.diminuta")), na.rm = TRUE),
           Caenis.latipennis = rowSums(across(c("Caenis.latipennis", "Caenis.cf.latipennis")), na.rm = TRUE),
           Haliplus.apicalis = rowSums(across(c("Haliplus.apicalis", "Haliplus.strigatus")), na.rm = TRUE),
           Gyrinus.sp. = rowSums(across(c("Gyrinus.sp.", "Gyrinus.wallisi.aeratus")), na.rm = TRUE),
           Ceratopogon.Culicoides.Stilobezzia = rowSums(across(c("Ceratopogon.Culicoides.Stilobezzia", "Culicoides")), na.rm = TRUE)
           ) %>%
    dplyr::select(!c(Caenis.diminuta, Caenis.cf.latipennis, Haliplus.strigatus, Gyrinus.wallisi.aeratus, Culicoides, Genus.9.RPH)) #note: deleting Genus.9.RPH b/c I have no idea what it is.
str(fine.data)

#### clean up UMOS data to match the lookup table ####
UMOS.data <- UMOS.data.2 %>%
    rename_with(~ gsub("/", ".", .x)) %>% #replace '/' with '.' to match names in the lookup table
    mutate(Lymnaea.stagnalis.jugularis = rowSums(across(c("Lymnaea.stagnalis.jugularis", "Lymnaea.stagnalis jugularis")), na.rm = TRUE)) %>%
    dplyr::select(!"Lymnaea.stagnalis jugularis") %>%
    rename_with(~ gsub(" ", ".", .x)) %>% #replace spaces with '.' to match names in lookup table
    rename_with(~ gsub("\\.+", ".", .x)) %>% #replace '..' with '.' to match names in lookup table
    mutate(Ephemeroptera.UID = rowSums(across(c("Ephemeroptera.UID", ".Ephemeroptera.UID")), na.rm = TRUE)) %>%
    dplyr::select(!.Ephemeroptera.UID) %>%
    rename_with(~ gsub("^\\.", "", .x)) %>% #remove '.' at the start of column names
    mutate(Cenocorixa.bifida = rowSums(across(c("Cenocorixa.bifida", "Cenocorixa.bifida.bifida")), na.rm = TRUE),
           Cordulia.shurtleffii = rowSums(across(c("Cordulia.shurtleffii", "Cordulia.shurtlefii")), na.rm = TRUE),
           Haliplus.cribrarius = rowSums(across(c("Haliplus.cribrarius", "Haliplus.cribarius")), na.rm = TRUE),
           Hesperocorixa.michiganensis = rowSums(across(c("Hesperocorixa.michiganensis", "Hesperocorixa.michigaenis")), na.rm = TRUE),
           Promenetus.exacuous = rowSums(across(c("Promenetus.exacuous", "Promenetus.exacuous.megas")), na.rm = TRUE),
           Lymnaea.stagnalis = rowSums(across(c("Lymnaea.stagnalis", "Lymnaea.stagnalis.jugularis")), na.rm = TRUE),
           Gyrinus.sp. = rowSums(across(c("Gyraulus.wallisi.aeratus", "Gyrinus.wallisi.aeratus")), na.rm = TRUE),
           Agrypnia.sp. = rowSums(across(c("Agrypnia.sp.", "Banksiola.sp.")), na.rm = TRUE),
           Limnoporus.dissortis = rowSums(across(c("Limnoporus.dissortis", "Limnoporus.sp.")), na.rm = TRUE),
           Coenagrionidae.UID = rowSums(across(c("Coenagrionidae.UID", "Coenagrion.Enallagma")), na.rm = TRUE),
           ) %>%
    dplyr::select(!c(Cenocorixa.bifida.bifida, Cordulia.shurtlefii, Haliplus.cribarius, Hesperocorixa.michigaenis, Promenetus.exacuous.megas,
                     Lymnaea.stagnalis.jugularis, Gyraulus.wallisi.aeratus, Gyrinus.wallisi.aeratus, Banksiola.sp., Limnoporus.sp., Coenagrion.Enallagma)) %>%
    rename(Mystacides.sp.other = Mystacides.sp.)
rownames(UMOS.data) <- UMOS.data$SiteYear

#### export ####
write_rds(coarse.data, "output/cleaned CABIN coarse data.rds")
write_rds(fine.data, "output/cleaned CABIN fine data.rds")
write_rds(UMOS.data, "output/cleaned CABIN UMOS data.rds")
