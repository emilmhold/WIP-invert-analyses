## Data setup and analysis final
## Author: Ermias T. Azeria
## Last edited by Emily: July 8, 2025
#Load raw data
#Look up table
#Site information

#install.packages("openxlsx")
#install.packages("abmi.themes")
#install.packages("tidyverse")
#install.packages("backports")

library(openxlsx)
library(abmi.themes)
library(tidyverse)
library(backports)

setwd("~/GitHub/New aquatic invertebrate analyses")

natcol <- scale_fill_abmi(palette = "main")$palette(6)

### this script requires functions from "Taxonomic ambiguous assignment abmi_clean.R"
source("scripts/Taxonomic ambiguous assignment abmi_clean.R")
    ##Em's note, there will be an error message:"Error in eval(ei, envir) :
        #object 'fix.sites' not found"
    ## that's ok - the fix.sites function was a fragment in Ermias' script that
        #I couldn't figure out where it came from.

#### Lookup table ####
#First step: USE "Analysis_Grouping" to aggregate "Analysis_Name";
#but keep taxon names with "Distribute to Children" in the Analysis_Grouping as is in raw data
#handle this in LOOK up table, add Taxon_Analysis

Invert_lookup <-  read.csv("Data/Look-up Table_Invert Taxonomy_2025-07-28.csv")
#Add a column with analysis name for ambiguous parents designated as "distribute to children" in analysis grouping
Invert_lookup$Taxon_Analysis <- Invert_lookup$Analysis_Grouping
parent.id <- which(Invert_lookup$Taxon_Analysis=="Distribute.to.Children")
Invert_lookup$Taxon_Analysis[parent.id ] <- Invert_lookup$Analysis_Name[parent.id]

#### load raw data ####
#Get fine species data, and use lookup table to first aggregate them to analysis-grouping
#use first step taxon grouping to aggregate data
SppData.raw <- read_rds("output/cleaned CABIN fine data.rds")
str(SppData.raw)

#read in UMOS data
umosData.raw <- read_rds("output/cleaned CABIN UMOS data.rds")
str(umosData.raw)

#read in coarse data (including MCN data)
SppData.coarse <- read_rds("output/cleaned CABIN coarse data.rds") %>%
    rename(Coleoptera.adults = "Coleoptera (adults)",
           Coleoptera.larvae = "Coleoptera (larvae)",
            Other.Diptera = "Other Diptera",
           Other.Laevicaudata = "Other: Laevicaudata")
str(SppData.coarse)

#site strata (aka climate info)
WetlandInfo <- get(load("~/GitHub/invertebrate-protocol-analyses_2021/Taxonomic assignment processing/Data/wetland-climate_2023 updated.RData"))

#### Add climate info for CABIN wetlands ####
cabin.df <- SppData.raw %>%
  filter(startsWith(SiteYear, "CABIN")) %>% ## select cabin sites
  rename(cabin_site = Site) %>% #rename to facilitate joining with WetlandInfo df
  mutate(Site = str_remove(cabin_site, "^CABIN-"), .before = "Year") %>% # removes CABIN at the start of the site columns
  dplyr::select(cabin_site, Site, Year)
str(cabin.df)

##subset duplicate observations
cabin.df.duplicates <- SppData.raw %>%
    filter(endsWith(Site, "-D")) %>% #identify duplicate sites
    mutate(base_site = Site %>%
               sub("-D$", "", .) %>%         # Remove -D at the end to facilitate joining with climate df
               sub("^CABIN-", "", .)) %>%    # Remove CABIN- at the start
dplyr::select(Site, base_site, Year) %>%
    rename(cabin_site = Site,
        Site = base_site) # rename to facilitate joining
str(cabin.df.duplicates)

## find rows in WetlandInfo where Site matches the Site column in cabin.df
cabin.WetlandInfo <- WetlandInfo %>%
  inner_join(cabin.df, by = c("Site","Year")) %>% #select rows where Site and Year aligns with cabin.site
  dplyr::select(!c("Site", "SiteYear")) %>%
  rename(Site = cabin_site) %>% #rename site to facilitate rbinding.
  relocate(Site, .before = Year)  %>% # Move Site before Year
  mutate(SiteYear = paste0(Site, "_", Year), .before = Site)
str(cabin.WetlandInfo) #only 93 sites are present here; missing climate data for 12 sites

# check that columns in cabin.WetlandInfo and WetlandInfo are in the same order
identical(names(cabin.WetlandInfo), names(WetlandInfo))

#find rows in WetlandInfo where Site matches the Site column in duplicates df
cabin.duplicate.WetlandInfo <- cabin.df.duplicates %>%
    inner_join(WetlandInfo, by = c("Site", "Year"), multiple = "all") %>% #select rows where Site aligns with cabin.site
    dplyr::select(!c("Site", "SiteYear")) %>%
    rename(Site = cabin_site) %>% #rename site to facilitate rbinding.
    relocate(Site, .before = Year)  %>% # Move Site before Year
    mutate(SiteYear = paste0(Site, "_", Year), .before = Site)
str(cabin.duplicate.WetlandInfo)
rownames(cabin.duplicate.WetlandInfo) <- cabin.duplicate.WetlandInfo$SiteYear

# check that columns are in the same order
identical(names(cabin.duplicate.WetlandInfo), names(WetlandInfo))

#append new dfs to WetlandInfo
new.WetlandInfo <- bind_rows(WetlandInfo, cabin.WetlandInfo) %>%
  bind_rows(cabin.duplicate.WetlandInfo) %>%
  dplyr::select(!SiteYear) %>%
  mutate(SiteYear = paste0(Site,"_",Year), .before = Site)  # make sure cabin data are reflected in SiteYear
str(new.WetlandInfo) ## looks good but there are duplicates
## pull out duplicates
dupes <- new.WetlandInfo[new.WetlandInfo$SiteYear %in% new.WetlandInfo$SiteYear[duplicated(new.WetlandInfo$SiteYear)], ]
    #duplicates are true duplicates. Can drop either the first or second occurrence
new.WetlandInfo <- new.WetlandInfo[!duplicated(new.WetlandInfo$SiteYear), ] #drop duplicates
rownames(new.WetlandInfo) <- new.WetlandInfo$SiteYear
##export df
write_rds(new.WetlandInfo, "output/cleaned wetland and climate info.rds")

#### prep fine species data ####
SppData <- SppData.raw %>%
    column_to_rownames("SiteYear") %>%
    dplyr::select(!c(Site,Year))
str(SppData)
#SppData <- adj.count
#SppData$Site <-  SppData$Year <- SppData$SiteYear <- NULL
#SppData$Site <- SppData$Protocol <- SppData$Marchant <- NULL
mnames <- match(colnames(SppData), Invert_lookup$Analysis_Name) #check
missing.names <- colnames(SppData)[is.na(mnames)] #check for missing names
missing.names
#find what year these species were observed in
missing.spp <- SppData.raw %>%
    dplyr::filter(if_any(all_of(missing.names), ~ .x > 0))
#export missing names
#write.csv(missing.spp, "output/fine data species names and observation years not in lookup table (June 19 2025).csv", row.names = FALSE)

#Match and aggregate to taxon group (analysis grouping)
mnames <- match(colnames(SppData), Invert_lookup$Analysis_Name)
taxon_gp <- Invert_lookup$Taxon_Analysis[mnames]
SppData.in  <- as.data.frame(as.matrix(mefa4::groupSums(as.matrix(SppData) ,2, by= as.character(taxon_gp), na.rm=TRUE)))
#site info
misssites <- setdiff(rownames(SppData.in), new.WetlandInfo$SiteYear)
print(misssites) #sites in new.WetlandInfo missing in SppData.in
##check for NAs
rownames(SppData.in)[apply(SppData.in, 1, anyNA)]


#### Auxiliary (UMOS) data ####
    ##Note: this is auxiliary data with 'resolved' children data to be used for resolving taxonomy. For ABMI this refers to UMOS (unique mature organisms search)
rownames(umosData.raw) <- NULL
umosData <- umosData.raw %>%
    column_to_rownames(var = "SiteYear") %>%
    dplyr::select(!c("Site", "Year"))
mnames <- match(colnames(umosData), Invert_lookup$Analysis_Name)
missing.UMOS.names <- colnames(umosData)[is.na(mnames)]
missing.UMOS.names #none
#identify missing names in data
#missing.UMOS.spp <- umosData %>%
 #   dplyr::filter(if_any(all_of(missing.UMOS.names), ~ .x > 0)) %>%
  #  dplyr::select(Site, Year, SiteYear, all_of(missing.UMOS.names))
##export these data
#write.csv(missing.UMOS.spp, "output/UMOS species names and observation years not in lookup table (June 20 2025).csv")

taxon_gp <- Invert_lookup$Taxon_Analysis[mnames]
umosData.in  <- as.data.frame(as.matrix(mefa4::groupSums(as.matrix(umosData) ,2, by= as.character(taxon_gp), na.rm=TRUE)))
setdiff(rownames(umosData.in ), new.WetlandInfo$SiteYear) ## Em's note: change new.WetlandInfo if you don't want to include CABIN data
umosData.in$Do.not.analyze<- NULL

#RUN ANALYSIS

Taxon_Lookup <- Invert_lookup

#Match main data and auxiliary (UMOS) data
mmdata <- fn_matchdata (data=SppData.in, aux.data=umosData.in)
SppData.in <- mmdata$data
umosData.in<- mmdata$aux.data



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### FIRST ADJUST FOR MARCHANT CELLS (balanced) ####
#Alternative is to run the analysis first then adjust for proportion of MCN after processing data

SppData.bal <- SppData.in
setdiff(rownames( SppData.bal) , rownames(SppData.coarse ))
setdiff(rownames(SppData.coarse ),rownames( SppData.bal)  )

commsites.b <- intersect(rownames(SppData.bal), rownames(SppData.coarse))
SppData.bal <- SppData.bal[ commsites.b ,  ]
SppData.coarse.bal <- SppData.coarse[ commsites.b ,  ]
SppData.bal  <- (SppData.bal /SppData.coarse.bal$MC_counted) *100 #adjust to 100 cells
SppData.bal  <- round(SppData.bal)

umosData.bal <- umosData.in[ commsites.b ,  ]

setdiff(commsites.b, rownames(new.WetlandInfo))
NR.site.b <- as.character(new.WetlandInfo[  commsites.b  , "NR"])
NSR.site.b <- as.character(new.WetlandInfo[  commsites.b , "NSR"])

##check if there are any NAs in dataframes (this will prevent the algorithm from running)
rownames(SppData.bal)[apply(SppData.bal, 1, anyNA)]
rownames(SppData.coarse.bal)[apply(SppData.coarse.bal, 1, anyNA)]
rownames(umosData.bal)[apply(umosData.bal, 1, anyNA)]
rownames(Taxon_Lookup)[apply(Taxon_Lookup, 1, anyNA)]

#### run algorithm ####
#WITH STRATA
#most frequent/most abundant
# APTC
#### Em's note: Jenet wants to go with APTC-f approach so I'm only running the algorithm for this approach! ####

resAPTC.frq.bal <- aptc_abmi (data=SppData.bal,aux.data=umosData.bal, taxalookup=Taxon_Lookup,
                                     strata.primary=NSR.site.b,strata.secondary=NR.site.b,  method.region="frequent", reg.prop=F)
resAPTC.abu.bal  <- aptc_abmi (data=SppData.bal, aux.data=umosData.bal,taxalookup=Taxon_Lookup,
                                      strata.primary=NSR.site.b,strata.secondary=NR.site.b,  method.region="abundant", reg.prop=F)



write.xlsx(resAPTC.frq.bal , rowNames=TRUE, file = "output/taxonomically adjusted data/CABIN Ambiguous parent resolved all taxa_APTC.frequent.xlsx")
write.xlsx(resAPTC.abu.bal, rowNames=TRUE, file = "output/taxonomically adjusted data/CABIN Ambiguous parent resolved all taxa_APTC.abundant.xlsx")


# #DPAC
# resDPAC.frq.bal <- dpac_abmi (data=SppData.bal,aux.data=umosData.bal, taxalookup=Taxon_Lookup,
#                                      strata.primary=NSR.site.b,strata.secondary=NR.site.b,  method.region="frequent", reg.prop=F)
# resDPAC.abu.bal  <- dpac_abmi (data=SppData.bal, aux.data=umosData.bal,taxalookup=Taxon_Lookup,
#                                       strata.primary=NSR.site.b,strata.secondary=NR.site.b,  method.region="abundant", reg.prop=F)
#
#
#
# write.xlsx(resDPAC.frq.bal , rowNames=TRUE, file = "~/GitHub/invertebrate-protocol-analyses_2021/Taxonomic assignment processing/Results/MCN adjusted abundance/CABIN Ambiguous parent resolved all taxa_DPAC.frequent.xlsx")
# write.xlsx(resDPAC.abu.bal, rowNames=TRUE, file = "~/GitHub/invertebrate-protocol-analyses_2021/Taxonomic assignment processing/Results/MCN adjusted abundance/CABIN Ambiguous parent resolved all taxa_DPAC.abundant.xlsx")
#
#
#
#PROPORTIONAL DATA
#APTC
resAPTC.frqp.bal  <- aptc_abmi (data=SppData.bal,aux.data=umosData.bal, taxalookup=Taxon_Lookup,
                                       strata.primary=NSR.site.b,strata.secondary=NR.site.b, method.region="frequent", reg.prop=T)
resAPTC.abup.bal  <- aptc_abmi (data=SppData.bal, aux.data=umosData.bal,taxalookup=Taxon_Lookup,
                                       strata.primary=NSR.site.b,strata.secondary=NR.site.b,  method.region="abundant", reg.prop=T)


write.xlsx(resAPTC.frqp.bal , rowNames=TRUE, file = "output/taxonomically adjusted data/Ambiguous parent resolved all taxa_APTC.Propfrequency.xlsx")
write.xlsx(resAPTC.abup.bal, rowNames=TRUE, file = "output/taxonomically adjusted data/Ambiguous parent resolved all taxa_APTC.Propabundance.xlsx")

#
#
# #DPAC
# resDPAC.frqp.bal  <- dpac_abmi (data=SppData.bal,aux.data=umosData.bal, taxalookup=Taxon_Lookup,
#                                        strata.primary=NSR.site.b,strata.secondary=NR.site.b, method.region="frequent", reg.prop=T)
# resDPAC.abup.bal  <- dpac_abmi (data=SppData.bal, aux.data=umosData.bal,taxalookup=Taxon_Lookup,
#                                        strata.primary=NSR.site.b,strata.secondary=NR.site.b,  method.region="abundant", reg.prop=T)
#
#
# write.xlsx(resDPAC.frqp.bal , rowNames=TRUE, file = "~/GitHub/invertebrate-protocol-analyses_2021/Taxonomic assignment processing/Results/MCN adjusted abundance/Ambiguous parent resolved all taxa_DPAC.Propfrequency.xlsx")
# write.xlsx(resDPAC.abup.bal, rowNames=TRUE, file = "~/GitHub/invertebrate-protocol-analyses_2021/Taxonomic assignment processing/Results/MCN adjusted abundance/Ambiguous parent resolved all taxa_DPAC.Propabundance.xlsx")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Biodiversity metrics
#Most frequent or most abundant in region
APTC.frq_bdiv.adj <- fn.commdiv (resAPTC.frq.bal)
APTC.abu_bdiv.adj <- fn.commdiv (resAPTC.abu.bal)
#DPAC.frq_bdiv.adj <- fn.commdiv (resDPAC.frq.bal)
#DPAC.abu_bdiv.adj <- fn.commdiv (resDPAC.abu.bal)

#proportional to regional abundance /frequency

APTC.Propfrq_bdiv.adj <- fn.commdiv (resAPTC.frqp.bal)
APTC.Propabu_bdiv.adj <- fn.commdiv (resAPTC.abup.bal)
#DPAC.Propfrq_bdiv.adj <- fn.commdiv (resDPAC.frqp.bal)
#DPAC.Propabu_bdiv.adj <- fn.commdiv (resDPAC.abup.bal)



#Export diversity data
BDIV.out_MCN.adj <- list(APTC.frq=APTC.frq_bdiv.adj , APTC.abu=APTC.abu_bdiv.adj ,DPAC.frq=DPAC.frq_bdiv.adj , DPAC.abu=DPAC.abu_bdiv.adj  ,
                         APTC.Propfrq=APTC.Propfrq_bdiv.adj , APTC.Propabu=APTC.Propabu_bdiv.adj ,DPAC.Propfrq=DPAC.Propfrq_bdiv.adj , DPAC.Propabu=DPAC.Propabu_bdiv.adj )
write.xlsx(BDIV.out_MCN.adj , rowNames=TRUE, file = "~/GitHub/invertebrate-protocol-analyses_2021/Taxonomic assignment processing/Results/Biodiversity metrics for resolved taxonomy all taxa_Marchant count adjusted_with CABIN.xlsx")

save(BDIV.out_MCN.adj ,  file = "~/GitHub/invertebrate-protocol-analyses_2021/Taxonomic assignment processing/Results/Biodiversity metrics for resolved taxonomy all taxa_Marchant count adjusted_with CABIN.RData")

#Most frequent or abundant
#~~~
#### Figures ####

for (j in 1:length( BDIV.out_MCN.adj )){
  datIN <- BDIV.out_MCN.adj[[j]]
  dataname <- names( BDIV.out_MCN.adj)[j]
  #  dataname <- deparse(substitute(DPAC.frq_bdiv.adj))
  # dataname <- gsub ("_bdiv", "",dataname, fixed=T )
  ###################


  pdf.out <- paste0("~/GitHub/invertebrate-protocol-analyses_2021/Taxonomic assignment processing/Results/Summary of wetland invertebrate biodiversity metrics by natural region resolved by ",dataname, "_MCN adjusted.pdf")
  pdf(pdf.out , width=6, height=5)
  for (i in 1:ncol(datIN)){
    varin<- 1
    agg <- aggregate(datIN[ , i] ,list(NR.site.b), mean)
    rownames(agg)<-agg$Group.1
    agg$Group.1 <- NULL
    agg<- t(agg)
    agg <- agg[ , c ("Rocky Mountain","Foothills" , "Boreal" ,"Canadian Shield" , "Parkland"  , "Grassland")]
    op<- par(oma=c(0,0,0,0),mar=c(6,4,2,2))
    x <- barplot(agg,srt=45, col=natcol[1:6] , xaxt="n", ylab="Average", main=colnames(datIN)[i])
    labs <- c ("Rocky Mountain","Foothills" , "Boreal" ,"Canadian Shield" , "Parkland"  , "Grassland")

    text(x, par("usr")[3]-(0.03*max(agg)),   srt = 45, adj = 1, xpd = TRUE,
         labels =  labs, cex = 1, col=natcol , font=2)
    box(bty="l",col="#2d2a2a")
  }
  dev.off()
  par(op)

}








