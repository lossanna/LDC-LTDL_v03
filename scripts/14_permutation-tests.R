# Created: 2026-05-04
# Updated: 2026-05-04

# Purpose: Run permutation tests (based on Ron's script).


library(tidyverse)

# Load data ---------------------------------------------------------------

load("RData/13.1_matched-data.RData")
geoindicators.raw <- read_csv("data/raw/downloaded/ldc-data-2026-03-11/geoindicators.csv")


# Data wrangling ----------------------------------------------------------

# Adjust column names of geoindicators
col_rename_map <- c(
  "Project Key" = "ProjectKey",
  "Primary Key" = "PrimaryKey",
  "Date Visited" = "DateVisited",
  "Ecological Site ID" = "EcoSiteID",
  "Latitude (decimal degrees, NAD83)" = "Latitude",
  "Longitude (decimal degrees, NAD83)" = "Longitude",
  "Location Status" = "LocationStatus",
  "Location Type" = "LocationType",
  "Latitude, Actual (decimal degrees, NAD83)" = "LatActual",
  "Longitude, Actual (decimal degrees, NAD83)" = "LonActual",
  "Bare Soil (% First Hit)" = "BareSoil_FH",
  "Total Foliar Cover (%)" = "TotalFoliarCover",
  "Annual Forb Cover (% Any Hit)" = "AnnForbCover_AH",
  "Annual Graminoid Cover (% Any Hit)" = "AnnGramCover_AH",
  "Forb Cover (% Any Hit)" = "ForbCover_AH",
  "Annual Forb and Graminoid Cover (% Any Hit)" = "AnnForbGramCover_AH",
  "Graminoid Cover (% Any Hit)" = "GramCover_AH",
  "Perennial Forb Cover (% Any Hit)" = "PerForbCover_AH",
  "Perennial Forb and Graminoid Cover (% Any Hit)" = "PerForbGramCover_AH",
  "Perennial Graminoid Cover (% Any Hit)" = "PerGramCover_AH",
  "Shrub Cover (% Any Hit)" = "ShrubCover_AH",
  "FH Cyanobacteria Cover (% First Hit)" = "CyanobacteriaCover_FH",
  "Deposited Soil Cover (% First Hit)" = "DepositedSoilCover_FH",
  "Duff Cover (% First Hit)" = "DuffCover_FH",
  "Embedded Litter Cover (% First Hit)" = "EmbeddedLitterCover_FH",
  "Herbaceous Litter Cover (% First Hit)" = "HerbLitterCover_FH",
  "Lichen Cover (% First Hit)" = "LichenCover_FH",
  "Moss Cover (% First Hit)" = "MossCover_FH",
  "Rock Cover (% First Hit)" = "RockCover_FH",
  "Total Litter Cover (% First Hit)" = "TotalLitterCover_FH",
  "Vagrant Lichen Cover (% First Hit)" = "VagrantLichenCover_FH",
  "Water Cover (% First Hit)" = "WaterCover_FH",
  "Woody Litter Cover (% First Hit)" = "WoodyLitterCover_FH",
  "Canopy Gaps 25 - 50 cm (%)" = "Gap25_50",
  "Canopy Gaps 51-100 cm (%)" = "Gap51_100",
  "Canopy Gaps 101 - 200 cm (%)" = "Gap101_200",
  "Canopy Gaps > 200 cm (%)" = "Gap200plus",
  "Canopy Gaps > 25 cm (%)" = "Gap25plus",
  "Mean Forb Height (cm)" = "MeanForbHgt",
  "Mean Graminoid Height (cm)" = "MeanGramHgt",
  "Mean Herbaceous Plant Height (cm)" = "MeanHerbHgt",
  "Mean Perennial Forb Height (cm)" = "MeanPerForbHgt",
  "Mean Perennial Forb Graminoid Height (cm)" = "MeanPFbGrHgt",
  "Mean Perennial Graminoid Height (cm)" = "MeanPerGramHgt",
  "Mean Woody Plant Height (cm)" = "MeanWoodyHgt",
  "Total Annual Production (Rangeland Health)" = "TotAnnualProduction_RH",
  "Bare Ground (Rangeland Health)" = "BareGround_RH",
  "Biotic Integrity (Rangeland Health)" = "BioticIntegrity_RH",
  "Comments: Biotic Integrity (Rangeland Health)" = "BioticIntegrity_comments",
  "Comments: Hydrologic Function (Rangeland Health)" = "HydrologicFunction_comments",
  "Comments: Soil and Site Stability (Rangeland Health)" = "SoilSiteStability_comments",
  "Compaction (Rangeland Health)" = "Compaction_RH",
  "Proportion of Dead or Dying Plant Parts (Rangeland Health)" = "PropDeadDyingPlants_RH",
  "Functional/Sructural Groups (Rangeland Health)" = "FunctionalStructuralGroups_RH",
  "Gullies (Rangeland Health)" = "Gullies_RH",
  "Hydrologic Function (Rangeland Health)" = "HydrologicFunction_RH",
  "Invasive Plants (Rangeland Health)" = "InvasivePlants_RH",
  "Litter Amount (Rangeland Health)" = "LitterAmount_RH",
  "Litter Movement (Rangeland Health)" = "LitterMovement_RH",
  "Pedestals/Terracettes (Rangeland Health)" = "Pedestals_RH",
  "Plant Community Composition (Rangeland Health)" = "PlantCommunityComposition_RH",
  "Perennial Reproductive Capability (Rangeland Health)" = "PerReproCapactiy_RH",
  "Rills (Rangeland Health)" = "Rills_RH",
  "Soil Site Stability (Rangeland Health)" = "SoilSiteStability_RH",
  "Soil Surface Loss/Degradation (Rangeland Health)" = "SoilSurfaceLoss_RH",
  "Soil Surface Erosion Resistance (Rangeland Health)" = "SoilErosionResistance_RH",
  "Water Flow Patterns (Rangeland Health)" = "WaterFlowPatterns_RH",
  "Wind Scoured Areas (Rangeland Health)" = "WindScouredAreas_RH",
  "Mean Soil Stability: Surface" = "MeanSoilStability_Surface",
  "Mean Soil Stability: Protected Samples" = "MeanSoilStability_Protected",
  "Mean Soil Stability: Unprotected Samples" = "MeanSoilStabilityUnprotected",
  "MLRA Description" = "MLRADesc_LDC",
  "MLRA Symbol" = "MLRASym_LDC",
  "Ecoregion Level I" = "EcoLvl1_LDC",
  "Ecoregion Level II" = "EcoLvl2_LDC",
  "Ecoregion Level III" = "EcoLvl3_LDC",
  "Ecoregion Level IV" = "EcoLvl4_LDC",
  "State" = "State",
  "MODIS IGBP Name" = "MODISName",
  "Database Key" = "DBKey",
  "Date Loaded in Database" = "DateLoad",
  "Total Horizontal Flux" = "TotalHorizontalFlux",
  "Total Vertical Flux" = "TotalVerticalFlux",
  "PM 2.5 Vertical Flux" = "PM25Flux",
  "PM 10 Vertical Flux" = "PM10Flux",
  "Long-Term Mean Precipitation" = "LongTermMeanPrecip",
  "Long-Term Mean Runoff" = "LongTermMeanRunoff",
  "Long-Term Mean Sediment Yield" = "LongTermMeanSedimentYield",
  "Long-Term Mean Soil Loss" = "LongTermMeanSoilLoss"
)

geoindicators <- geoindicators.raw |>
  rename(!!!setNames(names(col_rename_map), col_rename_map))

# Identify completely empty columns
empty_cols <- geoindicators |>
  summarise(across(everything(), ~ all(is.na(.)))) |>
  pivot_longer(everything(), names_to = "column", values_to = "is_empty") |>
  filter(is_empty) |>
  pull(column)

geoindicators <- geoindicators |>
  select(-all_of(empty_cols))


# geoindcator cols for joining (functional group cover)
geoindicators.join <- geoindicators |> 
  select(PrimaryKey, AnnForbCover_AH, AnnGramCover_AH, PerForbCover_AH, PerGramCover_AH,
         ShrubCover_AH)



# Arizona/New Mexico Mountains --------------------------------------------

## 1. Prescribed Burn -----------------------------------------------------

# Join cover cols
model01.matched <- anm.prescribed.burn.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join)

#   pivot_longer() for cover cols
model01.matched <- model01.matched |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "pct_cover"
    )

#   Rename functional group cover types
model01.matched <- model01.matched |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub"
           ))

# Boxplot
model01.bp <- model01.matched |> 
  ggplot(aes(x = indicators, y = pct_cover, fill = trt_control)) +
  geom_boxplot() +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "1. AZ/NM Mountains: Prescribed Burn") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model01.bp

# Calculate observed mean difference
model01.diff <- model01.matched |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(pct_cover),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model01.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model01.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model01.matched |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data %>%
      group_by(indicators, trt_control) %>%
      summarize(mean_cover = mean(pct_cover), .groups = "drop") %>%
      group_by(indicators) %>%
      summarize(mean_diff = diff(mean_cover), .groups = "drop") %>%
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each functional group
p_values <- model01.perm %>%
  inner_join(x = ., y = model01.diff, by = "indicators") %>%
  group_by(indicators) %>%
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values


# Plot frequency distribution
#   Annual forb
model01.annforb <- model01.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed") +
  labs(title = "1. AZ/NM Mts, Prescribed Burn: Annual forb",
       x = "Difference in means",
       y = "Frequency") +
  theme_bw()
model01.annforb

#   Annual grass
model01.anngrass <- model01.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed") +
  labs(title = "1. AZ/NM Mts, Prescribed Burn: Annual grass",
       x = "Difference in means",
       y = "Frequency") +
  theme_bw()
model01.anngrass

#   Perennial forb
model01.perforb <- model01.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed") +
  labs(title = "1. AZ/NM Mts, Prescribed Burn: Perennial forb",
       x = "Difference in means",
       y = "Frequency") +
  theme_bw()
model01.perforb

#   Perennial grass
model01.pergrass <- model01.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed") +
  labs(title = "1. AZ/NM Mts, Prescribed Burn: Perennial grass",
       x = "Difference in means",
       y = "Frequency") +
  theme_bw()
model01.pergrass

#   Shrub
model01.shrub <- model01.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed") +
  labs(title = "1. AZ/NM Mts, Prescribed Burn: Shrub",
       x = "Difference in means",
       y = "Frequency") +
  theme_bw()
model01.shrub
