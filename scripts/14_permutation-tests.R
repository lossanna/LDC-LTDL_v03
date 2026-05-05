# Created: 2026-05-04
# Updated: 2026-05-05

# Purpose: Run permutation tests (based on Ron's script).


library(tidyverse)
library(ggsignif)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)

# Load data ---------------------------------------------------------------

load("RData/13_matched-data.RData")
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

## 1. Prescribed burn -----------------------------------------------------

# Join cover cols
model01.matched2 <- model01.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join)

#   pivot_longer() for cover cols
model01.matched2 <- model01.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "pct_cover"
    )

#   Rename functional group cover types
model01.matched2 <- model01.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub"
           ))

# Calculate observed mean difference
model01.diff <- model01.matched2 |> 
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
    permuted_data <- model01.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(pct_cover), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each functional group
p_values <- model01.perm |>
  inner_join(model01.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values # p = 0.048 for perennial forb

# Boxplot
model01.bp <- model01.matched2 |> 
  ggplot(aes(x = indicators, y = pct_cover, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "1. AZ/NM Mountains: Prescribed burn") +
  theme(legend.title = element_blank()) +
  geom_signif(
    y_position = 25,
    xmin = 2.8,
    xmax = 3.2,
    annotations = c("*")
  ) +
  theme(plot.margin = margin(10, 10, 20, 10))
model01.bp


# Plot frequency distribution
#   Annual forb
model01.annforb <- model01.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.annforb

#   Annual grass
model01.anngrass <- model01.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.anngrass

#   Perennial forb
model01.perforb <- model01.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.perforb

#   Perennial grass
model01.pergrass <- model01.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.pergrass

#   Shrub
model01.shrub <- model01.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.shrub


# Combine plots
grid.arrange(
  model01.bp, model01.annforb, model01.anngrass,
  model01.perforb, model01.pergrass, model01.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)




# Arizona/New Mexico Plateau ----------------------------------------------

## 2. Herbicide -----------------------------------------------------------

# Join cover cols
model02.matched2 <- model02.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join)

#   pivot_longer() for cover cols
model02.matched2 <- model02.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "pct_cover"
  )

#   Rename functional group cover types
model02.matched2 <- model02.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub"
           ))

# Calculate observed mean difference
model02.diff <- model02.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(pct_cover),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model02.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model02.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model02.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(pct_cover), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each functional group
p_values <- model02.perm |>
  inner_join(model02.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values

# Boxplot
model02.bp <- model02.matched2 |> 
  ggplot(aes(x = indicators, y = pct_cover, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "2. AZ/NM Plateau: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(plot.margin = margin(10, 10, 20, 10))
model01.bp


# Plot frequency distribution
#   Annual forb
model02.annforb <- model02.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.annforb

#   Annual grass
model02.anngrass <- model02.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.anngrass

#   Perennial forb
model02.perforb <- model02.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.perforb

#   Perennial grass
model02.pergrass <- model02.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.pergrass

#   Shrub
model02.shrub <- model02.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.shrub


# Combine plots
grid.arrange(
  model02.bp, model02.annforb, model02.anngrass,
  model02.perforb, model02.pergrass, model02.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)



## 3. Prescribed Burn -----------------------------------------------------

# Join cover cols
model03.matched2 <- model03.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join)

#   pivot_longer() for cover cols
model03.matched2 <- model03.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "pct_cover"
  )

#   Rename functional group cover types
model03.matched2 <- model03.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub"
           ))

# Calculate observed mean difference
model03.diff <- model03.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(pct_cover),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model03.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model03.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model03.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(pct_cover), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each functional group
p_values <- model03.perm |>
  inner_join(model03.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values

# Boxplot
model03.bp <- model03.matched2 |> 
  ggplot(aes(x = indicators, y = pct_cover, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "3. AZ/NM Plateau: Prescribed burn") +
  theme(legend.title = element_blank()) +
  theme(plot.margin = margin(10, 10, 20, 10))
model03.bp


# Plot frequency distribution
#   Annual forb
model03.annforb <- model03.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.annforb

#   Annual grass
model03.anngrass <- model03.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.anngrass

#   Perennial forb
model03.perforb <- model03.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.perforb

#   Perennial grass
model03.pergrass <- model03.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.pergrass

#   Shrub
model03.shrub <- model03.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.shrub


# Combine plots
grid.arrange(
  model03.bp, model03.annforb, model03.anngrass,
  model03.perforb, model03.pergrass, model03.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)



## 4. Seeding -------------------------------------------------------------

# Join cover cols
model04.matched2 <- model04.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join)

#   pivot_longer() for cover cols
model04.matched2 <- model04.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "pct_cover"
  )

#   Rename functional group cover types
model04.matched2 <- model04.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub"
           ))

# Calculate observed mean difference
model04.diff <- model04.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(pct_cover),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model04.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model04.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model04.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(pct_cover), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each functional group
p_values <- model04.perm |>
  inner_join(model04.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values

# Boxplot
model04.bp <- model04.matched2 |> 
  ggplot(aes(x = indicators, y = pct_cover, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "4. AZ/NM Mountains: Seeding") +
  theme(legend.title = element_blank()) +
  theme(plot.margin = margin(10, 10, 20, 10))
model04.bp


# Plot frequency distribution
#   Annual forb
model04.annforb <- model04.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.annforb

#   Annual grass
model04.anngrass <- model04.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.anngrass

#   Perennial forb
model04.perforb <- model04.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.perforb

#   Perennial grass
model04.pergrass <- model04.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.pergrass

#   Shrub
model04.shrub <- model04.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.shrub


# Combine plots
grid.arrange(
  model04.bp, model04.annforb, model04.anngrass,
  model04.perforb, model04.pergrass, model04.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)



## 5. Soil Disturbance ----------------------------------------------------

# Join cover cols
model05.matched2 <- model05.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join)

#   pivot_longer() for cover cols
model05.matched2 <- model05.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "pct_cover"
  )

#   Rename functional group cover types
model05.matched2 <- model05.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub"
           ))

# Calculate observed mean difference
model05.diff <- model05.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(pct_cover),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model05.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model05.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model05.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(pct_cover), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each functional group
p_values <- model05.perm |>
  inner_join(model05.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values

# Boxplot
model05.bp <- model05.matched2 |> 
  ggplot(aes(x = indicators, y = pct_cover, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "5. AZ/NM Mountains: Soil disturbance") +
  theme(legend.title = element_blank()) +
  theme(plot.margin = margin(10, 10, 20, 10))
model05.bp


# Plot frequency distribution
#   Annual forb
model05.annforb <- model05.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.annforb

#   Annual grass
model05.anngrass <- model05.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.anngrass

#   Perennial forb
model05.perforb <- model05.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.perforb

#   Perennial grass
model05.pergrass <- model05.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.pergrass

#   Shrub
model05.shrub <- model05.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.shrub


# Combine plots
grid.arrange(
  model05.bp, model05.annforb, model05.anngrass,
  model05.perforb, model05.pergrass, model05.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)



# Blue Mountains ----------------------------------------------------------

## 6. Herbicide -----------------------------------------------------------

# Join cover cols
model06.matched2 <- model06.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join)

#   pivot_longer() for cover cols
model06.matched2 <- model06.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "pct_cover"
  )

#   Rename functional group cover types
model06.matched2 <- model06.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub"
           ))

# Calculate observed mean difference
model06.diff <- model06.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(pct_cover),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model06.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model06.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model06.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(pct_cover), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each functional group
p_values <- model06.perm |>
  inner_join(model06.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values

# Boxplot
model06.bp <- model06.matched2 |> 
  ggplot(aes(x = indicators, y = pct_cover, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "6. Blue Mountains: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(plot.margin = margin(10, 10, 20, 10))
model06.bp


# Plot frequency distribution
#   Annual forb
model06.annforb <- model06.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.annforb

#   Annual grass
model06.anngrass <- model06.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.anngrass

#   Perennial forb
model06.perforb <- model06.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.perforb

#   Perennial grass
model06.pergrass <- model06.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.pergrass

#   Shrub
model06.shrub <- model06.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.shrub


# Combine plots
grid.arrange(
  model06.bp, model06.annforb, model06.anngrass,
  model06.perforb, model06.pergrass, model06.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)



## 7. Vegetation Disturbance ----------------------------------------------

# Join cover cols
model07.matched2 <- model07.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join)

#   pivot_longer() for cover cols
model07.matched2 <- model07.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "pct_cover"
  )

#   Rename functional group cover types
model07.matched2 <- model07.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub"
           ))

# Calculate observed mean difference
model07.diff <- model07.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(pct_cover),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model07.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model07.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model07.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(pct_cover), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each functional group
p_values <- model07.perm |>
  inner_join(model07.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values

# Boxplot
model07.bp <- model07.matched2 |> 
  ggplot(aes(x = indicators, y = pct_cover, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "7. Blue Mountains: Vegetation disturbance") +
  theme(legend.title = element_blank()) +
  theme(plot.margin = margin(10, 10, 20, 10))
model07.bp


# Plot frequency distribution
#   Annual forb
model07.annforb <- model07.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.annforb

#   Annual grass
model07.anngrass <- model07.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.anngrass

#   Perennial forb
model07.perforb <- model07.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.perforb

#   Perennial grass
model07.pergrass <- model07.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.pergrass

#   Shrub
model07.shrub <- model07.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw() +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.shrub


# Combine plots
grid.arrange(
  model07.bp, model07.annforb, model07.anngrass,
  model07.perforb, model07.pergrass, model07.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)

save.image("RData/14_permutation-tests.RData")

