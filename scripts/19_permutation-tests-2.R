# Updated: 2026-06-01
# Updated: 2026-06-02

# Purpose: Run permutation tests for functional group cover, Shannon diversity,
#   and invasive species. (For now RHEM output is on hold due to all the missing values.)


library(tidyverse)
library(ggsignif)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)

# Load data ---------------------------------------------------------------

load("RData/13_matched-data.RData")
geoindicators.raw <- read_csv("data/raw/downloaded/ldc-data-2026-03-11/geoindicators.csv")
all_diversity <- read_csv("data/versions-from-R/16_shannon-diversity_all-models.csv")
load("RData/17_invasive-matched-data.RData")


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

# Join cover & shannon cols
model01.matched2 <- model01.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model01.invasive) |> 
  left_join(filter(all_diversity, Model == "model01")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model01.matched2 <- model01.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model01.matched2 <- model01.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus rubens"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus rubens")))

# Calculate observed mean difference
model01.diff <- model01.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
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
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values01 <- model01.perm |>
  inner_join(model01.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values01 # p = 0.04 for perennial forb; p = 0.02 for shannon

# Boxplot
model01.bp <- model01.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
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
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus rubens" = expression(italic("Bromus rubens")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model01.bp

# Plot frequency distribution
#   Annual forb
model01.annforb <- model01.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.annforb

#   Annual grass
model01.anngrass <- model01.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.anngrass

#   Perennial forb
model01.perforb <- model01.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.perforb

#   Perennial grass
model01.pergrass <- model01.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.pergrass

#   Shrub
model01.shrub <- model01.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.shrub

#   Shannon diversity
model01.shannon <- model01.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.shannon

#   Bromus rubens
model01.brru2 <- model01.perm |> 
  filter(indicators == "Bromus rubens") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model01.diff$obs_diff[model01.diff$indicators == "Bromus rubens"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus rubens"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model01.brru2

# Combine plots
grid.arrange(
  model01.bp, model01.annforb, model01.anngrass,
  model01.perforb, model01.pergrass, model01.shrub,
  model01.shannon, model01.brru2,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)



# Arizona/New Mexico Plateau ----------------------------------------------

## 2. Herbicide -----------------------------------------------------------

# Join cover & shannon cols
model02.matched2 <- model02.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model02.invasive) |> 
  left_join(filter(all_diversity, Model == "model02")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model02.matched2 <- model02.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
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
  summarise(mean_cover = mean(value),
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
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values02 <- model02.perm |>
  inner_join(model02.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values02

# Boxplot
model02.bp <- model02.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "2. AZ/NM Plateau: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("B" = expression(italic("B")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model01.bp

# Plot frequency distribution
#   Annual forb
model02.annforb <- model02.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.annforb

#   Annual grass
model02.anngrass <- model02.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.anngrass

#   Perennial forb
model02.perforb <- model02.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.perforb

#   Perennial grass
model02.pergrass <- model02.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model02.pergrass

#   Shrub
model02.shrub <- model02.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model02.diff$obs_diff[model02.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
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


## 3. Prescribed burn -----------------------------------------------------

# Join cover & shannon cols
model03.matched2 <- model03.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join)

#   pivot_longer() for cover & shannon cols
model03.matched2 <- model03.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
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
  summarise(mean_cover = mean(value),
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
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values03 <- model03.perm |>
  inner_join(model03.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values03

# Boxplot
model03.bp <- model03.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "3. AZ/NM Plateau: Prescribed burn") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model03.bp

# Plot frequency distribution
#   Annual forb
model03.annforb <- model03.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.annforb

#   Annual grass
model03.anngrass <- model03.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.anngrass

#   Perennial forb
model03.perforb <- model03.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.perforb

#   Perennial grass
model03.pergrass <- model03.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model03.pergrass

#   Shrub
model03.shrub <- model03.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model03.diff$obs_diff[model03.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
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

# Join cover & shannon cols
model04.matched2 <- model04.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join)

#   pivot_longer() for cover & shannon cols
model04.matched2 <- model04.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
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
  summarise(mean_cover = mean(value),
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
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values04 <- model04.perm |>
  inner_join(model04.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values04

# Boxplot
model04.bp <- model04.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "4. AZ/NM Mountains: Seeding") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model04.bp

# Plot frequency distribution
#   Annual forb
model04.annforb <- model04.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.annforb

#   Annual grass
model04.anngrass <- model04.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.anngrass

#   Perennial forb
model04.perforb <- model04.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.perforb

#   Perennial grass
model04.pergrass <- model04.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model04.pergrass

#   Shrub
model04.shrub <- model04.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model04.diff$obs_diff[model04.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
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


## 5. Soil disturbance ----------------------------------------------------

# Join cover & shannon cols
model05.matched2 <- model05.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join)

#   pivot_longer() for cover & shannon cols
model05.matched2 <- model05.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
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
  summarise(mean_cover = mean(value),
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
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values05 <- model05.perm |>
  inner_join(model05.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values05

# Boxplot
model05.bp <- model05.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "5. AZ/NM Mountains: Soil disturbance") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model05.bp

# Plot frequency distribution
#   Annual forb
model05.annforb <- model05.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.annforb

#   Annual grass
model05.anngrass <- model05.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.anngrass

#   Perennial forb
model05.perforb <- model05.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.perforb

#   Perennial grass
model05.pergrass <- model05.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model05.pergrass

#   Shrub
model05.shrub <- model05.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model05.diff$obs_diff[model05.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
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

# Join cover & shannon cols
model06.matched2 <- model06.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model06.invasive) |> 
  left_join(filter(all_diversity, Model == "model06")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model06.matched2 <- model06.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model06.matched2 <- model06.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model06.diff <- model06.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
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
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values06 <- model06.perm |>
  inner_join(model06.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values06

# Boxplot
model06.bp <- model06.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "6. Blue Mountains: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model06.bp

# Plot frequency distribution
#   Annual forb
model06.annforb <- model06.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.annforb

#   Annual grass
model06.anngrass <- model06.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.anngrass

#   Perennial forb
model06.perforb <- model06.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.perforb

#   Perennial grass
model06.pergrass <- model06.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.pergrass

#   Shrub
model06.shrub <- model06.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.shrub

#   Shannon diversity
model06.shannon <- model06.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.shannon

#   Bromus tectorum
model06.brte <- model06.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model06.diff$obs_diff[model06.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model06.brte

# Combine plots
grid.arrange(
  model06.bp, model06.annforb, model06.anngrass,
  model06.perforb, model06.pergrass, model06.shrub,
  model06.shannon, model06.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 7. Vegetation disturbance ----------------------------------------------

# Join cover & shannon cols
model07.matched2 <- model07.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model07.invasive) |> 
  left_join(filter(all_diversity, Model == "model07")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model07.matched2 <- model07.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model07.matched2 <- model07.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model07.diff <- model07.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
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
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values07 <- model07.perm |>
  inner_join(model07.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values07

# Boxplot
model07.bp <- model07.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "7. Blue Mountains: Vegetation disturbance") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model07.bp

# Plot frequency distribution
#   Annual forb
model07.annforb <- model07.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.annforb

#   Annual grass
model07.anngrass <- model07.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.anngrass

#   Perennial forb
model07.perforb <- model07.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.perforb

#   Perennial grass
model07.pergrass <- model07.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.pergrass

#   Shrub
model07.shrub <- model07.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.shrub

#   Shannon diversity
model07.shannon <- model07.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.shannon

#   Bromus tectorum
model07.brte <- model07.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model07.diff$obs_diff[model07.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model07.brte

# Combine plots
grid.arrange(
  model07.bp, model07.annforb, model07.anngrass,
  model07.perforb, model07.pergrass, model07.shrub,
  model07.shannon, model07.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 8. Post-burn herbicide -------------------------------------------------

# Join cover & shannon cols
model08.matched2 <- model08.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model08.invasive) |> 
  left_join(filter(all_diversity, Model == "model08")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model08.matched2 <- model08.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model08.matched2 <- model08.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model08.diff <- model08.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model08.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model08.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model08.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values08 <- model08.perm |>
  inner_join(model08.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values08

# Boxplot
model08.bp <- model08.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "8. Blue Mountains: Post-burn herbicide") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model08.bp

# Plot frequency distribution
#   Annual forb
model08.annforb <- model08.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model08.diff$obs_diff[model08.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model08.annforb

#   Annual grass
model08.anngrass <- model08.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model08.diff$obs_diff[model08.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model08.anngrass

#   Perennial forb
model08.perforb <- model08.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model08.diff$obs_diff[model08.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model08.perforb

#   Perennial grass
model08.pergrass <- model08.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model08.diff$obs_diff[model08.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model08.pergrass

#   Shrub
model08.shrub <- model08.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model08.diff$obs_diff[model08.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model08.shrub

#   Shannon diversity
model08.shannon <- model08.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model08.diff$obs_diff[model08.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model08.shannon

#   Bromus tectorum
model08.brte <- model08.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model08.diff$obs_diff[model08.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model08.brte

# Combine plots
grid.arrange(
  model08.bp, model08.annforb, model08.anngrass,
  model08.perforb, model08.pergrass, model08.shrub,
  model08.shannon, model08.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)



# Central Basin and Range -------------------------------------------------

## 9. Aerial seeding ------------------------------------------------------

# Join cover & shannon cols
model09.matched2 <- model09.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model09.invasive) |> 
  left_join(filter(all_diversity, Model == "model09")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model09.matched2 <- model09.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model09.matched2 <- model09.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model09.diff <- model09.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model09.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model09.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model09.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values09 <- model09.perm |>
  inner_join(model09.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values09 # p = 0.046 for perennial herb; p = 0.03 for shannon

# Boxplot
model09.bp <- model09.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "9. Central Basin and Range: Aerial seeding") +
  geom_signif(
    y_position = 25,
    xmin = 2.8,
    xmax = 3.2,
    annotations = c("*")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model09.bp

# Plot frequency distribution
#   Annual forb
model09.annforb <- model09.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model09.diff$obs_diff[model09.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model09.annforb

#   Annual grass
model09.anngrass <- model09.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model09.diff$obs_diff[model09.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model09.anngrass

#   Perennial forb
model09.perforb <- model09.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model09.diff$obs_diff[model09.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model09.perforb

#   Perennial grass
model09.pergrass <- model09.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model09.diff$obs_diff[model09.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model09.pergrass

#   Shrub
model09.shrub <- model09.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model09.diff$obs_diff[model09.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model09.shrub

#   Shannon diversity
model09.shannon <- model09.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model09.diff$obs_diff[model09.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model09.shannon

#   Bromus tectorum
model09.brte <- model09.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model09.diff$obs_diff[model09.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model09.brte

# Combine plots
grid.arrange(
  model09.bp, model09.annforb, model09.anngrass,
  model09.perforb, model09.pergrass, model09.shrub,
  model09.shannon, model09.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 10. Drill seeding & soil disturbance -----------------------------------

# Join cover & shannon cols
model10.matched2 <- model10.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model10.invasive) |> 
  left_join(filter(all_diversity, Model == "model10")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model10.matched2 <- model10.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model10.matched2 <- model10.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model10.diff <- model10.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model10.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model10.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model10.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values10 <- model10.perm |>
  inner_join(model10.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values10 # p = 0.0009 for perennial grass

# Boxplot
model10.bp <- model10.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "10. Central Basin and Range: Drill seeding & soil disturbance") +
  geom_signif(
    y_position = 82,
    xmin = 3.8,
    xmax = 4.2, 
    annotations = c("***")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  ylim(0, 90) +
  theme(plot.margin = margin(10, 10, 20, 10))
model10.bp

# Plot frequency distribution
#   Annual forb
model10.annforb <- model10.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model10.diff$obs_diff[model10.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model10.annforb

#   Annual grass
model10.anngrass <- model10.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model10.diff$obs_diff[model10.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model10.anngrass

#   Perennial forb
model10.perforb <- model10.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model10.diff$obs_diff[model10.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model10.perforb

#   Perennial grass
model10.pergrass <- model10.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model10.diff$obs_diff[model10.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model10.pergrass

#   Shrub
model10.shrub <- model10.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model10.diff$obs_diff[model10.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model10.shrub

#   Shannon diversity
model10.shannon <- model10.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model10.diff$obs_diff[model10.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model10.shannon

#   Bromus tectorum
model10.brte <- model10.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model10.diff$obs_diff[model10.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model10.brte

# Combine plots
grid.arrange(
  model10.bp, model10.annforb, model10.anngrass,
  model10.perforb, model10.pergrass, model10.shrub,
  model10.shannon, model10.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 11. Prescribed burn ----------------------------------------------------

# Join cover & shannon cols
model11.matched2 <- model11.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model11.invasive) |> 
  left_join(filter(all_diversity, Model == "model11")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model11.matched2 <- model11.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model11.matched2 <- model11.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model11.diff <- model11.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model11.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model11.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model11.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values11 <- model11.perm |>
  inner_join(model11.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values11 # p = 0.003 for perennial grass; p = 0.001 for shannon

# Boxplot
model11.bp <- model11.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "11. Central Basin and Range: Prescribed burn") +
  geom_signif(
    y_position = 65,
    xmin = 3.8,
    xmax = 4.2, 
    annotations = c("**")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model11.bp

# Plot frequency distribution
#   Annual forb
model11.annforb <- model11.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model11.diff$obs_diff[model11.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model11.annforb

#   Annual grass
model11.anngrass <- model11.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model11.diff$obs_diff[model11.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model11.anngrass

#   Perennial forb
model11.perforb <- model11.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model11.diff$obs_diff[model11.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model11.perforb

#   Perennial grass
model11.pergrass <- model11.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model11.diff$obs_diff[model11.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model11.pergrass

#   Shrub
model11.shrub <- model11.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model11.diff$obs_diff[model11.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model11.shrub

#   Shannon diversity
model11.shannon <- model11.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model11.diff$obs_diff[model11.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model11.shannon

#   Bromus tectorum
model11.brte <- model11.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model11.diff$obs_diff[model11.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model11.brte

# Combine plots
grid.arrange(
  model11.bp, model11.annforb, model11.anngrass,
  model11.perforb, model11.pergrass, model11.shrub,
  model11.shannon, model11.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 12. Vegetation disturbance ---------------------------------------------

# Join cover & shannon cols
model12.matched2 <- model12.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model12.invasive) |> 
  left_join(filter(all_diversity, Model == "model12")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model12.matched2 <- model12.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model12.matched2 <- model12.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model12.diff <- model12.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model12.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model12.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model12.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values12 <- model12.perm |>
  inner_join(model12.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values12 # p = 0.007 for shannon

# Boxplot
model12.bp <- model12.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "12. Central Basin and Range: Vegetation disturbance") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model12.bp

# Plot frequency distribution
#   Annual forb
model12.annforb <- model12.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model12.diff$obs_diff[model12.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model12.annforb

#   Annual grass
model12.anngrass <- model12.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model12.diff$obs_diff[model12.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model12.anngrass

#   Perennial forb
model12.perforb <- model12.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model12.diff$obs_diff[model12.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model12.perforb

#   Perennial grass
model12.pergrass <- model12.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model12.diff$obs_diff[model12.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model12.pergrass

#   Shrub
model12.shrub <- model12.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model12.diff$obs_diff[model12.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model12.shrub

#   Shannon diversity
model12.shannon <- model12.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model12.diff$obs_diff[model12.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model12.shannon

#   Bromus tectorum
model12.brte <- model12.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model12.diff$obs_diff[model12.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model12.brte

# Combine plots
grid.arrange(
  model12.bp, model12.annforb, model12.anngrass,
  model12.perforb, model12.pergrass, model12.shrub,
  model12.shannon, model12.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 13. Post-burn aerial seeding -------------------------------------------

# Join cover & shannon cols
model13.matched2 <- model13.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model13.invasive) |> 
  left_join(filter(all_diversity, Model == "model13")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n) |> 
  filter(!is.na(Shannon)) # one row must be removed because there were no cover 
#                             measurements for any species for this plot 
#               (Primary key: NV_NV-Elko-DO-ESR-2021_ELKO-ShafterComplex-01B_V12021-09-01)

#   pivot_longer() for cover & shannon cols
model13.matched2 <- model13.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model13.matched2 <- model13.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model13.diff <- model13.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model13.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model13.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model13.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values13 <- model13.perm |>
  inner_join(model13.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values13 # p = 0.003 for perennial forb, p < 0.001 for shannon

# Boxplot
model13.bp <- model13.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "13. Central Basin and Range: Post-burn aerial seeding") +
  geom_signif(
    y_position = 60,
    xmin = 2.8,
    xmax = 3.2, 
    annotations = c("**")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model13.bp

# Plot frequency distribution
#   Annual forb
model13.annforb <- model13.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model13.diff$obs_diff[model13.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model13.annforb

#   Annual grass
model13.anngrass <- model13.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model13.diff$obs_diff[model13.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model13.anngrass

#   Perennial forb
model13.perforb <- model13.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model13.diff$obs_diff[model13.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model13.perforb

#   Perennial grass
model13.pergrass <- model13.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model13.diff$obs_diff[model13.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model13.pergrass

#   Shrub
model13.shrub <- model13.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model13.diff$obs_diff[model13.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model13.shrub

#   Shannon diversity
model13.shannon <- model13.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model13.diff$obs_diff[model13.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model13.shannon

#   Bromus tectorum
model13.brte <- model13.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model13.diff$obs_diff[model13.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model13.brte

# Combine plots
grid.arrange(
  model13.bp, model13.annforb, model13.anngrass,
  model13.perforb, model13.pergrass, model13.shrub,
  model13.shannon, model13.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 14. Post-burn drill seeding --------------------------------------------

# Join cover & shannon cols
model14.matched2 <- model14.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model14.invasive) |> 
  left_join(filter(all_diversity, Model == "model14")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model14.matched2 <- model14.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model14.matched2 <- model14.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model14.diff <- model14.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model14.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model14.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model14.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values14 <- model14.perm |>
  inner_join(model14.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values14

# Boxplot
model14.bp <- model14.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "14. Central Basin and Range: Post-burn drill seeding") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model14.bp

# Plot frequency distribution
#   Annual forb
model14.annforb <- model14.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model14.diff$obs_diff[model14.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model14.annforb

#   Annual grass
model14.anngrass <- model14.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model14.diff$obs_diff[model14.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model14.anngrass

#   Perennial forb
model14.perforb <- model14.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model14.diff$obs_diff[model14.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model14.perforb

#   Perennial grass
model14.pergrass <- model14.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model14.diff$obs_diff[model14.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model14.pergrass

#   Shrub
model14.shrub <- model14.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model14.diff$obs_diff[model14.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model14.shrub

#   Shannon diversity
model14.shannon <- model14.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model14.diff$obs_diff[model14.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model14.shannon

#   Bromus tectorum
model14.brte <- model14.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model14.diff$obs_diff[model14.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model14.brte

# Combine plots
grid.arrange(
  model14.bp, model14.annforb, model14.anngrass,
  model14.perforb, model14.pergrass, model14.shrub,
  model14.shannon, model14.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 15. Post-burn ground seeding -------------------------------------------

# Join cover & shannon cols
model15.matched2 <- model15.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model15.invasive) |> 
  left_join(filter(all_diversity, Model == "model15")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model15.matched2 <- model15.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model15.matched2 <- model15.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model15.diff <- model15.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model15.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model15.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model15.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values15 <- model15.perm |>
  inner_join(model15.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values15

# Boxplot
model15.bp <- model15.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "15. Central Basin and Range: Post-burn ground seeding") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model15.bp

# Plot frequency distribution
#   Annual forb
model15.annforb <- model15.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model15.diff$obs_diff[model15.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model15.annforb

#   Annual grass
model15.anngrass <- model15.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model15.diff$obs_diff[model15.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model15.anngrass

#   Perennial forb
model15.perforb <- model15.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model15.diff$obs_diff[model15.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model15.perforb

#   Perennial grass
model15.pergrass <- model15.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model15.diff$obs_diff[model15.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model15.pergrass

#   Shrub
model15.shrub <- model15.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model15.diff$obs_diff[model15.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model15.shrub

#   Shannon diversity
model15.shannon <- model15.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model15.diff$obs_diff[model15.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model15.shannon

#   Bromus tectorum
model15.brte <- model15.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model15.diff$obs_diff[model15.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model15.brte

# Combine plots
grid.arrange(
  model15.bp, model15.annforb, model15.anngrass,
  model15.perforb, model15.pergrass, model15.shrub,
  model15.shannon, model15.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 16. Post-burn herbicide ------------------------------------------------

# Join cover & shannon cols
model16.matched2 <- model16.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model16.invasive) |> 
  left_join(filter(all_diversity, Model == "model16")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n) |> 
  filter(!is.na(Shannon)) # one row must be removed because there were no cover 
#                             measurements for any species for this plot 
#               (Primary key: 20184945202113B2)

#   pivot_longer() for cover & shannon cols
model16.matched2 <- model16.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model16.matched2 <- model16.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model16.diff <- model16.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model16.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model16.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model16.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values16 <- model16.perm |>
  inner_join(model16.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values16

# Boxplot
model16.bp <- model16.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "16. Central Basin and Range: Post-burn herbicide") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model16.bp

# Plot frequency distribution
#   Annual forb
model16.annforb <- model16.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model16.diff$obs_diff[model16.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model16.annforb

#   Annual grass
model16.anngrass <- model16.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model16.diff$obs_diff[model16.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model16.anngrass

#   Perennial forb
model16.perforb <- model16.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model16.diff$obs_diff[model16.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model16.perforb

#   Perennial grass
model16.pergrass <- model16.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model16.diff$obs_diff[model16.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model16.pergrass

#   Shrub
model16.shrub <- model16.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model16.diff$obs_diff[model16.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model16.shrub

#   Shannon diversity
model16.shannon <- model16.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model16.diff$obs_diff[model16.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model16.shannon

#   Bromus tectorum
model16.brte <- model16.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model16.diff$obs_diff[model16.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model16.brte

# Combine plots
grid.arrange(
  model16.bp, model16.annforb, model16.anngrass,
  model16.perforb, model16.pergrass, model16.shrub,
  model16.shannon, model16.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)



# Chihuahuan Desert -------------------------------------------------------

## 17. Herbicide ----------------------------------------------------------

# Join cover & shannon cols
model17.matched2 <- model17.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model17.invasive) |> 
  left_join(filter(all_diversity, Model == "model17")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model17.matched2 <- model17.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model17.matched2 <- model17.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Prosopis glandulosa"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Prosopis glandulosa")))

# Calculate observed mean difference
model17.diff <- model17.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model17.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model17.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model17.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values17 <- model17.perm |>
  inner_join(model17.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values17

# Boxplot
model17.bp <- model17.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "17. Chihuahuan Desert: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Prosopis glandulosa" = expression(italic("Prosopis glandulosa")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model17.bp

# Plot frequency distribution
#   Annual forb
model17.annforb <- model17.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model17.diff$obs_diff[model17.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model17.annforb

#   Annual grass
model17.anngrass <- model17.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model17.diff$obs_diff[model17.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model17.anngrass

#   Perennial forb
model17.perforb <- model17.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model17.diff$obs_diff[model17.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model17.perforb

#   Perennial grass
model17.pergrass <- model17.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model17.diff$obs_diff[model17.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model17.pergrass

#   Shrub
model17.shrub <- model17.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model17.diff$obs_diff[model17.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model17.shrub

#   Shannon diversity
model17.shannon <- model17.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model17.diff$obs_diff[model17.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model17.shannon

#   Prosopis glandulosa
model17.prgl2 <- model17.perm |> 
  filter(indicators == "Prosopis glandulosa") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model17.diff$obs_diff[model17.diff$indicators == "Prosopis glandulosa"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Prosopis glandulosa"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model17.prgl2

# Combine plots
grid.arrange(
  model17.bp, model17.annforb, model17.anngrass,
  model17.perforb, model17.pergrass, model17.shrub,
  model17.shannon, model17.prgl2,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)



# Colorado Plateaus -------------------------------------------------------

## 18. Aerial seeding & soil disturbance ----------------------------------

# Join cover & shannon cols
model18.matched2 <- model18.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model18.invasive) |> 
  left_join(filter(all_diversity, Model == "model18")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model18.matched2 <- model18.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model18.matched2 <- model18.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model18.diff <- model18.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model18.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model18.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model18.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values18 <- model18.perm |>
  inner_join(model18.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values18 # p = 0.011 for shannon

# Boxplot
model18.bp <- model18.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "18. Colorado Plateaus: Aerial seeding & soil disturbance") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model18.bp

# Plot frequency distribution
#   Annual forb
model18.annforb <- model18.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model18.diff$obs_diff[model18.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model18.annforb

#   Annual grass
model18.anngrass <- model18.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model18.diff$obs_diff[model18.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model18.anngrass

#   Perennial forb
model18.perforb <- model18.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model18.diff$obs_diff[model18.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model18.perforb

#   Perennial grass
model18.pergrass <- model18.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model18.diff$obs_diff[model18.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model18.pergrass

#   Shrub
model18.shrub <- model18.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model18.diff$obs_diff[model18.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model18.shrub

#   Shannon diversity
model18.shannon <- model18.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model18.diff$obs_diff[model18.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model18.shannon

#   Bromus tectorum
model18.brte <- model18.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model18.diff$obs_diff[model18.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model18.brte

# Combine plots
grid.arrange(
  model18.bp, model18.annforb, model18.anngrass,
  model18.perforb, model18.pergrass, model18.shrub,
  model18.shannon, model18.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 19. Herbicide ----------------------------------------------------------

# Join cover & shannon cols
model19.matched2 <- model19.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model19.invasive) |> 
  left_join(filter(all_diversity, Model == "model19")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model19.matched2 <- model19.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model19.matched2 <- model19.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model19.diff <- model19.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model19.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model19.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model19.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values19 <- model19.perm |>
  inner_join(model19.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values19

# Boxplot
model19.bp <- model19.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "19. Colorado Plateaus: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model19.bp

# Plot frequency distribution
#   Annual forb
model19.annforb <- model19.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model19.diff$obs_diff[model19.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model19.annforb

#   Annual grass
model19.anngrass <- model19.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model19.diff$obs_diff[model19.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model19.anngrass

#   Perennial forb
model19.perforb <- model19.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model19.diff$obs_diff[model19.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model19.perforb

#   Perennial grass
model19.pergrass <- model19.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model19.diff$obs_diff[model19.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model19.pergrass

#   Shrub
model19.shrub <- model19.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model19.diff$obs_diff[model19.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model19.shrub

#   Shannon diversity
model19.shannon <- model19.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model19.diff$obs_diff[model19.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model19.shannon

#   Bromus tectorum
model19.brte <- model19.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model19.diff$obs_diff[model19.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model19.brte

# Combine plots
grid.arrange(
  model19.bp, model19.annforb, model19.anngrass,
  model19.perforb, model19.pergrass, model19.shrub,
  model19.shannon, model19.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 20. Prescribed burn ----------------------------------------------------

# Join cover & shannon cols
model20.matched2 <- model20.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model20.invasive) |> 
  left_join(filter(all_diversity, Model == "model20")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model20.matched2 <- model20.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model20.matched2 <- model20.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model20.diff <- model20.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model20.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model20.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model20.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values20 <- model20.perm |>
  inner_join(model20.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values20

# Boxplot
model20.bp <- model20.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "20. Colorado Plateaus: Prescribed burn") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model20.bp

# Plot frequency distribution
#   Annual forb
model20.annforb <- model20.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model20.diff$obs_diff[model20.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model20.annforb

#   Annual grass
model20.anngrass <- model20.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model20.diff$obs_diff[model20.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model20.anngrass

#   Perennial forb
model20.perforb <- model20.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model20.diff$obs_diff[model20.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model20.perforb

#   Perennial grass
model20.pergrass <- model20.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model20.diff$obs_diff[model20.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model20.pergrass

#   Shrub
model20.shrub <- model20.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model20.diff$obs_diff[model20.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model20.shrub

#   Shannon diversity
model20.shannon <- model20.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model20.diff$obs_diff[model20.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model20.shannon

#   Bromus tectorum
model20.brte <- model20.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model20.diff$obs_diff[model20.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model20.brte

# Combine plots
grid.arrange(
  model20.bp, model20.annforb, model20.anngrass,
  model20.perforb, model20.pergrass, model20.shrub,
  model20.shannon, model20.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 21. Soil disturbance ---------------------------------------------------

# Join cover & shannon cols
model21.matched2 <- model21.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model21.invasive) |> 
  left_join(filter(all_diversity, Model == "model21")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model21.matched2 <- model21.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model21.matched2 <- model21.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model21.diff <- model21.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model21.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model21.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model21.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values21 <- model21.perm |>
  inner_join(model21.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values21

# Boxplot
model21.bp <- model21.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "21. Colorado Plateaus: Soil disturbance") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model21.bp

# Plot frequency distribution
#   Annual forb
model21.annforb <- model21.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model21.diff$obs_diff[model21.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model21.annforb

#   Annual grass
model21.anngrass <- model21.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model21.diff$obs_diff[model21.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model21.anngrass

#   Perennial forb
model21.perforb <- model21.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model21.diff$obs_diff[model21.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model21.perforb

#   Perennial grass
model21.pergrass <- model21.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model21.diff$obs_diff[model21.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model21.pergrass

#   Shrub
model21.shrub <- model21.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model21.diff$obs_diff[model21.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model21.shrub

#   Shannon diversity
model21.shannon <- model21.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model21.diff$obs_diff[model21.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model21.shannon

#   Bromus tectorum
model21.brte <- model21.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model21.diff$obs_diff[model21.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model21.brte

# Combine plots
grid.arrange(
  model21.bp, model21.annforb, model21.anngrass,
  model21.perforb, model21.pergrass, model21.shrub,
  model21.shannon, model21.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 22. Vegetation disturbance ---------------------------------------------

# Join cover & shannon cols
model22.matched2 <- model22.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model22.invasive) |> 
  left_join(filter(all_diversity, Model == "model22")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)


#   pivot_longer() for cover & shannon cols
model22.matched2 <- model22.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model22.matched2 <- model22.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model22.diff <- model22.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model22.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model22.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model22.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values22 <- model22.perm |>
  inner_join(model22.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values22 # p = 0.002 for annual grass; p = 0.007 for BRTE

# Boxplot
model22.bp <- model22.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "22. Colorado Plateaus: Vegetation disturbance") +
  geom_signif(
    y_position = 83,
    xmin = 1.8,
    xmax = 2.2, 
    annotations = c("**")
  ) +
  geom_signif(
    y_position = 65,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("**")
  ) +
  ylim(0, 90) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model22.bp

# Plot frequency distribution
#   Annual forb
model22.annforb <- model22.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model22.diff$obs_diff[model22.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model22.annforb

#   Annual grass
model22.anngrass <- model22.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model22.diff$obs_diff[model22.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model22.anngrass

#   Perennial forb
model22.perforb <- model22.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model22.diff$obs_diff[model22.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model22.perforb

#   Perennial grass
model22.pergrass <- model22.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model22.diff$obs_diff[model22.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model22.pergrass

#   Shrub
model22.shrub <- model22.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model22.diff$obs_diff[model22.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model22.shrub

#   Shannon diversity
model22.shannon <- model22.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model22.diff$obs_diff[model22.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model22.shannon

#   Bromus tectorum
model22.brte <- model22.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model22.diff$obs_diff[model22.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum") ~ "(**)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model22.brte

# Combine plots
grid.arrange(
  model22.bp, model22.annforb, model22.anngrass,
  model22.perforb, model22.pergrass, model22.shrub,
  model22.shannon, model22.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 23. Post-burn aerial seeding -------------------------------------------

# Join cover & shannon cols
model23.matched2 <- model23.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model23.invasive) |> 
  left_join(filter(all_diversity, Model == "model23")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model23.matched2 <- model23.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model23.matched2 <- model23.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model23.diff <- model23.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model23.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model23.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model23.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values23 <- model23.perm |>
  inner_join(model23.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values23 # p < 0.0001 for shannon 

# Boxplot
model23.bp <- model23.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "23. Colorado Plateaus: Post-burn aerial seeding") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model23.bp

# Plot frequency distribution
#   Annual forb
model23.annforb <- model23.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model23.diff$obs_diff[model23.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model23.annforb

#   Annual grass
model23.anngrass <- model23.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model23.diff$obs_diff[model23.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model23.anngrass

#   Perennial forb
model23.perforb <- model23.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model23.diff$obs_diff[model23.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model23.perforb

#   Perennial grass
model23.pergrass <- model23.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model23.diff$obs_diff[model23.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model23.pergrass

#   Shrub
model23.shrub <- model23.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model23.diff$obs_diff[model23.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model23.shrub

#   Shannon diversity
model23.shannon <- model23.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model23.diff$obs_diff[model23.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model23.shannon

#   Bromus tectorum
model23.brte <- model23.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model23.diff$obs_diff[model23.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model23.brte

# Combine plots
grid.arrange(
  model23.bp, model23.annforb, model23.anngrass,
  model23.perforb, model23.pergrass, model23.shrub,
  model23.shannon, model23.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)



# Middle Rockies ----------------------------------------------------------

## 24. Herbicide ----------------------------------------------------------

# Join cover & shannon cols
model24.matched2 <- model24.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join)

#   pivot_longer() for cover & shannon cols
model24.matched2 <- model24.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model24.matched2 <- model24.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub"
           ))

# Calculate observed mean difference
model24.diff <- model24.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model24.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model24.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model24.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values24 <- model24.perm |>
  inner_join(model24.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values24

# Boxplot
model24.bp <- model24.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "24. Middle Rockies: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model24.bp

# Plot frequency distribution
#   Annual forb
model24.annforb <- model24.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model24.diff$obs_diff[model24.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model24.annforb

#   Annual grass
model24.anngrass <- model24.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model24.diff$obs_diff[model24.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model24.anngrass

#   Perennial forb
model24.perforb <- model24.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model24.diff$obs_diff[model24.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model24.perforb

#   Perennial grass
model24.pergrass <- model24.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model24.diff$obs_diff[model24.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model24.pergrass

#   Shrub
model24.shrub <- model24.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model24.diff$obs_diff[model24.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model24.shrub

# Combine plots
grid.arrange(
  model24.bp, model24.annforb, model24.anngrass,
  model24.perforb, model24.pergrass, model24.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)



# Mojave Basin and Range --------------------------------------------------

## 25. Post-burn aerial seeding -------------------------------------------

# Join cover & shannon cols
model25.matched2 <- model25.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model25.invasive) |> 
  left_join(filter(all_diversity, Model == "model25")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model25.matched2 <- model25.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model25.matched2 <- model25.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus rubens"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus rubens")))

# Calculate observed mean difference
model25.diff <- model25.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model25.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model25.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model25.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values25 <- model25.perm |>
  inner_join(model25.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values25 # p = 0.004 for BRRU2

# Boxplot
model25.bp <- model25.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "25. Mojave Basin and Range: Post-burn aerial seeding") +
  geom_signif(
    y_position = 80,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("**")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus rubens" = expression(italic("Bromus rubens")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model25.bp

# Plot frequency distribution
#   Annual forb
model25.annforb <- model25.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model25.diff$obs_diff[model25.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model25.annforb

#   Annual grass
model25.anngrass <- model25.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model25.diff$obs_diff[model25.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model25.anngrass

#   Perennial forb
model25.perforb <- model25.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model25.diff$obs_diff[model25.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model25.perforb

#   Perennial grass
model25.pergrass <- model25.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model25.diff$obs_diff[model25.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model25.pergrass

#   Shrub
model25.shrub <- model25.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model25.diff$obs_diff[model25.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model25.shrub

#   Shannon diversity
model25.shannon <- model25.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model25.diff$obs_diff[model25.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model25.shannon

#   Bromus rubens
model25.brru2 <- model25.perm |> 
  filter(indicators == "Bromus rubens") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model25.diff$obs_diff[model25.diff$indicators == "Bromus rubens"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus rubens") ~ "(**)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model25.brru2

# Combine plots
grid.arrange(
  model25.bp, model25.annforb, model25.anngrass,
  model25.perforb, model25.pergrass, model25.shrub,
  model25.shannon, model25.brru2,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)



# Northern Basin and Range ------------------------------------------------

## 26. Drill seeding ------------------------------------------------------

# Join cover & shannon cols
model26.matched2 <- model26.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model26.invasive) |> 
  left_join(filter(all_diversity, Model == "model26")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model26.matched2 <- model26.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model26.matched2 <- model26.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model26.diff <- model26.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model26.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model26.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model26.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values26 <- model26.perm |>
  inner_join(model26.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values26 # p = 0.009 for perennial grass

# Boxplot
model26.bp <- model26.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "26. Northern Basin and Range: Drill seeding") +
  geom_signif(
    y_position = 85,
    xmin = 3.8,
    xmax = 4.2, 
    annotations = c("**")
  ) +
  ylim(0, 90) + 
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model26.bp

# Plot frequency distribution
#   Annual forb
model26.annforb <- model26.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model26.diff$obs_diff[model26.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model26.annforb

#   Annual grass
model26.anngrass <- model26.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model26.diff$obs_diff[model26.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model26.anngrass

#   Perennial forb
model26.perforb <- model26.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model26.diff$obs_diff[model26.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model26.perforb

#   Perennial grass
model26.pergrass <- model26.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model26.diff$obs_diff[model26.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model26.pergrass

#   Shrub
model26.shrub <- model26.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model26.diff$obs_diff[model26.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model26.shrub

#   Shannon diversity
model26.shannon <- model26.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model26.diff$obs_diff[model26.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model26.shannon

#   Bromus tectorum
model26.brte <- model26.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model26.diff$obs_diff[model26.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model26.brte

# Combine plots
grid.arrange(
  model26.bp, model26.annforb, model26.anngrass,
  model26.perforb, model26.pergrass, model26.shrub,
  model26.shannon, model26.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 27. Drill seeding & soil disturbance -----------------------------------

# Join cover & shannon cols
model27.matched2 <- model27.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model27.invasive) |> 
  left_join(filter(all_diversity, Model == "model27")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model27.matched2 <- model27.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model27.matched2 <- model27.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model27.diff <- model27.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model27.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model27.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model27.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values27 <- model27.perm |>
  inner_join(model27.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values27 # p = 0.04 for annual forb; p = 0.02 for perennial forb; p = 0.02 for BRTE

# Boxplot
model27.bp <- model27.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "27. Northern Basin and Range: Drill seeding & soil disturbance") +
  geom_signif(
    y_position = 43,
    xmin = 0.8,
    xmax = 1.2, 
    annotations = c("*")
  ) +
  geom_signif(
    y_position = 25,
    xmin = 2.8,
    xmax = 3.2, 
    annotations = c("*")
  ) +
  geom_signif(
    y_position = 95,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("*")
  ) +
  ylim(0, 103) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model27.bp

# Plot frequency distribution
#   Annual forb
model27.annforb <- model27.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model27.diff$obs_diff[model27.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model27.annforb

#   Annual grass
model27.anngrass <- model27.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model27.diff$obs_diff[model27.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model27.anngrass

#   Perennial forb
model27.perforb <- model27.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model27.diff$obs_diff[model27.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model27.perforb

#   Perennial grass
model27.pergrass <- model27.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model27.diff$obs_diff[model27.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model27.pergrass

#   Shrub
model27.shrub <- model27.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model27.diff$obs_diff[model27.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model27.shrub

#   Shannon diversity
model27.shannon <- model27.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model27.diff$obs_diff[model27.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model27.shannon

#   Bromus tectorum
model27.brte <- model27.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model27.diff$obs_diff[model27.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model27.brte

# Combine plots
grid.arrange(
  model27.bp, model27.annforb, model27.anngrass,
  model27.perforb, model27.pergrass, model27.shrub,
  model27.shannon, model27.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 28. Herbicide ----------------------------------------------------------

# Join cover & shannon cols
model28.matched2 <- model28.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model28.invasive) |> 
  left_join(filter(all_diversity, Model == "model28")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model28.matched2 <- model28.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model28.matched2 <- model28.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model28.diff <- model28.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model28.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model28.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model28.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values28 <- model28.perm |>
  inner_join(model28.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values28

# Boxplot
model28.bp <- model28.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "28. Northern Basin and Range: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model28.bp

# Plot frequency distribution
#   Annual forb
model28.annforb <- model28.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model28.diff$obs_diff[model28.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model28.annforb

#   Annual grass
model28.anngrass <- model28.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model28.diff$obs_diff[model28.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model28.anngrass

#   Perennial forb
model28.perforb <- model28.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model28.diff$obs_diff[model28.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model28.perforb

#   Perennial grass
model28.pergrass <- model28.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model28.diff$obs_diff[model28.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model28.pergrass

#   Shrub
model28.shrub <- model28.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model28.diff$obs_diff[model28.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model28.shrub

#   Shannon diversity
model28.shannon <- model28.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model28.diff$obs_diff[model28.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model28.shannon

#   Bromus tectorum
model28.brte <- model28.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model28.diff$obs_diff[model28.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model28.brte

# Combine plots
grid.arrange(
  model28.bp, model28.annforb, model28.anngrass,
  model28.perforb, model28.pergrass, model28.shrub,
  model28.shannon, model28.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 29. Prescribed burn ----------------------------------------------------

# Join cover & shannon cols
model29.matched2 <- model29.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model29.invasive) |> 
  left_join(filter(all_diversity, Model == "model29")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model29.matched2 <- model29.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model29.matched2 <- model29.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model29.diff <- model29.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model29.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model29.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model29.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values29 <- model29.perm |>
  inner_join(model29.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values29 # p = 0.03 for annual grass; p = 0.002 for perennial grass;
#             p = 0.049 for shannon; p = 0.04 for BRTE

# Boxplot
model29.bp <- model29.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "29. Northern Basin and Range: Prescribed burn") +
  geom_signif(
    y_position = 98,
    xmin = 1.8,
    xmax = 2.2, 
    annotations = c("*")
  ) +
  geom_signif(
    y_position = 95,
    xmin = 3.8,
    xmax = 4.2, 
    annotations = c("**")
  ) +
  geom_signif(
    y_position = 98,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("*")
  ) +
  ylim(0, 105) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model29.bp

# Plot frequency distribution
#   Annual forb
model29.annforb <- model29.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model29.diff$obs_diff[model29.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model29.annforb

#   Annual grass
model29.anngrass <- model29.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model29.diff$obs_diff[model29.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model29.anngrass

#   Perennial forb
model29.perforb <- model29.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model29.diff$obs_diff[model29.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model29.perforb

#   Perennial grass
model29.pergrass <- model29.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model29.diff$obs_diff[model29.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model29.pergrass

#   Shrub
model29.shrub <- model29.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model29.diff$obs_diff[model29.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model29.shrub

#   Shannon diversity
model29.shannon <- model29.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model29.diff$obs_diff[model29.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model29.shannon

#   Bromus tectorum
model29.brte <- model29.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model29.diff$obs_diff[model29.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model29.brte

# Combine plots
grid.arrange(
  model29.bp, model29.annforb, model29.anngrass,
  model29.perforb, model29.pergrass, model29.shrub,
  model29.shannon, model29.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 30. Vegetation disturbance ---------------------------------------------

# Join cover & shannon cols
model30.matched2 <- model30.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model30.invasive) |> 
  left_join(filter(all_diversity, Model == "model30")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model30.matched2 <- model30.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model30.matched2 <- model30.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model30.diff <- model30.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model30.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model30.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model30.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values30 <- model30.perm |>
  inner_join(model30.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values30

# Boxplot
model30.bp <- model30.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "30. Northern BR: Vegetation disturbance") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model30.bp

# Plot frequency distribution
#   Annual forb
model30.annforb <- model30.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model30.diff$obs_diff[model30.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model30.annforb

#   Annual grass
model30.anngrass <- model30.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model30.diff$obs_diff[model30.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model30.anngrass

#   Perennial forb
model30.perforb <- model30.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model30.diff$obs_diff[model30.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model30.perforb

#   Perennial grass
model30.pergrass <- model30.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model30.diff$obs_diff[model30.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model30.pergrass

#   Shrub
model30.shrub <- model30.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model30.diff$obs_diff[model30.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model30.shrub

#   Shannon diversity
model30.shannon <- model30.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model30.diff$obs_diff[model30.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model30.shannon

#   Bromus tectorum
model30.brte <- model30.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model30.diff$obs_diff[model30.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model30.brte

# Combine plots
grid.arrange(
  model30.bp, model30.annforb, model30.anngrass,
  model30.perforb, model30.pergrass, model30.shrub,
  model30.shannon, model30.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 31. Post-burn aerial seeding -------------------------------------------

# Join cover & shannon cols
model31.matched2 <- model31.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model31.invasive) |> 
  left_join(filter(all_diversity, Model == "model31")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model31.matched2 <- model31.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model31.matched2 <- model31.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model31.diff <- model31.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model31.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model31.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model31.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values31 <- model31.perm |>
  inner_join(model31.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values31 # p = 0.008 for shannon

# Boxplot
model31.bp <- model31.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "31. Northern Basin and Range: Post-burn aerial seeding") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model31.bp

# Plot frequency distribution
#   Annual forb
model31.annforb <- model31.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model31.diff$obs_diff[model31.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model31.annforb

#   Annual grass
model31.anngrass <- model31.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model31.diff$obs_diff[model31.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model31.anngrass

#   Perennial forb
model31.perforb <- model31.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model31.diff$obs_diff[model31.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model31.perforb

#   Perennial grass
model31.pergrass <- model31.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model31.diff$obs_diff[model31.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model31.pergrass

#   Shrub
model31.shrub <- model31.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model31.diff$obs_diff[model31.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model31.shrub

#   Shannon diversity
model31.shannon <- model31.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model31.diff$obs_diff[model31.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model31.shannon

#   Bromus tectorum
model31.brte <- model31.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model31.diff$obs_diff[model31.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model31.brte

# Combine plots
grid.arrange(
  model31.bp, model31.annforb, model31.anngrass,
  model31.perforb, model31.pergrass, model31.shrub,
  model31.shannon, model31.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 32. Post-burn aerial & drill seeding -----------------------------------

# Join cover & shannon cols
model32.matched2 <- model32.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model32.invasive) |> 
  left_join(filter(all_diversity, Model == "model32")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model32.matched2 <- model32.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model32.matched2 <- model32.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model32.diff <- model32.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model32.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model32.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model32.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values32 <- model32.perm |>
  inner_join(model32.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values32

# Boxplot
model32.bp <- model32.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "32. Northern Basin and Range: Post-burn aerial & drill seeding") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model32.bp

# Plot frequency distribution
#   Annual forb
model32.annforb <- model32.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model32.diff$obs_diff[model32.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model32.annforb

#   Annual grass
model32.anngrass <- model32.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model32.diff$obs_diff[model32.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model32.anngrass

#   Perennial forb
model32.perforb <- model32.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model32.diff$obs_diff[model32.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model32.perforb

#   Perennial grass
model32.pergrass <- model32.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model32.diff$obs_diff[model32.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model32.pergrass

#   Shrub
model32.shrub <- model32.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model32.diff$obs_diff[model32.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model32.shrub

#   Shannon diversity
model32.shannon <- model32.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model32.diff$obs_diff[model32.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model32.shannon

#   Bromus tectorum
model32.brte <- model32.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model32.diff$obs_diff[model32.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model32.brte

# Combine plots
grid.arrange(
  model32.bp, model32.annforb, model32.anngrass,
  model32.perforb, model32.pergrass, model32.shrub,
  model32.shannon, model32.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 33. Post-burn closure --------------------------------------------------

# Join cover & shannon cols
model33.matched2 <- model33.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model33.invasive) |> 
  left_join(filter(all_diversity, Model == "model33")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model33.matched2 <- model33.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model33.matched2 <- model33.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model33.diff <- model33.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model33.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model33.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model33.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values33 <- model33.perm |>
  inner_join(model33.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values33 # p < 0.001 for annual and perennial grass; p = 0.001 for BRTE

# Boxplot
model33.bp <- model33.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "33. Northern Basin and Range: Post-burn closure") +
  geom_signif(
    y_position = 107,
    xmin = 1.8,
    xmax = 2.2, 
    annotations = c("***")
  ) +
  geom_signif(
    y_position = 105,
    xmin = 3.8,
    xmax = 4.2, 
    annotations = c("***")
  ) +
  geom_signif(
    y_position = 95,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("**")
  ) +
  ylim(0, 120) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model33.bp 

# Plot frequency distribution
#   Annual forb
model33.annforb <- model33.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model33.diff$obs_diff[model33.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model33.annforb

#   Annual grass
model33.anngrass <- model33.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model33.diff$obs_diff[model33.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model33.anngrass

#   Perennial forb
model33.perforb <- model33.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model33.diff$obs_diff[model33.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model33.perforb

#   Perennial grass
model33.pergrass <- model33.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model33.diff$obs_diff[model33.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model33.pergrass

#   Shrub
model33.shrub <- model33.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model33.diff$obs_diff[model33.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model33.shrub

#   Shannon diversity
model33.shannon <- model33.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model33.diff$obs_diff[model33.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model33.shannon

#   Bromus tectorum
model33.brte <- model33.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model33.diff$obs_diff[model33.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum") ~ "(**)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model33.brte

# Combine plots
grid.arrange(
  model33.bp, model33.annforb, model33.anngrass,
  model33.perforb, model33.pergrass, model33.shrub,
  model33.shannon, model33.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 34. Post-burn drill seeding --------------------------------------------------

# Join cover & shannon cols
model34.matched2 <- model34.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model34.invasive) |> 
  left_join(filter(all_diversity, Model == "model34")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model34.matched2 <- model34.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model34.matched2 <- model34.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model34.diff <- model34.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model34.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model34.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model34.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values34 <- model34.perm |>
  inner_join(model34.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values34 # p = 0.0098 for perennial forb; p = 0.001 for shannon

# Boxplot
model34.bp <- model34.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "34. Northern Basin and Range: Post-burn drill seeding") +
  geom_signif(
    y_position = 50,
    xmin = 2.8,
    xmax = 3.2, 
    annotations = c("**")
  ) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model34.bp

# Plot frequency distribution
#   Annual forb
model34.annforb <- model34.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model34.diff$obs_diff[model34.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model34.annforb

#   Annual grass
model34.anngrass <- model34.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model34.diff$obs_diff[model34.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model34.anngrass

#   Perennial forb
model34.perforb <- model34.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model34.diff$obs_diff[model34.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model34.perforb

#   Perennial grass
model34.pergrass <- model34.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model34.diff$obs_diff[model34.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model34.pergrass

#   Shrub
model34.shrub <- model34.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model34.diff$obs_diff[model34.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model34.shrub

#   Shannon diversity
model34.shannon <- model34.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model34.diff$obs_diff[model34.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model34.shannon

#   Bromus tectorum
model34.brte <- model34.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model34.diff$obs_diff[model34.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model34.brte

# Combine plots
grid.arrange(
  model34.bp, model34.annforb, model34.anngrass,
  model34.perforb, model34.pergrass, model34.shrub,
  model34.shannon, model34.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 35. Post-burn herbicide ------------------------------------------------

# Join cover & shannon cols
model35.matched2 <- model35.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model35.invasive) |> 
  left_join(filter(all_diversity, Model == "model35")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model35.matched2 <- model35.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model35.matched2 <- model35.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model35.diff <- model35.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model35.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model35.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model35.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values35 <- model35.perm |>
  inner_join(model35.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values35 # p = 0.002 for perennial grass; p = 0.005 for shannon; p = 0.002 for BRTE

# Boxplot
model35.bp <- model35.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "35. Northern Basin and Range: Post-burn herbicide") +
  geom_signif(
    y_position = 95,
    xmin = 3.8,
    xmax = 4.2, 
    annotations = c("**")
  ) +
  geom_signif(
    y_position = 95,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("**")
  ) +
  ylim(0, 103) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model35.bp

# Plot frequency distribution
#   Annual forb
model35.annforb <- model35.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model35.diff$obs_diff[model35.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model35.annforb

#   Annual grass
model35.anngrass <- model35.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model35.diff$obs_diff[model35.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model35.anngrass

#   Perennial forb
model35.perforb <- model35.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model35.diff$obs_diff[model35.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model35.perforb

#   Perennial grass
model35.pergrass <- model35.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model35.diff$obs_diff[model35.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model35.pergrass

#   Shrub
model35.shrub <- model35.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model35.diff$obs_diff[model35.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model35.shrub

#   Shannon diversity
model35.shannon <- model35.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model35.diff$obs_diff[model35.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model35.shannon

#   Bromus tectorum
model35.brte <- model35.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model35.diff$obs_diff[model35.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum") ~ "(**)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model35.brte

# Combine plots
grid.arrange(
  model35.bp, model35.annforb, model35.anngrass,
  model35.perforb, model35.pergrass, model35.shrub,
  model35.shannon, model35.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 36. Post-burn seedling planting ----------------------------------------

# Join cover & shannon cols
model36.matched2 <- model36.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model36.invasive) |> 
  left_join(filter(all_diversity, Model == "model36")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model36.matched2 <- model36.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model36.matched2 <- model36.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model36.diff <- model36.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model36.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model36.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model36.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values36 <- model36.perm |>
  inner_join(model36.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values36 # p = 0.02 for annual grass; p = 0.02 for shannon; p = 0.006 for BRTE

# Boxplot
model36.bp <- model36.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "36. Northern Basin and Range: Post-burn seedling planting") +
  geom_signif(
    y_position = 105,
    xmin = 1.8,
    xmax = 2.2, 
    annotations = c("*")
  ) +
  geom_signif(
    y_position = 92,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("**")
  ) +
  ylim(0, 113) + 
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model36.bp

# Plot frequency distribution
#   Annual forb
model36.annforb <- model36.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model36.diff$obs_diff[model36.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model36.annforb

#   Annual grass
model36.anngrass <- model36.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model36.diff$obs_diff[model36.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model36.anngrass

#   Perennial forb
model36.perforb <- model36.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model36.diff$obs_diff[model36.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model36.perforb

#   Perennial grass
model36.pergrass <- model36.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model36.diff$obs_diff[model36.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model36.pergrass

#   Shrub
model36.shrub <- model36.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model36.diff$obs_diff[model36.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model36.shrub

#   Shannon diversity
model36.shannon <- model36.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model36.diff$obs_diff[model36.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model36.shannon

#   Bromus tectorum
model36.brte <- model36.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model36.diff$obs_diff[model36.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum") ~ "(**)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model36.brte

# Combine plots
grid.arrange(
  model36.bp, model36.annforb, model36.anngrass,
  model36.perforb, model36.pergrass, model36.shrub,
  model36.shannon, model36.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)



# Northwestern Great Plains -----------------------------------------------

## 37. Prescribed burn ----------------------------------------------------

# Join cover & shannon cols
model37.matched2 <- model37.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model37.invasive) |> 
  left_join(filter(all_diversity, Model == "model37")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model37.matched2 <- model37.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model37.matched2 <- model37.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model37.diff <- model37.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model37.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model37.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model37.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values37 <- model37.perm |>
  inner_join(model37.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values37 # p = 0.02 for shannon

# Boxplot
model37.bp <- model37.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "37. Northwestern Great Plains: Prescribed burn") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model37.bp

# Plot frequency distribution
#   Annual forb
model37.annforb <- model37.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model37.diff$obs_diff[model37.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model37.annforb

#   Annual grass
model37.anngrass <- model37.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model37.diff$obs_diff[model37.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model37.anngrass

#   Perennial forb
model37.perforb <- model37.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model37.diff$obs_diff[model37.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model37.perforb

#   Perennial grass
model37.pergrass <- model37.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model37.diff$obs_diff[model37.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model37.pergrass

#   Shrub
model37.shrub <- model37.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model37.diff$obs_diff[model37.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model37.shrub

#   Shannon diversity
model37.shannon <- model37.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model37.diff$obs_diff[model37.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model37.shannon

#   Bromus tectorum
model37.brte <- model37.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model37.diff$obs_diff[model37.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model37.brte

# Combine plots
grid.arrange(
  model37.bp, model37.annforb, model37.anngrass,
  model37.perforb, model37.pergrass, model37.shrub,
  model37.shannon, model37.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)



# Snake River Plain -------------------------------------------------------

## 38. Post-burn aerial seeding -------------------------------------------

# Join cover & shannon cols
model38.matched2 <- model38.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model38.invasive) |> 
  left_join(filter(all_diversity, Model == "model38")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model38.matched2 <- model38.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model38.matched2 <- model38.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model38.diff <- model38.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model38.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model38.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model38.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values38 <- model38.perm |>
  inner_join(model38.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values38 # p = 0.01 for perennial forb; p = 0.003 for perennial grass
#             p < 0.001 for shannon

# Boxplot
model38.bp <- model38.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "38. Snake River Plain: Post-burn aerial seeding") +
  geom_signif(
    y_position = 50,
    xmin = 2.8,
    xmax = 3.2, 
    annotations = c("*")
  ) +
  geom_signif(
    y_position = 95,
    xmin = 3.8,
    xmax = 4.2, 
    annotations = c("**")
  ) +
  ylim(0, 103) + 
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model38.bp

# Plot frequency distribution
#   Annual forb
model38.annforb <- model38.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model38.diff$obs_diff[model38.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model38.annforb

#   Annual grass
model38.anngrass <- model38.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model38.diff$obs_diff[model38.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model38.anngrass

#   Perennial forb
model38.perforb <- model38.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model38.diff$obs_diff[model38.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model38.perforb

#   Perennial grass
model38.pergrass <- model38.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model38.diff$obs_diff[model38.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model38.pergrass

#   Shrub
model38.shrub <- model38.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model38.diff$obs_diff[model38.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model38.shrub

#   Shannon diversity
model38.shannon <- model38.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model38.diff$obs_diff[model38.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model38.shannon

#   Bromus tectorum
model38.brte <- model38.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model38.diff$obs_diff[model38.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model38.brte

# Combine plots
grid.arrange(
  model38.bp, model38.annforb, model38.anngrass,
  model38.perforb, model38.pergrass, model38.shrub,
  model38.shannon, model38.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 39. Post-burn aerial & drill seeding -----------------------------------

# Join cover & shannon cols
model39.matched2 <- model39.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model39.invasive) |> 
  left_join(filter(all_diversity, Model == "model39")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model39.matched2 <- model39.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model39.matched2 <- model39.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model39.diff <- model39.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model39.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model39.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model39.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values39 <- model39.perm |>
  inner_join(model39.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values39 # p = 0.04 from BRTE

# Boxplot
model39.bp <- model39.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "39. Snake River Plain: Post-burn aerial & drill seeding") +
  geom_signif(
    y_position = 100,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("*")
  ) +
  ylim(0, 108) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model39.bp

# Plot frequency distribution
#   Annual forb
model39.annforb <- model39.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model39.diff$obs_diff[model39.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model39.annforb

#   Annual grass
model39.anngrass <- model39.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model39.diff$obs_diff[model39.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model39.anngrass

#   Perennial forb
model39.perforb <- model39.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model39.diff$obs_diff[model39.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model39.perforb

#   Perennial grass
model39.pergrass <- model39.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model39.diff$obs_diff[model39.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model39.pergrass

#   Shrub
model39.shrub <- model39.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model39.diff$obs_diff[model39.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model39.shrub

#   Shannon diversity
model39.shannon <- model39.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model39.diff$obs_diff[model39.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model39.shannon

#   Bromus tectorum
model39.brte <- model39.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model39.diff$obs_diff[model39.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model39.brte

# Combine plots
grid.arrange(
  model39.bp, model39.annforb, model39.anngrass,
  model39.perforb, model39.pergrass, model39.shrub,
  model39.shannon, model39.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 40. Post-burn closure --------------------------------------------------

# Join cover & shannon cols
model40.matched2 <- model40.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model40.invasive) |> 
  left_join(filter(all_diversity, Model == "model40")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model40.matched2 <- model40.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model40.matched2 <- model40.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model40.diff <- model40.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model40.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model40.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model40.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values40 <- model40.perm |>
  inner_join(model40.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values40 # p = 0.006 for annual grass; p = 0.0099 for perennial grass;
#             p < 0.001 for shannon; p = 0.2 for BRTE

# Boxplot
model40.bp <- model40.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "40. Snake River Plain: Post-burn closure") +
  geom_signif(
    y_position = 105,
    xmin = 1.8,
    xmax = 2.2, 
    annotations = c("**")
  ) +
  geom_signif(
    y_position = 90,
    xmin = 3.8,
    xmax = 4.2, 
    annotations = c("**")
  ) +
  geom_signif(
    y_position = 90,
    xmin = 5.8,
    xmax = 6.2, 
    annotations = c("*")
  ) +
  ylim(0, 112) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model40.bp

# Plot frequency distribution
#   Annual forb
model40.annforb <- model40.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model40.diff$obs_diff[model40.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model40.annforb

#   Annual grass
model40.anngrass <- model40.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model40.diff$obs_diff[model40.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model40.anngrass

#   Perennial forb
model40.perforb <- model40.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model40.diff$obs_diff[model40.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model40.perforb

#   Perennial grass
model40.pergrass <- model40.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model40.diff$obs_diff[model40.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass (**)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model40.pergrass

#   Shrub
model40.shrub <- model40.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model40.diff$obs_diff[model40.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model40.shrub

#   Shannon diversity
model40.shannon <- model40.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model40.diff$obs_diff[model40.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (***)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model40.shannon

#   Bromus tectorum
model40.brte <- model40.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model40.diff$obs_diff[model40.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum") ~ "(*)")) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model40.brte

# Combine plots
grid.arrange(
  model40.bp, model40.annforb, model40.anngrass,
  model40.perforb, model40.pergrass, model40.shrub,
  model40.shannon, model40.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 41. Post-burn drill seeding --------------------------------------------------

# Join cover & shannon cols
model41.matched2 <- model41.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model41.invasive) |> 
  left_join(filter(all_diversity, Model == "model41")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model41.matched2 <- model41.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model41.matched2 <- model41.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model41.diff <- model41.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model41.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model41.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model41.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values41 <- model41.perm |>
  inner_join(model41.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values41 # p = 0.02 for shannon

# Boxplot
model41.bp <- model41.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "41. Snake River Plain: Post-burn drill seeding") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model41.bp

# Plot frequency distribution
#   Annual forb
model41.annforb <- model41.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model41.diff$obs_diff[model41.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model41.annforb

#   Annual grass
model41.anngrass <- model41.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model41.diff$obs_diff[model41.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model41.anngrass

#   Perennial forb
model41.perforb <- model41.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model41.diff$obs_diff[model41.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model41.perforb

#   Perennial grass
model41.pergrass <- model41.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model41.diff$obs_diff[model41.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model41.pergrass

#   Shrub
model41.shrub <- model41.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model41.diff$obs_diff[model41.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model41.shrub

#   Shannon diversity
model41.shannon <- model41.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model41.diff$obs_diff[model41.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity (*)") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model41.shannon

#   Bromus tectorum
model41.brte <- model41.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model41.diff$obs_diff[model41.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model41.brte

# Combine plots
grid.arrange(
  model41.bp, model41.annforb, model41.anngrass,
  model41.perforb, model41.pergrass, model41.shrub,
  model41.shannon, model41.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)


## 42. Post-burn herbicide ------------------------------------------------

# Join cover & shannon cols
model42.matched2 <- model42.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model42.invasive) |> 
  left_join(filter(all_diversity, Model == "model42")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model42.matched2 <- model42.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model42.matched2 <- model42.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model42.diff <- model42.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model42.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model42.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model42.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values42 <- model42.perm |>
  inner_join(model42.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values42

# Boxplot
model42.bp <- model42.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "42. Snake River Plain: Post-burn herbicide") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model42.bp

# Plot frequency distribution
#   Annual forb
model42.annforb <- model42.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model42.diff$obs_diff[model42.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model42.annforb

#   Annual grass
model42.anngrass <- model42.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model42.diff$obs_diff[model42.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model42.anngrass

#   Perennial forb
model42.perforb <- model42.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model42.diff$obs_diff[model42.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model42.perforb

#   Perennial grass
model42.pergrass <- model42.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model42.diff$obs_diff[model42.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model42.pergrass

#   Shrub
model42.shrub <- model42.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model42.diff$obs_diff[model42.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model42.shrub

#   Shannon diversity
model42.shannon <- model42.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model42.diff$obs_diff[model42.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model42.shannon

#   Bromus tectorum
model42.brte <- model42.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model42.diff$obs_diff[model42.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model42.brte

# Combine plots
grid.arrange(
  model42.bp, model42.annforb, model42.anngrass,
  model42.perforb, model42.pergrass, model42.shrub,
  model42.shannon, model42.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)



# Southern Rockies --------------------------------------------------------

## 43. Herbicide ----------------------------------------------------------

# Join cover & shannon cols
model43.matched2 <- model43.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join)

#   pivot_longer() for cover & shannon cols
model43.matched2 <- model43.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model43.matched2 <- model43.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub"
           ))

# Calculate observed mean difference
model43.diff <- model43.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model43.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model43.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model43.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values43 <- model43.perm |>
  inner_join(model43.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values43

# Boxplot
model43.bp <- model43.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "43. Southern Rockies: Herbicide") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model43.bp

# Plot frequency distribution
#   Annual forb
model43.annforb <- model43.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model43.diff$obs_diff[model43.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model43.annforb

#   Annual grass
model43.anngrass <- model43.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model43.diff$obs_diff[model43.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model43.anngrass

#   Perennial forb
model43.perforb <- model43.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model43.diff$obs_diff[model43.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model43.perforb

#   Perennial grass
model43.pergrass <- model43.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model43.diff$obs_diff[model43.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model43.pergrass

#   Shrub
model43.shrub <- model43.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model43.diff$obs_diff[model43.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model43.shrub

# Combine plots
grid.arrange(
  model43.bp, model43.annforb, model43.anngrass,
  model43.perforb, model43.pergrass, model43.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)


## 44. Prescribed burn ----------------------------------------------------

# Join cover & shannon cols
model44.matched2 <- model44.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join)

#   pivot_longer() for cover & shannon cols
model44.matched2 <- model44.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model44.matched2 <- model44.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub"
           ))

# Calculate observed mean difference
model44.diff <- model44.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model44.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model44.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model44.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values44 <- model44.perm |>
  inner_join(model44.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values44

# Boxplot
model44.bp <- model44.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "44. Southern Rockies: Prescribed burn") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model44.bp

# Plot frequency distribution
#   Annual forb
model44.annforb <- model44.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model44.diff$obs_diff[model44.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model44.annforb

#   Annual grass
model44.anngrass <- model44.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model44.diff$obs_diff[model44.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model44.anngrass

#   Perennial forb
model44.perforb <- model44.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model44.diff$obs_diff[model44.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model44.perforb

#   Perennial grass
model44.pergrass <- model44.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model44.diff$obs_diff[model44.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model44.pergrass

#   Shrub
model44.shrub <- model44.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model44.diff$obs_diff[model44.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model44.shrub

# Combine plots
grid.arrange(
  model44.bp, model44.annforb, model44.anngrass,
  model44.perforb, model44.pergrass, model44.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)


## 45. Vegetation disturbance ---------------------------------------------

# Join cover & shannon cols
model45.matched2 <- model45.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join)

#   pivot_longer() for cover & shannon cols
model45.matched2 <- model45.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model45.matched2 <- model45.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub"
           ))

# Calculate observed mean difference
model45.diff <- model45.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model45.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model45.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model45.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values45 <- model45.perm |>
  inner_join(model45.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values45

# Boxplot
model45.bp <- model45.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "45. Southern Rockies: Vegetation disturbance") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model45.bp

# Plot frequency distribution
#   Annual forb
model45.annforb <- model45.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model45.diff$obs_diff[model45.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model45.annforb

#   Annual grass
model45.anngrass <- model45.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model45.diff$obs_diff[model45.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model45.anngrass

#   Perennial forb
model45.perforb <- model45.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model45.diff$obs_diff[model45.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model45.perforb

#   Perennial grass
model45.pergrass <- model45.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model45.diff$obs_diff[model45.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model45.pergrass

#   Shrub
model45.shrub <- model45.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model45.diff$obs_diff[model45.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model45.shrub

# Combine plots
grid.arrange(
  model45.bp, model45.annforb, model45.anngrass,
  model45.perforb, model45.pergrass, model45.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)



# Wyoming Basin -----------------------------------------------------------

## 46. Prescribed burn ----------------------------------------------------

# Join cover & shannon cols
model46.matched2 <- model46.matched |> 
  select(LDCpointID, PrimaryKey, trt_control) |> 
  left_join(geoindicators.join) |> 
  left_join(model46.invasive) |> 
  left_join(filter(all_diversity, Model == "model46")) |> 
  select(-Model, -EcoLvl3, -Species, -ScientificName, -Species_AH_n)

#   pivot_longer() for cover & shannon cols
model46.matched2 <- model46.matched2 |> 
  pivot_longer(
    cols = !c(LDCpointID, PrimaryKey, trt_control),
    names_to = "indicators",
    values_to = "value"
  )

#   Rename functional group cover types
model46.matched2 <- model46.matched2 |> 
  mutate(indicators = 
           case_when(
             indicators == "AnnForbCover_AH" ~ "Annual forb",
             indicators == "AnnGramCover_AH" ~ "Annual grass",
             indicators == "PerForbCover_AH" ~ "Perennial forb",
             indicators == "PerGramCover_AH" ~ "Perennial grass",
             indicators == "ShrubCover_AH" ~ "Shrub",
             indicators == "Shannon" ~ "Shannon diversity",
             indicators == "SpeciesCover_AH" ~ "Bromus tectorum"
           )) |> 
  mutate(indicators = factor(indicators,
                             levels = c("Annual forb",
                                        "Annual grass",
                                        "Perennial forb",
                                        "Perennial grass",
                                        "Shrub",
                                        "Shannon diversity",
                                        "Bromus tectorum")))

# Calculate observed mean difference
model46.diff <- model46.matched2 |> 
  group_by(indicators, trt_control) |> 
  summarise(mean_cover = mean(value),
            .groups = "drop") |> 
  group_by(indicators) |> 
  summarise(obs_diff = diff(mean_cover))
model46.diff

# Permutation test
n_perms <- 10000

set.seed(1)

model46.perm <- map_dfr(
  1:n_perms,
  ~ {
    # shuffle treatment labels
    permuted_data <- model46.matched2 |> 
      mutate(trt_control = sample(trt_control))
    
    # calculate mean differences for each functional group
    permuted_data |>
      group_by(indicators, trt_control) |>
      summarize(mean_cover = mean(value), .groups = "drop") |>
      group_by(indicators) |>
      summarize(mean_diff = diff(mean_cover), .groups = "drop") |>
      mutate(Iteration = .x)
  }
)

#   Calculate p-values for each variable
p_values46 <- model46.perm |>
  inner_join(model46.diff, by = "indicators") |>
  group_by(indicators) |>
  summarize(p_value = mean(abs(mean_diff) >= abs(obs_diff[1])))
p_values46

# Boxplot
model46.bp <- model46.matched2 |> 
  filter(indicators != "Shannon diversity") |> 
  ggplot(aes(x = indicators, y = value, fill = trt_control)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_bw() +
  labs(y = "Cover (%)",
       x = NULL,
       title = "46. Wyoming Basin: Prescribed burn") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) +
  scale_x_discrete(
    labels = c("Bromus tectorum" = expression(italic("Bromus tectorum")))) +
  theme(plot.margin = margin(10, 10, 20, 10))
model46.bp

# Plot frequency distribution
#   Annual forb
model46.annforb <- model46.perm |> 
  filter(indicators == "Annual forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model46.diff$obs_diff[model46.diff$indicators == "Annual forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model46.annforb

#   Annual grass
model46.anngrass <- model46.perm |> 
  filter(indicators == "Annual grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model46.diff$obs_diff[model46.diff$indicators == "Annual grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Annual grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model46.anngrass

#   Perennial forb
model46.perforb <- model46.perm |> 
  filter(indicators == "Perennial forb") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model46.diff$obs_diff[model46.diff$indicators == "Perennial forb"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial forb") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model46.perforb

#   Perennial grass
model46.pergrass <- model46.perm |> 
  filter(indicators == "Perennial grass") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model46.diff$obs_diff[model46.diff$indicators == "Perennial grass"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Perennial grass") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model46.pergrass

#   Shrub
model46.shrub <- model46.perm |> 
  filter(indicators == "Shrub") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model46.diff$obs_diff[model46.diff$indicators == "Shrub"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shrub") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model46.shrub

#   Shannon diversity
model46.shannon <- model46.perm |> 
  filter(indicators == "Shannon diversity") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model46.diff$obs_diff[model46.diff$indicators == "Shannon diversity"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = "Shannon diversity") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model46.shannon

#   Bromus tectorum
model46.brte <- model46.perm |> 
  filter(indicators == "Bromus tectorum") |> 
  ggplot(aes(x = mean_diff)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  geom_vline(xintercept = model46.diff$obs_diff[model46.diff$indicators == "Bromus tectorum"],
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Difference in means",
       y = "Frequency",
       title = expression(italic("Bromus tectorum"))) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 10, 10, 10))
model46.brte

# Combine plots
grid.arrange(
  model46.bp, model46.annforb, model46.anngrass,
  model46.perforb, model46.pergrass, model46.shrub,
  model46.shannon, model46.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)



# Write out figures -------------------------------------------------------

## AZ/NM Mountains --------------------------------------------------------

# 1. AZ/NM Mountains: Prescribed burn
tiff("figures/2026-06_permutation-tests-2/model01_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model01.bp, model01.annforb, model01.anngrass,
  model01.perforb, model01.pergrass, model01.shrub,
  model01.shannon, model01.brru2,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()


## AZ/NM Plateau ----------------------------------------------------------

# 2. AZ/NM Plateau: Herbicide
tiff("figures/2026-06_permutation-tests-2/model02_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model02.bp, model02.annforb, model02.anngrass,
  model02.perforb, model02.pergrass, model02.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)
dev.off()

# 3. AZ/NM Plateau: Prescribed burn
tiff("figures/2026-06_permutation-tests-2/model03_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model03.bp, model03.annforb, model03.anngrass,
  model03.perforb, model03.pergrass, model03.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)
dev.off()

# 4. AZ/NM Plateau: Seeding
tiff("figures/2026-06_permutation-tests-2/model04_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model04.bp, model04.annforb, model04.anngrass,
  model04.perforb, model04.pergrass, model04.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)
dev.off()

# 5. AZ/NM Plateau: Soil disturbance
tiff("figures/2026-06_permutation-tests-2/model05_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model05.bp, model05.annforb, model05.anngrass,
  model05.perforb, model05.pergrass, model05.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)
dev.off()


## Blue Mountains ---------------------------------------------------------

# 6. Blue Mountains: Herbicide
tiff("figures/2026-06_permutation-tests-2/model06_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model06.bp, model06.annforb, model06.anngrass,
  model06.perforb, model06.pergrass, model06.shrub,
  model06.shannon, model06.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 7. Blue Mountains: Vegetation disturbance
tiff("figures/2026-06_permutation-tests-2/model07_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model07.bp, model07.annforb, model07.anngrass,
  model07.perforb, model07.pergrass, model07.shrub,
  model07.shannon, model07.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 8. Blue Mountains: Post-burn herbicide
tiff("figures/2026-06_permutation-tests-2/model08_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model08.bp, model08.annforb, model08.anngrass,
  model08.perforb, model08.pergrass, model08.shrub,
  model08.shannon, model08.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()


## Central Basin and Range ------------------------------------------------

# 9. Central BR: Aerial seeding
tiff("figures/2026-06_permutation-tests-2/model09_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model09.bp, model09.annforb, model09.anngrass,
  model09.perforb, model09.pergrass, model09.shrub,
  model09.shannon, model09.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 10. Central BR: Drill seeding & soil disturbance
tiff("figures/2026-06_permutation-tests-2/model10_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model10.bp, model10.annforb, model10.anngrass,
  model10.perforb, model10.pergrass, model10.shrub,
  model10.shannon, model10.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 11. Central BR: Prescribed burn
tiff("figures/2026-06_permutation-tests-2/model11_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model11.bp, model11.annforb, model11.anngrass,
  model11.perforb, model11.pergrass, model11.shrub,
  model11.shannon, model11.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)

dev.off()

# 12. Central BR: Vegetation disturbance
tiff("figures/2026-06_permutation-tests-2/model12_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model12.bp, model12.annforb, model12.anngrass,
  model12.perforb, model12.pergrass, model12.shrub,
  model12.shannon, model12.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 13. Central BR: Post-burn aerial seeding
tiff("figures/2026-06_permutation-tests-2/model13_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model13.bp, model13.annforb, model13.anngrass,
  model13.perforb, model13.pergrass, model13.shrub,
  model13.shannon, model13.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 14. Central BR: Post-burn drill seeding
tiff("figures/2026-06_permutation-tests-2/model14_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model14.bp, model14.annforb, model14.anngrass,
  model14.perforb, model14.pergrass, model14.shrub,
  model14.shannon, model14.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 15. Central BR: Post-burn ground seeding
tiff("figures/2026-06_permutation-tests-2/model15_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model15.bp, model15.annforb, model15.anngrass,
  model15.perforb, model15.pergrass, model15.shrub,
  model15.shannon, model15.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 16. Central BR: Post-burn herbicide
tiff("figures/2026-06_permutation-tests-2/model16_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model16.bp, model16.annforb, model16.anngrass,
  model16.perforb, model16.pergrass, model16.shrub,
  model16.shannon, model16.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()


## Chihuahuan Desert ------------------------------------------------------

# 17. Chihuahuan Desert: Herbicide
tiff("figures/2026-06_permutation-tests-2/model17_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model17.bp, model17.annforb, model17.anngrass,
  model17.perforb, model17.pergrass, model17.shrub,
  model17.shannon, model17.prgl2,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()


## Colorado Plateaus ------------------------------------------------------

# 18. CO Plateaus: Aerial seeding & soil disturbance
tiff("figures/2026-06_permutation-tests-2/model18_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model18.bp, model18.annforb, model18.anngrass,
  model18.perforb, model18.pergrass, model18.shrub,
  model18.shannon, model18.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 19. CO Plateaus: Herbicide
tiff("figures/2026-06_permutation-tests-2/model19_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model19.bp, model19.annforb, model19.anngrass,
  model19.perforb, model19.pergrass, model19.shrub,
  model19.shannon, model19.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 20. CO Plateaus: Prescribed burn
tiff("figures/2026-06_permutation-tests-2/model20_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model20.bp, model20.annforb, model20.anngrass,
  model20.perforb, model20.pergrass, model20.shrub,
  model20.shannon, model20.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 21. CO Plateaus: Soil disturbance
tiff("figures/2026-06_permutation-tests-2/model21_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model21.bp, model21.annforb, model21.anngrass,
  model21.perforb, model21.pergrass, model21.shrub,
  model21.shannon, model21.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 22. CO Plateaus: Vegetation disturbance
tiff("figures/2026-06_permutation-tests-2/model22_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model22.bp, model22.annforb, model22.anngrass,
  model22.perforb, model22.pergrass, model22.shrub,
  model22.shannon, model22.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 23. CO Plateaus: Post-burn aerial seeding
tiff("figures/2026-06_permutation-tests-2/model23_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model23.bp, model23.annforb, model23.anngrass,
  model23.perforb, model23.pergrass, model23.shrub,
  model23.shannon, model23.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()


## Middle Rockies ---------------------------------------------------------

# 24. Middle Rockies: Herbicide
tiff("figures/2026-06_permutation-tests-2/model24_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model24.bp, model24.annforb, model24.anngrass,
  model24.perforb, model24.pergrass, model24.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)
dev.off()


## Mojave Basin and Range -------------------------------------------------

# 25. Mojave BR: Post-burn aerial seeding
tiff("figures/2026-06_permutation-tests-2/model25_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model25.bp, model25.annforb, model25.anngrass,
  model25.perforb, model25.pergrass, model25.shrub,
  model25.shannon, model25.brru2,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()


## Northern Basin and Range -----------------------------------------------

# 26. Northern BR: Drill seeding
tiff("figures/2026-06_permutation-tests-2/model26_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model26.bp, model26.annforb, model26.anngrass,
  model26.perforb, model26.pergrass, model26.shrub,
  model26.shannon, model26.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 27. Northern BR: Drill seeding & soil disturbance
tiff("figures/2026-06_permutation-tests-2/model27_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model27.bp, model27.annforb, model27.anngrass,
  model27.perforb, model27.pergrass, model27.shrub,
  model27.shannon, model27.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 28. Northern BR: Herbicide
tiff("figures/2026-06_permutation-tests-2/model28_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model28.bp, model28.annforb, model28.anngrass,
  model28.perforb, model28.pergrass, model28.shrub,
  model28.shannon, model28.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 29. Northern BR: Prescribed burn
tiff("figures/2026-06_permutation-tests-2/model29_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model29.bp, model29.annforb, model29.anngrass,
  model29.perforb, model29.pergrass, model29.shrub,
  model29.shannon, model29.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 30. Northern BR: Vegetation disturbance
tiff("figures/2026-06_permutation-tests-2/model30_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model30.bp, model30.annforb, model30.anngrass,
  model30.perforb, model30.pergrass, model30.shrub,
  model30.shannon, model30.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 31. Northern BR: Post-burn aerial seeding
tiff("figures/2026-06_permutation-tests-2/model31_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model31.bp, model31.annforb, model31.anngrass,
  model31.perforb, model31.pergrass, model31.shrub,
  model31.shannon, model31.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 32. Northern BR: Post-burn aerial and drill seeding
tiff("figures/2026-06_permutation-tests-2/model32_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model32.bp, model32.annforb, model32.anngrass,
  model32.perforb, model32.pergrass, model32.shrub,
  model32.shannon, model32.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 33. Northern BR: Post-burn closure
tiff("figures/2026-06_permutation-tests-2/model33_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model33.bp, model33.annforb, model33.anngrass,
  model33.perforb, model33.pergrass, model33.shrub,
  model33.shannon, model33.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 34. Northern BR: Post-burn drill seeding
tiff("figures/2026-06_permutation-tests-2/model34_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model34.bp, model34.annforb, model34.anngrass,
  model34.perforb, model34.pergrass, model34.shrub,
  model34.shannon, model34.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 35. Northern BR: Post-burn herbicide
tiff("figures/2026-06_permutation-tests-2/model35_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model35.bp, model35.annforb, model35.anngrass,
  model35.perforb, model35.pergrass, model35.shrub,
  model35.shannon, model35.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 36. Northern BR: Post-burn seedling planting
tiff("figures/2026-06_permutation-tests-2/model36_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model36.bp, model36.annforb, model36.anngrass,
  model36.perforb, model36.pergrass, model36.shrub,
  model36.shannon, model36.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()


## Northwestern Great Plains ----------------------------------------------

# 37. NW Great Plains: Prescribed burn
tiff("figures/2026-06_permutation-tests-2/model37_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model37.bp, model37.annforb, model37.anngrass,
  model37.perforb, model37.pergrass, model37.shrub,
  model37.shannon, model37.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()


## Snake River Plain ------------------------------------------------------

# 38. Snake River Plain: Post-burn aerial seeding
tiff("figures/2026-06_permutation-tests-2/model38_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model38.bp, model38.annforb, model38.anngrass,
  model38.perforb, model38.pergrass, model38.shrub,
  model38.shannon, model38.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 39. Snake River Plain: Post-burn aerial & drill seeding
tiff("figures/2026-06_permutation-tests-2/model39_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model39.bp, model39.annforb, model39.anngrass,
  model39.perforb, model39.pergrass, model39.shrub,
  model39.shannon, model39.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 40. Snake River Plain: Post-burn closure
tiff("figures/2026-06_permutation-tests-2/model40_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model40.bp, model40.annforb, model40.anngrass,
  model40.perforb, model40.pergrass, model40.shrub,
  model40.shannon, model40.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 41. Snake River Plain: Post-burn drill seeding
tiff("figures/2026-06_permutation-tests-2/model41_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model41.bp, model41.annforb, model41.anngrass,
  model41.perforb, model41.pergrass, model41.shrub,
  model41.shannon, model41.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()

# 42. Snake River Plain: Post-burn herbicide
tiff("figures/2026-06_permutation-tests-2/model42_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model42.bp, model42.annforb, model42.anngrass,
  model42.perforb, model42.pergrass, model42.shrub,
  model42.shannon, model42.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()


## Southern Rockies -------------------------------------------------------

# 43. Southern Rockies: Herbicide
tiff("figures/2026-06_permutation-tests-2/model43_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model43.bp, model43.annforb, model43.anngrass,
  model43.perforb, model43.pergrass, model43.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)
dev.off()

# 44. Southern Rockies: Prescribed burn
tiff("figures/2026-06_permutation-tests-2/model44_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model44.bp, model44.annforb, model44.anngrass,
  model44.perforb, model44.pergrass, model44.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)
dev.off()

# 45. Southern Rockies: Vegetation disturbance
tiff("figures/2026-06_permutation-tests-2/model45_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model45.bp, model45.annforb, model45.anngrass,
  model45.perforb, model45.pergrass, model45.shrub,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6)
  )
)
dev.off()


## Wyoming Basin ----------------------------------------------------------

# 46. Wyoming Basin: Prescribed burn
tiff("figures/2026-06_permutation-tests-2/model46_permutation-2.tiff",
     units = "in", width = 10, height = 9, res = 150)
grid.arrange(
  model46.bp, model46.annforb, model46.anngrass,
  model46.perforb, model46.pergrass, model46.shrub,
  model46.shannon, model46.brte,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 3, 3, NA),
    c(4, 4, 5, 5, 6, 6),
    c(NA, 7, 7, 8, 8, NA)
  )
)
dev.off()


save.image("RData/19_permutation-tests-2.RData")
