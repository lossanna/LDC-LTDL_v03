# Created: 2026-06-02
# Updated: 2026-06-02

# Purpose: Identify prominent invasive species by ecoregion. Create dataframes of matched
#   data with invasive species cover for each model.


library(tidyverse)

# Load data ---------------------------------------------------------------

geospecies.raw <- read_csv("data/raw/downloaded/ldc-data-2026-06-01/geospecies.csv")
ldc.007.raw <- read_csv("data/versions-from-R/12.3_LDC-points_v007.csv")
load("RData/13_matched-data.RData")


## Prepare geospecies -----------------------------------------------------

# Rename columns
col_rename_map <- c(
  "Project Key" = "ProjectKey",
  "Primary Key" = "PrimaryKey",
  "Date Visited" = "DateVisited",
  "Species" = "Species",
  "Scientific Name" = "ScientificName", 
  "AH Species Cover" = "SpeciesCover_AH",
  "AH Species Cover Count" = "Species_AH_n",
  "Mean Species Height (cm)" = "MeanSpeciesHgt",
  "Mean Species Height Count (n)" = "MeanSpeciesHgt_n",
  "Duration" = "Duration",
  "Growth Habitat" = "Woody",
  "Growth Habitat Subcategory" = "Lifeform",
  "Species Key" = "SpeciesKey",
  "Database Key" = "DatabaseKey",
  "Date Loaded in Database" = "DateLoaded"
)

geospecies <- geospecies.raw |>
  rename(!!!setNames(names(col_rename_map), col_rename_map))

# geospecies version for joining
geospecies.join <- geospecies |> 
  filter(!is.na(SpeciesCover_AH)) |> 
  select(PrimaryKey, Species, ScientificName, SpeciesCover_AH, Species_AH_n)


# Prepare LDC data --------------------------------------------------------

# Cols for joining
ldc.007 <- ldc.007.raw |> 
  select(EcoLvl3, Category, PrimaryKey) 



# Arizona/New Mexico Mountains --------------------------------------------

# Species found in ecoregion
az.nm.mts.species <- ldc.007 |> 
  filter(EcoLvl3 == "Arizona/New Mexico Mountains") |> 
  left_join(geospecies.join)

# Most abundant species
az.nm.mts.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRRU2


## 1. Prescribed burn -----------------------------------------------------

# Add invasive species cover
model01.invasive <- model01.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRRU2")

# Create df for plots with 0 invasive species
model01.add0 <- model01.matched |> 
  filter(!PrimaryKey %in% model01.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRRU2",
         ScientificName = "Bromus rubens L.",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model01.invasive <- bind_rows(model01.invasive, model01.add0)



# Arizona/New Mexico Plateau ----------------------------------------------

# Species found in ecoregion
az.nm.plat.species <- ldc.007 |> 
  filter(EcoLvl3 == "Arizona/New Mexico Plateau") |> 
  left_join(geospecies.join)

# Most abundant species
az.nm.plat.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # honestly I am not sure


## 2. Herbicide -----------------------------------------------------------

## 3. Prescribed burn -----------------------------------------------------

## 4. Seeding -------------------------------------------------------------

## 5. Soil disturbance ----------------------------------------------------



# Blue Mountains ----------------------------------------------------------

# Species found in ecoregion
blue.mts.species <- ldc.007 |> 
  filter(EcoLvl3 == "Blue Mountains") |> 
  left_join(geospecies.join)

# Most abundant species
blue.mts.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRTE


## 6. Herbicide -----------------------------------------------------------

# Add invasive species cover
model06.invasive <- model06.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model06.add0 <- model06.matched |> 
  filter(!PrimaryKey %in% model06.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model06.invasive <- bind_rows(model06.invasive, model06.add0)


## 7. Vegetation disturbance ----------------------------------------------

# Add invasive species cover
model07.invasive <- model07.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model07.add0 <- model07.matched |> 
  filter(!PrimaryKey %in% model07.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model07.invasive <- bind_rows(model07.invasive, model07.add0)


## 8. Post-burn herbicide -------------------------------------------------

# Add invasive species cover
model08.invasive <- model08.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model08.add0 <- model08.matched |> 
  filter(!PrimaryKey %in% model08.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model08.invasive <- bind_rows(model08.invasive, model08.add0)



# Central Basin and Range -------------------------------------------------

# Species found in ecoregion
cbr.species <- ldc.007 |> 
  filter(EcoLvl3 == "Central Basin and Range") |> 
  left_join(geospecies.join)

# Most abundant species
cbr.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRTE


## 9. Aerial seeding ------------------------------------------------------

# Add invasive species cover
model09.invasive <- model09.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model09.add0 <- model09.matched |> 
  filter(!PrimaryKey %in% model09.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model09.invasive <- bind_rows(model09.invasive, model09.add0)


## 10. Drill seeding & soil disturbance -----------------------------------

# Add invasive species cover
model10.invasive <- model10.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model10.add0 <- model10.matched |> 
  filter(!PrimaryKey %in% model10.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model10.invasive <- bind_rows(model10.invasive, model10.add0)


## 11. Prescribed burn ----------------------------------------------------

# Add invasive species cover
model11.invasive <- model11.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model11.add0 <- model11.matched |> 
  filter(!PrimaryKey %in% model11.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model11.invasive <- bind_rows(model11.invasive, model11.add0)


## 12. Vegetation disturbance ---------------------------------------------

# Add invasive species cover
model12.invasive <- model12.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model12.add0 <- model12.matched |> 
  filter(!PrimaryKey %in% model12.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model12.invasive <- bind_rows(model12.invasive, model12.add0)


## 13. Post-burn aerial seeding -------------------------------------------

# Add invasive species cover
model13.invasive <- model13.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model13.add0 <- model13.matched |> 
  filter(!PrimaryKey %in% model13.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model13.invasive <- bind_rows(model13.invasive, model13.add0)


# 14. Post-burn drill seeding ---------------------------------------------

# Add invasive species cover
model14.invasive <- model14.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model14.add0 <- model14.matched |> 
  filter(!PrimaryKey %in% model14.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model14.invasive <- bind_rows(model14.invasive, model14.add0)


## 15. Post-burn ground seeding -------------------------------------------

# Add invasive species cover
model15.invasive <- model15.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model15.add0 <- model15.matched |> 
  filter(!PrimaryKey %in% model15.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model15.invasive <- bind_rows(model15.invasive, model15.add0)


## 16. Post-burn herbicide ------------------------------------------------

# Add invasive species cover
model16.invasive <- model16.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model16.add0 <- model16.matched |> 
  filter(!PrimaryKey %in% model16.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model16.invasive <- bind_rows(model16.invasive, model16.add0)
