# Created: 2026-06-02
# Updated: 2026-06-04

# Purpose: Identify prominent invasive species by ecoregion. Create dataframes of matched
#   data with invasive species cover for each model.


library(tidyverse)

# Load data ---------------------------------------------------------------

geospecies.raw <- read_csv("data/raw/downloaded/ldc-data-2026-06-01/geospecies.csv")
ldc.007.raw <- read_csv("data/versions-from-R/12.3_LDC-points_v007.csv")
load("RData/13_matched-data-1.RData")


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
  arrange(desc(sum_cover)) |> 
  print(n = 20) # SATR12? honestly not sure


## 2. Herbicide -----------------------------------------------------------

# Add invasive species cover
model02.invasive <- model02.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "SATR12")

# Create df for plots with 0 invasive species
model02.add0 <- model02.matched |> 
  filter(!PrimaryKey %in% model02.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "SATR12",
         ScientificName = "Salsola tragus",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model02.invasive <- bind_rows(model02.invasive, model02.add0)


## 3. Prescribed burn -----------------------------------------------------

# Add invasive species cover
model03.invasive <- model03.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "SATR12")

# Create df for plots with 0 invasive species
model03.add0 <- model03.matched |> 
  filter(!PrimaryKey %in% model03.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "SATR12",
         ScientificName = "Salsola tragus",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model03.invasive <- bind_rows(model03.invasive, model03.add0)


## 4. Seeding -------------------------------------------------------------

# Add invasive species cover
model04.invasive <- model04.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "SATR12")

# Create df for plots with 0 invasive species
model04.add0 <- model04.matched |> 
  filter(!PrimaryKey %in% model04.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "SATR12",
         ScientificName = "Salsola tragus",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model04.invasive <- bind_rows(model04.invasive, model04.add0)


## 5. Soil disturbance ----------------------------------------------------

# Add invasive species cover
model05.invasive <- model05.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "SATR12")

# Create df for plots with 0 invasive species
model05.add0 <- model05.matched |> 
  filter(!PrimaryKey %in% model05.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "SATR12",
         ScientificName = "Salsola tragus",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model05.invasive <- bind_rows(model05.invasive, model05.add0)



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


## 14. Post-burn drill seeding --------------------------------------------

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



# Chihuahuan Desert -------------------------------------------------------

# Species found in ecoregion
chihuahuan.species <- ldc.007 |> 
  filter(EcoLvl3 == "Chihuahuan Desert") |> 
  left_join(geospecies.join)

# Most abundant species
chihuahuan.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # PRGL2 (technically native, but probably what the herbicide is for)


## 17. Herbicide ----------------------------------------------------------

# Add invasive species cover
model17.invasive <- model17.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "PRGL2")

# Create df for plots with 0 invasive species
model17.add0 <- model17.matched |> 
  filter(!PrimaryKey %in% model17.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "PRGL2",
         ScientificName = "Prosopis glandulosa",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model17.invasive <- bind_rows(model17.invasive, model17.add0)



# Colorado Plateaus -------------------------------------------------------

# Species found in ecoregion
cop.species <- ldc.007 |> 
  filter(EcoLvl3 == "Colorado Plateaus") |> 
  left_join(geospecies.join)

# Most abundant species
cop.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRTE


## 18. Aerial seeding & soil disturbance ----------------------------------

# Add invasive species cover
model18.invasive <- model18.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model18.add0 <- model18.matched |> 
  filter(!PrimaryKey %in% model18.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model18.invasive <- bind_rows(model18.invasive, model18.add0)


## 19. Herbicide ----------------------------------------------------------

# Add invasive species cover
model19.invasive <- model19.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model19.add0 <- model19.matched |> 
  filter(!PrimaryKey %in% model19.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model19.invasive <- bind_rows(model19.invasive, model19.add0)


## 20. Prescribed burn ----------------------------------------------------

# Add invasive species cover
model20.invasive <- model20.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model20.add0 <- model20.matched |> 
  filter(!PrimaryKey %in% model20.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model20.invasive <- bind_rows(model20.invasive, model20.add0)


## 21. Soil disturbance ---------------------------------------------------

# Add invasive species cover
model21.invasive <- model21.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model21.add0 <- model21.matched |> 
  filter(!PrimaryKey %in% model21.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model21.invasive <- bind_rows(model21.invasive, model21.add0)


## 22. Vegetation disturbance ---------------------------------------------

# Add invasive species cover
model22.invasive <- model22.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model22.add0 <- model22.matched |> 
  filter(!PrimaryKey %in% model22.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model22.invasive <- bind_rows(model22.invasive, model22.add0)


## 23. Post-burn aerial seeding -------------------------------------------

# Add invasive species cover
model23.invasive <- model23.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model23.add0 <- model23.matched |> 
  filter(!PrimaryKey %in% model23.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model23.invasive <- bind_rows(model23.invasive, model23.add0)



# Middle Rockies ----------------------------------------------------------

# Species found in ecoregion
m.rockies.species <- ldc.007 |> 
  filter(EcoLvl3 == "Middle Rockies") |> 
  left_join(geospecies.join)

# Most abundant species
m.rockies.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 20) # BRTE? honestly not sure


## 24. Herbicide ----------------------------------------------------------

# Add invasive species cover
model24.invasive <- model24.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model24.add0 <- model24.matched |> 
  filter(!PrimaryKey %in% model24.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model24.invasive <- bind_rows(model24.invasive, model24.add0)



# Mojave Basin and Range --------------------------------------------------

# Species found in ecoregion
mojave.species <- ldc.007 |> 
  filter(EcoLvl3 == "Mojave Basin and Range") |> 
  left_join(geospecies.join)

# Most abundant species
mojave.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRRU2, SCBA, BRTE, ERCI6



## 25. Post-burn aerial seeding -------------------------------------------

# Add invasive species cover
model25.invasive <- model25.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRRU2")

# Create df for plots with 0 invasive species
model25.add0 <- model25.matched |> 
  filter(!PrimaryKey %in% model25.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRRU2",
         ScientificName = "Bromus rubens",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model25.invasive <- bind_rows(model25.invasive, model25.add0)



# Northern Basin and Range ------------------------------------------------

# Species found in ecoregion
nbr.species <- ldc.007 |> 
  filter(EcoLvl3 == "Northern Basin and Range") |> 
  left_join(geospecies.join)

# Most abundant species
nbr.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRTE


## 26. Drill seeding ------------------------------------------------------

# Add invasive species cover
model26.invasive <- model26.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model26.add0 <- model26.matched |> 
  filter(!PrimaryKey %in% model26.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model26.invasive <- bind_rows(model26.invasive, model26.add0)


## 27. Drill seeding & soil disturbance -----------------------------------

# Add invasive species cover
model27.invasive <- model27.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model27.add0 <- model27.matched |> 
  filter(!PrimaryKey %in% model27.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model27.invasive <- bind_rows(model27.invasive, model27.add0)


## 28. Herbicide ----------------------------------------------------------

# Add invasive species cover
model28.invasive <- model28.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model28.add0 <- model28.matched |> 
  filter(!PrimaryKey %in% model28.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model28.invasive <- bind_rows(model28.invasive, model28.add0)


## 29. Prescribed burn ----------------------------------------------------

# Add invasive species cover
model29.invasive <- model29.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model29.add0 <- model29.matched |> 
  filter(!PrimaryKey %in% model29.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model29.invasive <- bind_rows(model29.invasive, model29.add0)


## 30. Vegetation disturbance ---------------------------------------------

# Add invasive species cover
model30.invasive <- model30.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model30.add0 <- model30.matched |> 
  filter(!PrimaryKey %in% model30.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model30.invasive <- bind_rows(model30.invasive, model30.add0)


## 31. Post-burn aerial seeding -------------------------------------------

# Add invasive species cover
model31.invasive <- model31.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model31.add0 <- model31.matched |> 
  filter(!PrimaryKey %in% model31.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model31.invasive <- bind_rows(model31.invasive, model31.add0)


## 32. Post-burn aerial & drill seeding -----------------------------------

# Add invasive species cover
model32.invasive <- model32.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model32.add0 <- model32.matched |> 
  filter(!PrimaryKey %in% model32.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model32.invasive <- bind_rows(model32.invasive, model32.add0)


## 33. Post-burn closure --------------------------------------------------

# Add invasive species cover
model33.invasive <- model33.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model33.add0 <- model33.matched |> 
  filter(!PrimaryKey %in% model33.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model33.invasive <- bind_rows(model33.invasive, model33.add0)


## 34. Post-burn drill seeding --------------------------------------------

# Add invasive species cover
model34.invasive <- model34.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model34.add0 <- model34.matched |> 
  filter(!PrimaryKey %in% model34.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model34.invasive <- bind_rows(model34.invasive, model34.add0)


## 35. Post-burn herbicide ------------------------------------------------

# Add invasive species cover
model35.invasive <- model35.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model35.add0 <- model35.matched |> 
  filter(!PrimaryKey %in% model35.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model35.invasive <- bind_rows(model35.invasive, model35.add0)


## 36. Post-burn seedling planting ----------------------------------------

# Add invasive species cover
model36.invasive <- model36.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model36.add0 <- model36.matched |> 
  filter(!PrimaryKey %in% model36.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model36.invasive <- bind_rows(model36.invasive, model36.add0)



# Northwestern Great Plains -----------------------------------------------

# Species found in ecoregion
ngp.species <- ldc.007 |> 
  filter(EcoLvl3 == "Northwestern Great Plains") |> 
  left_join(geospecies.join)

# Most abundant species
ngp.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRTE


## 37. Prescribed burn ----------------------------------------------------

# Add invasive species cover
model37.invasive <- model37.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model37.add0 <- model37.matched |> 
  filter(!PrimaryKey %in% model37.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model37.invasive <- bind_rows(model37.invasive, model37.add0)



# Snake River Plain -------------------------------------------------------

# Species found in ecoregion
srp.species <- ldc.007 |> 
  filter(EcoLvl3 == "Snake River Plain") |> 
  left_join(geospecies.join)

# Most abundant species
srp.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRTE


## 38. Post-burn aerial seeding -------------------------------------------

# Add invasive species cover
model38.invasive <- model38.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model38.add0 <- model38.matched |> 
  filter(!PrimaryKey %in% model38.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model38.invasive <- bind_rows(model38.invasive, model38.add0)


## 39. Post-burn aerial & drill seeding -----------------------------------

# Add invasive species cover
model39.invasive <- model39.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model39.add0 <- model39.matched |> 
  filter(!PrimaryKey %in% model39.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model39.invasive <- bind_rows(model39.invasive, model39.add0)


## 40. Post-burn closure --------------------------------------------------

# Add invasive species cover
model40.invasive <- model40.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model40.add0 <- model40.matched |> 
  filter(!PrimaryKey %in% model40.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model40.invasive <- bind_rows(model40.invasive, model40.add0)


## 41. Post-burn drill seeding --------------------------------------------

# Add invasive species cover
model41.invasive <- model41.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model41.add0 <- model41.matched |> 
  filter(!PrimaryKey %in% model41.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model41.invasive <- bind_rows(model41.invasive, model41.add0)


## 42. Post-burn herbicide ------------------------------------------------

# Add invasive species cover
model42.invasive <- model42.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model42.add0 <- model42.matched |> 
  filter(!PrimaryKey %in% model42.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model42.invasive <- bind_rows(model42.invasive, model42.add0)



# Southern Rockies --------------------------------------------------------

# Species found in ecoregion
s.rockies.species <- ldc.007 |> 
  filter(EcoLvl3 == "Southern Rockies") |> 
  left_join(geospecies.join)

# Most abundant species
s.rockies.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) |> 
  print(n = 30) # also not sure (all of the top 20 are native; doing BRTE for now)


## 43. Herbicide ----------------------------------------------------------

# Add invasive species cover
model43.invasive <- model43.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model43.add0 <- model43.matched |> 
  filter(!PrimaryKey %in% model43.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model43.invasive <- bind_rows(model43.invasive, model43.add0)


## 44. Prescribed burn ----------------------------------------------------

# Add invasive species cover
model44.invasive <- model44.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model44.add0 <- model44.matched |> 
  filter(!PrimaryKey %in% model44.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model44.invasive <- bind_rows(model44.invasive, model44.add0)


## 45. Vegetation disturbance ---------------------------------------------

# Add invasive species cover
model45.invasive <- model45.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model45.add0 <- model45.matched |> 
  filter(!PrimaryKey %in% model45.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model45.invasive <- bind_rows(model45.invasive, model45.add0)



# Wyoming Basin -----------------------------------------------------------

# Species found in ecoregion
wyb.species <- ldc.007 |> 
  filter(EcoLvl3 == "Wyoming Basin") |> 
  left_join(geospecies.join)

# Most abundant species
wyb.species |> 
  group_by(Species, ScientificName) |> 
  summarise(sum_cover = sum(SpeciesCover_AH)) |> 
  arrange(desc(sum_cover)) # BRTE


## 43. Prescribed burn ----------------------------------------------------

# Add invasive species cover
model43.invasive <- model43.matched |> 
  select(EcoLvl3, trt_control, PrimaryKey) |> 
  left_join(geospecies.join) |> 
  filter(Species == "BRTE")

# Create df for plots with 0 invasive species
model43.add0 <- model43.matched |> 
  filter(!PrimaryKey %in% model43.invasive$PrimaryKey) |> 
  select(EcoLvl3, trt_control, PrimaryKey) |>
  mutate(Species = "BRTE",
         ScientificName = "Bromus tectorum",
         SpeciesCover_AH = 0,
         Species_AH_n = 0)

# Combine
model43.invasive <- bind_rows(model43.invasive, model43.add0)



# Save --------------------------------------------------------------------

# Save matched dataframes with invasive species cover

save(list = ls(pattern = "\\.invasive$"), 
     file = "RData/17_invasive-matched-data-1.RData")


save.image("RData/17_calculate-invasive-species-1.RData")

