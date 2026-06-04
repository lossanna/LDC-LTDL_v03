# Created: 2026-06-01
# Updated: 2026-06-04

# Purpose: Calculate Shannon diversity for each LDC point within each ecoregion/treatment
#   data subset (for each of the 46 models).

library(tidyverse)
library(vegan)

# Load data ---------------------------------------------------------------

load("RData/13.1_matched-data-1.RData")
geospecies.raw <- read_csv("data/raw/downloaded/ldc-data-2026-06-01/geospecies.csv")


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

# geospecies for joining
geospecies.join <- geospecies |> 
  select(PrimaryKey, Species, Species_AH_n) |> 
  filter(!is.na(Species_AH_n))



# Write function to calculate diversity -----------------------------------

calc_shannon_diversity <- function(matched_df, geospecies.join) {
  
  # Narrow down species for model
  species_wide <- geospecies.join |>
    filter(PrimaryKey %in% matched_df$PrimaryKey) |>
    pivot_wider(
      names_from = Species,
      values_from = Species_AH_n,
      values_fill = 0
    )
  
  # Extract PrimaryKey as a separate column
  pk <- species_wide$PrimaryKey
  
  # Create species matrix
  species_matrix <- species_wide |>
    select(-PrimaryKey) |>
    as.data.frame()
  
  # Calculate Shannon diversity
  shannon_vals <- vegan::diversity(species_matrix, index = "shannon")
  
  # Write out diversity df
  tibble::tibble(
    PrimaryKey = pk,
    Shannon = as.numeric(shannon_vals)
  )
}



# Run function for all models ---------------------------------------------


model_list <- lapply(1:46, function(i) {
  matched_name <- sprintf("model%02d.matched", i)
  matched_df <- get(matched_name, envir = .GlobalEnv)
  
  calc_shannon_diversity(matched_df, geospecies.join) |>
    mutate(Model = sprintf("model%02d", i))
})

all_diversity <- bind_rows(model_list)


# Write to CSV ------------------------------------------------------------

write_csv(all_diversity,
          file = "data/versions-from-R/16_shannon-diversity_all-models.csv")
