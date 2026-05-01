# Created: 2026-03-12
# Updated: 2026-03-12

# Purpose: Create LDC points v001, which is the most recent monitoring instance for each
#   unique point in space.

# There are still the same issues with Project v02, where there are 22 points without a date visited,
#   and 62 points with multiple DateVisited values.

library(tidyverse)

# Load data ---------------------------------------------------------------

ldc.pts <- read_csv("data/GIS-exports/001_LDCpts_export.csv")
ldc.countoverlapping <- read_csv("data/GIS-exports/001_LDCpts-CountOverlapping_export.csv")
ldc.overlaptable <- read_csv("data/GIS-exports/001_LDCpts-OverlapTable_export.csv")
geoindicators.raw <- read_csv("data/raw/downloaded/ldc-data-2026-03-11/geoindicators.csv")

# Join tables -------------------------------------------------------------

# Join with OverlapTable
ldc.join <- ldc.pts |>
  left_join(ldc.overlaptable)

#   Look for NAs
apply(ldc.join, 2, anyNA)

# Join with CountOverlapping
ldc.join <- ldc.join |>
  left_join(ldc.countoverlapping)

#   Look for NAs
apply(ldc.join, 2, anyNA)

# Remove unnecessary cols and rename COUNT_
ldc.join <- ldc.join |>
  select(-COUNT_FC, -ORIG_NAME) |>
  rename(ldc_count = COUNT_)


# Add DateVisited column --------------------------------------------------

# Cols from geoindicators
geoindicators.join <- geoindicators.raw |>
  rename(
    ProjectKey = `Project Key`,
    PrimaryKey = `Primary Key`,
    DateVisited = `Date Visited`
  ) |>
  select(ProjectKey, PrimaryKey, DateVisited)

# Join with LDC points
ldc.join <- ldc.join |>
  left_join(geoindicators.join)


# Extract rows with no DateVisited ----------------------------------------

# NA for DateVisited
ldc.date.na <- ldc.join |>
  filter(is.na(DateVisited))
# For now, these rows should just be deleted. (still 22 points)


# Extract rows of most recent monitoring ----------------------------------

# Extract the most recent point for LDC plots that were monitored multiple times
most.recent <- ldc.join |>
  group_by(OVERLAP_OID) |>
  filter(DateVisited == max(DateVisited)) |>
  ungroup()
length(unique(most.recent$OVERLAP_OID)) == nrow(most.recent) # FALSE
#   this means that there are some points that have the same DateVisited, so multiple rows
#     for those cases are created

# Separate out points where there is only one most recent date for DateVisited
#   these ones are fine and don't need fixing
most.recent.single <- most.recent |>
  group_by(OVERLAP_OID) |>
  filter(n() == 1) |>
  ungroup()


## Multiple points/rows for DateVisited ------------------------------------

# Separate out points where DateVisited is the same for multiple rows
most.recent.multiple <- most.recent |>
  group_by(OVERLAP_OID) |>
  filter(n() > 1) |>
  ungroup() |>
  arrange(OVERLAP_OID)

# Unlike Project v02, I didn't retain any of the other columns from geoindicators, so I can't
#   tell what's different between the versions. I am again just going to take the first
#   instance of every duplicate, like how I did in Project v02.


# Retain only the first instance of duplicate rows
most.recent.multiple.fixed <- most.recent.multiple |>
  group_by(OVERLAP_OID) |>
  slice_head(n = 1) |>
  ungroup()


## Combine all with corrections -------------------------------------------

# Remove rows with NA for DateVisited and use correction when there are multiple most recent rows
most.recent.combined <- most.recent.single |>
  bind_rows(most.recent.multiple.fixed) |>
  filter(!is.na(DateVisited)) |>
  select(ORIG_OID, OVERLAP_OID, ProjectKey, PrimaryKey, DateVisited)


# Write LDC001 to CSV -----------------------------------------------------

# All columns
write_csv(most.recent.combined,
  file = "data/versions-from-R/03_LDC-points_v001.csv"
)


save.image("RData/03_LDC-points_v001.RData")
