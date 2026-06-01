# Created: 2026-06-01
# Updated: 2026-06-01

# Purpose: Collate Landscape Data Commons geospecies and geoindicators data from the   
#   four batches of downloads into single table, and write new CSV.

# Some of the primary keys are missing from this new version of geoindicators versus the 
#   previous one from 2026-03-11, but all of the primary keys needed for LDC points v007 
#   are present, so I'm not going to worry about the other missing primary keys.


library(tidyverse)

# Load data ---------------------------------------------------------------

# 2026-03-11 version of geoindicators
geoindicators.prev <- read_csv("data/raw/downloaded/ldc-data-2026-03-11/geoindicators.csv")

# LDC points v007
ldc.007 <- read_csv("data/versions-from-R/12.3_LDC-points_v007.csv")

# All
geoindicators.all <- read_csv("data/raw/downloaded/ldc-data-2026-06-01/all/geoindicators.csv")

# Batch 1
geoindicators1 <- read_csv("data/raw/downloaded/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180018/geoindicators.csv")
geospecies1 <- read_csv("data/raw/downloaded/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180018/geospecies.csv")

# Batch 2
geoindicators2 <- read_csv("data/raw/downloaded/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180224/geoindicators.csv")
geospecies2 <- read_csv("data/raw/downloaded/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180224/geospecies.csv")

# Batch 3
geoindicators3 <- read_csv("data/raw/downloaded/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180331/geoindicators.csv")
geospecies3 <- read_csv("data/raw/downloaded/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180331/geospecies.csv")

# Batch 4
geoindicators4 <- read_csv("data/raw/downloaded/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180523/geoindicators.csv")
geospecies4 <- read_csv("data/raw/downloaded/ldc-data-2026-06-01/ldc-lossanna-dot-nmsu-at-gmail-dot-com-20260601-180523/geospecies.csv")


# Data wrangling ----------------------------------------------------------

# Combine geoindicators to check for row count (1 row per plot)
geoindicators <- bind_rows(geoindicators1, geoindicators2, geoindicators3, geoindicators4) %>% 
  distinct(.keep_all = TRUE)
nrow(geoindicators) == nrow(geoindicators.all) # all 62,441 plots included

# Combine geospecies
geospecies <- bind_rows(geospecies1, geospecies2, geospecies3, geospecies4) %>% 
  distinct(.keep_all = TRUE)


# Check for missing primary keys from previous geoindicators version
setdiff(geoindicators.prev$`Primary Key`, geoindicators$`Primary Key`)

# # Check for missing primary keys from LDC points v007
setdiff(ldc.007$PrimaryKey, geoindicators$`Primary Key`)



# Write to CSV ------------------------------------------------------------

write_csv(geospecies,
          file = "data/raw/downloaded/ldc-data-2026-06-01/geospecies.csv")
