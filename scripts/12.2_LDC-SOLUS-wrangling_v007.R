# Created: 2026-04-27
# Updated: 2026-04-27

# Purpose: Drop rows/points with no SOLUS data.

library(tidyverse)

# Load data ---------------------------------------------------------------

ldc.cetwi.solus.raw <- read_csv("data/versions-from-R/12.1_LDC-CETWI_v007.csv")


# Data wrangling ----------------------------------------------------------

ldc.cetwi.solus <- ldc.cetwi.solus.raw |> 
  filter(!is.na(sandtotal_0_cm))


# Write to CSV ------------------------------------------------------------

write_csv(ldc.cetwi.solus,
          file = "data/versions-from-R/12.2_LDC-CETWI-SOLUS_v007.csv")
