# Created: 2026-03-30
# Updated: 2026-04-01

# Purpose: Format AERO data into single table that can be added to LDC points.

# AERO run output was read into Excel because there were a few rows where the
#   tab spacing was off and there were the wrong number of columns. Excel could automatically
#   fix that, but read_table() did not.

# horizontal_flux_total_MD is the column of interest; it is the median.
#   MN is mean, LPI is lower prediction interval, UPI is upper prediction interval, 
#   STD is standard deviation.

# This AERO data is the latest available, which is more up to date than what is currently
#   loaded into LDC and part of the downloaded geoindicators data.


library(tidyverse)
library(readxl)

# Load data ---------------------------------------------------------------

aero.aim.raw <- read_xlsx("data/raw/AERO/2026-03-12_AIM-and-LMF.xlsx",
                          sheet = "2026-03-12_summary_AIM_OUT")
aero.lmf.raw <- read_xlsx("data/raw/AERO/2026-03-12_AIM-and-LMF.xlsx",
                          sheet = "2026-03-12_summary_LMF_OUT")
ldc.003.sj <- read_csv("data/GIS-exports/004_LDC003-Eco3-MLRA-SpatialJoin_export.csv")


# Data wrangling ----------------------------------------------------------

# Combine AERO runs
aero <- bind_rows(aero.aim.raw, aero.lmf.raw)

# Create version for joining
aero.join <- aero %>% 
  select(PrimaryKey, horizontal_flux_total_MD)

# Join with LDC points to create v004
ldc.004 <- ldc.003.sj %>% 
  left_join(aero.join) %>% 
  filter(!is.na(horizontal_flux_total_MD))

# Arrange columns
ldc.004 <- ldc.004 %>% 
  select(-NA_L3CODE, -NA_L2CODE, -NA_L1CODE, -MLRA_ID, -LRRSYM, -LRR_NAME)


# Write to CSV ------------------------------------------------------------

write_csv(ldc.004,
          file = "data/versions-from-R/07_LDC-points_v004.csv",
          na = "")
