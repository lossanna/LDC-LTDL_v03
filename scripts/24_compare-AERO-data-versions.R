# Created: 2026-06-18
# Updated: 2026-06-18

# Purpose: Check current AERO data with new version that Ron sent on 2026-06-10.

library(tidyverse)
library(readxl)

# Load data ---------------------------------------------------------------

ldc.007.raw <- read_csv("data/versions-from-R/12.3_LDC-points_v007.csv")
aero.new <- read_csv("data/raw/AERO/all_aero_AIM_LMF_20260610.csv")

# Data wrangling ----------------------------------------------------------

# Narrow down new AERO data to LDC points v007
aero.new.ldc.007 <- aero.new |> 
  filter(PrimaryKey %in% ldc.007.raw$PrimaryKey)

# Retain only AERO cols and reorder rows
ldc.007 <- ldc.007.raw |> 
  select(PrimaryKey, horizontal_flux_total_MD) |> 
  arrange(PrimaryKey) |> 
  rename(h_flux_MD_old = horizontal_flux_total_MD)

aero.new.ldc.007 <- aero.new.ldc.007 |> 
  select(PrimaryKey, horizontal_flux_total_MD) |> 
  arrange(PrimaryKey) |> 
  rename(h_flux_MD_new = horizontal_flux_total_MD)

# Compare
identical(ldc.007$PrimaryKey, aero.new.ldc.007$PrimaryKey)
identical(ldc.007, aero.new.ldc.007)

# Combine into single dataframe
aero.compare <- aero.new.ldc.007 |> 
  left_join(ldc.007)

# Differences
aero.diff <- aero.compare |> 
  filter(h_flux_MD_new != h_flux_MD_old) |> 
  mutate(h_flux_diff = h_flux_MD_new - h_flux_MD_old)
summary(aero.diff$h_flux_diff)

# Differences in ln(Q)
aero.diff <- aero.diff |> 
  mutate(ln_q_new = if_else(h_flux_MD_new == 0, 0, log(h_flux_MD_new)),
         ln_q_old = if_else(h_flux_MD_old == 0, 0, log(h_flux_MD_old))) |> 
  mutate(ln_q_diff = ln_q_new - ln_q_old)
summary(aero.diff$ln_q_diff)
