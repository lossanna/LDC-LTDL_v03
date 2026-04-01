# Created: 2026-04-01
# Updated: 2026-04-01

# Purpose: Create LDC points v006, which only includes control points in MLRAs that also
#   have overlap with a treatment polygon.


library(tidyverse)

# Load data ---------------------------------------------------------------

trt.mlra.raw <- read_csv("data/GIS-exports/006_TrtPoly003-MLRA-PairwiseIntersect_export.csv")
ldc003.trtpolyid.raw <- read_csv("data/versions-from-R/06.3_LDC-points-with-TrtPolyID_v003.csv")
ldc.005.raw <- read_csv("data/versions-from-R/08_LDC-points_v005.csv")


# Data wrangling ----------------------------------------------------------

# LDCpointID & TrtPolyID cols for joining
trtid.join <- ldc003.trtpolyid.raw %>% 
  select(LDCpointID, TrtPolyID) 

# Join with LDC005 and filter for treated only
ldc.005.trt <- ldc.005.raw %>% 
  left_join(trtid.join) %>% 
  filter(!Category %in% c("Control, never burned", "Control, post-burn")) 


# Create list of MLRAs with treatment polygons based on LDC005
mlra.005 <- trt.mlra.raw %>% 
  filter(TrtPolyID %in% ldc.005.trt$TrtPolyID) %>% 
  select(MLRA_NAME) %>% 
  distinct(.keep_all = TRUE)

#   Inspect to ensure there are only control points
setdiff(ldc.005.raw$MLRA_name, mlra.005$MLRA_NAME)
mlra.ctrl.inspect <- ldc.005.raw %>% 
  filter(MLRA_name %in% setdiff(ldc.005.raw$MLRA_name, mlra.005$MLRA_NAME))
unique(mlra.ctrl.inspect$Category)


# Filter to include LDC points only in relevant MLRAs
ldc.006 <- ldc.005.raw %>% 
  filter(MLRA_name %in% mlra.005$MLRA_NAME)


# Write to CSV ------------------------------------------------------------

write_csv(ldc.006,
          file = "data/versions-from-R/09_LDC-points_v006.csv",
          na = "")
