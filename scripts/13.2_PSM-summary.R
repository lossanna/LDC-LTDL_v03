# Created: 2026-04-29
# Updated: 2026-04-29

# Purpose: Summarize PSM results (the 13.1.R script is very long). Grouped by type of result,
#   not by model.

library(tidyverse)
library(MatchIt)
library(cobalt)
library(marginaleffects)


# Load data ---------------------------------------------------------------

load("RData/13.1_PSM-calculations.RData")
level3.trt30 <- read_csv("data/versions-from-R/12.3_eco3-trt30_count-table.csv")


# Count table -------------------------------------------------------------

# Reformat level3.trt30 table
count.table0 <- level3.trt30 |> 
  select(EcoLvl3, Category, Trt_Type_Sub) |> 
  arrange(Trt_Type_Sub) |> 
  arrange(Category) |> 
  arrange(EcoLvl3) |> 
  mutate(Model = c(1:46), .before = EcoLvl3) 

# Add sample sizes from PSM
count.table.psm <- count.table0 |> 
  mutate(
    Treated_n = c(
      summary(anm.prescribed.burn.psm)[[2]][4, 2], summary(anp.herbicide.psm)[[2]][4, 2],
      summary(anp.prescribed.burn.psm)[[2]][4, 2], summary(anp.seeding.psm)[[2]][4, 2],
      summary(anp.soil.disturbance.psm)[[2]][4, 2], summary(bm.herbicide.psm)[[2]][4, 2],
      summary(bm.veg.disturbance.psm)[[2]][4, 2], summary(bm.pb.herbicide.psm)[[2]][4, 2],
      summary(cbr.aerial.seeding.psm)[[2]][4, 2], summary(cbr.drill.soil.psm)[[2]][4, 2],
      summary(cbr.prescribed.burn.psm)[[2]][4, 2], summary(cbr.veg.disturbance.psm)[[2]][4, 2],
      summary(cbr.pb.aerial.seeding.psm)[[2]][4, 2], summary(cbr.pb.drill.seeding.psm)[[2]][4, 2],
      summary(cbr.pb.ground.seeding.psm)[[2]][4, 2], summary(cbr.pb.herbicide.psm)[[2]][4, 2],
      summary(cd.herbicide.psm)[[2]][4, 2], summary(cp.aerial.soil.psm)[[2]][4, 2],
      summary(cp.herbicide.psm)[[2]][4, 2], summary(cp.prescribed.burn.psm)[[2]][4, 2],
      summary(cp.soil.disturbance.psm)[[2]][4, 2], summary(cp.veg.disturbance.psm)[[2]][4, 2],
      summary(cp.pb.aerial.seeding.psm)[[2]][4, 2], summary(mr.herbicide.psm)[[2]][4, 2],
      summary(mbr.pb.aerial.seeding.psm)[[2]][4, 2], summary(nbr.drill.seeding.psm)[[2]][4, 2],
      summary(nbr.drill.soil.psm)[[2]][4, 2], summary(nbr.herbicide.psm)[[2]][4, 2],
      summary(nbr.prescribed.burn.psm)[[2]][4, 2], summary(nbr.veg.disturbance.psm)[[2]][4, 2],
      summary(nbr.pb.aerial.seeding.psm)[[2]][4, 2], summary(nbr.pb.aerial.drill.psm)[[2]][4, 2],
      summary(nbr.pb.closure.psm)[[2]][4, 2], summary(nbr.pb.drill.seeding.psm)[[2]][4, 2],
      summary(nbr.pb.herbicide.psm)[[2]][4, 2], summary(nbr.pb.seedling.planting.psm)[[2]][4, 2],
      summary(ngp.prescribed.burn.psm)[[2]][4, 2], summary(srp.pb.aerial.seeding.psm)[[2]][4, 2],
      summary(srp.pb.aerial.drill.psm)[[2]][4, 2], summary(srp.pb.closure.psm)[[2]][4, 2],
      summary(srp.pb.drill.seeding.psm)[[2]][4, 2], summary(srp.pb.herbicide.psm)[[2]][4, 2],
      summary(sr.herbicide.psm)[[2]][4, 2], summary(sr.prescribed.burn.psm)[[2]][4, 2],
      summary(sr.veg.disturbance.psm)[[2]][4, 2], summary(wb.prescribed.burn.psm)[[2]][4, 2]
    )
  ) |> 
  mutate(Control_n = Treated_n * 2)



# Love plots --------------------------------------------------------------

love.plot(anm.prescribed.burn.psm, stars = "std",
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "AZ/NM Mt: Prescribed Burn") # 1
love.plot(anp.herbicide.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "AZ/NM Plateau: Herbicide") # 2
love.plot(anp.prescribed.burn.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "AZ/NM Plateau: Prescribed Burn") # 3
love.plot(anp.seeding.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "AZ/NM Plateau: Seeding") # 4
love.plot(anp.soil.disturbance.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "AZ/NM Plateau: Soil Disturbance") # 5
love.plot(bm.herbicide.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "Blue Mountains: Herbicide") # 6
love.plot(bm.veg.disturbance.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "Blue Mountains: Vegetation Disturbance") # 7
love.plot(bm.pb.herbicide.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "Blue Mountains: Herbicide") # 8
love.plot(cbr.aerial.seeding.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "CBR: Aerial seeding") # 9
love.plot(cbr.drill.soil.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "CBR: Drill Seeding, Soil Disturbance") # 10
love.plot(cbr.prescribed.burn.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "CBR: Prescribed Burn") # 11
love.plot(cbr.veg.disturbance.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "CBR: Vegetation Disturbance") # 12
love.plot(cbr.pb.aerial.seeding.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "CBR: Post-burn Aerial seeding") # 13
love.plot(cbr.pb.drill.seeding.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "CBR: Post-burn Drill seeding") # 14
love.plot(cbr.pb.ground.seeding.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "CBR: Post-burn Ground seeding") # 15
love.plot(cbr.pb.herbicide.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "CBR: Post-burn Herbicide") # 16
love.plot(cd.herbicide.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "Chihuahuan: Herbicide") # 17
love.plot(cp.aerial.soil.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "CO Plateaus: Aerial Seeding, Soil Disturbance") # 18
love.plot(cp.herbicide.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "CO Plateaus: Herbicide") # 19
love.plot(cp.prescribed.burn.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "CO Plateaus: Prescribed Burn") # 20
love.plot(cp.soil.disturbance.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "CO Plateaus: Soil Disturbance") # 21
love.plot(cp.veg.disturbance.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "CO Plateaus: Vegetation Disturbance") # 22
love.plot(cp.pb.aerial.seeding.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "CO Plateaus: Post-burn Aerial Seeding") # 23
love.plot(mr.herbicide.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "Middle Rockies: Herbicide") # 24
love.plot(mbr.pb.aerial.seeding.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "Mojave: Aerial Seeding") # 25
love.plot(nbr.drill.seeding.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "NBR: Drill Seeding") # 26
love.plot(nbr.drill.soil.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "NBR: Drill Seeding, Soil Disturbance") # 27
love.plot(nbr.herbicide.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "NBR: Herbicide") # 28 
love.plot(nbr.prescribed.burn.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "NBR: Prescribed Burn") # 29
love.plot(nbr.veg.disturbance.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "NBR: Vegetation Disturbance") # 30
love.plot(nbr.pb.aerial.seeding.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "NBR: Post-burn Aerial Seeding") # 31
love.plot(nbr.pb.aerial.drill.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "NBR: Post-burn Aerial & Drill") # 32
love.plot(nbr.pb.closure.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "NBR: Post-burn Closure") # 33
love.plot(nbr.pb.drill.seeding.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "NBR: Post-burn Drill Seeding") # 34
love.plot(nbr.pb.herbicide.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "NBR: Post-burn Herbicide") # 35
love.plot(nbr.pb.seedling.planting.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "NBR: Post-burn Seedling Planting") # 36
love.plot(ngp.prescribed.burn.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "NW Great Plains: Prescribed Burn") # 37
love.plot(srp.pb.aerial.seeding.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "Snake River: Post-burn Aerial Seeding") # 38
love.plot(srp.pb.aerial.drill.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "Snake River: Post-burn Aerial & Drill") # 39
love.plot(srp.pb.closure.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "Snake River: Post-burn Closure") # 40
love.plot(srp.pb.drill.seeding.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "Snake River: Post-burn Drill Seeding") # 41
love.plot(srp.pb.herbicide.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "Snake River: Post-burn Herbicide") # 42
love.plot(sr.herbicide.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "Southern Rockies: Herbicide") # 43
love.plot(sr.prescribed.burn.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "Southern Rockies: Prescribed Burn") # 44
love.plot(sr.veg.disturbance.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "Southern Rockies: Vegetation Disturbance") # 45
love.plot(wb.prescribed.burn.psm, stars = "std",           
          thresholds = c(m = 0.2, v = 2)) +
  labs(title = "WY Basin: Prescribed Burn") # 46

