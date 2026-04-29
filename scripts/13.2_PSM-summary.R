# Created: 2026-04-29
# Updated: 2026-04-29

# Purpose: Summarize PSM results (the 13.1.R script is very long).

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
  mutate(Model = c(1:46)) |> 
  select(Model, EcoLvl3, Category, Trt_Type_Sub)

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
  
