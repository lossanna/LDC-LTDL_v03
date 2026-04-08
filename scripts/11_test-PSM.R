# Created: 2026-04-07
# Updated: 2026-04-08

# Purpose: Trial run of propensity score matching based on Ron's code (restoration_psm_analysis_20250910.R)
#   with SPEI data (also from Ron) as substitute for CETWI (while we wait to recieve CETWI). 

library(tidyverse)
library(MatchIt)
library(cobalt)
library(marginaleffects)

# Load data ---------------------------------------------------------------

ldc.006.raw <- read_csv("data/versions-from-R/10_LDC-with-PSM-cols_v006.csv")
spei.raw <- read_csv("from-Ron/spei1yr_data.csv")

# Data wrangling ----------------------------------------------------------

# Create binary col for treatment and fire
ldc.006 <- ldc.006.raw %>% 
  mutate(trt_binary = if_else(is.na(Trt_Type_Sub), 0, 1))

# Join SPEI
ldc.006 <- ldc.006 %>% 
  left_join(spei.raw)

# Natural log transformation of horizontal flux (q) when greater than 0
ldc.006 <- ldc.006 %>% 
  mutate(ln_q = if_else(horizontal_flux_total_MD == 0, 0, log(horizontal_flux_total_MD)), 
         .after = horizontal_flux_total_MD)


# Check for NAs
apply(ldc.006, 2, anyNA)

#   Check which is missing gap data
ldc.006 %>% 
  filter(is.na(Gap100plus)) # only one treated missing gap data, in a category with 65, so dropping one is fine

#   Drop points without gap data
ldc.006 <- ldc.006 %>% 
  filter(!is.na(Gap100plus))


# Recalculate Treatment_count
treatment.count <- ldc.006 %>% 
  group_by(EcoLvl3, Category) %>% 
  count(Trt_Type_Sub) %>% 
  rename(Treatment_count = n)

#   Remove old Treatment_count col
ldc.006 <- ldc.006 %>% 
  select(-Treatment_count) 

#   Rejoin
ldc.006 <- treatment.count %>% 
  left_join(ldc.006)


# Treatments
sort(unique(ldc.006$Trt_Type_Sub))

# Ecoregions
sort(unique(ldc.006$EcoLvl3))



# Central Basin and Range -------------------------------------------------

## Aerial seeding, never burned -------------------------------------------

### PSM -------------------------------------------------------------------

# Filter data
cbr.aerial.seeding.dat <- ldc.006 %>% 
  filter(Trt_Type_Sub == "Aerial Seeding" | 
           is.na(Trt_Type_Sub)) %>% 
  filter(str_detect(Category, "never burned")) %>% 
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
cbr.aerial.seeding.psm <- matchit(
  data = cbr.aerial.seeding.dat, 
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + spei1yr,
  distance = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2)

# Summary
summary(cbr.aerial.seeding.psm)

# Diagnostic love plot
love.plot(cbr.aerial.seeding.psm, stars = "std") +
  labs(title = "Aerial seeding, never burned")


# Matched data
cbr.aerial.seeding.matched <- match_data(cbr.aerial.seeding.psm)

# Create trt_control variable
cbr.aerial.seeding.matched <- cbr.aerial.seeding.matched %>% 
  mutate(trt_control = if_else(trt_binary == 1, "Aerial Seeding", "Not Treated")) %>%
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Aerial Seeding")))

# Center and scale numeric variables
cbr.aerial.seeding.matched <- cbr.aerial.seeding.matched %>% 
  mutate(BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
         TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
         ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
         GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
         ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
         Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
         spei1yr_scaled = scale(spei1yr, center = TRUE, scale = TRUE)[, 1])

# Multiple linear regression with interactions
cbr.aerial.seeding.lm <- lm(ln_q ~ trt_control * (
  BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
    ShrubCover_AH_scaled + Gap100plus_scaled + spei1yr_scaled),
  data = cbr.aerial.seeding.matched,
  weights = weights)

# G computation
cbr.aerial.seeding.pred <- avg_predictions(
  model = cbr.aerial.seeding.lm,
  variables = "trt_control",
  vcov = ~ subclass,
  by = "trt_control")
cbr.aerial.seeding.pred

#   Table version
cbr.aerial.seeding.pred.df <- plot_predictions(
  cbr.aerial.seeding.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Aerial Seeding"),
    grid_type = "counterfactual"),
  draw = FALSE,
  vcov = ~ subclass)

# Plot
cbr.aerial.seeding.plot <- ggplot(data = cbr.aerial.seeding.pred.df, 
                                  aes(x = trt_control, y = estimate)) +
  geom_point(shape = 18,
             size = 4,
             color = "red") +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
cbr.aerial.seeding.plot

# Estimation of treatment effect
cbr.aerial.seeding.comp <- avg_comparisons(model = cbr.aerial.seeding.lm,
                                variables = "trt_control",
                                vcov = ~ subclass)
cbr.aerial.seeding.comp
