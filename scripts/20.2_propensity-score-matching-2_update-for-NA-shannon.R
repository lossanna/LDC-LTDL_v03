# Created: 2026-06-05
# Updated: 2026-06-05

# Purpose: Recalculate g computation and average treatment effect for models 29, 32,
#   and 43 after removing rows of plots with NA Shannon diversity values (NA because
#   no cover measurements for any species were recorded in the plots) and saved 
#   updated set of matched data, g computation, and average treatment effect 
#   (for all 46 models). Recreate control/treated figures for the three affected 
#   models (love plots stay the same).

# Plots with NA Shannon diversity values must be removed for permutation tests
#   (23_permutation-tests-2.R). The relevant primary keys/plots were identified
#   in 21_calculate-shannon-diversity-2.R.

# The affected models are:
#   29. Central Basin and Range post-burn aerial seeding
#   32. Central Basin and Range post-burn herbicide
#   43. AZ/NM Plateau soil disturbance


library(tidyverse)
library(MatchIt)
library(cobalt)
library(marginaleffects)
library(ggsignif)

# Load data ---------------------------------------------------------------

load("RData/20.1_matched-data-2.RData")
load("RData/20.1_PSM-2.RData")
load("RData/20.1_g-computation-2.RData")
load("RData/20.1_average-treatment-effect-2.RData")
ldc.007.raw <- read_csv("data/versions-from-R/12.3_LDC-points_v007.csv")
shannon.na <- read_csv("data/versions-from-R/21_shannon-diversity-2_NA.csv")


# Data wrangling ----------------------------------------------------------

# Create binary col for treatment and fire
ldc.007 <- ldc.007.raw |>
  mutate(trt_binary = if_else(is.na(Trt_Type_Sub), 0, 1))

# Natural log transformation of horizontal flux (q) when greater than 0
ldc.007 <- ldc.007 |>
  mutate(
    ln_q = if_else(horizontal_flux_total_MD == 0, 0, log(horizontal_flux_total_MD)),
    .after = horizontal_flux_total_MD
  )


# Check for NAs
apply(ldc.007, 2, anyNA)



# Central Basin and Range -------------------------------------------------

## 29. Post-burn aerial seeding -------------------------------------------

# Filter data
model29.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Aerial Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post-burn")) |>
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
model29.psm <- matchit(
  data = model29.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model29.psm
summary(model29.psm) # 345 treated matched

# Diagnostic love plot
model29.loveplot <- love.plot(model29.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "29. Central Basin and Range: Post-burn aerial seeding") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model29.loveplot

# eCDF plots
bal.plot(model29.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model29.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model29.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model29.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model29.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model29.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model29.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model29.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model29.psm, var = "BareSoil_FH", which = "both")
bal.plot(model29.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model29.psm, var = "ForbCover_AH", which = "both")
bal.plot(model29.psm, var = "GramCover_AH", which = "both")
bal.plot(model29.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model29.psm, var = "Gap100plus", which = "both")
bal.plot(model29.psm, var = "CETWI", which = "both")
bal.plot(model29.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model29.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model29.psm, type = "qq")

# Matched data
model29.matched <- match_data(model29.psm)

# Remove primary keys with NA Shannon diversity values
model29.matched <- model29.matched |> 
  filter(!PrimaryKey %in% shannon.na$PrimaryKey)

# Create trt_control variable
model29.matched <- model29.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn aerial seeding", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn aerial seeding")))

# Center and scale numeric variables
model29.matched <- model29.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model29.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model29.matched,
  weights = weights
)

# G computation to estimate marginal effects
model29.pred <- avg_predictions(
  model = model29.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 29, .before = trt_control)
model29.pred

# Estimation of average treatment effect
model29.comp <- avg_comparisons(
  model = model29.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 29, .before = term)
model29.comp

# Plot
model29.plot <- model29.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "29. Central Basin and Range: Post-burn aerial seeding",
       x = NULL)
model29.plot


## 32. Post-burn herbicide ------------------------------------------------

# Filter data
model32.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post-burn")) |>
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
model32.psm <- matchit(
  data = model32.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model32.psm
summary(model32.psm) # 80 treated matched

# Diagnostic love plot
model32.loveplot <- love.plot(model32.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "32. Central Basin and Range: Post-burn herbicide") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model32.loveplot

# eCDF plots
bal.plot(model32.psm, which = "both", type = "ecdf")
bal.plot(model32.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model32.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model32.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model32.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model32.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model32.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model32.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model32.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model32.psm, which = "both")
bal.plot(model32.psm, var = "BareSoil_FH", which = "both")
bal.plot(model32.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model32.psm, var = "ForbCover_AH", which = "both")
bal.plot(model32.psm, var = "GramCover_AH", which = "both")
bal.plot(model32.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model32.psm, var = "Gap100plus", which = "both")
bal.plot(model32.psm, var = "CETWI", which = "both")
bal.plot(model32.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model32.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model32.psm, type = "qq")

# Matched data
model32.matched <- match_data(model32.psm)

# Remove primary keys with NA Shannon diversity values
model32.matched <- model32.matched |> 
  filter(!PrimaryKey %in% shannon.na$PrimaryKey)

# Create trt_control variable
model32.matched <- model32.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Post-burn herbicide", "Post-burn control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Post-burn control", "Post-burn herbicide")))

# Center and scale numeric variables
model32.matched <- model32.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model32.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model32.matched,
  weights = weights
)

# G computation to estimate marginal effects
model32.pred <- avg_predictions(
  model = model32.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 32, .before = trt_control)
model32.pred # NS for both control & treated (not different from 0)

# Estimation of average treatment effect
model32.comp <- avg_comparisons(
  model = model32.lm,
  variables = "trt_control",
  vcov = ~subclass
) |>
  mutate(Model = 32, .before = term)
model32.comp

# Plot
model32.plot <- model32.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "32. Central Basin and Range: Post-burn herbicide",
       x = NULL)
model32.plot



# Arizona/New Mexico Plateau ----------------------------------------------

### 43. Soil Disturbance --------------------------------------------------

# Filter data
model43.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Soil Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Arizona/New Mexico Plateau")

# PSM
model43.psm <- matchit(
  data = model43.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
model43.psm
summary(model43.psm) # 41 treated matched

# Diagnostic love plot
model43.loveplot <- love.plot(model43.psm, stars = "std",           
                              thresholds = c(m = 0.2, v = 2)) +
  labs(title = "43. AZ/NM Plateau: Soil disturbance") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom")
model43.loveplot

# eCDF plots
bal.plot(model43.psm, which = "both", type = "ecdf")
bal.plot(model43.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(model43.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(model43.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(model43.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(model43.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(model43.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(model43.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(model43.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
bal.plot(model43.psm, which = "both")
bal.plot(model43.psm, var = "BareSoil_FH", which = "both")
bal.plot(model43.psm, var = "TotalFoliarCover", which = "both")
bal.plot(model43.psm, var = "ForbCover_AH", which = "both")
bal.plot(model43.psm, var = "GramCover_AH", which = "both")
bal.plot(model43.psm, var = "ShrubCover_AH", which = "both")
bal.plot(model43.psm, var = "Gap100plus", which = "both")
bal.plot(model43.psm, var = "CETWI", which = "both")
bal.plot(model43.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(model43.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(model43.psm, type = "qq")

# Matched data
model43.matched <- match_data(model43.psm)

# Remove primary keys with NA Shannon diversity values
model43.matched <- model43.matched |> 
  filter(!PrimaryKey %in% shannon.na$PrimaryKey)

# Create trt_control variable
model43.matched <- model43.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Soil disturbance", "Control")) |>
  mutate(trt_control = factor(trt_control, levels = c("Control", "Soil disturbance")))

# Center and scale numeric variables
model43.matched <- model43.matched |>
  mutate(
    BareSoil_FH_scaled = scale(BareSoil_FH, center = TRUE, scale = TRUE)[, 1],
    TotalFoliarCover_scaled = scale(TotalFoliarCover, center = TRUE, scale = TRUE)[, 1],
    ForbCover_AH_scaled = scale(ForbCover_AH, center = TRUE, scale = TRUE)[, 1],
    GramCover_AH_scaled = scale(GramCover_AH, center = TRUE, scale = TRUE)[, 1],
    ShrubCover_AH_scaled = scale(ShrubCover_AH, center = TRUE, scale = TRUE)[, 1],
    Gap100plus_scaled = scale(Gap100plus, center = TRUE, scale = TRUE)[, 1],
    CETWI_scaled = scale(CETWI, center = TRUE, scale = TRUE)[, 1],
    sandtotal_0_cm_scaled = scale(sandtotal_0_cm, center = TRUE, scale = TRUE)[, 1]
  )

# Linear model with covariates
model43.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = model43.matched,
  weights = weights
)

# G computation to estimate marginal effects
model43.pred <- avg_predictions(
  model = model43.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
) |> 
  mutate(Model = 43, .before = trt_control)
model43.pred

# Estimation of average treatment effect
model43.comp <- avg_comparisons(
  model = model43.lm,
  variables = "trt_control",
  vcov = ~subclass
) |> 
  mutate(Model = 43, .before = term)
model43.comp # p = 0.009

# Plot
model43.plot <- model43.pred |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black")) +
  labs(title = "43. AZ/NM Plateau: Soil disturbance",
       x = NULL) +
  geom_signif(
    comparisons = list(c("Control", "Soil disturbance")),
    annotations = c("**")
  )
model43.plot



# Save matched data -------------------------------------------------------

# Matched data
save(list = ls(pattern = "\\.matched$"), 
     file = "RData/20.2_matched-data-2_updated.RData")

# PSM
save(list = ls(pattern = "\\.psm$"), 
     file = "RData/20.2_PSM-2_updated.RData")

# G computation
save(list = ls(pattern = "\\.pred$"), 
     file = "RData/20.2_g-computation-2_updated.RData")

# Average treatment effect
save(list = ls(pattern = "\\.comp$"), 
     file = "RData/20.2_average-treatment-effect-2_updated.RData")



# Write out figures -------------------------------------------------------

## Central Basin and Range ------------------------------------------------

# 29. Central BR: Post-burn aerial seeding
#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests-2/model29_average-treatment-effect_updated.tiff",
     units = "in", width = 6, height = 4, res = 150)
model29.plot
dev.off()


# 32. Central BR: Post-burn herbicide
#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests-2/model32_average-treatment-effect_updated.tiff",
     units = "in", width = 6, height = 4, res = 150)
model32.plot
dev.off()


## AZ/NM Plateau ----------------------------------------------------------

# 43. AZ/NM Plateau: Soil disturbance
#   Treatment effect
tiff("figures/2026-06_PSM-and-permutation-tests-2/model43_average-treatment-effect_updated.tiff",
     units = "in", width = 6, height = 4, res = 150)
model43.plot
dev.off()

