# Created: 2026-04-29
# Updated: 2026-04-29

# Purpose: First attempt at propensity score matching

# https://kosukeimai.github.io/MatchIt/articles/MatchIt.html

library(tidyverse)
library(MatchIt)
library(cobalt)
library(marginaleffects)

# Load data ---------------------------------------------------------------

ldc.007.raw <- read_csv("data/versions-from-R/12.3_LDC-points_v007.csv")


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



# Arizona/New Mexico Mountains --------------------------------------------

## Prescribed Burn --------------------------------------------------------

# Filter data
anm.prescribed.burn.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Prescribed Burn" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Arizona/New Mexico Mountains")

# PSM
anm.prescribed.burn.psm <- matchit(
  data = anm.prescribed.burn.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
anm.prescribed.burn.psm
summary(anm.prescribed.burn.psm) # 30 treated matched


# Diagnostic love plot
love.plot(anm.prescribed.burn.psm, stars = "std") +
  labs(title = "AZ/NM Mt: Prescribed Burn")

# eCDF plots
plot(anm.prescribed.burn.psm, type = "ecdf")
bal.plot(anm.prescribed.burn.psm, which = "both", type = "ecdf")
bal.plot(anm.prescribed.burn.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(anm.prescribed.burn.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(anm.prescribed.burn.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(anm.prescribed.burn.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(anm.prescribed.burn.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(anm.prescribed.burn.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(anm.prescribed.burn.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(anm.prescribed.burn.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(anm.prescribed.burn.psm, type = "density")
bal.plot(anm.prescribed.burn.psm, which = "both")
bal.plot(anm.prescribed.burn.psm, var = "BareSoil_FH", which = "both")
bal.plot(anm.prescribed.burn.psm, var = "TotalFoliarCover", which = "both")
bal.plot(anm.prescribed.burn.psm, var = "ForbCover_AH", which = "both")
bal.plot(anm.prescribed.burn.psm, var = "GramCover_AH", which = "both")
bal.plot(anm.prescribed.burn.psm, var = "ShrubCover_AH", which = "both")
bal.plot(anm.prescribed.burn.psm, var = "Gap100plus", which = "both")
bal.plot(anm.prescribed.burn.psm, var = "CETWI", which = "both")
bal.plot(anm.prescribed.burn.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(anm.prescribed.burn.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(anm.prescribed.burn.psm, type = "qq")


# Matched data
anm.prescribed.burn.matched <- match_data(anm.prescribed.burn.psm)

# Create trt_control variable
anm.prescribed.burn.matched <- anm.prescribed.burn.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Prescribed Burn", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Prescribed Burn")))

# Center and scale numeric variables
anm.prescribed.burn.matched <- anm.prescribed.burn.matched |>
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
anm.prescribed.burn.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = anm.prescribed.burn.matched,
  weights = weights
)

# G computation to estimate marginal effects
anm.prescribed.burn.pred <- avg_predictions(
  model = anm.prescribed.burn.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
anm.prescribed.burn.pred

#   Table version
anm.prescribed.burn.pred.df <- plot_predictions(
  anm.prescribed.burn.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Prescribed Burn"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
anm.prescribed.burn.pred.df

# Plot
anm.prescribed.burn.plot <- anm.prescribed.burn.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
anm.prescribed.burn.plot

# Estimation of average treatment effect
anm.prescribed.burn.comp <- avg_comparisons(
  model = anm.prescribed.burn.lm,
  variables = "trt_control",
  vcov = ~subclass
)
anm.prescribed.burn.comp




# Arizona/New Mexico Plateau ----------------------------------------------

## Herbicide --------------------------------------------------------------

# Filter data
anp.herbicide.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Arizona/New Mexico Plateau")

# PSM
anp.herbicide.psm <- matchit(
  data = anp.herbicide.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
anp.herbicide.psm
summary(anp.herbicide.psm) # 49 treated matched


# Diagnostic love plot
love.plot(anp.herbicide.psm, stars = "std") +
  labs(title = "AZ/NM Plateau: Herbicide")

# eCDF plots
plot(anp.herbicide.psm, type = "ecdf")
bal.plot(anp.herbicide.psm, which = "both", type = "ecdf")
bal.plot(anp.herbicide.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(anp.herbicide.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(anp.herbicide.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(anp.herbicide.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(anp.herbicide.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(anp.herbicide.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(anp.herbicide.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(anp.herbicide.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(anp.herbicide.psm, type = "density")
bal.plot(anp.herbicide.psm, which = "both")
bal.plot(anp.herbicide.psm, var = "BareSoil_FH", which = "both")
bal.plot(anp.herbicide.psm, var = "TotalFoliarCover", which = "both")
bal.plot(anp.herbicide.psm, var = "ForbCover_AH", which = "both")
bal.plot(anp.herbicide.psm, var = "GramCover_AH", which = "both")
bal.plot(anp.herbicide.psm, var = "ShrubCover_AH", which = "both")
bal.plot(anp.herbicide.psm, var = "Gap100plus", which = "both")
bal.plot(anp.herbicide.psm, var = "CETWI", which = "both")
bal.plot(anp.herbicide.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(anp.herbicide.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(anp.herbicide.psm, type = "qq")


# Matched data
anp.herbicide.matched <- match_data(anp.herbicide.psm)

# Create trt_control variable
anp.herbicide.matched <- anp.herbicide.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Herbicide")))

# Center and scale numeric variables
anp.herbicide.matched <- anp.herbicide.matched |>
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
anp.herbicide.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = anp.herbicide.matched,
  weights = weights
)

# G computation to estimate marginal effects
anp.herbicide.pred <- avg_predictions(
  model = anp.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
anp.herbicide.pred

#   Table version
anp.herbicide.pred.df <- plot_predictions(
  anp.herbicide.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Herbicide"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
anp.herbicide.pred.df

# Plot
anp.herbicide.plot <- anp.herbicide.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
anp.herbicide.plot

# Estimation of average treatment effect
anp.herbicide.comp <- avg_comparisons(
  model = anp.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass
)
anp.herbicide.comp



## Prescribed Burn --------------------------------------------------------

# Filter data
anp.prescribed.burn.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Prescribed Burn" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Arizona/New Mexico Plateau")

# PSM
anp.prescribed.burn.psm <- matchit(
  data = anp.prescribed.burn.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
anp.prescribed.burn.psm
summary(anp.prescribed.burn.psm) # 28 treated matched


# Diagnostic love plot
love.plot(anp.prescribed.burn.psm, stars = "std") +
  labs(title = "AZ/NM Plateau: Prescribed Burn")

# eCDF plots
plot(anp.prescribed.burn.psm, type = "ecdf")
bal.plot(anp.prescribed.burn.psm, which = "both", type = "ecdf")
bal.plot(anp.prescribed.burn.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(anp.prescribed.burn.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(anp.prescribed.burn.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(anp.prescribed.burn.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(anp.prescribed.burn.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(anp.prescribed.burn.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(anp.prescribed.burn.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(anp.prescribed.burn.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(anp.prescribed.burn.psm, type = "density")
bal.plot(anp.prescribed.burn.psm, which = "both")
bal.plot(anp.prescribed.burn.psm, var = "BareSoil_FH", which = "both")
bal.plot(anp.prescribed.burn.psm, var = "TotalFoliarCover", which = "both")
bal.plot(anp.prescribed.burn.psm, var = "ForbCover_AH", which = "both")
bal.plot(anp.prescribed.burn.psm, var = "GramCover_AH", which = "both")
bal.plot(anp.prescribed.burn.psm, var = "ShrubCover_AH", which = "both")
bal.plot(anp.prescribed.burn.psm, var = "Gap100plus", which = "both")
bal.plot(anp.prescribed.burn.psm, var = "CETWI", which = "both")
bal.plot(anp.prescribed.burn.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(anp.prescribed.burn.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(anp.prescribed.burn.psm, type = "qq")


# Matched data
anp.prescribed.burn.matched <- match_data(anp.prescribed.burn.psm)

# Create trt_control variable
anp.prescribed.burn.matched <- anp.prescribed.burn.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Prescribed Burn", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Prescribed Burn")))

# Center and scale numeric variables
anp.prescribed.burn.matched <- anp.prescribed.burn.matched |>
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
anp.prescribed.burn.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = anp.prescribed.burn.matched,
  weights = weights
)

# G computation to estimate marginal effects
anp.prescribed.burn.pred <- avg_predictions(
  model = anp.prescribed.burn.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
anp.prescribed.burn.pred

#   Table version
anp.prescribed.burn.pred.df <- plot_predictions(
  anp.prescribed.burn.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Prescribed Burn"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
anp.prescribed.burn.pred.df

# Plot
anp.prescribed.burn.plot <- anp.prescribed.burn.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
anp.prescribed.burn.plot

# Estimation of average treatment effect
anp.prescribed.burn.comp <- avg_comparisons(
  model = anp.prescribed.burn.lm,
  variables = "trt_control",
  vcov = ~subclass
)
anp.prescribed.burn.comp



## Seeding ----------------------------------------------------------------

# Filter data
anp.seeding.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Arizona/New Mexico Plateau")

# PSM
anp.seeding.psm <- matchit(
  data = anp.seeding.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
anp.seeding.psm
summary(anp.seeding.psm) # 34 treated matched


# Diagnostic love plot
love.plot(anp.seeding.psm, stars = "std") +
  labs(title = "AZ/NM Plateau: Seeding")

# eCDF plots
plot(anp.seeding.psm, type = "ecdf")
bal.plot(anp.seeding.psm, which = "both", type = "ecdf")
bal.plot(anp.seeding.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(anp.seeding.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(anp.seeding.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(anp.seeding.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(anp.seeding.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(anp.seeding.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(anp.seeding.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(anp.seeding.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(anp.seeding.psm, type = "density")
bal.plot(anp.seeding.psm, which = "both")
bal.plot(anp.seeding.psm, var = "BareSoil_FH", which = "both")
bal.plot(anp.seeding.psm, var = "TotalFoliarCover", which = "both")
bal.plot(anp.seeding.psm, var = "ForbCover_AH", which = "both")
bal.plot(anp.seeding.psm, var = "GramCover_AH", which = "both")
bal.plot(anp.seeding.psm, var = "ShrubCover_AH", which = "both")
bal.plot(anp.seeding.psm, var = "Gap100plus", which = "both")
bal.plot(anp.seeding.psm, var = "CETWI", which = "both")
bal.plot(anp.seeding.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(anp.seeding.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(anp.seeding.psm, type = "qq")


# Matched data
anp.seeding.matched <- match_data(anp.seeding.psm)

# Create trt_control variable
anp.seeding.matched <- anp.seeding.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Seeding", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Seeding")))

# Center and scale numeric variables
anp.seeding.matched <- anp.seeding.matched |>
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
anp.seeding.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = anp.seeding.matched,
  weights = weights
)

# G computation to estimate marginal effects
anp.seeding.pred <- avg_predictions(
  model = anp.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
anp.seeding.pred

#   Table version
anp.seeding.pred.df <- plot_predictions(
  anp.seeding.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Seeding"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
anp.seeding.pred.df

# Plot
anp.seeding.plot <- anp.seeding.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
anp.seeding.plot

# Estimation of average treatment effect
anp.seeding.comp <- avg_comparisons(
  model = anp.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass
)
anp.seeding.comp



## Soil Disturbance -------------------------------------------------------

# Filter data
anp.soil.disturbance.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Soil Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Arizona/New Mexico Plateau")

# PSM
anp.soil.disturbance.psm <- matchit(
  data = anp.soil.disturbance.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
anp.soil.disturbance.psm
summary(anp.soil.disturbance.psm) # 41 treated matched


# Diagnostic love plot
love.plot(anp.soil.disturbance.psm, stars = "std") +
  labs(title = "AZ/NM Plateau: Soil Disturbance")

# eCDF plots
plot(anp.soil.disturbance.psm, type = "ecdf")
bal.plot(anp.soil.disturbance.psm, which = "both", type = "ecdf")
bal.plot(anp.soil.disturbance.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(anp.soil.disturbance.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(anp.soil.disturbance.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(anp.soil.disturbance.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(anp.soil.disturbance.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(anp.soil.disturbance.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(anp.soil.disturbance.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(anp.soil.disturbance.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(anp.soil.disturbance.psm, type = "density")
bal.plot(anp.soil.disturbance.psm, which = "both")
bal.plot(anp.soil.disturbance.psm, var = "BareSoil_FH", which = "both")
bal.plot(anp.soil.disturbance.psm, var = "TotalFoliarCover", which = "both")
bal.plot(anp.soil.disturbance.psm, var = "ForbCover_AH", which = "both")
bal.plot(anp.soil.disturbance.psm, var = "GramCover_AH", which = "both")
bal.plot(anp.soil.disturbance.psm, var = "ShrubCover_AH", which = "both")
bal.plot(anp.soil.disturbance.psm, var = "Gap100plus", which = "both")
bal.plot(anp.soil.disturbance.psm, var = "CETWI", which = "both")
bal.plot(anp.soil.disturbance.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(anp.soil.disturbance.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(anp.soil.disturbance.psm, type = "qq")


# Matched data
anp.soil.disturbance.matched <- match_data(anp.soil.disturbance.psm)

# Create trt_control variable
anp.soil.disturbance.matched <- anp.soil.disturbance.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Soil Disturbance", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Soil Disturbance")))

# Center and scale numeric variables
anp.soil.disturbance.matched <- anp.soil.disturbance.matched |>
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
anp.soil.disturbance.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = anp.soil.disturbance.matched,
  weights = weights
)

# G computation to estimate marginal effects
anp.soil.disturbance.pred <- avg_predictions(
  model = anp.soil.disturbance.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
anp.soil.disturbance.pred

#   Table version
anp.soil.disturbance.pred.df <- plot_predictions(
  anp.soil.disturbance.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Soil Disturbance"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
anp.soil.disturbance.pred.df

# Plot
anp.soil.disturbance.plot <- anp.soil.disturbance.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
anp.soil.disturbance.plot

# Estimation of average treatment effect
anp.soil.disturbance.comp <- avg_comparisons(
  model = anp.soil.disturbance.lm,
  variables = "trt_control",
  vcov = ~subclass
)
anp.soil.disturbance.comp




# Blue Mountains ----------------------------------------------------------

## Herbicide --------------------------------------------------------------

# Filter data
bm.herbicide.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Blue Mountains")

# PSM
bm.herbicide.psm <- matchit(
  data = bm.herbicide.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
bm.herbicide.psm
summary(bm.herbicide.psm) # 96 treated matched


# Diagnostic love plot
love.plot(bm.herbicide.psm, stars = "std") +
  labs(title = "Blue Mountains: Herbicide")

# eCDF plots
plot(bm.herbicide.psm, type = "ecdf")
bal.plot(bm.herbicide.psm, which = "both", type = "ecdf")
bal.plot(bm.herbicide.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(bm.herbicide.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(bm.herbicide.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(bm.herbicide.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(bm.herbicide.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(bm.herbicide.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(bm.herbicide.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(bm.herbicide.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(bm.herbicide.psm, type = "density")
bal.plot(bm.herbicide.psm, which = "both")
bal.plot(bm.herbicide.psm, var = "BareSoil_FH", which = "both")
bal.plot(bm.herbicide.psm, var = "TotalFoliarCover", which = "both")
bal.plot(bm.herbicide.psm, var = "ForbCover_AH", which = "both")
bal.plot(bm.herbicide.psm, var = "GramCover_AH", which = "both")
bal.plot(bm.herbicide.psm, var = "ShrubCover_AH", which = "both")
bal.plot(bm.herbicide.psm, var = "Gap100plus", which = "both")
bal.plot(bm.herbicide.psm, var = "CETWI", which = "both")
bal.plot(bm.herbicide.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(bm.herbicide.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(bm.herbicide.psm, type = "qq")


# Matched data
bm.herbicide.matched <- match_data(bm.herbicide.psm)

# Create trt_control variable
bm.herbicide.matched <- bm.herbicide.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Herbicide")))

# Center and scale numeric variables
bm.herbicide.matched <- bm.herbicide.matched |>
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
bm.herbicide.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = bm.herbicide.matched,
  weights = weights
)

# G computation to estimate marginal effects
bm.herbicide.pred <- avg_predictions(
  model = bm.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
bm.herbicide.pred

#   Table version
bm.herbicide.pred.df <- plot_predictions(
  bm.herbicide.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Herbicide"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
bm.herbicide.pred.df

# Plot
bm.herbicide.plot <- bm.herbicide.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
bm.herbicide.plot

# Estimation of average treatment effect
bm.herbicide.comp <- avg_comparisons(
  model = bm.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass
)
bm.herbicide.comp



## Vegetation Disturbance -------------------------------------------------

# Filter data
bm.veg.disturbance.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Vegetation Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Blue Mountains")

# PSM
bm.veg.disturbance.psm <- matchit(
  data = bm.veg.disturbance.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
bm.veg.disturbance.psm
summary(bm.veg.disturbance.psm) # 31 treated matched


# Diagnostic love plot
love.plot(bm.veg.disturbance.psm, stars = "std") +
  labs(title = "Blue Mountains: Vegetation Disturbance")

# eCDF plots
plot(bm.veg.disturbance.psm, type = "ecdf")
bal.plot(bm.veg.disturbance.psm, which = "both", type = "ecdf")
bal.plot(bm.veg.disturbance.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(bm.veg.disturbance.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(bm.veg.disturbance.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(bm.veg.disturbance.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(bm.veg.disturbance.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(bm.veg.disturbance.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(bm.veg.disturbance.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(bm.veg.disturbance.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(bm.veg.disturbance.psm, type = "density")
bal.plot(bm.veg.disturbance.psm, which = "both")
bal.plot(bm.veg.disturbance.psm, var = "BareSoil_FH", which = "both")
bal.plot(bm.veg.disturbance.psm, var = "TotalFoliarCover", which = "both")
bal.plot(bm.veg.disturbance.psm, var = "ForbCover_AH", which = "both")
bal.plot(bm.veg.disturbance.psm, var = "GramCover_AH", which = "both")
bal.plot(bm.veg.disturbance.psm, var = "ShrubCover_AH", which = "both")
bal.plot(bm.veg.disturbance.psm, var = "Gap100plus", which = "both")
bal.plot(bm.veg.disturbance.psm, var = "CETWI", which = "both")
bal.plot(bm.veg.disturbance.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(bm.veg.disturbance.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(bm.veg.disturbance.psm, type = "qq")


# Matched data
bm.veg.disturbance.matched <- match_data(bm.veg.disturbance.psm)

# Create trt_control variable
bm.veg.disturbance.matched <- bm.veg.disturbance.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Vegetation Disturbance", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Vegetation Disturbance")))

# Center and scale numeric variables
bm.veg.disturbance.matched <- bm.veg.disturbance.matched |>
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
bm.veg.disturbance.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = bm.veg.disturbance.matched,
  weights = weights
)

# G computation to estimate marginal effects
bm.veg.disturbance.pred <- avg_predictions(
  model = bm.veg.disturbance.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
bm.veg.disturbance.pred

#   Table version
bm.veg.disturbance.pred.df <- plot_predictions(
  bm.veg.disturbance.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Vegetation Disturbance"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
bm.veg.disturbance.pred.df

# Plot
bm.veg.disturbance.plot <- bm.veg.disturbance.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
bm.veg.disturbance.plot

# Estimation of average treatment effect
bm.veg.disturbance.comp <- avg_comparisons(
  model = bm.veg.disturbance.lm,
  variables = "trt_control",
  vcov = ~subclass
)
bm.veg.disturbance.comp



## Post-burn Herbicide ----------------------------------------------------

# Filter data
bm.pb.herbicide.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Blue Mountains")

# PSM
bm.pb.herbicide.psm <- matchit(
  data = bm.pb.herbicide.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
bm.pb.herbicide.psm
summary(bm.pb.herbicide.psm) # 42 treated matched


# Diagnostic love plot
love.plot(bm.pb.herbicide.psm, stars = "std") +
  labs(title = "Blue Mountains: Herbicide")

# eCDF plots
plot(bm.pb.herbicide.psm, type = "ecdf")
bal.plot(bm.pb.herbicide.psm, which = "both", type = "ecdf")
bal.plot(bm.pb.herbicide.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(bm.pb.herbicide.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(bm.pb.herbicide.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(bm.pb.herbicide.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(bm.pb.herbicide.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(bm.pb.herbicide.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(bm.pb.herbicide.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(bm.pb.herbicide.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(bm.pb.herbicide.psm, type = "density")
bal.plot(bm.pb.herbicide.psm, which = "both")
bal.plot(bm.pb.herbicide.psm, var = "BareSoil_FH", which = "both")
bal.plot(bm.pb.herbicide.psm, var = "TotalFoliarCover", which = "both")
bal.plot(bm.pb.herbicide.psm, var = "ForbCover_AH", which = "both")
bal.plot(bm.pb.herbicide.psm, var = "GramCover_AH", which = "both")
bal.plot(bm.pb.herbicide.psm, var = "ShrubCover_AH", which = "both")
bal.plot(bm.pb.herbicide.psm, var = "Gap100plus", which = "both")
bal.plot(bm.pb.herbicide.psm, var = "CETWI", which = "both")
bal.plot(bm.pb.herbicide.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(bm.pb.herbicide.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(bm.pb.herbicide.psm, type = "qq")


# Matched data
bm.pb.herbicide.matched <- match_data(bm.pb.herbicide.psm)

# Create trt_control variable
bm.pb.herbicide.matched <- bm.pb.herbicide.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Herbicide")))

# Center and scale numeric variables
bm.pb.herbicide.matched <- bm.pb.herbicide.matched |>
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
bm.pb.herbicide.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = bm.pb.herbicide.matched,
  weights = weights
)

# G computation to estimate marginal effects
bm.pb.herbicide.pred <- avg_predictions(
  model = bm.pb.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
bm.pb.herbicide.pred

#   Table version
bm.pb.herbicide.pred.df <- plot_predictions(
  bm.pb.herbicide.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Herbicide"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
bm.pb.herbicide.pred.df

# Plot
bm.pb.herbicide.plot <- bm.pb.herbicide.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
bm.pb.herbicide.plot

# Estimation of average treatment effect
bm.pb.herbicide.comp <- avg_comparisons(
  model = bm.pb.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass
)
bm.pb.herbicide.comp




# Central Basin and Range -------------------------------------------------

## Aerial seeding ---------------------------------------------------------

# Filter data
cbr.aerial.seeding.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Aerial Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never burned")) |>
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
cbr.aerial.seeding.psm <- matchit(
  data = cbr.aerial.seeding.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
cbr.aerial.seeding.psm
summary(cbr.aerial.seeding.psm) # 31 treated matched


# Diagnostic love plot
love.plot(cbr.aerial.seeding.psm, stars = "std") +
  labs(title = "CBR: Aerial seeding")

# eCDF plots
plot(cbr.aerial.seeding.psm, type = "ecdf")
bal.plot(cbr.aerial.seeding.psm, which = "both", type = "ecdf")
bal.plot(cbr.aerial.seeding.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(cbr.aerial.seeding.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(cbr.aerial.seeding.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.aerial.seeding.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.aerial.seeding.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.aerial.seeding.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(cbr.aerial.seeding.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(cbr.aerial.seeding.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(cbr.aerial.seeding.psm, type = "density")
bal.plot(cbr.aerial.seeding.psm, which = "both")
bal.plot(cbr.aerial.seeding.psm, var = "BareSoil_FH", which = "both")
bal.plot(cbr.aerial.seeding.psm, var = "TotalFoliarCover", which = "both")
bal.plot(cbr.aerial.seeding.psm, var = "ForbCover_AH", which = "both")
bal.plot(cbr.aerial.seeding.psm, var = "GramCover_AH", which = "both")
bal.plot(cbr.aerial.seeding.psm, var = "ShrubCover_AH", which = "both")
bal.plot(cbr.aerial.seeding.psm, var = "Gap100plus", which = "both")
bal.plot(cbr.aerial.seeding.psm, var = "CETWI", which = "both")
bal.plot(cbr.aerial.seeding.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(cbr.aerial.seeding.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(cbr.aerial.seeding.psm, type = "qq")


# Matched data
cbr.aerial.seeding.matched <- match_data(cbr.aerial.seeding.psm)

# Create trt_control variable
cbr.aerial.seeding.matched <- cbr.aerial.seeding.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Aerial Seeding", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Aerial Seeding")))

# Center and scale numeric variables
cbr.aerial.seeding.matched <- cbr.aerial.seeding.matched |>
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
cbr.aerial.seeding.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = cbr.aerial.seeding.matched,
  weights = weights
)

# G computation to estimate marginal effects
cbr.aerial.seeding.pred <- avg_predictions(
  model = cbr.aerial.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
cbr.aerial.seeding.pred

#   Table version
cbr.aerial.seeding.pred.df <- plot_predictions(
  cbr.aerial.seeding.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Aerial Seeding"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
cbr.aerial.seeding.pred.df

# Plot
cbr.aerial.seeding.plot <- cbr.aerial.seeding.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
cbr.aerial.seeding.plot

# Estimation of average treatment effect
cbr.aerial.seeding.comp <- avg_comparisons(
  model = cbr.aerial.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass
)
cbr.aerial.seeding.comp



## Drill Seeding, Soil Disturbance ----------------------------------------

# Filter data
cbr.drill.soil.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Drill Seeding, Soil Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never burned")) |>
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
cbr.drill.soil.psm <- matchit(
  data = cbr.drill.soil.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
cbr.drill.soil.psm
summary(cbr.drill.soil.psm) # 36 treated matched


# Diagnostic love plot
love.plot(cbr.drill.soil.psm, stars = "std") +
  labs(title = "CBR: Drill Seeding, Soil Disturbance")

# eCDF plots
plot(cbr.drill.soil.psm, type = "ecdf")
bal.plot(cbr.drill.soil.psm, which = "both", type = "ecdf")
bal.plot(cbr.drill.soil.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(cbr.drill.soil.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(cbr.drill.soil.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.drill.soil.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.drill.soil.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.drill.soil.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(cbr.drill.soil.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(cbr.drill.soil.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(cbr.drill.soil.psm, type = "density")
bal.plot(cbr.drill.soil.psm, which = "both")
bal.plot(cbr.drill.soil.psm, var = "BareSoil_FH", which = "both")
bal.plot(cbr.drill.soil.psm, var = "TotalFoliarCover", which = "both")
bal.plot(cbr.drill.soil.psm, var = "ForbCover_AH", which = "both")
bal.plot(cbr.drill.soil.psm, var = "GramCover_AH", which = "both")
bal.plot(cbr.drill.soil.psm, var = "ShrubCover_AH", which = "both")
bal.plot(cbr.drill.soil.psm, var = "Gap100plus", which = "both")
bal.plot(cbr.drill.soil.psm, var = "CETWI", which = "both")
bal.plot(cbr.drill.soil.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(cbr.drill.soil.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(cbr.drill.soil.psm, type = "qq")


# Matched data
cbr.drill.soil.matched <- match_data(cbr.drill.soil.psm)

# Create trt_control variable
cbr.drill.soil.matched <- cbr.drill.soil.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Drill Seeding, Soil Disturbance", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Drill Seeding, Soil Disturbance")))

# Center and scale numeric variables
cbr.drill.soil.matched <- cbr.drill.soil.matched |>
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
cbr.drill.soil.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = cbr.drill.soil.matched,
  weights = weights
)

# G computation to estimate marginal effects
cbr.drill.soil.pred <- avg_predictions(
  model = cbr.drill.soil.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
cbr.drill.soil.pred

#   Table version
cbr.drill.soil.pred.df <- plot_predictions(
  cbr.drill.soil.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Drill Seeding, Soil Disturbance"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
cbr.drill.soil.pred.df

# Plot
cbr.drill.soil.plot <- cbr.drill.soil.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
cbr.drill.soil.plot

# Estimation of average treatment effect
cbr.drill.soil.comp <- avg_comparisons(
  model = cbr.drill.soil.lm,
  variables = "trt_control",
  vcov = ~subclass
)
cbr.drill.soil.comp



## Prescribed Burn --------------------------------------------------------

# Filter data
cbr.prescribed.burn.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Prescribed Burn" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never burned")) |>
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
cbr.prescribed.burn.psm <- matchit(
  data = cbr.prescribed.burn.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
cbr.prescribed.burn.psm
summary(cbr.prescribed.burn.psm) # 31 treated matched


# Diagnostic love plot
love.plot(cbr.prescribed.burn.psm, stars = "std") +
  labs(title = "CBR: Prescribed Burn")

# eCDF plots
plot(cbr.prescribed.burn.psm, type = "ecdf")
bal.plot(cbr.prescribed.burn.psm, which = "both", type = "ecdf")
bal.plot(cbr.prescribed.burn.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(cbr.prescribed.burn.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(cbr.prescribed.burn.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.prescribed.burn.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.prescribed.burn.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.prescribed.burn.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(cbr.prescribed.burn.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(cbr.prescribed.burn.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(cbr.prescribed.burn.psm, type = "density")
bal.plot(cbr.prescribed.burn.psm, which = "both")
bal.plot(cbr.prescribed.burn.psm, var = "BareSoil_FH", which = "both")
bal.plot(cbr.prescribed.burn.psm, var = "TotalFoliarCover", which = "both")
bal.plot(cbr.prescribed.burn.psm, var = "ForbCover_AH", which = "both")
bal.plot(cbr.prescribed.burn.psm, var = "GramCover_AH", which = "both")
bal.plot(cbr.prescribed.burn.psm, var = "ShrubCover_AH", which = "both")
bal.plot(cbr.prescribed.burn.psm, var = "Gap100plus", which = "both")
bal.plot(cbr.prescribed.burn.psm, var = "CETWI", which = "both")
bal.plot(cbr.prescribed.burn.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(cbr.prescribed.burn.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(cbr.prescribed.burn.psm, type = "qq")


# Matched data
cbr.prescribed.burn.matched <- match_data(cbr.prescribed.burn.psm)

# Create trt_control variable
cbr.prescribed.burn.matched <- cbr.prescribed.burn.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Prescribed Burn", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Prescribed Burn")))

# Center and scale numeric variables
cbr.prescribed.burn.matched <- cbr.prescribed.burn.matched |>
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
cbr.prescribed.burn.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = cbr.prescribed.burn.matched,
  weights = weights
)

# G computation to estimate marginal effects
cbr.prescribed.burn.pred <- avg_predictions(
  model = cbr.prescribed.burn.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
cbr.prescribed.burn.pred

#   Table version
cbr.prescribed.burn.pred.df <- plot_predictions(
  cbr.prescribed.burn.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Prescribed Burn"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
cbr.prescribed.burn.pred.df

# Plot
cbr.prescribed.burn.plot <- cbr.prescribed.burn.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
cbr.prescribed.burn.plot

# Estimation of average treatment effect
cbr.prescribed.burn.comp <- avg_comparisons(
  model = cbr.prescribed.burn.lm,
  variables = "trt_control",
  vcov = ~subclass
)
cbr.prescribed.burn.comp



## Vegetation Disturbance -------------------------------------------------

# Filter data
cbr.veg.disturbance.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Vegetation Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never burned")) |>
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
cbr.veg.disturbance.psm <- matchit(
  data = cbr.veg.disturbance.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
cbr.veg.disturbance.psm
summary(cbr.veg.disturbance.psm) # 99 treated matched


# Diagnostic love plot
love.plot(cbr.veg.disturbance.psm, stars = "std") +
  labs(title = "CBR: Vegetation Disturbance")

# eCDF plots
plot(cbr.veg.disturbance.psm, type = "ecdf")
bal.plot(cbr.veg.disturbance.psm, which = "both", type = "ecdf")
bal.plot(cbr.veg.disturbance.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(cbr.veg.disturbance.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(cbr.veg.disturbance.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.veg.disturbance.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.veg.disturbance.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.veg.disturbance.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(cbr.veg.disturbance.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(cbr.veg.disturbance.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(cbr.veg.disturbance.psm, type = "density")
bal.plot(cbr.veg.disturbance.psm, which = "both")
bal.plot(cbr.veg.disturbance.psm, var = "BareSoil_FH", which = "both")
bal.plot(cbr.veg.disturbance.psm, var = "TotalFoliarCover", which = "both")
bal.plot(cbr.veg.disturbance.psm, var = "ForbCover_AH", which = "both")
bal.plot(cbr.veg.disturbance.psm, var = "GramCover_AH", which = "both")
bal.plot(cbr.veg.disturbance.psm, var = "ShrubCover_AH", which = "both")
bal.plot(cbr.veg.disturbance.psm, var = "Gap100plus", which = "both")
bal.plot(cbr.veg.disturbance.psm, var = "CETWI", which = "both")
bal.plot(cbr.veg.disturbance.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(cbr.veg.disturbance.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(cbr.veg.disturbance.psm, type = "qq")


# Matched data
cbr.veg.disturbance.matched <- match_data(cbr.veg.disturbance.psm)

# Create trt_control variable
cbr.veg.disturbance.matched <- cbr.veg.disturbance.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Vegetation Disturbance", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Vegetation Disturbance")))

# Center and scale numeric variables
cbr.veg.disturbance.matched <- cbr.veg.disturbance.matched |>
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
cbr.veg.disturbance.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = cbr.veg.disturbance.matched,
  weights = weights
)

# G computation to estimate marginal effects
cbr.veg.disturbance.pred <- avg_predictions(
  model = cbr.veg.disturbance.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
cbr.veg.disturbance.pred

#   Table version
cbr.veg.disturbance.pred.df <- plot_predictions(
  cbr.veg.disturbance.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Vegetation Disturbance"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
cbr.veg.disturbance.pred.df

# Plot
cbr.veg.disturbance.plot <- cbr.veg.disturbance.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
cbr.veg.disturbance.plot

# Estimation of average treatment effect
cbr.veg.disturbance.comp <- avg_comparisons(
  model = cbr.veg.disturbance.lm,
  variables = "trt_control",
  vcov = ~subclass
)
cbr.veg.disturbance.comp



## Post-burn Aerial Seeding -----------------------------------------------

# Filter data
cbr.pb.aerial.seeding.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Aerial Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post-burn")) |>
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
cbr.pb.aerial.seeding.psm <- matchit(
  data = cbr.pb.aerial.seeding.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
cbr.pb.aerial.seeding.psm
summary(cbr.pb.aerial.seeding.psm) # 345 treated matched


# Diagnostic love plot
love.plot(cbr.pb.aerial.seeding.psm, stars = "std") +
  labs(title = "CBR: Post-burn Aerial seeding")

# eCDF plots
plot(cbr.pb.aerial.seeding.psm, type = "ecdf")
bal.plot(cbr.pb.aerial.seeding.psm, which = "both", type = "ecdf")
bal.plot(cbr.pb.aerial.seeding.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(cbr.pb.aerial.seeding.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(cbr.pb.aerial.seeding.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.pb.aerial.seeding.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.pb.aerial.seeding.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.pb.aerial.seeding.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(cbr.pb.aerial.seeding.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(cbr.pb.aerial.seeding.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(cbr.pb.aerial.seeding.psm, type = "density")
bal.plot(cbr.pb.aerial.seeding.psm, which = "both")
bal.plot(cbr.pb.aerial.seeding.psm, var = "BareSoil_FH", which = "both")
bal.plot(cbr.pb.aerial.seeding.psm, var = "TotalFoliarCover", which = "both")
bal.plot(cbr.pb.aerial.seeding.psm, var = "ForbCover_AH", which = "both")
bal.plot(cbr.pb.aerial.seeding.psm, var = "GramCover_AH", which = "both")
bal.plot(cbr.pb.aerial.seeding.psm, var = "ShrubCover_AH", which = "both")
bal.plot(cbr.pb.aerial.seeding.psm, var = "Gap100plus", which = "both")
bal.plot(cbr.pb.aerial.seeding.psm, var = "CETWI", which = "both")
bal.plot(cbr.pb.aerial.seeding.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(cbr.pb.aerial.seeding.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(cbr.pb.aerial.seeding.psm, type = "qq")


# Matched data
cbr.pb.aerial.seeding.matched <- match_data(cbr.pb.aerial.seeding.psm)

# Create trt_control variable
cbr.pb.aerial.seeding.matched <- cbr.pb.aerial.seeding.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Aerial Seeding", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Aerial Seeding")))

# Center and scale numeric variables
cbr.pb.aerial.seeding.matched <- cbr.pb.aerial.seeding.matched |>
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
cbr.pb.aerial.seeding.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = cbr.pb.aerial.seeding.matched,
  weights = weights
)

# G computation to estimate marginal effects
cbr.pb.aerial.seeding.pred <- avg_predictions(
  model = cbr.pb.aerial.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
cbr.pb.aerial.seeding.pred

#   Table version
cbr.pb.aerial.seeding.pred.df <- plot_predictions(
  cbr.pb.aerial.seeding.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Aerial Seeding"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
cbr.pb.aerial.seeding.pred.df

# Plot
cbr.pb.aerial.seeding.plot <- cbr.pb.aerial.seeding.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
cbr.pb.aerial.seeding.plot

# Estimation of average treatment effect
cbr.pb.aerial.seeding.comp <- avg_comparisons(
  model = cbr.pb.aerial.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass
)
cbr.pb.aerial.seeding.comp



## Post-burn Drill Seeding ------------------------------------------------

# Filter data
cbr.pb.drill.seeding.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Drill Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post-burn")) |>
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
cbr.pb.drill.seeding.psm <- matchit(
  data = cbr.pb.drill.seeding.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
cbr.pb.drill.seeding.psm
summary(cbr.pb.drill.seeding.psm) # 88 treated matched


# Diagnostic love plot
love.plot(cbr.pb.drill.seeding.psm, stars = "std") +
  labs(title = "CBR: Post-burn Drill seeding")

# eCDF plots
plot(cbr.pb.drill.seeding.psm, type = "ecdf")
bal.plot(cbr.pb.drill.seeding.psm, which = "both", type = "ecdf")
bal.plot(cbr.pb.drill.seeding.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(cbr.pb.drill.seeding.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(cbr.pb.drill.seeding.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.pb.drill.seeding.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.pb.drill.seeding.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.pb.drill.seeding.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(cbr.pb.drill.seeding.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(cbr.pb.drill.seeding.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(cbr.pb.drill.seeding.psm, type = "density")
bal.plot(cbr.pb.drill.seeding.psm, which = "both")
bal.plot(cbr.pb.drill.seeding.psm, var = "BareSoil_FH", which = "both")
bal.plot(cbr.pb.drill.seeding.psm, var = "TotalFoliarCover", which = "both")
bal.plot(cbr.pb.drill.seeding.psm, var = "ForbCover_AH", which = "both")
bal.plot(cbr.pb.drill.seeding.psm, var = "GramCover_AH", which = "both")
bal.plot(cbr.pb.drill.seeding.psm, var = "ShrubCover_AH", which = "both")
bal.plot(cbr.pb.drill.seeding.psm, var = "Gap100plus", which = "both")
bal.plot(cbr.pb.drill.seeding.psm, var = "CETWI", which = "both")
bal.plot(cbr.pb.drill.seeding.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(cbr.pb.drill.seeding.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(cbr.pb.drill.seeding.psm, type = "qq")


# Matched data
cbr.pb.drill.seeding.matched <- match_data(cbr.pb.drill.seeding.psm)

# Create trt_control variable
cbr.pb.drill.seeding.matched <- cbr.pb.drill.seeding.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Drill Seeding", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Drill Seeding")))

# Center and scale numeric variables
cbr.pb.drill.seeding.matched <- cbr.pb.drill.seeding.matched |>
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
cbr.pb.drill.seeding.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = cbr.pb.drill.seeding.matched,
  weights = weights
)

# G computation to estimate marginal effects
cbr.pb.drill.seeding.pred <- avg_predictions(
  model = cbr.pb.drill.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
cbr.pb.drill.seeding.pred

#   Table version
cbr.pb.drill.seeding.pred.df <- plot_predictions(
  cbr.pb.drill.seeding.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Drill Seeding"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
cbr.pb.drill.seeding.pred.df

# Plot
cbr.pb.drill.seeding.plot <- cbr.pb.drill.seeding.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
cbr.pb.drill.seeding.plot

# Estimation of average treatment effect
cbr.pb.drill.seeding.comp <- avg_comparisons(
  model = cbr.pb.drill.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass
)
cbr.pb.drill.seeding.comp



## Post-burn Ground Seeding -----------------------------------------------

# Filter data
cbr.pb.ground.seeding.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Ground Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post-burn")) |>
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
cbr.pb.ground.seeding.psm <- matchit(
  data = cbr.pb.ground.seeding.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
cbr.pb.ground.seeding.psm
summary(cbr.pb.ground.seeding.psm) # 38 treated matched


# Diagnostic love plot
love.plot(cbr.pb.ground.seeding.psm, stars = "std") +
  labs(title = "CBR: Post-burn Ground seeding")

# eCDF plots
plot(cbr.pb.ground.seeding.psm, type = "ecdf")
bal.plot(cbr.pb.ground.seeding.psm, which = "both", type = "ecdf")
bal.plot(cbr.pb.ground.seeding.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(cbr.pb.ground.seeding.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(cbr.pb.ground.seeding.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.pb.ground.seeding.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.pb.ground.seeding.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.pb.ground.seeding.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(cbr.pb.ground.seeding.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(cbr.pb.ground.seeding.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(cbr.pb.ground.seeding.psm, type = "density")
bal.plot(cbr.pb.ground.seeding.psm, which = "both")
bal.plot(cbr.pb.ground.seeding.psm, var = "BareSoil_FH", which = "both")
bal.plot(cbr.pb.ground.seeding.psm, var = "TotalFoliarCover", which = "both")
bal.plot(cbr.pb.ground.seeding.psm, var = "ForbCover_AH", which = "both")
bal.plot(cbr.pb.ground.seeding.psm, var = "GramCover_AH", which = "both")
bal.plot(cbr.pb.ground.seeding.psm, var = "ShrubCover_AH", which = "both")
bal.plot(cbr.pb.ground.seeding.psm, var = "Gap100plus", which = "both")
bal.plot(cbr.pb.ground.seeding.psm, var = "CETWI", which = "both")
bal.plot(cbr.pb.ground.seeding.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(cbr.pb.ground.seeding.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(cbr.pb.ground.seeding.psm, type = "qq")


# Matched data
cbr.pb.ground.seeding.matched <- match_data(cbr.pb.ground.seeding.psm)

# Create trt_control variable
cbr.pb.ground.seeding.matched <- cbr.pb.ground.seeding.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Ground Seeding", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Ground Seeding")))

# Center and scale numeric variables
cbr.pb.ground.seeding.matched <- cbr.pb.ground.seeding.matched |>
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
cbr.pb.ground.seeding.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = cbr.pb.ground.seeding.matched,
  weights = weights
)

# G computation to estimate marginal effects
cbr.pb.ground.seeding.pred <- avg_predictions(
  model = cbr.pb.ground.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
cbr.pb.ground.seeding.pred

#   Table version
cbr.pb.ground.seeding.pred.df <- plot_predictions(
  cbr.pb.ground.seeding.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Ground Seeding"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
cbr.pb.ground.seeding.pred.df

# Plot
cbr.pb.ground.seeding.plot <- cbr.pb.ground.seeding.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
cbr.pb.ground.seeding.plot

# Estimation of average treatment effect
cbr.pb.ground.seeding.comp <- avg_comparisons(
  model = cbr.pb.ground.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass
)
cbr.pb.ground.seeding.comp



## Post-burn Herbicide ----------------------------------------------------

# Filter data
cbr.pb.herbicide.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post-burn")) |>
  filter(EcoLvl3 == "Central Basin and Range")

# PSM
cbr.pb.herbicide.psm <- matchit(
  data = cbr.pb.herbicide.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
cbr.pb.herbicide.psm
summary(cbr.pb.herbicide.psm) # 80 treated matched


# Diagnostic love plot
love.plot(cbr.pb.herbicide.psm, stars = "std") +
  labs(title = "CBR: Post-burn Herbicide")

# eCDF plots
plot(cbr.pb.herbicide.psm, type = "ecdf")
bal.plot(cbr.pb.herbicide.psm, which = "both", type = "ecdf")
bal.plot(cbr.pb.herbicide.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(cbr.pb.herbicide.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(cbr.pb.herbicide.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.pb.herbicide.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.pb.herbicide.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(cbr.pb.herbicide.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(cbr.pb.herbicide.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(cbr.pb.herbicide.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(cbr.pb.herbicide.psm, type = "density")
bal.plot(cbr.pb.herbicide.psm, which = "both")
bal.plot(cbr.pb.herbicide.psm, var = "BareSoil_FH", which = "both")
bal.plot(cbr.pb.herbicide.psm, var = "TotalFoliarCover", which = "both")
bal.plot(cbr.pb.herbicide.psm, var = "ForbCover_AH", which = "both")
bal.plot(cbr.pb.herbicide.psm, var = "GramCover_AH", which = "both")
bal.plot(cbr.pb.herbicide.psm, var = "ShrubCover_AH", which = "both")
bal.plot(cbr.pb.herbicide.psm, var = "Gap100plus", which = "both")
bal.plot(cbr.pb.herbicide.psm, var = "CETWI", which = "both")
bal.plot(cbr.pb.herbicide.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(cbr.pb.herbicide.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(cbr.pb.herbicide.psm, type = "qq")


# Matched data
cbr.pb.herbicide.matched <- match_data(cbr.pb.herbicide.psm)

# Create trt_control variable
cbr.pb.herbicide.matched <- cbr.pb.herbicide.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Herbicide")))

# Center and scale numeric variables
cbr.pb.herbicide.matched <- cbr.pb.herbicide.matched |>
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
cbr.pb.herbicide.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = cbr.pb.herbicide.matched,
  weights = weights
)

# G computation to estimate marginal effects
cbr.pb.herbicide.pred <- avg_predictions(
  model = cbr.pb.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
cbr.pb.herbicide.pred

#   Table version
cbr.pb.herbicide.pred.df <- plot_predictions(
  cbr.pb.herbicide.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Herbicide"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
cbr.pb.herbicide.pred.df

# Plot
cbr.pb.herbicide.plot <- cbr.pb.herbicide.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
cbr.pb.herbicide.plot

# Estimation of average treatment effect
cbr.pb.herbicide.comp <- avg_comparisons(
  model = cbr.pb.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass
)
cbr.pb.herbicide.comp




# Chihuahuan Desert -------------------------------------------------------

## Herbicide --------------------------------------------------------------

# Filter data
cd.herbicide.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Chihuahuan Desert")

# PSM
cd.herbicide.psm <- matchit(
  data = cd.herbicide.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
cd.herbicide.psm
summary(cd.herbicide.psm) # 64 treated matched


# Diagnostic love plot
love.plot(cd.herbicide.psm, stars = "std") +
  labs(title = "Chihuahuan: Herbicide")

# eCDF plots
plot(cd.herbicide.psm, type = "ecdf")
bal.plot(cd.herbicide.psm, which = "both", type = "ecdf")
bal.plot(cd.herbicide.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(cd.herbicide.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(cd.herbicide.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(cd.herbicide.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(cd.herbicide.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(cd.herbicide.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(cd.herbicide.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(cd.herbicide.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(cd.herbicide.psm, type = "density")
bal.plot(cd.herbicide.psm, which = "both")
bal.plot(cd.herbicide.psm, var = "BareSoil_FH", which = "both")
bal.plot(cd.herbicide.psm, var = "TotalFoliarCover", which = "both")
bal.plot(cd.herbicide.psm, var = "ForbCover_AH", which = "both")
bal.plot(cd.herbicide.psm, var = "GramCover_AH", which = "both")
bal.plot(cd.herbicide.psm, var = "ShrubCover_AH", which = "both")
bal.plot(cd.herbicide.psm, var = "Gap100plus", which = "both")
bal.plot(cd.herbicide.psm, var = "CETWI", which = "both")
bal.plot(cd.herbicide.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(cd.herbicide.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(cd.herbicide.psm, type = "qq")


# Matched data
cd.herbicide.matched <- match_data(cd.herbicide.psm)

# Create trt_control variable
cd.herbicide.matched <- cd.herbicide.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Herbicide")))

# Center and scale numeric variables
cd.herbicide.matched <- cd.herbicide.matched |>
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
cd.herbicide.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = cd.herbicide.matched,
  weights = weights
)

# G computation to estimate marginal effects
cd.herbicide.pred <- avg_predictions(
  model = cd.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
cd.herbicide.pred

#   Table version
cd.herbicide.pred.df <- plot_predictions(
  cd.herbicide.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Herbicide"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
cd.herbicide.pred.df

# Plot
cd.herbicide.plot <- cd.herbicide.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
cd.herbicide.plot

# Estimation of average treatment effect
cd.herbicide.comp <- avg_comparisons(
  model = cd.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass
)
cd.herbicide.comp




# Colorado Plateaus -------------------------------------------------------

## Aerial Seeding, Soil Disturbance ---------------------------------------

# Filter data
cp.aerial.soil.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Aerial Seeding, Soil Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never burned")) |>
  filter(EcoLvl3 == "Colorado Plateaus")

# PSM
cp.aerial.soil.psm <- matchit(
  data = cp.aerial.soil.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
cp.aerial.soil.psm
summary(cp.aerial.soil.psm) # 60 treated matched


# Diagnostic love plot
love.plot(cp.aerial.soil.psm, stars = "std") +
  labs(title = "CO Plateaus: Aerial Seeding, Soil Disturbance")

# eCDF plots
plot(cp.aerial.soil.psm, type = "ecdf")
bal.plot(cp.aerial.soil.psm, which = "both", type = "ecdf")
bal.plot(cp.aerial.soil.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(cp.aerial.soil.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(cp.aerial.soil.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(cp.aerial.soil.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(cp.aerial.soil.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(cp.aerial.soil.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(cp.aerial.soil.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(cp.aerial.soil.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(cp.aerial.soil.psm, type = "density")
bal.plot(cp.aerial.soil.psm, which = "both")
bal.plot(cp.aerial.soil.psm, var = "BareSoil_FH", which = "both")
bal.plot(cp.aerial.soil.psm, var = "TotalFoliarCover", which = "both")
bal.plot(cp.aerial.soil.psm, var = "ForbCover_AH", which = "both")
bal.plot(cp.aerial.soil.psm, var = "GramCover_AH", which = "both")
bal.plot(cp.aerial.soil.psm, var = "ShrubCover_AH", which = "both")
bal.plot(cp.aerial.soil.psm, var = "Gap100plus", which = "both")
bal.plot(cp.aerial.soil.psm, var = "CETWI", which = "both")
bal.plot(cp.aerial.soil.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(cp.aerial.soil.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(cp.aerial.soil.psm, type = "qq")


# Matched data
cp.aerial.soil.matched <- match_data(cp.aerial.soil.psm)

# Create trt_control variable
cp.aerial.soil.matched <- cp.aerial.soil.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Aerial Seeding, Soil Disturbance", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Aerial Seeding, Soil Disturbance")))

# Center and scale numeric variables
cp.aerial.soil.matched <- cp.aerial.soil.matched |>
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
cp.aerial.soil.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = cp.aerial.soil.matched,
  weights = weights
)

# G computation to estimate marginal effects
cp.aerial.soil.pred <- avg_predictions(
  model = cp.aerial.soil.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
cp.aerial.soil.pred

#   Table version
cp.aerial.soil.pred.df <- plot_predictions(
  cp.aerial.soil.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Aerial Seeding, Soil Disturbance"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
cp.aerial.soil.pred.df

# Plot
cp.aerial.soil.plot <- cp.aerial.soil.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
cp.aerial.soil.plot

# Estimation of average treatment effect
cp.aerial.soil.comp <- avg_comparisons(
  model = cp.aerial.soil.lm,
  variables = "trt_control",
  vcov = ~subclass
)
cp.aerial.soil.comp



## Herbicide --------------------------------------------------------------

# Filter data
cp.herbicide.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Colorado Plateaus")

# PSM
cp.herbicide.psm <- matchit(
  data = cp.herbicide.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
cp.herbicide.psm
summary(cp.herbicide.psm) # 47 treated matched


# Diagnostic love plot
love.plot(cp.herbicide.psm, stars = "std") +
  labs(title = "CO Plateaus: Herbicide")

# eCDF plots
plot(cp.herbicide.psm, type = "ecdf")
bal.plot(cp.herbicide.psm, which = "both", type = "ecdf")
bal.plot(cp.herbicide.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(cp.herbicide.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(cp.herbicide.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(cp.herbicide.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(cp.herbicide.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(cp.herbicide.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(cp.herbicide.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(cp.herbicide.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(cp.herbicide.psm, type = "density")
bal.plot(cp.herbicide.psm, which = "both")
bal.plot(cp.herbicide.psm, var = "BareSoil_FH", which = "both")
bal.plot(cp.herbicide.psm, var = "TotalFoliarCover", which = "both")
bal.plot(cp.herbicide.psm, var = "ForbCover_AH", which = "both")
bal.plot(cp.herbicide.psm, var = "GramCover_AH", which = "both")
bal.plot(cp.herbicide.psm, var = "ShrubCover_AH", which = "both")
bal.plot(cp.herbicide.psm, var = "Gap100plus", which = "both")
bal.plot(cp.herbicide.psm, var = "CETWI", which = "both")
bal.plot(cp.herbicide.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(cp.herbicide.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(cp.herbicide.psm, type = "qq")


# Matched data
cp.herbicide.matched <- match_data(cp.herbicide.psm)

# Create trt_control variable
cp.herbicide.matched <- cp.herbicide.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Herbicide")))

# Center and scale numeric variables
cp.herbicide.matched <- cp.herbicide.matched |>
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
cp.herbicide.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = cp.herbicide.matched,
  weights = weights
)

# G computation to estimate marginal effects
cp.herbicide.pred <- avg_predictions(
  model = cp.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
cp.herbicide.pred

#   Table version
cp.herbicide.pred.df <- plot_predictions(
  cp.herbicide.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Herbicide"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
cp.herbicide.pred.df

# Plot
cp.herbicide.plot <- cp.herbicide.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
cp.herbicide.plot

# Estimation of average treatment effect
cp.herbicide.comp <- avg_comparisons(
  model = cp.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass
)
cp.herbicide.comp



## Prescribed Burn --------------------------------------------------------

# Filter data
cp.prescribed.burn.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Prescribed Burn" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Colorado Plateaus")

# PSM
cp.prescribed.burn.psm <- matchit(
  data = cp.prescribed.burn.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
cp.prescribed.burn.psm
summary(cp.prescribed.burn.psm) # 76 treated matched


# Diagnostic love plot
love.plot(cp.prescribed.burn.psm, stars = "std") +
  labs(title = "CO Plateaus: Prescribed Burn")

# eCDF plots
plot(cp.prescribed.burn.psm, type = "ecdf")
bal.plot(cp.prescribed.burn.psm, which = "both", type = "ecdf")
bal.plot(cp.prescribed.burn.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(cp.prescribed.burn.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(cp.prescribed.burn.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(cp.prescribed.burn.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(cp.prescribed.burn.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(cp.prescribed.burn.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(cp.prescribed.burn.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(cp.prescribed.burn.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(cp.prescribed.burn.psm, type = "density")
bal.plot(cp.prescribed.burn.psm, which = "both")
bal.plot(cp.prescribed.burn.psm, var = "BareSoil_FH", which = "both")
bal.plot(cp.prescribed.burn.psm, var = "TotalFoliarCover", which = "both")
bal.plot(cp.prescribed.burn.psm, var = "ForbCover_AH", which = "both")
bal.plot(cp.prescribed.burn.psm, var = "GramCover_AH", which = "both")
bal.plot(cp.prescribed.burn.psm, var = "ShrubCover_AH", which = "both")
bal.plot(cp.prescribed.burn.psm, var = "Gap100plus", which = "both")
bal.plot(cp.prescribed.burn.psm, var = "CETWI", which = "both")
bal.plot(cp.prescribed.burn.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(cp.prescribed.burn.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(cp.prescribed.burn.psm, type = "qq")


# Matched data
cp.prescribed.burn.matched <- match_data(cp.prescribed.burn.psm)

# Create trt_control variable
cp.prescribed.burn.matched <- cp.prescribed.burn.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Prescribed Burn", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Prescribed Burn")))

# Center and scale numeric variables
cp.prescribed.burn.matched <- cp.prescribed.burn.matched |>
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
cp.prescribed.burn.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = cp.prescribed.burn.matched,
  weights = weights
)

# G computation to estimate marginal effects
cp.prescribed.burn.pred <- avg_predictions(
  model = cp.prescribed.burn.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
cp.prescribed.burn.pred

#   Table version
cp.prescribed.burn.pred.df <- plot_predictions(
  cp.prescribed.burn.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Prescribed Burn"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
cp.prescribed.burn.pred.df

# Plot
cp.prescribed.burn.plot <- cp.prescribed.burn.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
cp.prescribed.burn.plot

# Estimation of average treatment effect
cp.prescribed.burn.comp <- avg_comparisons(
  model = cp.prescribed.burn.lm,
  variables = "trt_control",
  vcov = ~subclass
)
cp.prescribed.burn.comp



## Soil Disturbance -------------------------------------------------------

# Filter data
cp.soil.disturbance.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Soil Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Colorado Plateaus")

# PSM
cp.soil.disturbance.psm <- matchit(
  data = cp.soil.disturbance.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
cp.soil.disturbance.psm
summary(cp.soil.disturbance.psm) # 40 treated matched


# Diagnostic love plot
love.plot(cp.soil.disturbance.psm, stars = "std") +
  labs(title = "CO Plateaus: Soil Disturbance")

# eCDF plots
plot(cp.soil.disturbance.psm, type = "ecdf")
bal.plot(cp.soil.disturbance.psm, which = "both", type = "ecdf")
bal.plot(cp.soil.disturbance.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(cp.soil.disturbance.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(cp.soil.disturbance.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(cp.soil.disturbance.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(cp.soil.disturbance.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(cp.soil.disturbance.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(cp.soil.disturbance.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(cp.soil.disturbance.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(cp.soil.disturbance.psm, type = "density")
bal.plot(cp.soil.disturbance.psm, which = "both")
bal.plot(cp.soil.disturbance.psm, var = "BareSoil_FH", which = "both")
bal.plot(cp.soil.disturbance.psm, var = "TotalFoliarCover", which = "both")
bal.plot(cp.soil.disturbance.psm, var = "ForbCover_AH", which = "both")
bal.plot(cp.soil.disturbance.psm, var = "GramCover_AH", which = "both")
bal.plot(cp.soil.disturbance.psm, var = "ShrubCover_AH", which = "both")
bal.plot(cp.soil.disturbance.psm, var = "Gap100plus", which = "both")
bal.plot(cp.soil.disturbance.psm, var = "CETWI", which = "both")
bal.plot(cp.soil.disturbance.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(cp.soil.disturbance.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(cp.soil.disturbance.psm, type = "qq")


# Matched data
cp.soil.disturbance.matched <- match_data(cp.soil.disturbance.psm)

# Create trt_control variable
cp.soil.disturbance.matched <- cp.soil.disturbance.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Soil Disturbance", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Soil Disturbance")))

# Center and scale numeric variables
cp.soil.disturbance.matched <- cp.soil.disturbance.matched |>
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
cp.soil.disturbance.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = cp.soil.disturbance.matched,
  weights = weights
)

# G computation to estimate marginal effects
cp.soil.disturbance.pred <- avg_predictions(
  model = cp.soil.disturbance.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
cp.soil.disturbance.pred

#   Table version
cp.soil.disturbance.pred.df <- plot_predictions(
  cp.soil.disturbance.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Soil Disturbance"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
cp.soil.disturbance.pred.df

# Plot
cp.soil.disturbance.plot <- cp.soil.disturbance.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
cp.soil.disturbance.plot

# Estimation of average treatment effect
cp.soil.disturbance.comp <- avg_comparisons(
  model = cp.soil.disturbance.lm,
  variables = "trt_control",
  vcov = ~subclass
)
cp.soil.disturbance.comp



## Vegetation Disturbance -------------------------------------------------

# Filter data
cp.veg.disturbance.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Vegetation Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Colorado Plateaus")

# PSM
cp.veg.disturbance.psm <- matchit(
  data = cp.veg.disturbance.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
cp.veg.disturbance.psm
summary(cp.veg.disturbance.psm) # 34 treated matched


# Diagnostic love plot
love.plot(cp.veg.disturbance.psm, stars = "std") +
  labs(title = "CO Plateaus: Vegetation Disturbance")

# eCDF plots
plot(cp.veg.disturbance.psm, type = "ecdf")
bal.plot(cp.veg.disturbance.psm, which = "both", type = "ecdf")
bal.plot(cp.veg.disturbance.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(cp.veg.disturbance.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(cp.veg.disturbance.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(cp.veg.disturbance.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(cp.veg.disturbance.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(cp.veg.disturbance.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(cp.veg.disturbance.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(cp.veg.disturbance.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(cp.veg.disturbance.psm, type = "density")
bal.plot(cp.veg.disturbance.psm, which = "both")
bal.plot(cp.veg.disturbance.psm, var = "BareSoil_FH", which = "both")
bal.plot(cp.veg.disturbance.psm, var = "TotalFoliarCover", which = "both")
bal.plot(cp.veg.disturbance.psm, var = "ForbCover_AH", which = "both")
bal.plot(cp.veg.disturbance.psm, var = "GramCover_AH", which = "both")
bal.plot(cp.veg.disturbance.psm, var = "ShrubCover_AH", which = "both")
bal.plot(cp.veg.disturbance.psm, var = "Gap100plus", which = "both")
bal.plot(cp.veg.disturbance.psm, var = "CETWI", which = "both")
bal.plot(cp.veg.disturbance.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(cp.veg.disturbance.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(cp.veg.disturbance.psm, type = "qq")


# Matched data
cp.veg.disturbance.matched <- match_data(cp.veg.disturbance.psm)

# Create trt_control variable
cp.veg.disturbance.matched <- cp.veg.disturbance.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Vegetation Disturbance", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Vegetation Disturbance")))

# Center and scale numeric variables
cp.veg.disturbance.matched <- cp.veg.disturbance.matched |>
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
cp.veg.disturbance.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = cp.veg.disturbance.matched,
  weights = weights
)

# G computation to estimate marginal effects
cp.veg.disturbance.pred <- avg_predictions(
  model = cp.veg.disturbance.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
cp.veg.disturbance.pred

#   Table version
cp.veg.disturbance.pred.df <- plot_predictions(
  cp.veg.disturbance.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Vegetation Disturbance"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
cp.veg.disturbance.pred.df

# Plot
cp.veg.disturbance.plot <- cp.veg.disturbance.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
cp.veg.disturbance.plot

# Estimation of average treatment effect
cp.veg.disturbance.comp <- avg_comparisons(
  model = cp.veg.disturbance.lm,
  variables = "trt_control",
  vcov = ~subclass
)
cp.veg.disturbance.comp



## Post-burn Aerial Seeding -----------------------------------------------

# Filter data
cp.pb.aerial.seeding.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Aerial Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Colorado Plateaus")

# PSM
cp.pb.aerial.seeding.psm <- matchit(
  data = cp.pb.aerial.seeding.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
cp.pb.aerial.seeding.psm
summary(cp.pb.aerial.seeding.psm) # 58 treated matched


# Diagnostic love plot
love.plot(cp.pb.aerial.seeding.psm, stars = "std") +
  labs(title = "CO Plateaus: Post-burn Aerial Seeding")

# eCDF plots
plot(cp.pb.aerial.seeding.psm, type = "ecdf")
bal.plot(cp.pb.aerial.seeding.psm, which = "both", type = "ecdf")
bal.plot(cp.pb.aerial.seeding.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(cp.pb.aerial.seeding.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(cp.pb.aerial.seeding.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(cp.pb.aerial.seeding.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(cp.pb.aerial.seeding.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(cp.pb.aerial.seeding.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(cp.pb.aerial.seeding.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(cp.pb.aerial.seeding.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(cp.pb.aerial.seeding.psm, type = "density")
bal.plot(cp.pb.aerial.seeding.psm, which = "both")
bal.plot(cp.pb.aerial.seeding.psm, var = "BareSoil_FH", which = "both")
bal.plot(cp.pb.aerial.seeding.psm, var = "TotalFoliarCover", which = "both")
bal.plot(cp.pb.aerial.seeding.psm, var = "ForbCover_AH", which = "both")
bal.plot(cp.pb.aerial.seeding.psm, var = "GramCover_AH", which = "both")
bal.plot(cp.pb.aerial.seeding.psm, var = "ShrubCover_AH", which = "both")
bal.plot(cp.pb.aerial.seeding.psm, var = "Gap100plus", which = "both")
bal.plot(cp.pb.aerial.seeding.psm, var = "CETWI", which = "both")
bal.plot(cp.pb.aerial.seeding.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(cp.pb.aerial.seeding.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(cp.pb.aerial.seeding.psm, type = "qq")


# Matched data
cp.pb.aerial.seeding.matched <- match_data(cp.pb.aerial.seeding.psm)

# Create trt_control variable
cp.pb.aerial.seeding.matched <- cp.pb.aerial.seeding.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Aerial Seeding", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Aerial Seeding")))

# Center and scale numeric variables
cp.pb.aerial.seeding.matched <- cp.pb.aerial.seeding.matched |>
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
cp.pb.aerial.seeding.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = cp.pb.aerial.seeding.matched,
  weights = weights
)

# G computation to estimate marginal effects
cp.pb.aerial.seeding.pred <- avg_predictions(
  model = cp.pb.aerial.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
cp.pb.aerial.seeding.pred

#   Table version
cp.pb.aerial.seeding.pred.df <- plot_predictions(
  cp.pb.aerial.seeding.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Aerial Seeding"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
cp.pb.aerial.seeding.pred.df

# Plot
cp.pb.aerial.seeding.plot <- cp.pb.aerial.seeding.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
cp.pb.aerial.seeding.plot

# Estimation of average treatment effect
cp.pb.aerial.seeding.comp <- avg_comparisons(
  model = cp.pb.aerial.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass
)
cp.pb.aerial.seeding.comp




# Middle Rockies ----------------------------------------------------------

## Herbicide --------------------------------------------------------------

# Filter data
mr.herbicide.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Middle Rockies")

# PSM
mr.herbicide.psm <- matchit(
  data = mr.herbicide.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
mr.herbicide.psm
summary(mr.herbicide.psm) # 33 treated matched


# Diagnostic love plot
love.plot(mr.herbicide.psm, stars = "std") +
  labs(title = "Middle Rockies: Herbicide")

# eCDF plots
plot(mr.herbicide.psm, type = "ecdf")
bal.plot(mr.herbicide.psm, which = "both", type = "ecdf")
bal.plot(mr.herbicide.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(mr.herbicide.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(mr.herbicide.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(mr.herbicide.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(mr.herbicide.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(mr.herbicide.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(mr.herbicide.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(mr.herbicide.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(mr.herbicide.psm, type = "density")
bal.plot(mr.herbicide.psm, which = "both")
bal.plot(mr.herbicide.psm, var = "BareSoil_FH", which = "both")
bal.plot(mr.herbicide.psm, var = "TotalFoliarCover", which = "both")
bal.plot(mr.herbicide.psm, var = "ForbCover_AH", which = "both")
bal.plot(mr.herbicide.psm, var = "GramCover_AH", which = "both")
bal.plot(mr.herbicide.psm, var = "ShrubCover_AH", which = "both")
bal.plot(mr.herbicide.psm, var = "Gap100plus", which = "both")
bal.plot(mr.herbicide.psm, var = "CETWI", which = "both")
bal.plot(mr.herbicide.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(mr.herbicide.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(mr.herbicide.psm, type = "qq")


# Matched data
mr.herbicide.matched <- match_data(mr.herbicide.psm)

# Create trt_control variable
mr.herbicide.matched <- mr.herbicide.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Herbicide")))

# Center and scale numeric variables
mr.herbicide.matched <- mr.herbicide.matched |>
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
mr.herbicide.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = mr.herbicide.matched,
  weights = weights
)

# G computation to estimate marginal effects
mr.herbicide.pred <- avg_predictions(
  model = mr.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
mr.herbicide.pred

#   Table version
mr.herbicide.pred.df <- plot_predictions(
  mr.herbicide.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Herbicide"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
mr.herbicide.pred.df

# Plot
mr.herbicide.plot <- mr.herbicide.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
mr.herbicide.plot

# Estimation of average treatment effect
mr.herbicide.comp <- avg_comparisons(
  model = mr.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass
)
mr.herbicide.comp




# Mojave Basin and Range --------------------------------------------------

## Post-burn Aerial Seeding -----------------------------------------------

# Filter data
mbr.pb.aerial.seeding.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Aerial Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Mojave Basin and Range")

# PSM
mbr.pb.aerial.seeding.psm <- matchit(
  data = mbr.pb.aerial.seeding.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
mbr.pb.aerial.seeding.psm
summary(mbr.pb.aerial.seeding.psm) # 65 treated matched


# Diagnostic love plot
love.plot(mbr.pb.aerial.seeding.psm, stars = "std") +
  labs(title = "Mojave: Aerial Seeding")

# eCDF plots
plot(mbr.pb.aerial.seeding.psm, type = "ecdf")
bal.plot(mbr.pb.aerial.seeding.psm, which = "both", type = "ecdf")
bal.plot(mbr.pb.aerial.seeding.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(mbr.pb.aerial.seeding.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(mbr.pb.aerial.seeding.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(mbr.pb.aerial.seeding.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(mbr.pb.aerial.seeding.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(mbr.pb.aerial.seeding.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(mbr.pb.aerial.seeding.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(mbr.pb.aerial.seeding.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(mbr.pb.aerial.seeding.psm, type = "density")
bal.plot(mbr.pb.aerial.seeding.psm, which = "both")
bal.plot(mbr.pb.aerial.seeding.psm, var = "BareSoil_FH", which = "both")
bal.plot(mbr.pb.aerial.seeding.psm, var = "TotalFoliarCover", which = "both")
bal.plot(mbr.pb.aerial.seeding.psm, var = "ForbCover_AH", which = "both")
bal.plot(mbr.pb.aerial.seeding.psm, var = "GramCover_AH", which = "both")
bal.plot(mbr.pb.aerial.seeding.psm, var = "ShrubCover_AH", which = "both")
bal.plot(mbr.pb.aerial.seeding.psm, var = "Gap100plus", which = "both")
bal.plot(mbr.pb.aerial.seeding.psm, var = "CETWI", which = "both")
bal.plot(mbr.pb.aerial.seeding.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(mbr.pb.aerial.seeding.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(mbr.pb.aerial.seeding.psm, type = "qq")


# Matched data
mbr.pb.aerial.seeding.matched <- match_data(mbr.pb.aerial.seeding.psm)

# Create trt_control variable
mbr.pb.aerial.seeding.matched <- mbr.pb.aerial.seeding.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Aerial Seeding", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Aerial Seeding")))

# Center and scale numeric variables
mbr.pb.aerial.seeding.matched <- mbr.pb.aerial.seeding.matched |>
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
mbr.pb.aerial.seeding.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = mbr.pb.aerial.seeding.matched,
  weights = weights
)

# G computation to estimate marginal effects
mbr.pb.aerial.seeding.pred <- avg_predictions(
  model = mbr.pb.aerial.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
mbr.pb.aerial.seeding.pred

#   Table version
mbr.pb.aerial.seeding.pred.df <- plot_predictions(
  mbr.pb.aerial.seeding.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Aerial Seeding"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
mbr.pb.aerial.seeding.pred.df

# Plot
mbr.pb.aerial.seeding.plot <- mbr.pb.aerial.seeding.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
mbr.pb.aerial.seeding.plot

# Estimation of average treatment effect
mbr.pb.aerial.seeding.comp <- avg_comparisons(
  model = mbr.pb.aerial.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass
)
mbr.pb.aerial.seeding.comp




# Northern Basin and Range ------------------------------------------------

## Drill Seeding ----------------------------------------------------------

# Filter data
nbr.drill.seeding.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Drill Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
nbr.drill.seeding.psm <- matchit(
  data = nbr.drill.seeding.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
nbr.drill.seeding.psm
summary(nbr.drill.seeding.psm) # 82 treated matched


# Diagnostic love plot
love.plot(nbr.drill.seeding.psm, stars = "std") +
  labs(title = "NBR: Drill Seeding")

# eCDF plots
plot(nbr.drill.seeding.psm, type = "ecdf")
bal.plot(nbr.drill.seeding.psm, which = "both", type = "ecdf")
bal.plot(nbr.drill.seeding.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(nbr.drill.seeding.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(nbr.drill.seeding.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.drill.seeding.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.drill.seeding.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.drill.seeding.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(nbr.drill.seeding.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(nbr.drill.seeding.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(nbr.drill.seeding.psm, type = "density")
bal.plot(nbr.drill.seeding.psm, which = "both")
bal.plot(nbr.drill.seeding.psm, var = "BareSoil_FH", which = "both")
bal.plot(nbr.drill.seeding.psm, var = "TotalFoliarCover", which = "both")
bal.plot(nbr.drill.seeding.psm, var = "ForbCover_AH", which = "both")
bal.plot(nbr.drill.seeding.psm, var = "GramCover_AH", which = "both")
bal.plot(nbr.drill.seeding.psm, var = "ShrubCover_AH", which = "both")
bal.plot(nbr.drill.seeding.psm, var = "Gap100plus", which = "both")
bal.plot(nbr.drill.seeding.psm, var = "CETWI", which = "both")
bal.plot(nbr.drill.seeding.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(nbr.drill.seeding.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(nbr.drill.seeding.psm, type = "qq")


# Matched data
nbr.drill.seeding.matched <- match_data(nbr.drill.seeding.psm)

# Create trt_control variable
nbr.drill.seeding.matched <- nbr.drill.seeding.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Drill Seeding", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Drill Seeding")))

# Center and scale numeric variables
nbr.drill.seeding.matched <- nbr.drill.seeding.matched |>
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
nbr.drill.seeding.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = nbr.drill.seeding.matched,
  weights = weights
)

# G computation to estimate marginal effects
nbr.drill.seeding.pred <- avg_predictions(
  model = nbr.drill.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
nbr.drill.seeding.pred

#   Table version
nbr.drill.seeding.pred.df <- plot_predictions(
  nbr.drill.seeding.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Drill Seeding"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
nbr.drill.seeding.pred.df

# Plot
nbr.drill.seeding.plot <- nbr.drill.seeding.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
nbr.drill.seeding.plot

# Estimation of average treatment effect
nbr.drill.seeding.comp <- avg_comparisons(
  model = nbr.drill.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass
)
nbr.drill.seeding.comp



## Drill Seeding, Soil Disturbance ----------------------------------------

# Filter data
nbr.drill.soil.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Drill Seeding, Soil Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
nbr.drill.soil.psm <- matchit(
  data = nbr.drill.soil.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
nbr.drill.soil.psm
summary(nbr.drill.soil.psm) # 61 treated matched


# Diagnostic love plot
love.plot(nbr.drill.soil.psm, stars = "std") +
  labs(title = "NBR: Drill Seeding, Soil Disturbance")

# eCDF plots
plot(nbr.drill.soil.psm, type = "ecdf")
bal.plot(nbr.drill.soil.psm, which = "both", type = "ecdf")
bal.plot(nbr.drill.soil.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(nbr.drill.soil.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(nbr.drill.soil.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.drill.soil.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.drill.soil.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.drill.soil.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(nbr.drill.soil.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(nbr.drill.soil.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(nbr.drill.soil.psm, type = "density")
bal.plot(nbr.drill.soil.psm, which = "both")
bal.plot(nbr.drill.soil.psm, var = "BareSoil_FH", which = "both")
bal.plot(nbr.drill.soil.psm, var = "TotalFoliarCover", which = "both")
bal.plot(nbr.drill.soil.psm, var = "ForbCover_AH", which = "both")
bal.plot(nbr.drill.soil.psm, var = "GramCover_AH", which = "both")
bal.plot(nbr.drill.soil.psm, var = "ShrubCover_AH", which = "both")
bal.plot(nbr.drill.soil.psm, var = "Gap100plus", which = "both")
bal.plot(nbr.drill.soil.psm, var = "CETWI", which = "both")
bal.plot(nbr.drill.soil.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(nbr.drill.soil.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(nbr.drill.soil.psm, type = "qq")


# Matched data
nbr.drill.soil.matched <- match_data(nbr.drill.soil.psm)

# Create trt_control variable
nbr.drill.soil.matched <- nbr.drill.soil.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Drill Seeding, Soil Disturbance", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Drill Seeding, Soil Disturbance")))

# Center and scale numeric variables
nbr.drill.soil.matched <- nbr.drill.soil.matched |>
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
nbr.drill.soil.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = nbr.drill.soil.matched,
  weights = weights
)

# G computation to estimate marginal effects
nbr.drill.soil.pred <- avg_predictions(
  model = nbr.drill.soil.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
nbr.drill.soil.pred

#   Table version
nbr.drill.soil.pred.df <- plot_predictions(
  nbr.drill.soil.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Drill Seeding, Soil Disturbance"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
nbr.drill.soil.pred.df

# Plot
nbr.drill.soil.plot <- nbr.drill.soil.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
nbr.drill.soil.plot

# Estimation of average treatment effect
nbr.drill.soil.comp <- avg_comparisons(
  model = nbr.drill.soil.lm,
  variables = "trt_control",
  vcov = ~subclass
)
nbr.drill.soil.comp



## Herbicide --------------------------------------------------------------

# Filter data
nbr.herbicide.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
nbr.herbicide.psm <- matchit(
  data = nbr.herbicide.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
nbr.herbicide.psm
summary(nbr.herbicide.psm) # 259 treated matched


# Diagnostic love plot
love.plot(nbr.herbicide.psm, stars = "std") +
  labs(title = "NBR: Herbicide")

# eCDF plots
plot(nbr.herbicide.psm, type = "ecdf")
bal.plot(nbr.herbicide.psm, which = "both", type = "ecdf")
bal.plot(nbr.herbicide.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(nbr.herbicide.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(nbr.herbicide.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.herbicide.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.herbicide.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.herbicide.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(nbr.herbicide.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(nbr.herbicide.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(nbr.herbicide.psm, type = "density")
bal.plot(nbr.herbicide.psm, which = "both")
bal.plot(nbr.herbicide.psm, var = "BareSoil_FH", which = "both")
bal.plot(nbr.herbicide.psm, var = "TotalFoliarCover", which = "both")
bal.plot(nbr.herbicide.psm, var = "ForbCover_AH", which = "both")
bal.plot(nbr.herbicide.psm, var = "GramCover_AH", which = "both")
bal.plot(nbr.herbicide.psm, var = "ShrubCover_AH", which = "both")
bal.plot(nbr.herbicide.psm, var = "Gap100plus", which = "both")
bal.plot(nbr.herbicide.psm, var = "CETWI", which = "both")
bal.plot(nbr.herbicide.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(nbr.herbicide.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(nbr.herbicide.psm, type = "qq")


# Matched data
nbr.herbicide.matched <- match_data(nbr.herbicide.psm)

# Create trt_control variable
nbr.herbicide.matched <- nbr.herbicide.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Herbicide")))

# Center and scale numeric variables
nbr.herbicide.matched <- nbr.herbicide.matched |>
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
nbr.herbicide.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = nbr.herbicide.matched,
  weights = weights
)

# G computation to estimate marginal effects
nbr.herbicide.pred <- avg_predictions(
  model = nbr.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
nbr.herbicide.pred

#   Table version
nbr.herbicide.pred.df <- plot_predictions(
  nbr.herbicide.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Herbicide"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
nbr.herbicide.pred.df

# Plot
nbr.herbicide.plot <- nbr.herbicide.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
nbr.herbicide.plot

# Estimation of average treatment effect
nbr.herbicide.comp <- avg_comparisons(
  model = nbr.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass
)
nbr.herbicide.comp



## Prescribed Burn --------------------------------------------------------

# Filter data
nbr.prescribed.burn.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Prescribed Burn" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
nbr.prescribed.burn.psm <- matchit(
  data = nbr.prescribed.burn.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
nbr.prescribed.burn.psm
summary(nbr.prescribed.burn.psm) # 158 treated matched


# Diagnostic love plot
love.plot(nbr.prescribed.burn.psm, stars = "std") +
  labs(title = "NBR: Prescribed Burn")

# eCDF plots
plot(nbr.prescribed.burn.psm, type = "ecdf")
bal.plot(nbr.prescribed.burn.psm, which = "both", type = "ecdf")
bal.plot(nbr.prescribed.burn.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(nbr.prescribed.burn.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(nbr.prescribed.burn.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.prescribed.burn.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.prescribed.burn.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.prescribed.burn.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(nbr.prescribed.burn.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(nbr.prescribed.burn.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(nbr.prescribed.burn.psm, type = "density")
bal.plot(nbr.prescribed.burn.psm, which = "both")
bal.plot(nbr.prescribed.burn.psm, var = "BareSoil_FH", which = "both")
bal.plot(nbr.prescribed.burn.psm, var = "TotalFoliarCover", which = "both")
bal.plot(nbr.prescribed.burn.psm, var = "ForbCover_AH", which = "both")
bal.plot(nbr.prescribed.burn.psm, var = "GramCover_AH", which = "both")
bal.plot(nbr.prescribed.burn.psm, var = "ShrubCover_AH", which = "both")
bal.plot(nbr.prescribed.burn.psm, var = "Gap100plus", which = "both")
bal.plot(nbr.prescribed.burn.psm, var = "CETWI", which = "both")
bal.plot(nbr.prescribed.burn.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(nbr.prescribed.burn.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(nbr.prescribed.burn.psm, type = "qq")


# Matched data
nbr.prescribed.burn.matched <- match_data(nbr.prescribed.burn.psm)

# Create trt_control variable
nbr.prescribed.burn.matched <- nbr.prescribed.burn.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Prescribed Burn", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Prescribed Burn")))

# Center and scale numeric variables
nbr.prescribed.burn.matched <- nbr.prescribed.burn.matched |>
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
nbr.prescribed.burn.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = nbr.prescribed.burn.matched,
  weights = weights
)

# G computation to estimate marginal effects
nbr.prescribed.burn.pred <- avg_predictions(
  model = nbr.prescribed.burn.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
nbr.prescribed.burn.pred

#   Table version
nbr.prescribed.burn.pred.df <- plot_predictions(
  nbr.prescribed.burn.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Prescribed Burn"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
nbr.prescribed.burn.pred.df

# Plot
nbr.prescribed.burn.plot <- nbr.prescribed.burn.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
nbr.prescribed.burn.plot

# Estimation of average treatment effect
nbr.prescribed.burn.comp <- avg_comparisons(
  model = nbr.prescribed.burn.lm,
  variables = "trt_control",
  vcov = ~subclass
)
nbr.prescribed.burn.comp



## Vegetation Disturbance -------------------------------------------------

# Filter data
nbr.veg.disturbance.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Vegetation Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
nbr.veg.disturbance.psm <- matchit(
  data = nbr.veg.disturbance.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
nbr.veg.disturbance.psm
summary(nbr.veg.disturbance.psm) # 105 treated matched


# Diagnostic love plot
love.plot(nbr.veg.disturbance.psm, stars = "std") +
  labs(title = "NBR: Vegetation Disturbance")

# eCDF plots
plot(nbr.veg.disturbance.psm, type = "ecdf")
bal.plot(nbr.veg.disturbance.psm, which = "both", type = "ecdf")
bal.plot(nbr.veg.disturbance.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(nbr.veg.disturbance.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(nbr.veg.disturbance.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.veg.disturbance.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.veg.disturbance.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.veg.disturbance.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(nbr.veg.disturbance.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(nbr.veg.disturbance.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(nbr.veg.disturbance.psm, type = "density")
bal.plot(nbr.veg.disturbance.psm, which = "both")
bal.plot(nbr.veg.disturbance.psm, var = "BareSoil_FH", which = "both")
bal.plot(nbr.veg.disturbance.psm, var = "TotalFoliarCover", which = "both")
bal.plot(nbr.veg.disturbance.psm, var = "ForbCover_AH", which = "both")
bal.plot(nbr.veg.disturbance.psm, var = "GramCover_AH", which = "both")
bal.plot(nbr.veg.disturbance.psm, var = "ShrubCover_AH", which = "both")
bal.plot(nbr.veg.disturbance.psm, var = "Gap100plus", which = "both")
bal.plot(nbr.veg.disturbance.psm, var = "CETWI", which = "both")
bal.plot(nbr.veg.disturbance.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(nbr.veg.disturbance.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(nbr.veg.disturbance.psm, type = "qq")


# Matched data
nbr.veg.disturbance.matched <- match_data(nbr.veg.disturbance.psm)

# Create trt_control variable
nbr.veg.disturbance.matched <- nbr.veg.disturbance.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Vegetation Disturbance", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Vegetation Disturbance")))

# Center and scale numeric variables
nbr.veg.disturbance.matched <- nbr.veg.disturbance.matched |>
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
nbr.veg.disturbance.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = nbr.veg.disturbance.matched,
  weights = weights
)

# G computation to estimate marginal effects
nbr.veg.disturbance.pred <- avg_predictions(
  model = nbr.veg.disturbance.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
nbr.veg.disturbance.pred

#   Table version
nbr.veg.disturbance.pred.df <- plot_predictions(
  nbr.veg.disturbance.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Vegetation Disturbance"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
nbr.veg.disturbance.pred.df

# Plot
nbr.veg.disturbance.plot <- nbr.veg.disturbance.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
nbr.veg.disturbance.plot

# Estimation of average treatment effect
nbr.veg.disturbance.comp <- avg_comparisons(
  model = nbr.veg.disturbance.lm,
  variables = "trt_control",
  vcov = ~subclass
)
nbr.veg.disturbance.comp



## Post-burn Aerial Seeding -----------------------------------------------

# Filter data
nbr.pb.aerial.seeding.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Aerial Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
nbr.pb.aerial.seeding.psm <- matchit(
  data = nbr.pb.aerial.seeding.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
nbr.pb.aerial.seeding.psm
summary(nbr.pb.aerial.seeding.psm) # 542 treated matched


# Diagnostic love plot
love.plot(nbr.pb.aerial.seeding.psm, stars = "std") +
  labs(title = "NBR: Post-burn Aerial Seeding")

# eCDF plots
plot(nbr.pb.aerial.seeding.psm, type = "ecdf")
bal.plot(nbr.pb.aerial.seeding.psm, which = "both", type = "ecdf")
bal.plot(nbr.pb.aerial.seeding.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(nbr.pb.aerial.seeding.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(nbr.pb.aerial.seeding.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.aerial.seeding.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.aerial.seeding.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.aerial.seeding.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(nbr.pb.aerial.seeding.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(nbr.pb.aerial.seeding.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(nbr.pb.aerial.seeding.psm, type = "density")
bal.plot(nbr.pb.aerial.seeding.psm, which = "both")
bal.plot(nbr.pb.aerial.seeding.psm, var = "BareSoil_FH", which = "both")
bal.plot(nbr.pb.aerial.seeding.psm, var = "TotalFoliarCover", which = "both")
bal.plot(nbr.pb.aerial.seeding.psm, var = "ForbCover_AH", which = "both")
bal.plot(nbr.pb.aerial.seeding.psm, var = "GramCover_AH", which = "both")
bal.plot(nbr.pb.aerial.seeding.psm, var = "ShrubCover_AH", which = "both")
bal.plot(nbr.pb.aerial.seeding.psm, var = "Gap100plus", which = "both")
bal.plot(nbr.pb.aerial.seeding.psm, var = "CETWI", which = "both")
bal.plot(nbr.pb.aerial.seeding.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(nbr.pb.aerial.seeding.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(nbr.pb.aerial.seeding.psm, type = "qq")


# Matched data
nbr.pb.aerial.seeding.matched <- match_data(nbr.pb.aerial.seeding.psm)

# Create trt_control variable
nbr.pb.aerial.seeding.matched <- nbr.pb.aerial.seeding.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Aerial Seeding", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Aerial Seeding")))

# Center and scale numeric variables
nbr.pb.aerial.seeding.matched <- nbr.pb.aerial.seeding.matched |>
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
nbr.pb.aerial.seeding.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = nbr.pb.aerial.seeding.matched,
  weights = weights
)

# G computation to estimate marginal effects
nbr.pb.aerial.seeding.pred <- avg_predictions(
  model = nbr.pb.aerial.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
nbr.pb.aerial.seeding.pred

#   Table version
nbr.pb.aerial.seeding.pred.df <- plot_predictions(
  nbr.pb.aerial.seeding.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Aerial Seeding"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
nbr.pb.aerial.seeding.pred.df

# Plot
nbr.pb.aerial.seeding.plot <- nbr.pb.aerial.seeding.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
nbr.pb.aerial.seeding.plot

# Estimation of average treatment effect
nbr.pb.aerial.seeding.comp <- avg_comparisons(
  model = nbr.pb.aerial.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass
)
nbr.pb.aerial.seeding.comp



## Post-burn Aerial Seeding, Drill Seeding --------------------------------

# Filter data
nbr.pb.aerial.drill.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Aerial Seeding, Drill Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
nbr.pb.aerial.drill.psm <- matchit(
  data = nbr.pb.aerial.drill.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
nbr.pb.aerial.drill.psm
summary(nbr.pb.aerial.drill.psm) # 88 treated matched


# Diagnostic love plot
love.plot(nbr.pb.aerial.drill.psm, stars = "std") +
  labs(title = "NBR: Post-burn Aerial & Drill")

# eCDF plots
plot(nbr.pb.aerial.drill.psm, type = "ecdf")
bal.plot(nbr.pb.aerial.drill.psm, which = "both", type = "ecdf")
bal.plot(nbr.pb.aerial.drill.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(nbr.pb.aerial.drill.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(nbr.pb.aerial.drill.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.aerial.drill.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.aerial.drill.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.aerial.drill.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(nbr.pb.aerial.drill.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(nbr.pb.aerial.drill.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(nbr.pb.aerial.drill.psm, type = "density")
bal.plot(nbr.pb.aerial.drill.psm, which = "both")
bal.plot(nbr.pb.aerial.drill.psm, var = "BareSoil_FH", which = "both")
bal.plot(nbr.pb.aerial.drill.psm, var = "TotalFoliarCover", which = "both")
bal.plot(nbr.pb.aerial.drill.psm, var = "ForbCover_AH", which = "both")
bal.plot(nbr.pb.aerial.drill.psm, var = "GramCover_AH", which = "both")
bal.plot(nbr.pb.aerial.drill.psm, var = "ShrubCover_AH", which = "both")
bal.plot(nbr.pb.aerial.drill.psm, var = "Gap100plus", which = "both")
bal.plot(nbr.pb.aerial.drill.psm, var = "CETWI", which = "both")
bal.plot(nbr.pb.aerial.drill.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(nbr.pb.aerial.drill.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(nbr.pb.aerial.drill.psm, type = "qq")


# Matched data
nbr.pb.aerial.drill.matched <- match_data(nbr.pb.aerial.drill.psm)

# Create trt_control variable
nbr.pb.aerial.drill.matched <- nbr.pb.aerial.drill.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Aerial Seeding, Drill Seeding", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Aerial Seeding, Drill Seeding")))

# Center and scale numeric variables
nbr.pb.aerial.drill.matched <- nbr.pb.aerial.drill.matched |>
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
nbr.pb.aerial.drill.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = nbr.pb.aerial.drill.matched,
  weights = weights
)

# G computation to estimate marginal effects
nbr.pb.aerial.drill.pred <- avg_predictions(
  model = nbr.pb.aerial.drill.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
nbr.pb.aerial.drill.pred

#   Table version
nbr.pb.aerial.drill.pred.df <- plot_predictions(
  nbr.pb.aerial.drill.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Aerial Seeding, Drill Seeding"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
nbr.pb.aerial.drill.pred.df

# Plot
nbr.pb.aerial.drill.plot <- nbr.pb.aerial.drill.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
nbr.pb.aerial.drill.plot

# Estimation of average treatment effect
nbr.pb.aerial.drill.comp <- avg_comparisons(
  model = nbr.pb.aerial.drill.lm,
  variables = "trt_control",
  vcov = ~subclass
)
nbr.pb.aerial.drill.comp



## Post-burn Closure ------------------------------------------------------

# Filter data
nbr.pb.closure.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Closure" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
nbr.pb.closure.psm <- matchit(
  data = nbr.pb.closure.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
nbr.pb.closure.psm
summary(nbr.pb.closure.psm) # 151 treated matched


# Diagnostic love plot
love.plot(nbr.pb.closure.psm, stars = "std") +
  labs(title = "NBR: Post-burn Closure")

# eCDF plots
plot(nbr.pb.closure.psm, type = "ecdf")
bal.plot(nbr.pb.closure.psm, which = "both", type = "ecdf")
bal.plot(nbr.pb.closure.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(nbr.pb.closure.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(nbr.pb.closure.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.closure.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.closure.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.closure.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(nbr.pb.closure.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(nbr.pb.closure.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(nbr.pb.closure.psm, type = "density")
bal.plot(nbr.pb.closure.psm, which = "both")
bal.plot(nbr.pb.closure.psm, var = "BareSoil_FH", which = "both")
bal.plot(nbr.pb.closure.psm, var = "TotalFoliarCover", which = "both")
bal.plot(nbr.pb.closure.psm, var = "ForbCover_AH", which = "both")
bal.plot(nbr.pb.closure.psm, var = "GramCover_AH", which = "both")
bal.plot(nbr.pb.closure.psm, var = "ShrubCover_AH", which = "both")
bal.plot(nbr.pb.closure.psm, var = "Gap100plus", which = "both")
bal.plot(nbr.pb.closure.psm, var = "CETWI", which = "both")
bal.plot(nbr.pb.closure.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(nbr.pb.closure.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(nbr.pb.closure.psm, type = "qq")


# Matched data
nbr.pb.closure.matched <- match_data(nbr.pb.closure.psm)

# Create trt_control variable
nbr.pb.closure.matched <- nbr.pb.closure.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Closure", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Closure")))

# Center and scale numeric variables
nbr.pb.closure.matched <- nbr.pb.closure.matched |>
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
nbr.pb.closure.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = nbr.pb.closure.matched,
  weights = weights
)

# G computation to estimate marginal effects
nbr.pb.closure.pred <- avg_predictions(
  model = nbr.pb.closure.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
nbr.pb.closure.pred

#   Table version
nbr.pb.closure.pred.df <- plot_predictions(
  nbr.pb.closure.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Closure"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
nbr.pb.closure.pred.df

# Plot
nbr.pb.closure.plot <- nbr.pb.closure.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
nbr.pb.closure.plot

# Estimation of average treatment effect
nbr.pb.closure.comp <- avg_comparisons(
  model = nbr.pb.closure.lm,
  variables = "trt_control",
  vcov = ~subclass
)
nbr.pb.closure.comp



## Post-burn Drill Seeding ------------------------------------------------------

# Filter data
nbr.pb.drill.seeding.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Drill Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
nbr.pb.drill.seeding.psm <- matchit(
  data = nbr.pb.drill.seeding.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
nbr.pb.drill.seeding.psm
summary(nbr.pb.drill.seeding.psm) # 217 treated matched


# Diagnostic love plot
love.plot(nbr.pb.drill.seeding.psm, stars = "std") +
  labs(title = "NBR: Post-burn Drill Seeding")

# eCDF plots
plot(nbr.pb.drill.seeding.psm, type = "ecdf")
bal.plot(nbr.pb.drill.seeding.psm, which = "both", type = "ecdf")
bal.plot(nbr.pb.drill.seeding.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(nbr.pb.drill.seeding.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(nbr.pb.drill.seeding.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.drill.seeding.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.drill.seeding.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.drill.seeding.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(nbr.pb.drill.seeding.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(nbr.pb.drill.seeding.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(nbr.pb.drill.seeding.psm, type = "density")
bal.plot(nbr.pb.drill.seeding.psm, which = "both")
bal.plot(nbr.pb.drill.seeding.psm, var = "BareSoil_FH", which = "both")
bal.plot(nbr.pb.drill.seeding.psm, var = "TotalFoliarCover", which = "both")
bal.plot(nbr.pb.drill.seeding.psm, var = "ForbCover_AH", which = "both")
bal.plot(nbr.pb.drill.seeding.psm, var = "GramCover_AH", which = "both")
bal.plot(nbr.pb.drill.seeding.psm, var = "ShrubCover_AH", which = "both")
bal.plot(nbr.pb.drill.seeding.psm, var = "Gap100plus", which = "both")
bal.plot(nbr.pb.drill.seeding.psm, var = "CETWI", which = "both")
bal.plot(nbr.pb.drill.seeding.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(nbr.pb.drill.seeding.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(nbr.pb.drill.seeding.psm, type = "qq")


# Matched data
nbr.pb.drill.seeding.matched <- match_data(nbr.pb.drill.seeding.psm)

# Create trt_control variable
nbr.pb.drill.seeding.matched <- nbr.pb.drill.seeding.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Drill Seeding", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Drill Seeding")))

# Center and scale numeric variables
nbr.pb.drill.seeding.matched <- nbr.pb.drill.seeding.matched |>
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
nbr.pb.drill.seeding.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = nbr.pb.drill.seeding.matched,
  weights = weights
)

# G computation to estimate marginal effects
nbr.pb.drill.seeding.pred <- avg_predictions(
  model = nbr.pb.drill.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
nbr.pb.drill.seeding.pred

#   Table version
nbr.pb.drill.seeding.pred.df <- plot_predictions(
  nbr.pb.drill.seeding.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Drill Seeding"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
nbr.pb.drill.seeding.pred.df

# Plot
nbr.pb.drill.seeding.plot <- nbr.pb.drill.seeding.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
nbr.pb.drill.seeding.plot

# Estimation of average treatment effect
nbr.pb.drill.seeding.comp <- avg_comparisons(
  model = nbr.pb.drill.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass
)
nbr.pb.drill.seeding.comp



## Post-burn Herbicide ----------------------------------------------------

# Filter data
nbr.pb.herbicide.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
nbr.pb.herbicide.psm <- matchit(
  data = nbr.pb.herbicide.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
nbr.pb.herbicide.psm
summary(nbr.pb.herbicide.psm) # 370 treated matched


# Diagnostic love plot
love.plot(nbr.pb.herbicide.psm, stars = "std") +
  labs(title = "NBR: Post-burn Herbicide")

# eCDF plots
plot(nbr.pb.herbicide.psm, type = "ecdf")
bal.plot(nbr.pb.herbicide.psm, which = "both", type = "ecdf")
bal.plot(nbr.pb.herbicide.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(nbr.pb.herbicide.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(nbr.pb.herbicide.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.herbicide.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.herbicide.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.herbicide.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(nbr.pb.herbicide.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(nbr.pb.herbicide.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(nbr.pb.herbicide.psm, type = "density")
bal.plot(nbr.pb.herbicide.psm, which = "both")
bal.plot(nbr.pb.herbicide.psm, var = "BareSoil_FH", which = "both")
bal.plot(nbr.pb.herbicide.psm, var = "TotalFoliarCover", which = "both")
bal.plot(nbr.pb.herbicide.psm, var = "ForbCover_AH", which = "both")
bal.plot(nbr.pb.herbicide.psm, var = "GramCover_AH", which = "both")
bal.plot(nbr.pb.herbicide.psm, var = "ShrubCover_AH", which = "both")
bal.plot(nbr.pb.herbicide.psm, var = "Gap100plus", which = "both")
bal.plot(nbr.pb.herbicide.psm, var = "CETWI", which = "both")
bal.plot(nbr.pb.herbicide.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(nbr.pb.herbicide.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(nbr.pb.herbicide.psm, type = "qq")


# Matched data
nbr.pb.herbicide.matched <- match_data(nbr.pb.herbicide.psm)

# Create trt_control variable
nbr.pb.herbicide.matched <- nbr.pb.herbicide.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Herbicide")))

# Center and scale numeric variables
nbr.pb.herbicide.matched <- nbr.pb.herbicide.matched |>
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
nbr.pb.herbicide.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = nbr.pb.herbicide.matched,
  weights = weights
)

# G computation to estimate marginal effects
nbr.pb.herbicide.pred <- avg_predictions(
  model = nbr.pb.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
nbr.pb.herbicide.pred

#   Table version
nbr.pb.herbicide.pred.df <- plot_predictions(
  nbr.pb.herbicide.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Herbicide"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
nbr.pb.herbicide.pred.df

# Plot
nbr.pb.herbicide.plot <- nbr.pb.herbicide.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
nbr.pb.herbicide.plot

# Estimation of average treatment effect
nbr.pb.herbicide.comp <- avg_comparisons(
  model = nbr.pb.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass
)
nbr.pb.herbicide.comp



## Post-burn Seedling Planting --------------------------------------------

# Filter data
nbr.pb.seedling.planting.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Seedling Planting" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Northern Basin and Range")

# PSM
nbr.pb.seedling.planting.psm <- matchit(
  data = nbr.pb.seedling.planting.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
nbr.pb.seedling.planting.psm
summary(nbr.pb.seedling.planting.psm) # 56 treated matched


# Diagnostic love plot
love.plot(nbr.pb.seedling.planting.psm, stars = "std") +
  labs(title = "NBR: Post-burn Seedling Planting")

# eCDF plots
plot(nbr.pb.seedling.planting.psm, type = "ecdf")
bal.plot(nbr.pb.seedling.planting.psm, which = "both", type = "ecdf")
bal.plot(nbr.pb.seedling.planting.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(nbr.pb.seedling.planting.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(nbr.pb.seedling.planting.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.seedling.planting.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.seedling.planting.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(nbr.pb.seedling.planting.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(nbr.pb.seedling.planting.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(nbr.pb.seedling.planting.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(nbr.pb.seedling.planting.psm, type = "density")
bal.plot(nbr.pb.seedling.planting.psm, which = "both")
bal.plot(nbr.pb.seedling.planting.psm, var = "BareSoil_FH", which = "both")
bal.plot(nbr.pb.seedling.planting.psm, var = "TotalFoliarCover", which = "both")
bal.plot(nbr.pb.seedling.planting.psm, var = "ForbCover_AH", which = "both")
bal.plot(nbr.pb.seedling.planting.psm, var = "GramCover_AH", which = "both")
bal.plot(nbr.pb.seedling.planting.psm, var = "ShrubCover_AH", which = "both")
bal.plot(nbr.pb.seedling.planting.psm, var = "Gap100plus", which = "both")
bal.plot(nbr.pb.seedling.planting.psm, var = "CETWI", which = "both")
bal.plot(nbr.pb.seedling.planting.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(nbr.pb.seedling.planting.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(nbr.pb.seedling.planting.psm, type = "qq")


# Matched data
nbr.pb.seedling.planting.matched <- match_data(nbr.pb.seedling.planting.psm)

# Create trt_control variable
nbr.pb.seedling.planting.matched <- nbr.pb.seedling.planting.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Seedling Planting", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Seedling Planting")))

# Center and scale numeric variables
nbr.pb.seedling.planting.matched <- nbr.pb.seedling.planting.matched |>
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
nbr.pb.seedling.planting.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = nbr.pb.seedling.planting.matched,
  weights = weights
)

# G computation to estimate marginal effects
nbr.pb.seedling.planting.pred <- avg_predictions(
  model = nbr.pb.seedling.planting.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
nbr.pb.seedling.planting.pred

#   Table version
nbr.pb.seedling.planting.pred.df <- plot_predictions(
  nbr.pb.seedling.planting.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Seedling Planting"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
nbr.pb.seedling.planting.pred.df

# Plot
nbr.pb.seedling.planting.plot <- nbr.pb.seedling.planting.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
nbr.pb.seedling.planting.plot

# Estimation of average treatment effect
nbr.pb.seedling.planting.comp <- avg_comparisons(
  model = nbr.pb.seedling.planting.lm,
  variables = "trt_control",
  vcov = ~subclass
)
nbr.pb.seedling.planting.comp




# Northwestern Great Plains -----------------------------------------------

## Prescribed Burn --------------------------------------------------------

# Filter data
ngp.prescribed.burn.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Prescribed Burn" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Northwestern Great Plains")

# PSM
ngp.prescribed.burn.psm <- matchit(
  data = ngp.prescribed.burn.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
ngp.prescribed.burn.psm
summary(ngp.prescribed.burn.psm) # 51 treated matched


# Diagnostic love plot
love.plot(ngp.prescribed.burn.psm, stars = "std") +
  labs(title = "NW Great Plains: Prescribed Burn")

# eCDF plots
plot(ngp.prescribed.burn.psm, type = "ecdf")
bal.plot(ngp.prescribed.burn.psm, which = "both", type = "ecdf")
bal.plot(ngp.prescribed.burn.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(ngp.prescribed.burn.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(ngp.prescribed.burn.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(ngp.prescribed.burn.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(ngp.prescribed.burn.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(ngp.prescribed.burn.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(ngp.prescribed.burn.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(ngp.prescribed.burn.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(ngp.prescribed.burn.psm, type = "density")
bal.plot(ngp.prescribed.burn.psm, which = "both")
bal.plot(ngp.prescribed.burn.psm, var = "BareSoil_FH", which = "both")
bal.plot(ngp.prescribed.burn.psm, var = "TotalFoliarCover", which = "both")
bal.plot(ngp.prescribed.burn.psm, var = "ForbCover_AH", which = "both")
bal.plot(ngp.prescribed.burn.psm, var = "GramCover_AH", which = "both")
bal.plot(ngp.prescribed.burn.psm, var = "ShrubCover_AH", which = "both")
bal.plot(ngp.prescribed.burn.psm, var = "Gap100plus", which = "both")
bal.plot(ngp.prescribed.burn.psm, var = "CETWI", which = "both")
bal.plot(ngp.prescribed.burn.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(ngp.prescribed.burn.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(ngp.prescribed.burn.psm, type = "qq")


# Matched data
ngp.prescribed.burn.matched <- match_data(ngp.prescribed.burn.psm)

# Create trt_control variable
ngp.prescribed.burn.matched <- ngp.prescribed.burn.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Prescribed Burn", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Prescribed Burn")))

# Center and scale numeric variables
ngp.prescribed.burn.matched <- ngp.prescribed.burn.matched |>
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
ngp.prescribed.burn.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = ngp.prescribed.burn.matched,
  weights = weights
)

# G computation to estimate marginal effects
ngp.prescribed.burn.pred <- avg_predictions(
  model = ngp.prescribed.burn.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
ngp.prescribed.burn.pred

#   Table version
ngp.prescribed.burn.pred.df <- plot_predictions(
  ngp.prescribed.burn.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Prescribed Burn"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
ngp.prescribed.burn.pred.df

# Plot
ngp.prescribed.burn.plot <- ngp.prescribed.burn.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
ngp.prescribed.burn.plot

# Estimation of average treatment effect
ngp.prescribed.burn.comp <- avg_comparisons(
  model = ngp.prescribed.burn.lm,
  variables = "trt_control",
  vcov = ~subclass
)
ngp.prescribed.burn.comp




# Snake River Plain -------------------------------------------------------

## Post-burn Aerial Seeding -----------------------------------------------

# Filter data
srp.pb.aerial.seeding.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Aerial Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Snake River Plain")

# PSM
srp.pb.aerial.seeding.psm <- matchit(
  data = srp.pb.aerial.seeding.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
srp.pb.aerial.seeding.psm
summary(srp.pb.aerial.seeding.psm) # 153 treated matched


# Diagnostic love plot
love.plot(srp.pb.aerial.seeding.psm, stars = "std") +
  labs(title = "Snake River: Post-burn Aerial Seeding")

# eCDF plots
plot(srp.pb.aerial.seeding.psm, type = "ecdf")
bal.plot(srp.pb.aerial.seeding.psm, which = "both", type = "ecdf")
bal.plot(srp.pb.aerial.seeding.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(srp.pb.aerial.seeding.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(srp.pb.aerial.seeding.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(srp.pb.aerial.seeding.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(srp.pb.aerial.seeding.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(srp.pb.aerial.seeding.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(srp.pb.aerial.seeding.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(srp.pb.aerial.seeding.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(srp.pb.aerial.seeding.psm, type = "density")
bal.plot(srp.pb.aerial.seeding.psm, which = "both")
bal.plot(srp.pb.aerial.seeding.psm, var = "BareSoil_FH", which = "both")
bal.plot(srp.pb.aerial.seeding.psm, var = "TotalFoliarCover", which = "both")
bal.plot(srp.pb.aerial.seeding.psm, var = "ForbCover_AH", which = "both")
bal.plot(srp.pb.aerial.seeding.psm, var = "GramCover_AH", which = "both")
bal.plot(srp.pb.aerial.seeding.psm, var = "ShrubCover_AH", which = "both")
bal.plot(srp.pb.aerial.seeding.psm, var = "Gap100plus", which = "both")
bal.plot(srp.pb.aerial.seeding.psm, var = "CETWI", which = "both")
bal.plot(srp.pb.aerial.seeding.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(srp.pb.aerial.seeding.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(srp.pb.aerial.seeding.psm, type = "qq")


# Matched data
srp.pb.aerial.seeding.matched <- match_data(srp.pb.aerial.seeding.psm)

# Create trt_control variable
srp.pb.aerial.seeding.matched <- srp.pb.aerial.seeding.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Aerial Seeding", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Aerial Seeding")))

# Center and scale numeric variables
srp.pb.aerial.seeding.matched <- srp.pb.aerial.seeding.matched |>
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
srp.pb.aerial.seeding.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = srp.pb.aerial.seeding.matched,
  weights = weights
)

# G computation to estimate marginal effects
srp.pb.aerial.seeding.pred <- avg_predictions(
  model = srp.pb.aerial.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
srp.pb.aerial.seeding.pred

#   Table version
srp.pb.aerial.seeding.pred.df <- plot_predictions(
  srp.pb.aerial.seeding.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Aerial Seeding"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
srp.pb.aerial.seeding.pred.df

# Plot
srp.pb.aerial.seeding.plot <- srp.pb.aerial.seeding.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
srp.pb.aerial.seeding.plot

# Estimation of average treatment effect
srp.pb.aerial.seeding.comp <- avg_comparisons(
  model = srp.pb.aerial.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass
)
srp.pb.aerial.seeding.comp



## Post-burn Aerial Seeding, Drill Seeding --------------------------------

# Filter data
srp.pb.aerial.drill.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Aerial Seeding, Drill Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Snake River Plain")

# PSM
srp.pb.aerial.drill.psm <- matchit(
  data = srp.pb.aerial.drill.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
srp.pb.aerial.drill.psm
summary(srp.pb.aerial.drill.psm) # 77 treated matched


# Diagnostic love plot
love.plot(srp.pb.aerial.drill.psm, stars = "std") +
  labs(title = "Snake River: Post-burn Aerial & Drill")

# eCDF plots
plot(srp.pb.aerial.drill.psm, type = "ecdf")
bal.plot(srp.pb.aerial.drill.psm, which = "both", type = "ecdf")
bal.plot(srp.pb.aerial.drill.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(srp.pb.aerial.drill.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(srp.pb.aerial.drill.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(srp.pb.aerial.drill.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(srp.pb.aerial.drill.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(srp.pb.aerial.drill.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(srp.pb.aerial.drill.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(srp.pb.aerial.drill.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(srp.pb.aerial.drill.psm, type = "density")
bal.plot(srp.pb.aerial.drill.psm, which = "both")
bal.plot(srp.pb.aerial.drill.psm, var = "BareSoil_FH", which = "both")
bal.plot(srp.pb.aerial.drill.psm, var = "TotalFoliarCover", which = "both")
bal.plot(srp.pb.aerial.drill.psm, var = "ForbCover_AH", which = "both")
bal.plot(srp.pb.aerial.drill.psm, var = "GramCover_AH", which = "both")
bal.plot(srp.pb.aerial.drill.psm, var = "ShrubCover_AH", which = "both")
bal.plot(srp.pb.aerial.drill.psm, var = "Gap100plus", which = "both")
bal.plot(srp.pb.aerial.drill.psm, var = "CETWI", which = "both")
bal.plot(srp.pb.aerial.drill.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(srp.pb.aerial.drill.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(srp.pb.aerial.drill.psm, type = "qq")


# Matched data
srp.pb.aerial.drill.matched <- match_data(srp.pb.aerial.drill.psm)

# Create trt_control variable
srp.pb.aerial.drill.matched <- srp.pb.aerial.drill.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Aerial Seeding, Drill Seeding", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Aerial Seeding, Drill Seeding")))

# Center and scale numeric variables
srp.pb.aerial.drill.matched <- srp.pb.aerial.drill.matched |>
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
srp.pb.aerial.drill.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = srp.pb.aerial.drill.matched,
  weights = weights
)

# G computation to estimate marginal effects
srp.pb.aerial.drill.pred <- avg_predictions(
  model = srp.pb.aerial.drill.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
srp.pb.aerial.drill.pred

#   Table version
srp.pb.aerial.drill.pred.df <- plot_predictions(
  srp.pb.aerial.drill.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Aerial Seeding, Drill Seeding"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
srp.pb.aerial.drill.pred.df

# Plot
srp.pb.aerial.drill.plot <- srp.pb.aerial.drill.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
srp.pb.aerial.drill.plot

# Estimation of average treatment effect
srp.pb.aerial.drill.comp <- avg_comparisons(
  model = srp.pb.aerial.drill.lm,
  variables = "trt_control",
  vcov = ~subclass
)
srp.pb.aerial.drill.comp



## Post-burn Closure ------------------------------------------------------

# Filter data
srp.pb.closure.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Closure" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Snake River Plain")

# PSM
srp.pb.closure.psm <- matchit(
  data = srp.pb.closure.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
srp.pb.closure.psm
summary(srp.pb.closure.psm) # 87 treated matched


# Diagnostic love plot
love.plot(srp.pb.closure.psm, stars = "std") +
  labs(title = "Snake River: Post-burn Closure")

# eCDF plots
plot(srp.pb.closure.psm, type = "ecdf")
bal.plot(srp.pb.closure.psm, which = "both", type = "ecdf")
bal.plot(srp.pb.closure.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(srp.pb.closure.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(srp.pb.closure.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(srp.pb.closure.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(srp.pb.closure.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(srp.pb.closure.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(srp.pb.closure.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(srp.pb.closure.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(srp.pb.closure.psm, type = "density")
bal.plot(srp.pb.closure.psm, which = "both")
bal.plot(srp.pb.closure.psm, var = "BareSoil_FH", which = "both")
bal.plot(srp.pb.closure.psm, var = "TotalFoliarCover", which = "both")
bal.plot(srp.pb.closure.psm, var = "ForbCover_AH", which = "both")
bal.plot(srp.pb.closure.psm, var = "GramCover_AH", which = "both")
bal.plot(srp.pb.closure.psm, var = "ShrubCover_AH", which = "both")
bal.plot(srp.pb.closure.psm, var = "Gap100plus", which = "both")
bal.plot(srp.pb.closure.psm, var = "CETWI", which = "both")
bal.plot(srp.pb.closure.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(srp.pb.closure.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(srp.pb.closure.psm, type = "qq")


# Matched data
srp.pb.closure.matched <- match_data(srp.pb.closure.psm)

# Create trt_control variable
srp.pb.closure.matched <- srp.pb.closure.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Closure", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Closure")))

# Center and scale numeric variables
srp.pb.closure.matched <- srp.pb.closure.matched |>
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
srp.pb.closure.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = srp.pb.closure.matched,
  weights = weights
)

# G computation to estimate marginal effects
srp.pb.closure.pred <- avg_predictions(
  model = srp.pb.closure.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
srp.pb.closure.pred

#   Table version
srp.pb.closure.pred.df <- plot_predictions(
  srp.pb.closure.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Closure"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
srp.pb.closure.pred.df

# Plot
srp.pb.closure.plot <- srp.pb.closure.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
srp.pb.closure.plot

# Estimation of average treatment effect
srp.pb.closure.comp <- avg_comparisons(
  model = srp.pb.closure.lm,
  variables = "trt_control",
  vcov = ~subclass
)
srp.pb.closure.comp



## Post-burn Drill Seeding ------------------------------------------------------

# Filter data
srp.pb.drill.seeding.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Drill Seeding" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Snake River Plain")

# PSM
srp.pb.drill.seeding.psm <- matchit(
  data = srp.pb.drill.seeding.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
srp.pb.drill.seeding.psm
summary(srp.pb.drill.seeding.psm) # 36 treated matched


# Diagnostic love plot
love.plot(srp.pb.drill.seeding.psm, stars = "std") +
  labs(title = "Snake River: Post-burn Drill Seeding")

# eCDF plots
plot(srp.pb.drill.seeding.psm, type = "ecdf")
bal.plot(srp.pb.drill.seeding.psm, which = "both", type = "ecdf")
bal.plot(srp.pb.drill.seeding.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(srp.pb.drill.seeding.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(srp.pb.drill.seeding.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(srp.pb.drill.seeding.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(srp.pb.drill.seeding.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(srp.pb.drill.seeding.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(srp.pb.drill.seeding.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(srp.pb.drill.seeding.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(srp.pb.drill.seeding.psm, type = "density")
bal.plot(srp.pb.drill.seeding.psm, which = "both")
bal.plot(srp.pb.drill.seeding.psm, var = "BareSoil_FH", which = "both")
bal.plot(srp.pb.drill.seeding.psm, var = "TotalFoliarCover", which = "both")
bal.plot(srp.pb.drill.seeding.psm, var = "ForbCover_AH", which = "both")
bal.plot(srp.pb.drill.seeding.psm, var = "GramCover_AH", which = "both")
bal.plot(srp.pb.drill.seeding.psm, var = "ShrubCover_AH", which = "both")
bal.plot(srp.pb.drill.seeding.psm, var = "Gap100plus", which = "both")
bal.plot(srp.pb.drill.seeding.psm, var = "CETWI", which = "both")
bal.plot(srp.pb.drill.seeding.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(srp.pb.drill.seeding.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(srp.pb.drill.seeding.psm, type = "qq")


# Matched data
srp.pb.drill.seeding.matched <- match_data(srp.pb.drill.seeding.psm)

# Create trt_control variable
srp.pb.drill.seeding.matched <- srp.pb.drill.seeding.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Drill Seeding", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Drill Seeding")))

# Center and scale numeric variables
srp.pb.drill.seeding.matched <- srp.pb.drill.seeding.matched |>
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
srp.pb.drill.seeding.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = srp.pb.drill.seeding.matched,
  weights = weights
)

# G computation to estimate marginal effects
srp.pb.drill.seeding.pred <- avg_predictions(
  model = srp.pb.drill.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
srp.pb.drill.seeding.pred

#   Table version
srp.pb.drill.seeding.pred.df <- plot_predictions(
  srp.pb.drill.seeding.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Drill Seeding"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
srp.pb.drill.seeding.pred.df

# Plot
srp.pb.drill.seeding.plot <- srp.pb.drill.seeding.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
srp.pb.drill.seeding.plot

# Estimation of average treatment effect
srp.pb.drill.seeding.comp <- avg_comparisons(
  model = srp.pb.drill.seeding.lm,
  variables = "trt_control",
  vcov = ~subclass
)
srp.pb.drill.seeding.comp



## Post-burn Herbicide ----------------------------------------------------

# Filter data
srp.pb.herbicide.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "post")) |>
  filter(EcoLvl3 == "Snake River Plain")

# PSM
srp.pb.herbicide.psm <- matchit(
  data = srp.pb.herbicide.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
srp.pb.herbicide.psm
summary(srp.pb.herbicide.psm) # 28 treated matched


# Diagnostic love plot
love.plot(srp.pb.herbicide.psm, stars = "std") +
  labs(title = "Snake River: Post-burn Herbicide")

# eCDF plots
plot(srp.pb.herbicide.psm, type = "ecdf")
bal.plot(srp.pb.herbicide.psm, which = "both", type = "ecdf")
bal.plot(srp.pb.herbicide.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(srp.pb.herbicide.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(srp.pb.herbicide.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(srp.pb.herbicide.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(srp.pb.herbicide.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(srp.pb.herbicide.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(srp.pb.herbicide.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(srp.pb.herbicide.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(srp.pb.herbicide.psm, type = "density")
bal.plot(srp.pb.herbicide.psm, which = "both")
bal.plot(srp.pb.herbicide.psm, var = "BareSoil_FH", which = "both")
bal.plot(srp.pb.herbicide.psm, var = "TotalFoliarCover", which = "both")
bal.plot(srp.pb.herbicide.psm, var = "ForbCover_AH", which = "both")
bal.plot(srp.pb.herbicide.psm, var = "GramCover_AH", which = "both")
bal.plot(srp.pb.herbicide.psm, var = "ShrubCover_AH", which = "both")
bal.plot(srp.pb.herbicide.psm, var = "Gap100plus", which = "both")
bal.plot(srp.pb.herbicide.psm, var = "CETWI", which = "both")
bal.plot(srp.pb.herbicide.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(srp.pb.herbicide.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(srp.pb.herbicide.psm, type = "qq")


# Matched data
srp.pb.herbicide.matched <- match_data(srp.pb.herbicide.psm)

# Create trt_control variable
srp.pb.herbicide.matched <- srp.pb.herbicide.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Herbicide")))

# Center and scale numeric variables
srp.pb.herbicide.matched <- srp.pb.herbicide.matched |>
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
srp.pb.herbicide.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = srp.pb.herbicide.matched,
  weights = weights
)

# G computation to estimate marginal effects
srp.pb.herbicide.pred <- avg_predictions(
  model = srp.pb.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
srp.pb.herbicide.pred

#   Table version
srp.pb.herbicide.pred.df <- plot_predictions(
  srp.pb.herbicide.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Herbicide"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
srp.pb.herbicide.pred.df

# Plot
srp.pb.herbicide.plot <- srp.pb.herbicide.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
srp.pb.herbicide.plot

# Estimation of average treatment effect
srp.pb.herbicide.comp <- avg_comparisons(
  model = srp.pb.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass
)
srp.pb.herbicide.comp




# Southern Rockies --------------------------------------------------------

## Herbicide --------------------------------------------------------------

# Filter data
sr.herbicide.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Herbicide" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Southern Rockies")

# PSM
sr.herbicide.psm <- matchit(
  data = sr.herbicide.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
sr.herbicide.psm
summary(sr.herbicide.psm) # 76 treated matched


# Diagnostic love plot
love.plot(sr.herbicide.psm, stars = "std") +
  labs(title = "Southern Rockies: Herbicide")

# eCDF plots
plot(sr.herbicide.psm, type = "ecdf")
bal.plot(sr.herbicide.psm, which = "both", type = "ecdf")
bal.plot(sr.herbicide.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(sr.herbicide.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(sr.herbicide.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(sr.herbicide.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(sr.herbicide.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(sr.herbicide.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(sr.herbicide.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(sr.herbicide.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(sr.herbicide.psm, type = "density")
bal.plot(sr.herbicide.psm, which = "both")
bal.plot(sr.herbicide.psm, var = "BareSoil_FH", which = "both")
bal.plot(sr.herbicide.psm, var = "TotalFoliarCover", which = "both")
bal.plot(sr.herbicide.psm, var = "ForbCover_AH", which = "both")
bal.plot(sr.herbicide.psm, var = "GramCover_AH", which = "both")
bal.plot(sr.herbicide.psm, var = "ShrubCover_AH", which = "both")
bal.plot(sr.herbicide.psm, var = "Gap100plus", which = "both")
bal.plot(sr.herbicide.psm, var = "CETWI", which = "both")
bal.plot(sr.herbicide.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(sr.herbicide.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(sr.herbicide.psm, type = "qq")


# Matched data
sr.herbicide.matched <- match_data(sr.herbicide.psm)

# Create trt_control variable
sr.herbicide.matched <- sr.herbicide.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Herbicide", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Herbicide")))

# Center and scale numeric variables
sr.herbicide.matched <- sr.herbicide.matched |>
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
sr.herbicide.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = sr.herbicide.matched,
  weights = weights
)

# G computation to estimate marginal effects
sr.herbicide.pred <- avg_predictions(
  model = sr.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
sr.herbicide.pred

#   Table version
sr.herbicide.pred.df <- plot_predictions(
  sr.herbicide.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Herbicide"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
sr.herbicide.pred.df

# Plot
sr.herbicide.plot <- sr.herbicide.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
sr.herbicide.plot

# Estimation of average treatment effect
sr.herbicide.comp <- avg_comparisons(
  model = sr.herbicide.lm,
  variables = "trt_control",
  vcov = ~subclass
)
sr.herbicide.comp



## Prescribed Burn --------------------------------------------------------

# Filter data
sr.prescribed.burn.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Prescribed Burn" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Southern Rockies")

# PSM
sr.prescribed.burn.psm <- matchit(
  data = sr.prescribed.burn.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
sr.prescribed.burn.psm
summary(sr.prescribed.burn.psm) # 37 treated matched


# Diagnostic love plot
love.plot(sr.prescribed.burn.psm, stars = "std") +
  labs(title = "Southern Rockies: Prescribed Burn")

# eCDF plots
plot(sr.prescribed.burn.psm, type = "ecdf")
bal.plot(sr.prescribed.burn.psm, which = "both", type = "ecdf")
bal.plot(sr.prescribed.burn.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(sr.prescribed.burn.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(sr.prescribed.burn.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(sr.prescribed.burn.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(sr.prescribed.burn.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(sr.prescribed.burn.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(sr.prescribed.burn.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(sr.prescribed.burn.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(sr.prescribed.burn.psm, type = "density")
bal.plot(sr.prescribed.burn.psm, which = "both")
bal.plot(sr.prescribed.burn.psm, var = "BareSoil_FH", which = "both")
bal.plot(sr.prescribed.burn.psm, var = "TotalFoliarCover", which = "both")
bal.plot(sr.prescribed.burn.psm, var = "ForbCover_AH", which = "both")
bal.plot(sr.prescribed.burn.psm, var = "GramCover_AH", which = "both")
bal.plot(sr.prescribed.burn.psm, var = "ShrubCover_AH", which = "both")
bal.plot(sr.prescribed.burn.psm, var = "Gap100plus", which = "both")
bal.plot(sr.prescribed.burn.psm, var = "CETWI", which = "both")
bal.plot(sr.prescribed.burn.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(sr.prescribed.burn.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(sr.prescribed.burn.psm, type = "qq")


# Matched data
sr.prescribed.burn.matched <- match_data(sr.prescribed.burn.psm)

# Create trt_control variable
sr.prescribed.burn.matched <- sr.prescribed.burn.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Prescribed Burn", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Prescribed Burn")))

# Center and scale numeric variables
sr.prescribed.burn.matched <- sr.prescribed.burn.matched |>
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
sr.prescribed.burn.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = sr.prescribed.burn.matched,
  weights = weights
)

# G computation to estimate marginal effects
sr.prescribed.burn.pred <- avg_predictions(
  model = sr.prescribed.burn.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
sr.prescribed.burn.pred

#   Table version
sr.prescribed.burn.pred.df <- plot_predictions(
  sr.prescribed.burn.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Prescribed Burn"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
sr.prescribed.burn.pred.df

# Plot
sr.prescribed.burn.plot <- sr.prescribed.burn.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
sr.prescribed.burn.plot

# Estimation of average treatment effect
sr.prescribed.burn.comp <- avg_comparisons(
  model = sr.prescribed.burn.lm,
  variables = "trt_control",
  vcov = ~subclass
)
sr.prescribed.burn.comp



## Vegetation Disturbance -------------------------------------------------

# Filter data
sr.veg.disturbance.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Vegetation Disturbance" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Southern Rockies")

# PSM
sr.veg.disturbance.psm <- matchit(
  data = sr.veg.disturbance.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
sr.veg.disturbance.psm
summary(sr.veg.disturbance.psm) # 28 treated matched


# Diagnostic love plot
love.plot(sr.veg.disturbance.psm, stars = "std") +
  labs(title = "Southern Rockies: Vegetation Disturbance")

# eCDF plots
plot(sr.veg.disturbance.psm, type = "ecdf")
bal.plot(sr.veg.disturbance.psm, which = "both", type = "ecdf")
bal.plot(sr.veg.disturbance.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(sr.veg.disturbance.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(sr.veg.disturbance.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(sr.veg.disturbance.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(sr.veg.disturbance.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(sr.veg.disturbance.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(sr.veg.disturbance.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(sr.veg.disturbance.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(sr.veg.disturbance.psm, type = "density")
bal.plot(sr.veg.disturbance.psm, which = "both")
bal.plot(sr.veg.disturbance.psm, var = "BareSoil_FH", which = "both")
bal.plot(sr.veg.disturbance.psm, var = "TotalFoliarCover", which = "both")
bal.plot(sr.veg.disturbance.psm, var = "ForbCover_AH", which = "both")
bal.plot(sr.veg.disturbance.psm, var = "GramCover_AH", which = "both")
bal.plot(sr.veg.disturbance.psm, var = "ShrubCover_AH", which = "both")
bal.plot(sr.veg.disturbance.psm, var = "Gap100plus", which = "both")
bal.plot(sr.veg.disturbance.psm, var = "CETWI", which = "both")
bal.plot(sr.veg.disturbance.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(sr.veg.disturbance.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(sr.veg.disturbance.psm, type = "qq")


# Matched data
sr.veg.disturbance.matched <- match_data(sr.veg.disturbance.psm)

# Create trt_control variable
sr.veg.disturbance.matched <- sr.veg.disturbance.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Vegetation Disturbance", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Vegetation Disturbance")))

# Center and scale numeric variables
sr.veg.disturbance.matched <- sr.veg.disturbance.matched |>
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
sr.veg.disturbance.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = sr.veg.disturbance.matched,
  weights = weights
)

# G computation to estimate marginal effects
sr.veg.disturbance.pred <- avg_predictions(
  model = sr.veg.disturbance.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
sr.veg.disturbance.pred

#   Table version
sr.veg.disturbance.pred.df <- plot_predictions(
  sr.veg.disturbance.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Vegetation Disturbance"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
sr.veg.disturbance.pred.df

# Plot
sr.veg.disturbance.plot <- sr.veg.disturbance.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
sr.veg.disturbance.plot

# Estimation of average treatment effect
sr.veg.disturbance.comp <- avg_comparisons(
  model = sr.veg.disturbance.lm,
  variables = "trt_control",
  vcov = ~subclass
)
sr.veg.disturbance.comp




# Wyoming Basin -----------------------------------------------------------

## Prescribed Burn --------------------------------------------------------

# Filter data
wb.prescribed.burn.dat <- ldc.007 |>
  filter(Trt_Type_Sub == "Prescribed Burn" |
           is.na(Trt_Type_Sub)) |>
  filter(str_detect(Category, "never")) |>
  filter(EcoLvl3 == "Wyoming Basin")

# PSM
wb.prescribed.burn.psm <- matchit(
  data = wb.prescribed.burn.dat,
  formula = trt_binary ~ MLRARSYM + BareSoil_FH + TotalFoliarCover + ForbCover_AH +
    GramCover_AH + ShrubCover_AH + Gap100plus + CETWI + sandtotal_0_cm,
  distance = "glm",
  link = "logit",
  method = "nearest",
  caliper = 0.2,
  ratio = 2
)
wb.prescribed.burn.psm
summary(wb.prescribed.burn.psm) # 37 treated matched


# Diagnostic love plot
love.plot(wb.prescribed.burn.psm, stars = "std") +
  labs(title = "WY Basin: Prescribed Burn")

# eCDF plots
plot(wb.prescribed.burn.psm, type = "ecdf")
bal.plot(wb.prescribed.burn.psm, which = "both", type = "ecdf")
bal.plot(wb.prescribed.burn.psm, var = "BareSoil_FH", which = "both", type = "ecdf")
bal.plot(wb.prescribed.burn.psm, var = "TotalFoliarCover", which = "both", type = "ecdf")
bal.plot(wb.prescribed.burn.psm, var = "ForbCover_AH", which = "both", type = "ecdf")
bal.plot(wb.prescribed.burn.psm, var = "GramCover_AH", which = "both", type = "ecdf")
bal.plot(wb.prescribed.burn.psm, var = "ShrubCover_AH", which = "both", type = "ecdf")
bal.plot(wb.prescribed.burn.psm, var = "Gap100plus", which = "both", type = "ecdf")
bal.plot(wb.prescribed.burn.psm, var = "CETWI", which = "both", type = "ecdf")
bal.plot(wb.prescribed.burn.psm, var = "sandtotal_0_cm", which = "both", type = "ecdf")

# Density plots
plot(wb.prescribed.burn.psm, type = "density")
bal.plot(wb.prescribed.burn.psm, which = "both")
bal.plot(wb.prescribed.burn.psm, var = "BareSoil_FH", which = "both")
bal.plot(wb.prescribed.burn.psm, var = "TotalFoliarCover", which = "both")
bal.plot(wb.prescribed.burn.psm, var = "ForbCover_AH", which = "both")
bal.plot(wb.prescribed.burn.psm, var = "GramCover_AH", which = "both")
bal.plot(wb.prescribed.burn.psm, var = "ShrubCover_AH", which = "both")
bal.plot(wb.prescribed.burn.psm, var = "Gap100plus", which = "both")
bal.plot(wb.prescribed.burn.psm, var = "CETWI", which = "both")
bal.plot(wb.prescribed.burn.psm, var = "sandtotal_0_cm", which = "both")

# Bar plot
bal.plot(wb.prescribed.burn.psm, var = "MLRARSYM", which = "both")

# eQQ plots
plot(wb.prescribed.burn.psm, type = "qq")


# Matched data
wb.prescribed.burn.matched <- match_data(wb.prescribed.burn.psm)

# Create trt_control variable
wb.prescribed.burn.matched <- wb.prescribed.burn.matched |>
  mutate(trt_control = if_else(trt_binary == 1, "Prescribed Burn", "Not Treated")) |>
  mutate(trt_control = factor(trt_control, levels = c("Not Treated", "Prescribed Burn")))

# Center and scale numeric variables
wb.prescribed.burn.matched <- wb.prescribed.burn.matched |>
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
wb.prescribed.burn.lm <- lm(
  ln_q ~ trt_control * (
    BareSoil_FH_scaled + TotalFoliarCover_scaled + ForbCover_AH_scaled + GramCover_AH_scaled +
      ShrubCover_AH_scaled + Gap100plus_scaled + CETWI_scaled + sandtotal_0_cm_scaled),
  data = wb.prescribed.burn.matched,
  weights = weights
)

# G computation to estimate marginal effects
wb.prescribed.burn.pred <- avg_predictions(
  model = wb.prescribed.burn.lm,
  variables = "trt_control",
  vcov = ~subclass,
  by = "trt_control"
)
wb.prescribed.burn.pred

#   Table version
wb.prescribed.burn.pred.df <- plot_predictions(
  wb.prescribed.burn.lm,
  by = "trt_control",
  newdata = datagrid(
    trt_control = c("Not Treated", "Prescribed Burn"),
    grid_type = "counterfactual"
  ),
  draw = FALSE,
  vcov = ~subclass
)
wb.prescribed.burn.pred.df

# Plot
wb.prescribed.burn.plot <- wb.prescribed.burn.pred.df |>
  ggplot(aes(x = trt_control, y = estimate)) +
  geom_point(
    shape = 18,
    size = 4,
    color = "red"
  ) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw()
wb.prescribed.burn.plot

# Estimation of average treatment effect
wb.prescribed.burn.comp <- avg_comparisons(
  model = wb.prescribed.burn.lm,
  variables = "trt_control",
  vcov = ~subclass
)
wb.prescribed.burn.comp



save.image("RData/13_PSM-calculations.RData")
