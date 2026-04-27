# Created: 2026-04-27
# Updated: 2026-04-27

# Purpose: Recalculate groups to find treatments per ecoregion with 30+ points, now that
#   CETWI and SOLUS data have been added (and some LDC points therefore removed).

# Result: Nothing with the treatment groups actually changed.


# Load data ---------------------------------------------------------------

ldc.cetwi.solus.raw <- read_csv("data/versions-from-R/12.2_LDC-CETWI-SOLUS_v007.csv")

# Convert date columns ----------------------------------------------------

# Convert relevant columns to date format
ldc.cetwi.solus <- ldc.cetwi.solus.raw |>
  mutate(
    DateVisited = as.Date(DateVisited, format = "%m/%d/%Y"),
    MR_trt_comp = as.Date(MR_trt_comp, format = "%m/%d/%Y"),
    MR_wildfire = as.Date(MR_wildfire, format = "%m/%d/%Y")
  )


# LDC points by Ecoregion 3 -----------------------------------------------

# Count table of at least 30 points per treatment group
level3.trt30 <- ldc.cetwi.solus |>
  group_by(EcoLvl1, EcoLvl2, EcoLvl3, Category) |>
  count(Trt_Type_Sub) |>
  filter(
    n >= 30,
    !is.na(Trt_Type_Sub)
  ) |>
  ungroup()
length(unique(level3.trt30$EcoLvl3)) # 13
unique(level3.trt30$Trt_Type_Sub)


# LDC points with at least 30 points per treatment group
eco3.trt30 <- level3.trt30 |>
  left_join(ldc.cetwi.solus)

#   Control equivalents
eco3.trt30.ctrl.count <- ldc.cetwi.solus |>
  filter(EcoLvl3 %in% eco3.trt30$EcoLvl3 & str_detect(Category, "Control")) |>
  group_by(EcoLvl3, Category) |>
  summarise(
    n = n(),
    .groups = "keep"
  ) |>
  ungroup()

eco3.trt30.ctrl <- ldc.cetwi.solus |>
  filter(EcoLvl3 %in% eco3.trt30$EcoLvl3 & str_detect(Category, "Control")) |>
  left_join(eco3.trt30.ctrl.count)

#   Combine & order cols
eco3.trt30.all <- eco3.trt30 |>
  bind_rows(eco3.trt30.ctrl) |>
  select(
    EcoLvl3, Category, Trt_Type_Sub, n, MR_trt_comp, LDCpointID, DateVisited,
    ProjectKey, PrimaryKey,
    Trt_Type_Major, recent_trt_count, FirePolyID, USGS_Assigned_ID, MR_wildfire,
    Fire_freq, Fire_freq_post_trt,
    EcoLvl1, NA_L1KEY, EcoLvl2, NA_L2KEY, NA_L3KEY, MLRARSYM, MLRA_name,
    horizontal_flux_total_MD, n_CETWI, CETWI, sandtotal_0_cm
  ) |>
  rename(Treatment_count = n) |>
  arrange(LDCpointID)
nrow(eco3.trt30.all) / nrow(ldc.cetwi.solus.raw) # all points used; none dropped


# Write to CSV ------------------------------------------------------------

write_csv(eco3.trt30.all,
          file = "data/versions-from-R/12.3_LDC-points_v007.csv",
          na = ""
)
