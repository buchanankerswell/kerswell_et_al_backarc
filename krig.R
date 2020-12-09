rm(list = ls())
source('functions.R')
load('data/hf.RData')
rm.list <- ls()

# Kriging
# Filter data
dat <- shp.hf %>% filter(Heat_Flow > 10 & Heat_Flow < 120 & minD < maxD & Gradient > 0 & Gradient < 500)

# Calculate experimental variogram, fit experimental variogram, krig, and interpolate
k100.w <- purrr::map(
  unique(shp.hf$segment),
  ~ krg(
    data = dat %>% filter(segment == .x),
    lags = 100,
    lag.cutoff = 3,
    param = 'Heat_Flow',
    krg.shp = shp.sa.segs.robin.pacific.buffer %>% filter(segment == .x),
    seg = shp.sa.segs.robin.pacific %>% filter(segment == .x),
    contours = shp.sa.countours.robin.pacific %>% filter(segment == .x),
    crs = proj4.robin.pacific,
    ngrid = 3.2e5,
    grid.method = 'hexagonal',
    v.mod = 'Sph',
    plot = T
  )
) %>% purrr::set_names(nm = unique(shp.hf$segment))

# Calculate experimental variogram, fit experimental variogram, krig, and interpolate
k100.n <- purrr::map(
  unique(shp.hf$segment),
  ~ krg(
    data = dat %>% filter(segment == .x),
    lags = 100,
    lag.cutoff = 6,
    param = 'Heat_Flow',
    krg.shp = shp.sa.segs.robin.pacific.buffer %>% filter(segment == .x),
    seg = shp.sa.segs.robin.pacific %>% filter(segment == .x),
    contours = shp.sa.countours.robin.pacific %>% filter(segment == .x),
    crs = proj4.robin.pacific,
    ngrid = 3.2e5,
    grid.method = 'hexagonal',
    v.mod = 'Sph',
    plot = T
  )
) %>% purrr::set_names(nm = unique(shp.hf$segment))

# Check number of point pairs
seg.np <- bind_rows(purrr::map(k100.n, ~.x$variogram), .id = 'segment') %>%
  group_by(segment) %>%
  summarise(np.max = max(np), np.min = min(np), np.mean = mean(np), np.median = median(np)) %>%
  arrange(desc(np.median))

# Filter out segments with median point pairs less than 30
seg.30np <- np %>% filter(np.median > 30)

# Clean up environment
rm(list = rm.list)
rm(dat, rm.list)

# Save
save.image('data/krig.RData')
