rm(list = ls())
source('functions.R')
load('data/hf.RData')
rm.list <- ls()

# Kriging

# Whole segment
# Calculate experimental variogram, fit experimental variogram, krig, and interpolate
k <- purrr::map(
  unique(shp.hf$segment),
  ~ krg(
    data = shp.hf %>% filter(segment == .x) %>% filter(Heat_Flow > 10 & Heat_Flow < 120 & minD < maxD & Gradient > 0 & Gradient < 500),
    lags = 100,
    lag.cutoff = 3,
    param = 'Heat_Flow',
    krg.shp = shp.sa.segs.robin.pacific.buffer %>% filter(segment == .x),
    seg = shp.sa.segs.robin.pacific %>% filter(segment == .x),
    contours = shp.sa.countours.robin.pacific %>% filter(segment == .x),
    crs = proj4.robin.pacific,
    ngrid = 3e5,
    grid.method = 'hexagonal',
    v.mod <- c('Sph', 'Nug', 'Gau', 'Mat', 'Lin', 'Cir', 'Per', 'Wav', 'Log', 'Pow'),
    plot = T
  )
) %>% purrr::set_names(nm = unique(shp.hf$segment))

# Subsegment
k.subseg <- purrr::map(seg.names, ~{
  d <- shp.hf.subseg %>% filter(segment == .x) %>% filter(Heat_Flow > 10 & Heat_Flow < 120 & minD < maxD & Gradient > 0 & Gradient < 500)
  sg <- shp.sa.segs.robin.pacific.subseg.buffer %>% filter(segment == .x)
  buf <- shp.sa.segs.robin.pacific.subseg.buffer %>% filter(segment == .x)
  cntr <- shp.sa.countours.robin.pacific %>% filter(segment == .x)
  purrr::map(unique(d$subseg) %>% sort(), ~{
    krg(
      data = d %>% filter(subseg == .x),
      lags = 100,
      lag.cutoff = 3,
      param = 'Heat_Flow',
      krg.shp = buf %>% filter(subseg == .x),
      seg = sg %>% filter(subseg == .x),
      contours = cntr,
      crs = proj4.robin.pacific,
      ngrid = 3e5,
      grid.method = 'hexagonal',
      v.mod <- c('Sph', 'Nug', 'Gau', 'Mat', 'Lin', 'Cir', 'Per', 'Wav', 'Log', 'Pow'),
      plot = T
    )
  })
}) %>% purrr::set_names(nm = seg.names)

# Summarise experimental variograms (wide variograms)
nugs <- bind_rows(purrr::map(k, ~.x$v.model), .id = 'segment') %>%
  filter(model == 'Nug') %>%
  select(c(segment, psill)) %>%
  as_tibble() %>%
  rename(nug = psill)
params <- bind_rows(purrr::map(k, ~.x$v.model), .id = 'segment') %>%
  filter(model == 'Sph') %>%
  select(c(segment, psill, range)) %>%
  as_tibble() %>%
  rename(sill = psill)
v.mods <- left_join(nugs, params) %>% mutate('sill' = sill+nug)

# Check number of point pairs
seg.np <- bind_rows(purrr::map(k, ~.x$variogram), .id = 'segment') %>%
  group_by(segment) %>%
  summarise(np.max = max(np), np.min = min(np), np.mean = mean(np), np.median = median(np)) %>%
  arrange(desc(np.median))

# Filter out segments with median point pairs less than 30
seg.30np <- seg.np %>% filter(np.median > 30)

# Clean up environment
rm(list = rm.list)
rm(dat, rm.list, nugs.n, nugs, params.n, params)

# Save
save.image('data/krig.RData')
# save(k.subseg, file = 'data/krig_subseg.RData')