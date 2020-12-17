rm(list = ls())
source('functions.R')
load('data/hf.RData')
rm.list <- ls()

# Kriging
# Filter data
dat <- shp.hf %>% filter(Heat_Flow > 10 & Heat_Flow < 120 & minD < maxD & Gradient > 0 & Gradient < 500)
dat.subseg <- shp.hf.subseg %>% filter(Heat_Flow > 10 & Heat_Flow < 120 & minD < maxD & Gradient > 0 & Gradient < 500)

andes <- shp.hf.subseg %>% filter(segment == 'Andes') %>% filter(Heat_Flow > 10 & Heat_Flow < 120 & minD < maxD & Gradient > 0 & Gradient < 500)
andes.seg <- shp.sa.segs.robin.pacific.subseg %>% filter(segment == 'Andes')
andes.buf <- shp.sa.segs.robin.pacific.subseg.buffer %>% filter(segment == 'Andes')

cairo_pdf('figs/andes.sub.krig/plots.pdf', onefile = T)
purrr::map(
  unique(andes$subseg) %>% sort(),
  ~ krg(
    data = andes %>% filter(subseg == .x),
    lags = 20,
    lag.cutoff = 3,
    param = 'Heat_Flow',
    krg.shp = andes.buf %>% filter(subseg == .x),
    seg = andes.seg %>% filter(subseg == .x),
    contours = NULL,
    crs = proj4.robin.pacific,
    ngrid = 3e5,
    grid.method = 'hexagonal',
    v.mod = 'Sph',
    plot = T
  )
)
dev.off()

# Calculate experimental variogram, fit experimental variogram, krig, and interpolate
k100 <- purrr::map(
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

# Summarise experimental variograms (wide variograms)
nugs <- bind_rows(purrr::map(k100, ~.x$v.model), .id = 'segment') %>%
  filter(model == 'Nug') %>%
  select(c(segment, psill)) %>%
  as_tibble() %>%
  rename(nug = psill)
params <- bind_rows(purrr::map(k100, ~.x$v.model), .id = 'segment') %>%
  filter(model == 'Sph') %>%
  select(c(segment, psill, range)) %>%
  as_tibble() %>%
  rename(sill = psill)
v.mods <- left_join(nugs, params) %>% mutate('sill' = sill+nug)

# Check number of point pairs
seg.np <- bind_rows(purrr::map(k100, ~.x$variogram), .id = 'segment') %>%
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
