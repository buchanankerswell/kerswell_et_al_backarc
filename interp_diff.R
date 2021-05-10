# Load functions and libraries
source('functions.R')
load('data/hf.RData')

# Bounding Boxes
cat('\nDefining bounding boxes')
purrr::map(seg.names,
  ~shp.sa.segs.robin.pacific.buffer %>%
  filter(segment == .x) %>%
  st_bbox() %>%
  bbox_widen(crs = proj4.robin.pacific,
  borders = c('top' = 0.1,
              'bottom' = 0.1,
              'left' = 0.1,
              'right' = 0.1))) %>%
  purrr::set_names(nm = seg.names) -> shp.box

# Crop data
cat('\nCropping heat flow data')
purrr::map(shp.box,
  ~shp.hf %>%
  rename(hf = `heat-flow (mW/m2)`) %>%
  st_crop(.x)) %>%
  purrr::set_names(nm = seg.names) -> shp.hf.crop

# Cropped Grid
cat('\nDefining interpolation grids')
purrr::map(shp.box, ~shp.grid %>% st_crop(.x)) %>%
  purrr::set_names(nm = seg.names) -> shp.grid.crop

# Define variograms manually
# Define variogram cutoffs for experimental variograms
cat('\nDefining variogram cutoffs for experimental variograms')
cutoff <- list(5, 6, 4, 5, 3, 8, 6, 3, 5, 15, 7, 3, 5) %>%
  purrr::set_names(seg.names)
lags <- list(15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15) %>%
  purrr::set_names(seg.names)
lag.start <- list(3, 1, 1, 2, 5, 1, 3, 1, 2, 1, 1, 1, 3)

# Calculate experimental variogram
cat('\nCalculating experimental variograms')
purrr::pmap(list(
  shp.hf.crop,
  cutoff,
  lags,
  lag.start,
  seg.names), ~{
  cutoff <- max(st_distance(..1))/..2
  lags <- ..3
  width <- as.vector(cutoff)/lags
  shift.cutoff <- width*(..3+..4)
  v <- variogram(hf~1,
  locations = .x,
  cutoff = shift.cutoff,
  width = width)
  v[..4:nrow(v),]
}) %>%
purrr::set_names(nm = seg.names) -> v.grms

# Arbitrary variogram models
cat('\nDefining variogram models')
purrr::map(1:13, ~vgm(model = c('Sph', 'Exp'))) %>%
purrr::set_names(nm = seg.names)-> v.mods

# Fit experimental variograms
cat('\nFitting variograms')
purrr::map2(v.grms, v.mods, ~{fit.variogram(.x, model = .y)}) %>%
purrr::set_names(seg.names) -> v.fits

cat('\nPlotting variogram models and saving to figs/vgrms.png')
purrr::pmap(list(v.grms, v.fits, seg.names), ~plot_vgrm(..1, ..2, ..3)) %>%
wrap_plots() -> p

# Save plot
ggsave(
  'figs/vgrms.png',
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 7,
  height = 7
)

# Interpolation difference
cat('\nCalculating interpolation differences')
purrr::pwalk(list(
  seg.names,
  shp.hf.crop,
  v.grms,
  v.mods,
  shp.grid.crop),
  Krige_diff,
  param = 'hf',
  data.compare = shp.hf.pred,
  path = 'data/diff/')