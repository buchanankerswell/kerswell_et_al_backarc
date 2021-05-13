# Load functions and libraries
source('functions.R')
load('data/hf.RData')

# Load genetic algorithm results and decode chromosome into a variogram model
cat('\nLoading genes from data/ga')
flist <- list.files(path = 'data/ga', full.name = T)
fnames <- list.files(path = 'data/ga') %>%
  stringr::str_replace_all('\\_opt.RData', '')

# Bounding Boxes
cat('\nDefining bounding boxes')
purrr::map(fnames,
  ~shp.sa.segs.robin.pacific.buffer %>%
  filter(segment == .x) %>%
  st_bbox() %>%
  bbox_widen(crs = proj4.robin.pacific,
  borders = c('top' = 0.1,
              'bottom' = 0.1,
              'left' = 0.1,
              'right' = 0.1))) %>%
purrr::set_names(nm = fnames) -> shp.box

# Crop data
cat('\nCropping heat flow data')
purrr::map(shp.box,
  ~shp.hf %>%
  rename(hf = `heat-flow (mW/m2)`) %>%
  st_crop(.x)) %>%
purrr::set_names(nm = fnames) -> shp.hf.crop

# Load ga results
for(i in flist) load(i)

# Decode chromosomes
# Variogram models
cat('\nDecoding genes')
purrr::map(fnames,
  ~{
  # Get ga result
   o <- get(paste0(.x, '_opt'))@solution
  # Variogram model discritization formula
   if(o[3] >= 0 && o[3] < 1) {
     v.mod <- 'Sph'
   } else if(o[3] >= 1 && o[3] <= 2) {
     v.mod <- 'Exp'
   }
  # Construct variogram model
  cat('\nConstructing variogram model for', .x)
   vgm(psill = o[4], model = v.mod, range = o[5], nugget = o[6])
 }) %>%
purrr::set_names(fnames) -> v.mods

# Calculate experimental variogram
purrr::map(fnames, ~{
  d <- shp.hf.crop[[.x]]
  o <- get(paste0(.x, '_opt'))@solution
  cutoff <- max(st_distance(d))/o[1]
  lags <- 15
  lag.start <- round(o[2])
  width <- as.vector(cutoff)/lags
  cutoff.shifted <- width*(lags + lag.start)
  cat('\nCalculating experimental variogram for', .x)
  v <- variogram(hf~1,
    locations = d,
    cutoff = cutoff.shifted,
    width = width)
    v[lag.start:nrow(v),]
 }) %>%
purrr::set_names(nm = fnames) -> v.grms

# Fit experimental variograms
cat('\nFitting variograms')
purrr::map2(v.grms, v.mods,
  ~{fit.variogram(.x, model = .y)}) %>%
purrr::set_names(fnames) -> v.fits

cat('\nPlotting variogram models and saving to figs/vgrms_decoded.png')
purrr::pmap(
  list(v.grms, v.fits, fnames),
  ~plot_vgrm(..1, ..2, ..3 %>% stringr::str_replace_all('_', ' '))) %>%
wrap_plots() -> p

# Save plot
ggsave(
  'figs/vgrms_decoded.png',
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 7,
  height = 7
)

# Save
cat('\nSaving variograms to data/genes_decoded.RData')
save(v.mods, v.grms, v.fits, file = 'data/genes_decoded.RData')