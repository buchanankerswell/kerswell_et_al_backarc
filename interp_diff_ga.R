# Load functions and libraries
source('functions.R')
load('data/hf.RData')

# Load decoded genes
cat('\nLoading decoded chromosomes from data/genes_decoded.RData')
load('data/genes_decoded.RData')

list.files(path = 'data/ga') %>%
stringr::str_replace_all('\\_opt.RData', '') -> fnames

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

# Cropped Grid
cat('\nDefining interpolation grids')
purrr::map(shp.box, ~shp.grid %>% st_crop(.x)) %>%
  purrr::set_names(nm = fnames) -> shp.grid.crop

# Interpolation difference
cat('\nCalculating interpolation differences\n')
pb <- progress::progress_bar$new(
  format = "Interpolating [:bar] :percent in :elapsed",
  total = length(seg.names), clear = FALSE, width= 60)
purrr::pwalk(list(
  fnames,
  shp.hf.crop,
  v.grms,
  v.mods,
  shp.grid.crop),
  ~{pb$tick()
    cat('\n', ..1, '\n')
    Krige_diff(..1, ..2, ..3, ..4, ..5,
               param = 'hf',
               data.compare = shp.hf.pred,
               path = 'data/diff_ga/')})

cat('\nDone\n')
