# Load packages
# Quiet loading
sshhh <- function(p){
  suppressWarnings(
    suppressPackageStartupMessages(
      library(p, quietly = T, character.only=TRUE)))}

# Package list
c('magrittr', 'ggplot2', 'tidyr', 'readr', 'purrr',
  'gstat', 'ggsflabel', 'sf', 'ggrepel', 'patchwork',
  'cowplot', 'dplyr') -> p.list

#cat('Loading libraries:', p.list, sep = '\n')

# auto-load quietly
sapply(p.list, sshhh)

#cat('Loading functions\n')

# Draw a widened box from a st_bbox object
bbox_widen <- function(
  bbox,
  crs,
  borders = c('left' = 0.5,
  'right' = 0.5,
  'top' = 0,
  'bottom' = 0)) {
  b <- bbox # current bounding box
  xrange <- b$xmax - b$xmin # range of x values
  yrange <- b$ymax - b$ymin # range of y values
  b[1] <- b[1] - (borders['left'] * xrange) # xmin - left
  b[3] <- b[3] + (borders['right'] * xrange) # xmax - right
  b[2] <- b[2] - (borders['bottom'] * yrange) # ymin - bottom
  b[4] <- b[4] + (borders['top'] * yrange) # ymax - top
  box <- st_polygon(list(matrix(c(
    b$xmin,
    b$ymax,
    b$xmin,
    b$ymin,
    b$xmax,
    b$ymin,
    b$xmax,
    b$ymax,
    b$xmin,
    b$ymax),
    ncol = 2,
    byrow = TRUE))) %>%
  st_sfc(crs = crs)
  return(box)
}

# Read gmt files, wrap the dateline to avoid plotting horizontal lines on map,
# make into tibble, add segment names, transform projection, and bind into one sf object
read_latlong <- function(file, fname, crs) {
  st_read(
    file,
    crs = '+proj=longlat +lon_wrap=180 +ellps=WGS84 +datum=WGS84 +no_defs',
    quiet = TRUE) %>%
  st_transform(crs) %>%
  tibble::as_tibble() %>%
  st_as_sf() %>%
  mutate(segment = fname, .before = geometry)
}

# Krige optimization algorithmcross-validation used for optimizing variogram model
f.obj <- function(
  data,
  param,
  cutoff = 3,
  lag.start = 1,
  v.mod = 0,
  v.sill = NA,
  v.range = NA,
  v.nug = NA,
  maxdist = Inf,
  fold = 10
) {
  # Variogram model discritization formula
  if(v.mod >= 0 && v.mod < 1) {
    v.mod <- 'Sph'
  } else if(v.mod >= 1 && v.mod <= 2) {
    v.mod <- 'Exp'
  }
  # Experimental variogram
  cutoff <- max(st_distance(data))/cutoff
  lags <- 15
  width <- as.vector(cutoff)/lags
  shift.cutoff <- width*(lags + lag.start)
  v <- variogram(
    get(param)~1,
    locations = data,
    cutoff = shift.cutoff,
    width = width)
  v.grm <- v[lag.start:nrow(v),]
  # Model variogram
  fit.variogram(
    v.grm,
    vgm(psill = v.sill,
        model = v.mod,
        range = v.range,
        nugget = v.nug),
    fit.method = 7) -> f
  # Kriging with n-fold cross validation
  krige.cv(
    formula = get(param)~1,
    locations = data,
    model = f,
    maxdist = maxdist,
    nfold = fold,
    verbose = T) %>%
  drop_na() %>%
  try() -> k
  if(class(k) == 'try-error') {
    return(Inf)
  } else {
  # Calculating cost function after Li et al 2018
  # Simultaneously minimizes misfit on variogram model and kriged interpolation errors
  # Weights
  wi <- 0.7 # interpolation error
  wf <- 0.3 # variogram fit error
  # Calculate variogram fit error
  # Root mean weighted (Nj/hj^2) sum of squared errors * (1-wi)/v.sigma
  sqrt((attr(f, "SSErr"))/nrow(v)) * ((1-wi)/sd(v$gamma)) -> v.s
  # Calculate interpolation error
  # RMSE * wi/i.sigma
  sqrt(sum((k$residual^2))/nrow(k)) * wi/sd(k$var1.pred) -> i.s
  # Cost
  return(v.s + i.s)
  }
}

# Kriging
Krige <- function(
  data,
  v.mod,
  param,
  grid,
  cv = FALSE) {
  # Kriging
  krige(
    formula = get(param)~1,
    locations = data,
    newdata = grid,
    model = v.mod) %>%
  as_tibble() %>%
  st_as_sf() -> k
  if(cv == TRUE) {
    cat('\nCalculating cross-validation error')
    krige.cv(
        formula = get(param)~1,
        locations = data,
        model = v.mod,
        nfold = 10,
        verbose = F) %>%
    drop_na() -> k.cv
  }
  attr(k, 'cv.rmse') <- sqrt(sum(k.cv$residual^2)/nrow(k.cv))
  return(k)
}

# Optimize krige results
Krige_opt <- function(
  seg.name,
  data,
  param,
  n.init = 50,
  maxitr = 200,
  run = 50,
  nfold = 10){
  # Cost function (to minimize)
  cat(
    nfold, 'fold cross-validation over',
    nrow(data), 'grid points\n')
  cat('Defining cost function\n')
  v.opt <- function(x){
    f.obj(
      data = data,
      param = param,
      cutoff = x[1],
      lag.start = x[2],
      v.mod = x[3],
      v.sill = x[4],
      v.range = x[5],
      v.nug = x[6],
      maxdist = x[7],
      fold = nfold)}
  # Suggested chromosomes for initial population
  cat('Initializing', n.init, 'chromosomes\n')
  tibble(
    cutoff = runif(n.init, 3, 15),
    lag.start = runif(n.init, 1, 5),
    v.mod = runif(n.init, 0, 2),
    v.sill = runif(n.init, 1, 2000000),
    v.range = runif(n.init, 1, 1000000),
    v.nug = runif(n.init, 0, 2000000),
    maxdist = runif(n.init, 1, 10000000)
  ) %>%
  as.matrix() -> suggestions
  # Genetic algorithm (GA) optimization
  cat('Optimizing using genetic algorithm with:\n',
    'Population:', n.init, '\n',
    'Max generations:', maxitr, '\n',
    'Run cutoff:', run, '\n',
    'Cross-fold validation:', nfold, '\n')
  GA::ga(
    type = "real-valued", 
    fitness = function(x) -v.opt(x),
    lower = c(3, 1, 0, 1, 1, 0, 1),
    upper = c(15, 5, 2, 2000000, 1000000, 2000000, 10000000),
    names = c(
      'cutoff',
      'lag.start',
      'v.mod',
      'v.sill',
      'v.range',
      'v.nug',
      'maxdist'),
   suggestions = suggestions,
   popSize = n.init,
   maxiter = maxitr,
   run = run,
   monitor = T,
   parallel = T) -> opt
  # Save
  dir.create('data/ga', showWarnings = FALSE)
  fname <- paste0(
    seg.name %>%
    stringr::str_replace_all(' ', '_') %>%
    stringr::str_replace_all('\\.', ''), '_opt')
  cat('Saving results to:', paste0(fname, '.RData'), '\n')
  assign(fname, opt)
  save(list = fname, file = paste0('data/ga/', fname, '.RData'))
}

# Krige, take difference, and visualize
Krige_diff <- function(
  seg.name,
  data,
  v.grm,
  v.mod,
  grid,
  param,
  data.compare,
  path){
  v <- fit.variogram(v.grm, model = v.mod)
  k <- Krige(data, v, param, grid, cv = T)
  # Difference
  shp.hf.pred <-
  data.compare %>%
  rename(
    hf.pred.luca = HF_pred,
    sigma.luca = sHF_pred,
    hf.obs.luca = Hf_obs) %>%
  st_crop(k) %>%
  mutate(
    hf.pred.krige = round(k$var1.pred, 1),
    sigma.krige = round(sqrt(k$var1.var), 1),
    hf.diff = round(hf.pred.luca - hf.pred.krige, 1),
    .before = geometry)
  # Save
  dir.create('data/diff', showWarnings = FALSE)
  fname <- seg.name %>%
  stringr::str_replace_all(' ', '_') %>%
  stringr::str_replace_all('\\.', '')
  assign(fname, list('k' = k, 'v.grm' = v.grm, 'v.mod' = v, 'diff' = shp.hf.pred))
  save(list = fname, file = paste0(path, fname, '.RData'))
}

# Plot variogram
plot_vgrm <- function(v.grm, v.mod, seg.name){
  ggplot() +
  geom_line(data = variogramLine(v.mod, maxdist = max(v.grm$dist)),
    aes(x = dist/1000, y = gamma)) +
  geom_point(data = v.grm, aes(x = dist/1000, y = gamma), size = 0.6) +
  labs(x = NULL, y = 'Semivariance', title = seg.name) +
  theme_classic() +
  theme(axis.text = element_text(color = 'black'),
  plot.title = element_text(hjust = 0.5),
  axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.y = element_blank())
}
