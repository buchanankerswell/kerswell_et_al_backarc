# Load functions
source('functions.R')
load('data/hf.RData')
library(doParallel)

# Bounding Boxes
purrr::map(
  seg.names,
  ~shp.sa.segs.robin.pacific.buffer %>%
  filter(segment == .x) %>%
  st_bbox() %>%
  bbox_widen(
    crs = proj4.robin.pacific,
    borders = c('top' = 0.1,
    'bottom' = 0.1,
    'left' = 0.1,
    'right' = 0.1))) %>%
  purrr::set_names(nm = seg.names) -> shp.box

# Crop data
purrr::map(
  shp.box,
  ~shp.hf %>%
  rename(hf = `heat-flow (mW/m2)`) %>%
  st_crop(.x)) -> shp.hf.crop

# Find optimal kriging parameters by genetic algorithm
args <- commandArgs(trailingOnly = TRUE)
cat('Running Genetic Algorithm for', args[1], 'segment\n')

seg.name <- args[1]
data <- shp.hf.crop[[args[1]]]
nfold <- as.numeric(args[2])

# Cost function (to minimize)
cat(
  nfold, 'fold cross-validation over',
  nrow(data), 'grid points\n')
cat('Defining cost function\n')
v.opt <- function(x, data, nfold){
  -f.obj(
    data = data,
    param = 'hf',
    fold = nfold,
    cutoff = x[1],
    lag.start = x[2],
    v.mod = x[3],
    v.sill = x[4],
    v.range = x[5],
    v.nug = x[6],
    maxdist = x[7])}

# Suggested chromosomes for initial population
tibble(
 cutoff = runif(50, 3, 15),
  lag.start = runif(50, 1, 5),
  v.mod = runif(50, 0, 2),
  v.sill = runif(50, 1, 2000000),
  v.range = runif(50, 1, 1000000),
  v.nug = runif(50, 0, 2000000),
  maxdist = runif(50, 1, 10000000)
) %>%
as.matrix() -> suggestions

# Cluster setup
nodelist <- Sys.getenv("NODES")
nodelist <- strsplit(nodelist, "\n", fixed=TRUE)[[1]]
nodelist <- nodelist[!(nodelist %in% nodelist[1])]
workers <- rep(nodelist, each = 40)

cl <- makeCluster(workers, type = "PSOCK")
registerDoParallel(cl)

clusterCall(cl, library, package = 'magrittr', character.only = TRUE, quietly = TRUE)
clusterCall(cl, library, package = 'dplyr', character.only = TRUE, quietly = TRUE)
clusterCall(cl, library, package = 'tidyr', character.only = TRUE, quietly = TRUE)
clusterCall(cl, library, package = 'purrr', character.only = TRUE, quietly = TRUE)
clusterCall(cl, library, package = 'GA', character.only = TRUE, quietly = TRUE)
clusterCall(cl, library, package = 'gstat', character.only = TRUE, quietly = TRUE)
clusterCall(cl, library, package = 'sf', character.only = TRUE, quietly = TRUE)

# Export variables to nodes
clusterExport(
  cl,
  varlist = c(
    'seg.name',
    'data',
    'nfold',
    'suggestions',
    'f.obj',
    'v.opt'))

# Print info
cat('Running Genetic Algorithm for', args[1], 'segment\n')
cat(
  nfold, 'fold cross-validation over',
  nrow(data), 'grid points\n')

# Genetic algorithm (GA) optimization
print(nodelist)
GA::ga(
  type = "real-valued", 
  fitness = function(x) v.opt(x = x, data = data, nfold = nfold),
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
 popSize = 50,
 maxiter = 250,
 run = 50,
 monitor = T,
 parallel = cl) -> opt

# Save
dir.create('data/ga', showWarnings = FALSE)
fname <- paste0(seg.name, '_opt')
cat('Saving results to:', paste0(fname, '.RData'), '\n')
assign(fname, opt)
save(list = fname, file = paste0('data/ga/', fname, '.RData'))

stopCluster(cl)
cat('\nDone')
