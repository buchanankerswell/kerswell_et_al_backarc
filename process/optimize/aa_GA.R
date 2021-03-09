# Load functions and packages
cat('Loading packages')
source('../../functions.R')

# Load data
cat('\nLoading data')
load('../../data/hf.RData')

# Buffer
shp.buf.aa <-
  shp.sa.segs.robin.pacific.buffer %>%
  filter(segment == 'Alaska Aleutians')

# Box
shp.box.aa <-
  bbox_widen(st_bbox(shp.buf.aa),
             crs = st_crs(shp.buf.aa),
             borders = c(
               'left' = 0.1,
               'right' = 0.1,
               'top' = 0.1,
               'bottom' = 0.1)) %>%
  st_as_sf()

# Crop luca hf
cat('\nCropping heat flow data\n')
shp.hf.aa <-
  shp.hf %>%
  filter(`heat-flow (mW/m2)` < 250) %>%
  st_crop(shp.box.aa)

# shp.hf.aa <- shp.hf.aa %>% slice_tail(n = 100)

# Clean up environment
rm(list = rm.list)

# Cost function (to minimize)
v.opt <- function(x){
  f.obj(data = shp.hf.aa %>% rename(hf = `heat-flow (mW/m2)`),
            param = 'hf',
            v.mod = x[1],
            v.sill = x[2],
            v.range = x[3],
            maxdist = x[4])}

# Suggested chromosomes for initial population
tibble(
  v.mod = runif(50, 0, 3),
  v.sill = runif(50, 1, 5000),
  v.range = runif(50, 1, 1000000),
	v.nug = runif(50, 0, 5000),
  maxdist = runif(50, 1, 1000000)
) %>%
  as.matrix() -> suggestions

# Start the clock!
ptm <- proc.time()

# Genetic algorithm (GA) optimization
GA <- ga(type = "real-valued", 
         fitness = function(x) -v.opt(x),
         lower = c(1, 1, 1, 0, 1),
         upper = c(3, 5000, 1000000, 5000, 1000000),
         names = c('v.mod', 'v.sill', 'v.range', 'v.nug', 'maxdist'),
         suggestions = suggestions,
         popSize = 50,
         maxiter = 200,
         run = 30,
         parallel = T)

# Stop the clock
print(proc.time() - ptm)

# Save
save(GA, file = 'aa_GA.RData')
