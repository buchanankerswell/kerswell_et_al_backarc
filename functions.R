# Load packages
# Quiet loading
function(p){
	library(p, character.only=TRUE) %>%
	suppressPackageStartupMessages() %>%
	suppressWarnings()} -> sshhh

# Package list
c('magrittr', 'ggplot2', 'tidyr', 'readr', 'purrr',
	'gstat', 'ggsflabel', 'sf', 'ggrepel', 'patchwork',
	'cowplot', 'dplyr') -> p.list

cat('Loading libraries', p.list, sep = '\n')

# auto-load quietly
sapply(p.list, sshhh)

cat('Loading functions\n')

# Draw a widened box from a st_bbox object
bbox_widen <- function(bbox,
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
  box <- st_polygon(list(matrix(c(b$xmin,
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
read_latlong <- function(files, filenames = NULL, crs) {
		f <-
		purrr::map2(files,
                filenames,
                ~ st_read(.x,
                          crs = '+proj=longlat +lon_wrap=180 +ellps=WGS84 +datum=WGS84 +no_defs',
                          quiet = TRUE) %>%
                  st_transform(crs) %>%
                  tibble::as_tibble() %>%
                  tibble::add_column('segment' = .y) %>%
                  dplyr::group_by(segment) %>%
                  st_as_sf()) %>%
		     bind_rows() %>%
				 purrr::set_names(nm = filenames) %>%
				 dplyr::bind_rows()
}

# Krige optimization algorithmcross-validation used for optimizing variogram model
f.obj <- function(
  data,
  param,
  lags = 15,
  lag.cutoff = 3,
  v.mod = 0,
  v.sill = NA,
  v.range = NA,
	v.nug = NA,
  maxdist = Inf,
  fold = 10
) {
  # Variogram cutoff
  max(st_distance(data))/lag.cutoff -> cutoff
  # Experimental variogram
  variogram(get(param)~1,
            locations = data,
            cutoff = as.vector(cutoff),
            width = as.vector(cutoff/lags)) -> v
  # Variogram model discritization formula
  if(v.mod >= 0 && v.mod < 1) {
    v.mod <- 'Sph'
  } else if(v.mod >= 1 && v.mod < 2) {
    v.mod <- 'Exp'
  } else if(v.mod >= 2 && v.mod <= 3) {
    v.mod <- 'Gau'
  }
  # Model variogram
  fit.variogram(v, 
		vgm(psill = v.sill, 
		    model = v.mod, 
		    range = v.range,
				nugget = v.nug), 
		fit.method = 7) -> f
  # Kriging with n-fold cross validation
  krige.cv(formula = get(param)~1,
           locations = data,
           model = f,
           maxdist = maxdist,
           verbose = T) %>%
	drop_na()-> k
  # Calculating cost function after Li et al 2018
  # Simultaneously minimizes misfit on variogram model and kriged interpolation errors
  # Weights
  wi <- 0.3 # interpolation error
  wf <- 0.7 # variogram fit error
  # Calculate variogram fit error
  # Root mean weighted (Nj/hj^2) sum of squared errors * (1-wi)/v.sigma
  sqrt((attr(f, "SSErr"))/nrow(v)) * ((1-wi)/sd(v$gamma)) -> v.s
  # Calculate interpolation error
  # RMSE * wi/i.sigma
  sqrt(sum((k$residual^2))/nrow(k)) * wi/sd(k$var1.pred) -> i.s
  # Cost
  v.s + i.s
}

# Kriging
Krige <- function(data,
									v.mod,
									param,
									grid) {
  # Kriging
  krige(formula = get(param)~1,
        locations = data,
        newdata = grid,
        model = v.mod) %>%
    as_tibble() %>%
    st_as_sf()
}

# Optimize krige results
Krige_opt <- function(seg.name,
											data,
											domain,
											param,
											n.init = 50,
											maxitr = 200,
											run = 50){
				# Cost function (to minimize)
				cat('Defining cost function\n')
				v.opt <- function(x){
				  f.obj(data = data,
				            param = param,
				            v.mod = x[1],
				            v.sill = x[2],
				            v.range = x[3],
				            maxdist = x[4])}
				# Suggested chromosomes for initial population
				cat('Initializing', n.init, 'chromosomes\n')
				tibble(
				  v.mod = runif(n.init, 0, 3),
				  v.sill = runif(n.init, 1, 5000),
				  v.range = runif(n.init, 1, 1000000),
					v.nug = runif(n.init, 0, 5000),
				  maxdist = runif(n.init, 1, 1000000)
				) %>%
 				as.matrix() -> suggestions
				# Genetic algorithm (GA) optimization
				cat('Optimizing using genetic algorithm with:\n',
						'Population:', n.init, '\n',
						'Max iterations:', maxitr, '\n',
						'Run cutoff:', run, '\n')
				opt <- try(
					GA::ga(type = "real-valued", 
				         fitness = function(x) -v.opt(x),
				         lower = c(1, 1, 1, 0, 1),
				         upper = c(3, 8000, 1000000, 8000, 1000000),
				         names = c('v.mod', 'v.sill', 'v.range', 'v.nug', 'maxdist'),
				         suggestions = suggestions,
				         popSize = n.init,
				         maxiter = maxitr,
				         run = run,
								 monitor = T,
				         parallel = T))
				if(class(opt) == 'try-error') opt <- NULL
				# Save
				fname <- paste0(seg.name %>% stringr::str_replace(' ', '_'), '_opt')
				cat('Saving results to:', paste0(fname, '.RData'), '\n')
				assign(fname, opt)
				save(list = fname, file = paste0('data/ga/', fname, '.RData'))
}

# Krige, take difference, and visualize
Krige_diff <- function(seg.name,
											 data,
											 v.mod,
											 param,
											 grid,
											 domain,
											 data.compare){
				k <-
				  Krige(data = data,
								v.mod = v.mod,
								param = param,
								grid = grid)
				# Difference
				shp.hf.pred <-
				  data.compare %>%
					rename(hf.pred.luca = HF_pred,
								 sigma.luca = sHF_pred,
								 hf.obs.luca = Hf_obs) %>%
				  st_crop(k$krige.results) %>%
				  mutate(hf.pred.krige = k$krige.results$hf,
								 sigma.krige = round(sqrt(k$krige.results$variance)),
				         hf.diff = hf.pred.luca - hf.pred.krige,
				         .before = geometry)
			# Save
				fname <- seg.name %>% stringr::str_replace(' ', '_')
				assign(fname, list('k' = k, 'diff' = shp.hf.pred))
				save(list = fname, file = paste0('data/diff/', fname, '_diff.RData'))
}
