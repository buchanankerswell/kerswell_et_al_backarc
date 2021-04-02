# Load functions and libraries
source('functions.R')
load('data/hf.RData')

# Load genetic algorithm results and decode chromosome into a variogram model
# flist <- list.files(pattern = 'opt.RData', full.name = T)
# fnames <- list.files(pattern = 'opt.RData') %>%
#   purrr::map(~.x %>% stringr::str_replace('_opt.RData', '_vmod'))

# Save environment
# init.list <- ls()

# Load ga results
# for(i in flist) load(i)

# Decode chromosomes
# v.mods <- purrr::map(ls()[!ls() %in% c(init.list, 'i', 'init.list')],
#            ~{
# Get ga result
#   o <- get(.x)@solution
  # Variogram model discritization formula
#   if(o[1] >= 0 && o[1] < 1) {
#     m <- 'Sph'
#   } else if(o[1] >= 1 && o[1] < 2) {
#     m <- 'Exp'
#   } else if(o[1] >= 2 && o[1] <= 3) {
#     m <- 'Gau'
#   }
	# Construct variogram model
#   vgm(psill = o[2], model = m, range = o[3], nugget = o[4])
# }) %>% purrr::set_names(fnames)

# Krige segments and compare with Lucazeau (2019) predicted heat flow
# s.names <- purrr::map(names(v.mods), ~.x %>%
#                       stringr::str_replace('_vmod', '') %>%
#                       stringr::str_replace('_', ' '))

# Bounding Boxes
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
purrr::map(shp.box,
					 ~shp.hf %>%
					 rename(hf = `heat-flow (mW/m2)`) %>%
					 st_crop(.x)) %>%
purrr::set_names(nm = seg.names) -> shp.hf.crop

# Cropped Grid
purrr::map(shp.box, ~shp.grid %>% st_crop(.x)) %>%
purrr::set_names(nm = seg.names) -> shp.grid.crop

# Define variogram cutoffs for experimental variograms
cutoff <- list(5, 6, 4, 5, 3, 8, 6, 3, 5, 15, 7, 3, 5) %>%
		purrr::set_names(seg.names)
lags <- list(15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15) %>%
		purrr::set_names(seg.names)
lag.start <- list(3, 1, 1, 2, 5, 1, 3, 1, 2, 1, 1, 1, 3)

# Calculate experimental variogram
purrr::pmap(list(shp.hf.crop, cutoff, lags, lag.start, seg.names), ~{
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
purrr::map(1:13, ~vgm(model = c('Sph', 'Exp'))) %>%
purrr::set_names(nm = seg.names)-> v.mods

# Fit experimental variograms
purrr::map2(v.grms, v.mods, ~{fit.variogram(.x, model = .y)}) %>%
purrr::set_names(seg.names) -> v.fits

purrr::pmap(list(v.grms, v.fits, seg.names), ~plot_vgrm(..1, ..2, ..3)) %>%
wrap_plots()

# Kriging
purrr::pwalk(list(seg.names,
									shp.hf.crop,
									v.grms,
									v.mods,
									shp.grid.crop),
						 Krige_diff,
						 param = 'hf',
						 data.compare = shp.hf.pred)