# Load functions and libraries
source('functions.R')
load('data/hf.RData')

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
					 st_crop(.x)) -> shp.hf.crop

# Find optimal kriging parameters by genetic algorithm
purrr::walk(shp.hf.crop,
						Krige_opt,
						param = 'hf',
						n.init = 50,
						maxitr = 250,
						run = 50)
