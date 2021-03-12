# Load functions and libraries
source('functions.R')
load('data/hf.RData')

# Kriging domain
cat('Drawing Kriging domain boxes\n')
shp.box <-
	purrr::map(seg.names, ~{
		shp.buf <-
		  shp.sa.segs.robin.pacific.buffer %>%
		  filter(segment == .x)
		bbox_widen(st_bbox(shp.buf),
							 crs = st_crs(shp.buf),
							 borders = c('left' = 0.1,
													 'right' = 0.1,
													 'top' = 0.1,
													 'bottom' = 0.1)) %>%
		st_as_sf()
}) %>% purrr::set_names(nm = seg.names)

# Crop luca hf
cat('Cropping heat flow data\n')
shp.hf.crop <-
	purrr::map(shp.box, ~{
	  shp.hf %>%
		rename(hf = `heat-flow (mW/m2)`) %>%
	  st_crop(.x)
})

# Find optimal kriging parameters by genetic algorithm
purrr::pwalk(list(seg.names[[1]],
									shp.hf.crop[[1]],
									shp.box[[1]]),
						 ~{Krige_opt(seg.name = ..1,
												 data = ..2,
												 domain = ..3,
												 param = 'hf',
												 n.init = 15,
												 maxitr = 15,
												 run = 15)})