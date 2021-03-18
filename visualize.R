# Load functions and libraries
source('functions.R')
load('data/hf.Rdata')

# Define paths and names
fpaths <- list.files('data/diff', pattern = '.RData', full.names = T)
fnames <- purrr::map_chr(list.files('data/diff', pattern = '.RData'),
					 ~.x %>%
					 stringr::str_replace('_diff.RData', ''))

# Load data
for (i in fpaths) load(i)

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

# Color scale
viridis.scale.white <- scale_color_viridis_c(option = 'magma', limits = c(0, 250), na.value = 'white')
viridis.scale.grey <- scale_color_viridis_c(option = 'magma', limits = c(0, 250), na.value = 'grey50')

# Plot sizes
p.widths <- list(11, 11, 8, 8, 8, 8, 8, 8, 8, 8, 11, 8, 8)
p.heights <- list(8, 9, 8, 12, 9, 11, 12, 8.5, 13, 9, 11, 14, 11)

# Draw plots
purrr::pwalk(list(seg.names, fnames, p.heights, p.widths), ~{
	cat('Plotting', ..1, '\n')
	shp.hf <- shp.hf.crop[[..1]]
	shp.cont <- shp.sa.countours.robin.pacific %>% filter(segment == ..1)
	shp.seg <- shp.sa.segs.robin.pacific %>% filter(segment == ..1)
	shp.diff <- get(..2)$diff
	shp.krige <- get(..2)$k
	v.grm <- get(..2)$v.grm
	v.mod <- get(..2)$v.mod
	v.line <- variogramLine(v.mod, maxdist = max(v.grm$dist))
	p.luca.pred <- ggplot() +
									geom_sf(data = shp.diff, aes(color = hf.pred.luca),
													size = 2, shape = 15) +
				    			geom_sf(data = shp.seg, size = 1.1, color = 'white') +
								  geom_sf(data = shp.cont, size = 0.15, color = 'white') +
									labs(color = bquote(mWm^-2),
											 title = 'Similarity prediction',
											 caption = 'data from Lucazeau (2019)') +
						viridis.scale.white +
				    coord_sf(expand = F,
				              xlim = c(st_bbox(shp.hf)$xmin, st_bbox(shp.hf)$xmax),
				              ylim = c(st_bbox(shp.hf)$ymin, st_bbox(shp.hf)$ymax)) +
				    theme_map(font_size = 11) +
				    theme(axis.text = element_text(),
				          panel.grid = element_line(size = 0.1, color = rgb(1,1,1,0.5)),
				          panel.ontop = T)
	p.krige.pred <- ggplot() +
										geom_sf(data = shp.diff, aes(color = hf.pred.krige),
														size = 2, shape = 15) +
					    			geom_sf(data = shp.seg, size = 1.1, color = 'white') +
									  geom_sf(data = shp.cont, size = 0.15, color = 'white') +
										labs(color = bquote(mWm^-2),
												 title = 'Krige prediction') +
							viridis.scale.white +
					    coord_sf(expand = F,
					              xlim = c(st_bbox(shp.hf)$xmin, st_bbox(shp.hf)$xmax),
					              ylim = c(st_bbox(shp.hf)$ymin, st_bbox(shp.hf)$ymax)) +
					    theme_map(font_size = 11) +
					    theme(axis.text = element_text(),
					          panel.grid = element_line(size = 0.1, color = rgb(1,1,1,0.5)),
					          panel.ontop = T)
	p.diff.pred <- ggplot() +
										geom_sf(data = shp.diff, aes(color = abs(hf.diff)),
														size = 2, shape = 15) +
					    			geom_sf(data = shp.seg, size = 1.1, color = 'white') +
									  geom_sf(data = shp.cont, size = 0.15, color = 'white') +
										labs(color = bquote(mWm^-2),
												 title = 'Absolute difference') +
							viridis.scale.white +
					    coord_sf(expand = F,
					              xlim = c(st_bbox(shp.hf)$xmin, st_bbox(shp.hf)$xmax),
					              ylim = c(st_bbox(shp.hf)$ymin, st_bbox(shp.hf)$ymax)) +
					    theme_map(font_size = 11) +
					    theme(axis.text = element_text(),
					          panel.grid = element_line(size = 0.1, color = rgb(1,1,1,0.5)),
					          panel.ontop = T)
	p.pts <- ggplot() +
						geom_sf(data = shp.world.robin.pacific, color = NA) +
				    geom_sf(data = shp.seg, size = 1.1) +
				    geom_sf(data = shp.cont, size = 0.15) +
				    geom_sf(data = shp.hf, aes(color = hf), shape = 20, size = 0.8) +
						guides(color = guide_colorbar(barheight = 10)) +
				    labs(color = bquote(mWm^-2),
								 title = bquote(bold('Observations')~n==.(nrow(shp.hf))),
								 caption = 'data from Lucazeau (2019)') +
						viridis.scale.grey +
				    coord_sf(expand = F,
				              xlim = c(st_bbox(shp.hf)$xmin, st_bbox(shp.hf)$xmax),
				              ylim = c(st_bbox(shp.hf)$ymin, st_bbox(shp.hf)$ymax)) +
				    theme_map(font_size = 11) +
				    theme(axis.text = element_text(),
				          panel.grid = element_line(size = 0.1, color = rgb(0,0,0,0.2)),
				          panel.ontop = T)
	p.hist <- ggplot() +
							geom_histogram(data = shp.diff, aes(x = hf.diff), binwidth = 2) +
							scale_x_continuous(limits = c(median(shp.diff$hf.diff)-(2*IQR(shp.diff$hf.diff)),
																						median(shp.diff$hf.diff)+(2*IQR(shp.diff$hf.diff)))) +
							labs(x = bquote(mWm^-2),
									 y = 'Frequency',
									 title = 'Prediction difference') +
							theme_classic() +
							theme(plot.title = element_text(face = 'bold'))
	p.vgrm <- ggplot() +
							geom_point(data = v.grm, aes(x = dist/1000, y = gamma), size = 0.8) +
							geom_line(data = v.line, aes(x = dist/1000, y = gamma)) +
							labs(x = 'Lag (km)', y = 'Semivariance',title = 'Variogram') +
							theme_classic() +
							theme(plot.title = element_text(face = 'bold'))
	# Composition
	p <-
	(p.luca.pred + theme(legend.position = 'none')) +
	(p.krige.pred + theme(legend.position = 'none')) +
	(p.diff.pred + theme(legend.position = 'none')) +
	p.pts +
	p.hist +
	p.vgrm +
	plot_annotation(title = ..1,
									tag_levels = 'a',
									theme = theme(legend.position = 'right',
																legend.direction = 'vertical'))
	if(..1 != 'Andes') {
		p.comp <- p + plot_layout(nrow = 3, ncol = 2, widths = 1, heights = c(1, 1, 0.5), guides = 'collect')
	} else {
		design <- '
			11223344
			55556666
		'
		p.comp <- p + plot_layout(design = design, heights = c(1, 0.5))
	}
	# Save
	cat('Saving plot to', paste0('figs/diff/', ..2, '_diff.png'), '\n')
	ggsave(file = paste0('figs/diff/', ..2, '_diff.png'),
				 device = 'png',
				 type = 'cairo',
				 plot = p.comp,
				 height = ..3,
				 width = ..4)
})
