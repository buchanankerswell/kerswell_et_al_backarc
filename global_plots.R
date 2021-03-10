source('functions.R')
cat('Loading segment and heat flow data\n')
load('data/hf.RData')

# Define bounding boxes
# Intividual segments bounding boxes
shp.segs.box <- 
  shp.sa.segs.robin.pacific.buffer$geometry %>% 
  purrr::map(~st_bbox(.x)) %>% 
  purrr::set_names(shp.sa.segs.robin.pacific$segment) %>%
	purrr::map(~.x %>% bbox_widen(proj4.robin.pacific,
											borders = c('top' = 0.1,
											'bottom' = 0.1,
											'left' = 0.1,
											'right' = 0.1)))

cat('\nPlotting ...', 'Saving figures to figs/base/')
cat('\nPlotting world map with segments ...')

# Color scale
viridis.scale.white <- scale_color_viridis_c(option = 'magma', limits = c(0, 250), na.value = 'white')
viridis.scale.grey <- scale_color_viridis_c(option = 'magma', limits = c(0, 250), na.value = 'grey50')

# World map
p <- ggplot() +
  geom_sf(data = shp.world.robin.pacific, color = NA) +
  geom_sf(data = shp.sa.segs.robin.pacific, size = 1.1) +
  geom_sf(data = shp.sa.countours.robin.pacific, size = 0.1) +
  geom_sf_label_repel(data = shp.sa.segs.robin.pacific,
                      aes(label = segment),
                      size = 3,
                      alpha = 0.8,
                      max.overlaps = 30) +
  labs(title = 'Subduction Zone Segments',
			 caption = 'from Syracuse & Abers (2006)') +
  coord_sf(expand = F) +
  theme_map(font_size = 11) +
  theme(axis.text = element_text(),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        panel.grid = element_line(size = 0.2,
																	color = rgb(0,0,0,0.2)),
				plot.margin = margin(),
        panel.ontop = T)

# Save plot
ggsave(filename = paste0('figs/base/world.png'),
       plot = p,
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.world.robin.pacific)$xmax/2.5e6 -
									 st_bbox(shp.world.robin.pacific)$xmin/2.5e6),
       height = abs(st_bbox(shp.world.robin.pacific)$ymax/2.5e6 -
										st_bbox(shp.world.robin.pacific)$ymin/2.5e6))

cat('\nPlotting world heat flow map ...')

# All Heat Flow
p.nghf <- ggplot() +
  geom_sf(data = shp.world.robin.pacific, color = NA) +
  geom_sf(data = shp.hf,
					aes(color = `heat-flow (mW/m2)`), size = 0.1) +
	viridis.scale.grey +
  labs(title = 'New Global Heat Flow',
			 caption = 'from Lucazeau (2019)',
			 color = bquote(mWm^-2)) +
  guides(color = guide_colorbar(title.vjust = 1,
																barwidth = unit(abs(st_bbox(shp.world.robin.pacific)$xmax/2.5e6 -
																										st_bbox(shp.world.robin.pacific)$xmin/2.5e6)/5, 'in'))) +
  coord_sf(expand = F) +
  theme_map() +
  theme(axis.text = element_text(),
        legend.position = 'top',
        legend.direction = 'horizontal',
        legend.justification = 'right',
				plot.title = element_text(vjust = -6),
        panel.grid = element_line(size = 0.2,
																	color = rgb(0,0,0,0.2)),
				plot.margin = margin(),
        panel.ontop = T)

# Save plot
ggsave(filename = paste0('figs/base/hf.png'),
       plot = p.nghf,
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.world.robin.pacific)$xmax/2.5e6 -
									 st_bbox(shp.world.robin.pacific)$xmin/2.5e6),
       height = abs(st_bbox(shp.world.robin.pacific)$ymax/2.5e6 -
										st_bbox(shp.world.robin.pacific)$ymin/2.5e6))

cat('\nPlotting Lucazeau (2019) predicted heat flow ...')

# Lucazeau (2019) predicted
p.pred <- shp.hf.pred %>%
  ggplot() +
  geom_sf(aes(color = HF_pred)) +
  geom_sf(data = shp.world.robin.pacific,
					color = rgb(0, 0, 0, 0.7),
					fill = NA,
					size = 0.2) +
  labs(title = 'Predicted Heat Flow',
			 caption = 'from Lucazeau (2019)',
			 subtitle = 'similarity method',
			 color = bquote(mWm^-2)) +
	viridis.scale.white +
  guides(color = guide_colorbar(title.vjust = 1,
																barwidth = unit(abs(st_bbox(shp.world.robin.pacific)$xmax/2.5e6
																										- st_bbox(shp.world.robin.pacific)$xmin/2.5e6)/5, 'in'))) +
  coord_sf(expand = F) +
  theme_map() +
  theme(axis.text = element_text(),
        legend.position = 'top',
        legend.direction = 'horizontal',
        legend.justification = 'right',
				plot.title = element_text(vjust = -12),
				plot.subtitle = element_text(vjust = -15),
				plot.margin = margin(-12, 0, 0, 0),
        panel.grid = element_line(size = 0.2, color = rgb(1,1,1,0.7)),
        panel.ontop = T)

# Save plot
ggsave(filename = paste0('figs/base/hf_pred.png'),
       device = 'png',
       plot = p.pred,
       type = 'cairo',
       width = abs(st_bbox(shp.world.robin.pacific)$xmax/2.5e6 -
									 st_bbox(shp.world.robin.pacific)$xmin/2.5e6),
       height = abs(st_bbox(shp.world.robin.pacific)$ymax/2.5e6 -
										st_bbox(shp.world.robin.pacific)$ymin/2.5e6))

cat('\nPlotting individual segments ...')

# Lucazeau composition

p <- (p.nghf + theme(plot.caption = element_blank())) /
				p.pred + plot_annotation(tag_levels = 'a')

cat('\nPlotting Lucazeau (2019) composition ...')

# Save plot
ggsave(filename = paste0('figs/base/hf_luca.png'),
       device = 'png',
       plot = p,
       type = 'cairo',
       width = abs(st_bbox(shp.world.robin.pacific)$xmax/2.5e6 -
									 st_bbox(shp.world.robin.pacific)$xmin/2.5e6),
       height = 2*abs(st_bbox(shp.world.robin.pacific)$ymax/2.5e6 -
										st_bbox(shp.world.robin.pacific)$ymin/2.5e6))

# Individual Segments

# Small segments
shp.segs.small <- shp.sa.segs.robin.pacific %>%
				filter(segment %in% c('N. Philippines', 'New Britain Solomon'))

shp.segs.box.small <-
				shp.segs.box[names(shp.segs.box) %in% c('N. Philippines', 'New Britain Solomon')]

# Med to large segments
shp.segs.large <- shp.sa.segs.robin.pacific %>%
				filter(!segment %in% c('N. Philippines', 'New Britain Solomon'))

shp.segs.box.large <-
				shp.segs.box[!names(shp.segs.box) %in% c('N. Philippines', 'New Britain Solomon')]

# Plot and save large segments
purrr::walk2(shp.segs.box.large, shp.segs.large$segment, ~{
  cat('\nPlotting', .y, ' ...')

  hf <- shp.hf %>% st_crop(.x)

  p <- ggplot() +
    geom_sf(data = shp.world.robin.pacific, color = NA) +
    geom_sf(data = shp.sa.segs.robin.pacific %>% filter(segment == .y), size = 1.1) +
    geom_sf(data = shp.sa.countours.robin.pacific %>% filter(segment == .y), size = 0.1) +
    geom_sf(data = hf,
            aes(color = `heat-flow (mW/m2)`),
            size = 1) +
    labs(color = bquote(mWm^-2),
				 title = .y,
				 subtitle = bquote(n==.(nrow(hf))),
				 caption = 'from Lucazeau (2019)') +
		viridis.scale.grey +
    guides(color = guide_colorbar(title.vjust = 1,
																	barwidth = unit(abs(st_bbox(.x)$xmax/1e6 - st_bbox(.x)$xmin/1e6)/1.5, 'in'))) +
    coord_sf(expand = F,
             xlim = c(st_bbox(.x)$xmin, st_bbox(.x)$xmax),
             ylim = c(st_bbox(.x)$ymin, st_bbox(.x)$ymax)) +
    theme_map(font_size = 11) +
    theme(axis.text = element_text(),
          legend.position = 'top',
          legend.justification = 'right',
          legend.direction = 'horizontal',
					plot.title = element_text(vjust = -12),
					plot.subtitle = element_text(vjust = -15),
					plot.margin = margin(-12, 0, 0, 0),
          panel.grid = element_line(size = 0.1, color = rgb(0,0,0,0.2)),
          panel.ontop = T)

	# Save plot
  ggsave(filename = paste0('figs/base/', .y, '.png'),
         plot = p,
         device = 'png',
         type = 'cairo',
         width = abs(st_bbox(.x)$xmax/5e5 - st_bbox(.x)$xmin/5e5),
         height = abs(st_bbox(.x)$ymax/5e5 - st_bbox(.x)$ymin/5e5))

})

# Plot and save small segments
purrr::walk2(shp.segs.box.small, shp.segs.small$segment, ~{
  cat('\nPlotting', .y, ' ...')

  hf <- shp.hf %>% st_crop(.x)

	viridis.scale.grey <- scale_color_viridis_c(option = 'magma', limits = c(0, 250))

  p <- ggplot() +
    geom_sf(data = shp.world.robin.pacific, color = NA) +
    geom_sf(data = shp.sa.segs.robin.pacific %>% filter(segment == .y), size = 1.1) +
    geom_sf(data = shp.sa.countours.robin.pacific %>% filter(segment == .y), size = 0.1) +
    geom_sf(data = hf,
            aes(color = `heat-flow (mW/m2)`),
            size = 1) +
    labs(color = bquote(mWm^-2),
				 title = .y,
				 subtitle = bquote(n==.(nrow(hf))),
				 caption = 'from Lucazeau (2019)') +
		viridis.scale.grey +
    guides(color = guide_colorbar(title.vjust = 1,
																	barwidth = unit(abs(st_bbox(.x)$xmax/1e6 -
																											st_bbox(.x)$xmin/1e6)/1.6, 'in'))) +
    coord_sf(expand = F,
             xlim = c(st_bbox(.x)$xmin, st_bbox(.x)$xmax),
             ylim = c(st_bbox(.x)$ymin, st_bbox(.x)$ymax)) +
    theme_map(font_size = 11) +
    theme(axis.text = element_text(),
          legend.position = 'top',
          legend.justification = 'right',
          legend.direction = 'horizontal',
					plot.title = element_text(vjust = -12),
					plot.subtitle = element_text(vjust = -15),
					plot.margin = margin(-12, 0, 0, 0),
          panel.grid = element_line(size = 0.1, color = rgb(0,0,0,0.2)),
          panel.ontop = T)

	# Save plot
  ggsave(filename = paste0('figs/base/', .y, '.png'),
         plot = p,
         device = 'png',
         type = 'cairo',
         width = abs(st_bbox(.x)$xmax/5e5 - st_bbox(.x)$xmin/5e5),
         height = abs(st_bbox(.x)$ymax/5e5 - st_bbox(.x)$ymin/5e5))

})

# Segments comp
cat('\nPlotting segments composition ...')

p.segs <- ggplot() +
  geom_sf(data = shp.world.robin.pacific, color = NA) +
  geom_sf(data = shp.sa.segs.robin.pacific, size = 1.1) +
  geom_sf(data = shp.sa.countours.robin.pacific, size = 0.1) +
  geom_sf_label_repel(data = shp.sa.segs.robin.pacific,
                      aes(label = segment),
                      size = 5,
                      alpha = 0.8,
                      max.overlaps = 30) +
  labs(title = 'Subduction Zone Segments',
			 caption = 'from Syracuse & Abers (2006)') +
  coord_sf(expand = F) +
  theme_map(font_size = 18) +
  theme(axis.text = element_text(),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        panel.grid = element_line(size = 0.2,
																	color = rgb(0,0,0,0.2)),
				plot.margin = margin(),
        panel.ontop = T)

# Define Vanuatu segment and bounding boxes
shp.van <- shp.sa.countours.robin.pacific %>%
				filter(segment == 'Vanuatu') %>%
				slice(1)
shp.van.buf <- shp.sa.segs.robin.pacific.buffer %>%
				filter(segment == 'Vanuatu')
shp.van.box <- shp.sa.segs.robin.pacific.buffer %>%
			 filter(segment == 'Vanuatu') %>%
			 st_bbox() %>%
			 bbox_widen(proj4.robin.pacific,
									borders = c('top' = 0,
									'bottom' = 0,
									'left' = 0,
									'right' = 0))
shp.van.box.wide <- shp.segs.box$Vanuatu
shp.van.hf <- shp.hf %>% st_crop(shp.van.box.wide)

# Segment with buffer
p.van.buf <- ggplot() +
    geom_sf(data = shp.world.robin.pacific, color = NA) +
		geom_sf(data = shp.van.buf, fill = 'cornflowerblue', alpha = 0.1) +
		geom_sf(data = shp.van.box %>% st_difference(shp.van.buf), fill = 'cornflowerblue', alpha = 0.3) +
		geom_sf(data = shp.van.box.wide %>% st_difference(shp.van.box), fill = 'cornflowerblue', alpha = 0.5) +
    geom_sf(data = shp.van, size = 1.1) +
    geom_sf(data = shp.sa.countours.robin.pacific %>% filter(segment == 'Vanuatu'), size = 0.2) +
    labs(color = bquote(mWm^-2),
				 title = 'Vanuatu',
				 subtitle = 'Interpolation Domain') +
		viridis.scale.grey +
    guides(color = guide_colorbar(title.vjust = 1,
																	barwidth = unit(abs(st_bbox(shp.van.buf)$xmax/1e6 - st_bbox(shp.van.buf)$xmin/1e6)/1.5, 'in'))) +
    coord_sf(expand = F,
             xlim = c(st_bbox(shp.van.box.wide)$xmin, st_bbox(shp.van.box.wide)$xmax),
             ylim = c(st_bbox(shp.van.box.wide)$ymin, st_bbox(shp.van.box.wide)$ymax)) +
    theme_map(font_size = 18) +
    theme(axis.text = element_text(),
          panel.grid = element_line(size = 0.1, color = rgb(0,0,0,0.2)),
					plot.margin = margin(0, 0, 20, 0),
          panel.ontop = T)

# With hf
p.van <- ggplot() +
    geom_sf(data = shp.world.robin.pacific, color = NA) +
		geom_sf(data = shp.van.buf, fill = 'transparent', size = 0.1) +
		geom_sf(data = shp.van.box %>% st_difference(shp.van.buf), size = 0.1, fill = 'transparent') +
		geom_sf(data = shp.van.box.wide %>% st_difference(shp.van.box), size = 0.1, fill = 'transparent') +
    geom_sf(data = shp.van, size = 1.1) +
    geom_sf(data = shp.sa.countours.robin.pacific %>% filter(segment == 'Vanuatu'), size = 0.2) +
    geom_sf(data = shp.van.hf,
            aes(color = `heat-flow (mW/m2)`),
            size = 1) +
    labs(color = bquote(mWm^-2),
				 title = 'Vanuatu',
				 subtitle = bquote(n==.(nrow(shp.van.hf))),
				 caption = 'from Lucazeau (2019)') +
		viridis.scale.grey +
    coord_sf(expand = F,
             xlim = c(st_bbox(shp.van.box.wide)$xmin, st_bbox(shp.van.box.wide)$xmax),
             ylim = c(st_bbox(shp.van.box.wide)$ymin, st_bbox(shp.van.box.wide)$ymax)) +
    theme_map(font_size = 18) +
    theme(axis.text = element_text(),
          panel.grid = element_line(size = 0.1, color = rgb(0,0,0,0.2)),
          panel.ontop = T)

# Plot composition
p <- p.segs / (p.van.buf + p.van) +
				plot_annotation(tag_levels = 'a', theme = theme(plot.margin = margin())) +
				plot_layout(heights = 1)

# Save plot
ggsave(filename = paste0('figs/base/segs_comp.png'),
       plot = p,
       device = 'png',
       type = 'cairo',
       width = 0.9 * abs(st_bbox(shp.world.robin.pacific)$xmax/2.5e6 -
									 st_bbox(shp.world.robin.pacific)$xmin/2.5e6),
       height = 1 * (abs(st_bbox(shp.world.robin.pacific)$ymax/2.5e6 -
										st_bbox(shp.world.robin.pacific)$ymin/2.5e6) +
								abs(st_bbox(shp.van.box.wide)$ymax/5e5 -
										st_bbox(shp.van.box.wide)$ymin/5e5)))

cat('\nDone')