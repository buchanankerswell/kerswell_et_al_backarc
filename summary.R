# Load functions and libraries
source('functions.R')
load('data/hf.Rdata')

# Define paths and names
fpath <- list.files('data/diff', pattern = '.RData', full.names = T)
fname <- purrr::map_chr(list.files('data/diff', pattern = '.RData'),
					 ~.x %>%
					 stringr::str_replace('_diff.RData', ''))

# Load data
for (i in fpath) load(i)

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
purrr::map2_df(shp.box, seg.names,
					 ~shp.hf %>%
					 rename(hf = `heat-flow (mW/m2)`) %>%
					 st_crop(.x) %>%
					 mutate(segment = .y, .before = country)) -> shp.hf.crop

# Summarize heat flow data
hf.summary <-
				shp.hf.crop %>%
				st_set_geometry(NULL) %>%
				group_by(segment) %>%
				rename(Segment = segment) %>%
				summarise(n = n(),
									Min = round(min(hf)),
									Max = round(max(hf)),
									Median =round(median(hf)),
									IQR = round(IQR(hf)),
									Mean = round(mean(hf)),
									Sigma = round(sd(hf)))

cat('Heat flow summary:\n')
print(hf.summary)

# Visualize
p <- shp.hf.crop %>%
				group_by(segment) %>%
				ggplot() +
				geom_boxplot(aes(x = hf, y = segment, group = segment),
										 width = 0.5,
										 outlier.size = 0.5,
										 outlier.color = rgb(0.5, 0.5, 0.5, 0.1)) +
				labs(x = bquote('Surface heat flow'~(mWm^-2)),
						 y = NULL,
						 title = 'Heat flow observations by segment',
						 caption = 'data from Lucazeau (2019)') +
				scale_x_continuous(limits = c(0, 250)) +
				scale_y_discrete(limits = rev(levels(as.factor(seg.names)))) +
				theme_classic(base_size = 14) +
				theme(strip.background = element_blank())

# Save plot
cat('Saving heat flow summary plot to:\nfigs/hf_summary.png\n')

ggsave(file = 'figs/hf_summary.png',
			 plot = p,
			 device = 'png',
			 type = 'cairo',
			 width = 7,
			 height = 7)

# Summarise variogram models (add rmse's from f.obj)
purrr::map_df(fname, ~get(.x)$v.mod, .id = 'segment') %>%
	as_tibble() %>%
	mutate('segment' = fname %>% stringr::str_replace('_', ' ')) %>%
	dplyr::select(segment, model, psill, range) %>%
	rename(Segment = segment, Model = model,
				 Sill = psill, Range = range) %>%
	mutate('Sill' = round(sqrt(Sill)), 'Range' = round(Range/1000, 1)) -> variogram.summary
cat('Variogram summary:\n')
print(variogram.summary)

# Visualize
p1 <- variogram.summary %>%
		ggplot() +
		geom_bar(aes(x = Range, y = Segment), stat = 'identity') +
		labs(x = 'Range (km)', y = NULL) +
		scale_y_discrete(limits = rev(levels(as.factor(seg.names)))) +
		theme_classic(base_size = 14)

p2 <- variogram.summary %>%
		mutate(seg.length = shp.sa.segs.robin.pacific %>%
					 							st_length() %>%
												as.vector()) %>%
		ggplot() +
		geom_point(aes(x = Range, y = seg.length/1000)) +
		geom_text_repel(aes(x = Range, y = seg.length/1000, label = Segment), size = 3) +
		labs(x = 'Range (km)', y = 'Segment Length (km)') +
		theme_classic(base_size = 14)

p3 <- variogram.summary %>%
		mutate(n = hf.summary$n) %>%
		ggplot() +
		geom_point(aes(x = Range, y = n)) +
		geom_text_repel(aes(x = Range, y = n, label = Segment), size = 3) +
		labs(x = 'Range (km)', y = 'Number of observations') +
		theme_classic(base_size = 14)

purrr::map_df(shp.box, st_bbox) %>%
		mutate(area = (xmax-xmin)*(ymax-ymin)) -> shp.area

p4 <- variogram.summary %>%
		mutate(area = as.vector(shp.area$area)) %>%
		ggplot() +
		geom_point(aes(x = Range, y = area/10^6)) +
		geom_text_repel(aes(x = Range, y = area/10^6, label = Segment), size = 3) +
		labs(x = 'Range (km)', y = bquote('Domain area'~(km^2))) +
		theme_classic(base_size = 14)

# Composition
p <- p1 + p2 + p4 + p3 +
		plot_annotation(tag_levels = 'a',
										title = 'Variogram model range summary and correlations')

# Save plot
cat('Saving variogram model summary plot to:\nfigs/variogram_summary.png\n')

ggsave(file = 'figs/variogram_summary.png',
			 plot = p,
			 device = 'png',
			 type = 'cairo',
			 width = 10,
			 height = 10)

# Interpolation difference
purrr::map(fname, ~get(.x)$diff %>% st_set_geometry(NULL)) %>%
purrr::set_names(nm = fname %>% stringr::str_replace('_', ' ')) -> hf.diff

hf.diff %>%
		bind_rows(.id = 'Segment') %>%
		group_by(Segment) %>%
		summarise(Min = round(min(hf.diff)),
		Max = round(max(hf.diff)),
		Median = round(median(hf.diff)),
		IQR = round(IQR(hf.diff)),
		Mean = round(mean(hf.diff)),
		Sigma = round(sd(hf.diff))) -> hf.diff.summary

cat('Heat flow difference summary:\n')
print(hf.diff.summary)

# Visualize
p <- hf.diff %>%
				bind_rows(.id = 'segment') %>%
				group_by(segment) %>%
				ggplot() +
				geom_boxplot(aes(x = hf.diff, y = segment, group = segment),
										 width = 0.5,
										 outlier.size = 0.5,
										 outlier.color = rgb(0.5, 0.5, 0.5, 0.1)) +
				labs(x = bquote('Heat flow difference'~(mWm^-2)),
						 y = NULL,
						 title = 'Similarity vs. Kriging difference by segment') +
				scale_x_continuous(limits = c(-2*max(hf.diff.summary$IQR), 2*max(hf.diff.summary$IQR))) +
				scale_y_discrete(limits = rev(levels(as.factor(seg.names)))) +
				theme_classic(base_size = 14) +
				theme(strip.background = element_blank())

# Save plot
cat('Saving heat flow difference summary plot to:\nfigs/hf_diff_summary.png')

ggsave(file = 'figs/hf_diff_summary.png',
			 plot = p,
			 device = 'png',
			 type = 'cairo',
			 width = 7,
			 height = 7)