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
										 outlier.size = 0.2,
										 outlier.color = rgb(0.5, 0.5, 0.5, 0.1)) +
				labs(x = bquote(~mWm^-2),
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
										 outlier.size = 0.2,
										 outlier.color = rgb(0.5, 0.5, 0.5, 0.1)) +
				labs(x = bquote(mWm^-2),
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