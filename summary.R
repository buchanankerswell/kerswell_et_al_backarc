source('functions.R')
load('data/hf.Rdata')

# Define paths and names
fname <- c('Alaska Aleutians', 'Andes', 'Central America', 'Kamchatka Marianas',
					 'Kyushu Ryukyu', 'Lesser Antilles', 'New Britain Solomon', 'N. Philippines',
					 'Sumatra Banda Sea', 'Scotia', 'S. Philippines', 'Tonga New Zealand', 'Vanuatu')
fpath <- list.files('process/', pattern = '.RData', full.names = T)

# Function to read RData
load_data <- function(x){
#loads an RData file, and returns it
    load(x)
		l <- ls()[!ls() %in% x] %>%
			purrr::map(~get(.x))
		l[1:length(l)-1]
}

# Read data
d <- purrr::map2(fpath, fname, ~{
	d <- .x %>% load_data()
	d %>% set_names(c('k', 'shp.hf.pred'))
	}) %>% set_names(fname)

# Draw bounding boxes and crop heat flow for each segment
shp.hf.crop <-
	purrr::map_df(shp.sa.segs.robin.pacific$segment,
		~ {box <- shp.sa.segs.robin.pacific.buffer %>%
				filter(segment == .x) %>%
				st_bbox() %>%
				bbox_widen(proj4.robin.pacific,
									 borders = c('top' = 0.1,
															'bottom' = 0.1,
															'left' = 0.1,
															'right' = 0.1))
				shp.hf %>%
				st_crop(box) %>%
				mutate('segment' = .x, .before = country)})

# Summarize heat flow data
hf.summary <-
				shp.hf.crop %>%
				st_set_geometry(NULL) %>%
				group_by(segment) %>%
				rename(hf = `heat-flow (mW/m2)`,
							 Segment = segment) %>%
				summarise(n = n(),
									Min = round(min(hf)),
									Max = round(max(hf)),
									Median =round(median(hf)),
									IQR = round(IQR(hf)),
									Mean = round(mean(hf)),
									Sigma = round(sd(hf)))

# Visualize
p <- shp.hf.crop %>%
				group_by(segment) %>%
				ggplot() +
				geom_histogram(aes(x = `heat-flow (mW/m2)`), binwidth = 5) +
				labs(x = bquote('Heat Flow'~mWm^-2), y = 'Frequency') +
				scale_x_continuous(limits = c(0, 250)) +
				facet_wrap(~segment, scales = 'free_y', ncol = 2) +
				theme_classic() +
				theme(
					strip.background = element_blank())

# Save plot
ggsave(file = 'figs/summary/hf_summary.png',
			 plot = p,
			 device = 'png',
			 type = 'cairo',
			 width = 6,
			 height = 11)

# Summarise variogram models (add rmse's from f.obj)
variogram.summary <-
				purrr::map_df(d, ~.x$k$model.variogram, .id = 'segment') %>%
				as_tibble() %>%
				dplyr::select(segment, model, psill, range) %>%
				rename(Segment = segment, Model = model,
							 Sill = psill, Range = range)

# Visualize
p <- purrr::map2(d, names(d),
					 ~{
		cutoff <- .x$k$variogram$dist %>% max()
		v.line <- .x$k$model.variogram %>% variogramLine(cutoff)
		v <- .x$k$variogram
    p <- ggplot() +
      geom_path(data = v.line, aes(x = dist/1000, y = gamma)) +
      geom_point(data = v, aes(x = dist/1000, y = semivar), size = 0.6) +
      labs(x = NULL, y = 'Semivariance', title = .y) +
      theme_classic() +
      theme(axis.text = element_text(color = 'black'),
						axis.text.y = element_blank(),
						axis.title.y = element_blank(),
						axis.ticks.y = element_blank(),
						axis.line.y = element_blank(),
						plot.title = element_text(hjust = 0.5, size = 9))
			if(.y %in% c('Vanuatu', 'Tonga New Zealand')) {
			p <- p + labs(x = 'Lag (km)')}
			p
					 }) %>%
wrap_plots(ncol = 2)

# Save plot
ggsave(file = 'figs/summary/variogram_summary.png',
			 plot = p,
			 device = 'png',
			 type = 'cairo',
			 width = 6,
			 height = 11)

# Interpolation difference
hf.diff.summary <-
				purrr::map_df(d,
											~.x$shp.hf.pred %>%
											st_set_geometry(NULL) %>%
											summarise(Min = round(min(HF_diff)),
																Max = round(max(HF_diff)),
																Median = round(median(HF_diff)),
																IQR = round(IQR(HF_diff)),
																Mean = round(mean(HF_diff)),
																Sigma = round(sd(HF_diff))),
											.id = 'Segment')