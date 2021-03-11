# Load functions and packages
cat('Loading packages')
source('../functions.R')

# Load data
cat('\nLoading data')
load('../data/hf.RData')

# Segment
cat('\nDefining segment')
shp.seg.sb <-
  shp.sa.segs.robin.pacific %>%
  filter(segment == 'Sumatra Banda Sea')

# Contours
shp.con.sb <-
  shp.sa.countours.robin.pacific %>%
  filter(segment == 'Sumatra Banda Sea')

# Buffer
shp.buf.sb <-
  shp.sa.segs.robin.pacific.buffer %>%
  filter(segment == 'Sumatra Banda Sea')

# Box
shp.box.sb <-
  bbox_widen(st_bbox(shp.buf.sb),
             crs = st_crs(shp.buf.sb),
             borders = c(
               'left' = 0.1,
               'right' = 0.1,
               'top' = 0.1,
               'bottom' = 0.1)) %>%
  st_as_sf()

# Crop luca hf
cat('\nCropping heat flow data\n')
shp.hf.pred.sb <-
  shp.hf.pred %>%
  st_crop(shp.box.sb)
shp.hf.sb <-
  shp.hf %>%
  filter(`heat-flow (mW/m2)` < 250) %>%
  st_crop(shp.box.sb)

# Kriging options
lags <- 15
lag.cutoff <- 3
rotate.angle <- 0
v.mod <- c('Sph', 'Exp', 'Gau')

# Krige entire segment
cat('\nKriging with',
    '\nLags:',
    lags,
    '\nLag cutoff at:',
    paste0('1/',lag.cutoff),
    'max lag distance',
    '\nTrying variogram models:',
    v.mod,
    '\n')

k <-
  krige_interp(
  data = shp.hf.sb %>% rename(hf = `heat-flow (mW/m2)`),
  lags = lags,
  lag.cutoff = lag.cutoff,
  param = 'hf',
  grid = shp.grid %>% st_crop(shp.box.sb),
  grid.rotate = rotate.angle,
  v.mod = v.mod,
  krige = T,
	plot = T
)

# Variogram
cat('\nPlotting variogram ... ')
ggsave('../figs/sb/sb_variogram.png',
      plot = k$variogram.plot,
      device = 'png',
      type = 'cairo',
      width = 7,
      height = 4)

# Compare
# Visualize difference bw Lucazeau pred and krige results
# Lucazeau predicted
cat('\nVisualizing results ... ',
    '\nSaving plots to figs/sb/\n')
v.magma.scale <-
  scale_color_viridis_c(
    name = bquote(mWm^-2),
    option = 'magma',
    limits = c(0, 250))
p1 <-
shp.hf.pred.sb %>%
  ggplot() +
  geom_sf(aes(color = HF_pred), shape = 15, size = 10) +
  geom_sf(data = shp.seg.sb, color = 'white', alpha = 0.6) +
  geom_sf(data = shp.con.sb, size = 0.2, color = 'white', alpha = 0.3) +
  labs(title = 'Predicted Heat Flow', subtitle = 'from Lucazeau (2019)') +
  v.magma.scale +
  coord_sf(expand = F) +
  theme(panel.ontop = T,
        panel.background = element_rect(fill = 'transparent', color = NA),
        panel.grid = element_line(color = rgb(1, 1, 1, 0.3), size = 0.3),
        axis.ticks = element_blank())

# Visualize krige results
p2 <-
  k$krige.results %>%
  ggplot() +
  geom_sf(aes(color = hf), shape = 15, size = 10) +
  geom_sf(data = shp.seg.sb, color = 'white', alpha = 0.6) +
  geom_sf(data = shp.con.sb, size = 0.2, color = 'white', alpha = 0.3) +
  labs(title = 'Kriged Heat Flow') +
  v.magma.scale +
  coord_sf(expand = F) +
  theme(panel.ontop = T,
        panel.background = element_rect(fill = 'transparent', color = NA),
        panel.grid = element_line(color = rgb(1, 1, 1, 0.3), size = 0.3),
        axis.ticks = element_blank())

# Composition
p <- p1 / p2 +
  plot_annotation(title = 'Sumatra Banda Sea', tag_levels = 'a')
ggsave('../figs/sb/sb_comp.png',
       plot = p,
       device = 'png',
       type = 'cairo',
       width = 7,
       height = 8)

# Difference
shp.hf.pred.sb <-
  shp.hf.pred.sb %>%
  st_crop(k$krige.results) %>%
  mutate(HF_krige = k$krige.results$hf,
         HF_diff = HF_pred - HF_krige,
         .before = geometry)
p3 <-
  shp.hf.pred.sb %>%
  ggplot() +
  geom_sf(aes(color = abs(HF_diff)), shape = 15, size = 10) +
  geom_sf(data = shp.seg.sb, color = 'white', alpha = 0.6) +
  geom_sf(data = shp.con.sb, size = 0.2, color = 'white', alpha = 0.3) +
  labs(title = 'Absolute Difference') +
  v.magma.scale +
  coord_sf(expand = F) +
  theme(panel.ontop = T,
        panel.background = element_rect(fill = 'transparent', color = NA),
        panel.grid = element_line(color = rgb(1, 1, 1, 0.3), size = 0.3),
        axis.ticks = element_blank())

# Histogram
hf.diff.summary <-
  shp.hf.pred.sb %>%
  summarise(min = min(HF_diff),
            max = max(HF_diff),
            med = median(HF_diff),
            mean = mean(HF_diff),
            sd = sd(HF_diff),
            iqr = IQR(HF_diff)) %>%
  st_set_geometry(NULL)
p4 <-
  shp.hf.pred.sb %>%
  ggplot() +
  geom_histogram(aes(x = HF_diff), binwidth = 2) +
  scale_x_continuous(limits = c(median(shp.hf.pred.sb$HF_diff)-(3*IQR(shp.hf.pred.sb$HF_diff)), median(shp.hf.pred.sb$HF_diff)+(3*IQR(shp.hf.pred.sb$HF_diff)))) +
  annotate('text',
           x = (max(shp.hf.pred.sb$HF_diff) - min(shp.hf.pred.sb$HF_diff))/4,
           y = -Inf,
           label = paste0('Summary\n',
                          'median: ',
                          round(hf.diff.summary$med, 1),
                          '\niqr: ',
                          round(hf.diff.summary$iqr, 1),
                          '\nmean: ',
                          round(hf.diff.summary$mean, 1),
                          '\nsd: ',
                          round(hf.diff.summary$sd, 1)),
           hjust = 0,
           vjust = -0.5) +
  labs(x = bquote('Heat Flow Difference'~mWm^-2), y = 'Frequency') +
  theme_classic()

# Composition
p <- p1 + p2 + p3 + p4 +
  plot_layout(nrow = 2, ncol = 2) +
  plot_annotation(title = 'Sumatra Banda Sea', tag_levels = 'a')
ggsave('../figs/sb/sb_diff_comp.png',
       plot = p,
       device = 'png',
       type = 'cairo',
       width = 11,
       height = 8)

# Save results
cat('Saving results to process/sb.RData')
save(k, shp.hf.pred.sb, file = 'sb.RData')

cat('\nDone\n')
