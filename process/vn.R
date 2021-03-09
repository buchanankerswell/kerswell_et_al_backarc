# Load functions and packages
cat('Loading packages')
source('../functions.R')

# Load data
cat('\nLoading data')
load('../data/hf.RData')

# Segment
cat('\nDefining segment')
shp.seg.vn <-
  shp.sa.segs.robin.pacific %>%
  filter(segment == 'Vanuatu')

# Contours
shp.con.vn <-
  shp.sa.countours.robin.pacific %>%
  filter(segment == 'Vanuatu')

# Buffer
shp.buf.vn <-
  shp.sa.segs.robin.pacific.buffer %>%
  filter(segment == 'Vanuatu')

# Box
shp.box.vn <-
  bbox_widen(st_bbox(shp.buf.vn),
             crs = st_crs(shp.buf.vn),
             borders = c(
               'left' = 0.1,
               'right' = 0.1,
               'top' = 0.1,
               'bottom' = 0.1)) %>%
  st_as_sf()

# Crop luca hf
cat('\nCropping heat flow data\n')
shp.hf.pred.vn <-
  shp.hf.pred %>%
  st_crop(shp.box.vn)
shp.hf.vn <-
  shp.hf %>%
  filter(`heat-flow (mW/m2)` < 250) %>%
  st_crop(shp.box.vn)

# Kriging options
lags <- 12
lag.cutoff <- 14
rotate.vngle <- 0
v.mod <- c('Sph', 'Exp', 'Gau')

# # Rotated projection
# proj.rot <- rotate_proj(shp.hf.vn, rotate.vngle)

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
  model_variogram(
  data = shp.hf.vn %>% rename(hf = `heat-flow (mW/m2)`),
  lags = lags,
  lag.cutoff = lag.cutoff,
  param = 'hf',
  grid = shp.grid %>% st_crop(shp.box.vn),
  grid.rotate = rotate.vngle,
  v.mod = v.mod,
  krige = T
)

# Variogram
cat('\nPlotting variogram ... ')
ggsave('../figs/vn/vn_variogram.png',
      plot = k$variogram.plot,
      device = 'png',
      type = 'cairo',
      width = 7,
      height = 4)

# Compare
# Visualize difference bw Lucazeau pred and krige results
# Lucazeau predicted
cat('\nVisualizing results ... ',
    '\nSaving plots to figs/vn/\n')
v.magma.scale <-
  scale_color_viridis_c(
    name = bquote(mWm^-2),
    option = 'magma',
    limits = c(0, 250))
p1 <-
shp.hf.pred.vn %>%
  ggplot() +
  geom_sf(aes(color = HF_pred), shape = 15) +

  geom_sf(data = shp.seg.vn, color = 'white', alpha = 0.6) +
  geom_sf(data = shp.con.vn, size = 0.2, color = 'white', alpha = 0.3) +
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
  geom_sf(aes(color = hf), shape = 15) +

  geom_sf(data = shp.seg.vn, color = 'white', alpha = 0.6) +
  geom_sf(data = shp.con.vn, size = 0.2, color = 'white', alpha = 0.3) +
  labs(title = 'Kriged Heat Flow') +
  v.magma.scale +
  coord_sf(expand = F) +
  theme(panel.ontop = T,
        panel.background = element_rect(fill = 'transparent', color = NA),
        panel.grid = element_line(color = rgb(1, 1, 1, 0.3), size = 0.3),
        axis.ticks = element_blank())

# Composition
p <- p1 / p2 +
  plot_annotation(title = 'Vanuatu', tag_levels = 'a')
ggsave('../figs/vn/vn_comp.png',
       plot = p,
       device = 'png',
       type = 'cairo',
       width = 7,
       height = 8)

# Difference
shp.hf.pred.vn <-
  shp.hf.pred.vn %>%
  st_crop(k$krige.results) %>%
  mutate(HF_krige = k$krige.results$hf,
         HF_diff = HF_pred - HF_krige,
         .before = geometry)
p3 <-
  shp.hf.pred.vn %>%
  ggplot() +
  geom_sf(aes(color = abs(HF_diff)), shape = 15) +

  geom_sf(data = shp.seg.vn, color = 'white', alpha = 0.6) +
  geom_sf(data = shp.con.vn, size = 0.2, color = 'white', alpha = 0.3) +
  labs(title = 'Absolute Difference') +
  v.magma.scale +
  coord_sf(expand = F) +
  theme(panel.ontop = T,
        panel.background = element_rect(fill = 'transparent', color = NA),
        panel.grid = element_line(color = rgb(1, 1, 1, 0.3), size = 0.3),
        axis.ticks = element_blank())

# Histogram
hf.diff.summary <-
  shp.hf.pred.vn %>%
  summarise(min = min(HF_diff),
            max = max(HF_diff),
            med = median(HF_diff),
            mean = mean(HF_diff),
            sd = sd(HF_diff),
            iqr = IQR(HF_diff)) %>%
  st_set_geometry(NULL)
p4 <-
  shp.hf.pred.vn %>%
  ggplot() +
  geom_histogram(aes(x = HF_diff), binwidth = 2) +
  scale_x_continuous(limits = c(median(shp.hf.pred.vn$HF_diff)-(3*IQR(shp.hf.pred.vn$HF_diff)), median(shp.hf.pred.vn$HF_diff)+(3*IQR(shp.hf.pred.vn$HF_diff)))) +
  annotate('text',
           x = (max(shp.hf.pred.vn$HF_diff) - min(shp.hf.pred.vn$HF_diff))/4,
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
  plot_annotation(title = 'Vanuatu', tag_levels = 'a')
ggsave('../figs/vn/vn_diff_comp.png',
       plot = p,
       device = 'png',
       type = 'cairo',
       width = 11,
       height = 8)

# Save results
cat('Saving results to process/vn.RData')
save(k, shp.hf.pred.vn, file = 'vn.RData')

cat('\nDone\n')
