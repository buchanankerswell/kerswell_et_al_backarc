rm(list = ls())
suppressMessages(source('../functions.R'))
load('../data/hf.RData')

suppressWarnings({
# Compile and transform data__________________________________________________
# Lesser Antilles

# Segment
shp.seg.la <- shp.sa.segs.robin.pacific %>%
  filter(segment == 'Lesser Antilles')

# Buffer
shp.buf.la <- shp.sa.segs.robin.pacific.buffer$km.1500 %>%
  filter(segment == 'Lesser Antilles')

# Contours
shp.con.la <- shp.sa.countours.robin.pacific %>%
  filter(segment == 'Lesser Antilles')

# Hf data
shp.hf.la <- shp.hf %>% st_join(shp.buf.la, left = F) %>%
  relocate(segment, .before = geometry)

# Hf data without simple feature geometry
d.hf.la <- shp.hf.la %>% st_set_geometry(NULL)

# Determine hf position and distance relative to trench
shp.hf.pos.la <- shp.hf.la %>%
  mutate(pos = pts_position(shp.hf.la, shp.seg.la, 5e6, 5e6, 'leftright', 'left'),
         t.dist = trench_distance(shp.hf.la, pos, shp.seg.la),
         .before = geometry)

# Volcanoes
shp.volc.la <- rbind(
  shp.volc %>% st_intersection(shp.buf.la),
  shp.volc.no.spread %>% st_intersection(shp.buf.la))

# Split segment
shp.seg.la.splt <- shp.seg.la %>%
  splt(cut.prop = 20,
      buffer = F,
      buffer.dist = 1500000,
      cap.style = 'SQUARE') %>%
  mutate(segment = 'Lesser Antilles', .before = geometry)

# Sample lines
shp.seg.la.splt <- shp.seg.la.splt[seq(1, nrow(shp.seg.la.splt), 2),]

# Split buffer
shp.buf.la.splt <- shp.seg.la %>%
  st_zm() %>%
  splt(cut.prop = 20,
       buffer = T,
       buffer.dist = 1500000,
       cap.style = 'SQUARE') %>%
  mutate(segment = 'Lesser Antilles', .before = geometry)

# Sample buffers
shp.buf.la.splt <- shp.buf.la.splt[seq(1, nrow(shp.buf.la.splt), 2),]

# Kriging ____________________________________________________________________
# Krige entire segment
cat('Kriging entire segment: Lesser Antilles ...', sep = '\n\n')
k.la <- suppressWarnings(
  krg(
    data = shp.hf.pos.la %>% select(-segment),
    lags = 100,
    lag.cutoff = 2,
    param = 'hf',
    buf = shp.buf.la,
    seg = shp.seg.la,
    contours = shp.con.la,
    volc = shp.volc.la,
    crs = proj4.robin.pacific,
    ngrid = 3e1,
    grid.method = 'hexagonal',
    v.mod <- c('Sph', 'Nug', 'Gau', 'Mat', 'Lin', 'Cir', 'Per', 'Wav', 'Log'),
    plot = T
  )
)

# Krige individual segments
cat('Kriging individual segments: Lesser Antilles ...', sep = '\n\n')
k.la.subseg <- purrr::map(unique(shp.buf.la.splt$subseg) %>% sort(), ~{
  suppressWarnings(
    krg(
      data = shp.hf.la %>% select(-segment),
      lags = 100,
      lag.cutoff = 2,
      param = 'hf',
      buf = shp.buf.la.splt %>% filter(subseg == .x),
      seg = shp.seg.la.splt %>% filter(subseg == .x),
      contours = shp.con.la,
      volc = shp.volc.la,
      crs = proj4.robin.pacific,
      ngrid = 3e1,
      grid.method = 'hexagonal',
      v.mod <- c('Sph', 'Nug', 'Gau', 'Mat', 'Lin', 'Cir', 'Per', 'Wav', 'Log'),
      plot = T
    )
  )
})

# Summarise and visualize data________________________________________________
# Print data
cat('All data: Lesser Antilles', sep = '\n\n')
print(d.hf.la %>% select(c(hf, grad, con, ele)) %>% arrange(hf),
      n = nrow(d.hf.la))

# Summarise data
cat('Data summary: Lesser Antilles', sep = '\n\n')
print(d.hf.la %>%
  select(c(hf, con, ele, grad)) %>%
  pivot_longer(everything(), 'variable') %>%
  mutate(variable = factor(variable,
         levels = c('hf', 'grad', 'con', 'ele', 'tempNo'))) %>%
  group_by(variable) %>%
  summarise(min = min(value, na.rm = T),
            max = max(value, na.rm = T),
            median = median(value, na.rm = T),
            mean = mean(value, na.rm = T),
            .groups = 'keep'))

# Summary Plot
d.hf.la %>%
  select(c(hf, con, ele, grad)) %>%
  rename('Heat Flow' = hf,
         'Conductivity' = con,
         'Elevation' = ele,
         'Gradient' = grad) %>%
  pivot_longer(everything(), 'variable') %>%
  group_by(variable) %>%
  ggplot() +
  geom_boxplot(aes(x = value,
                   group = variable,
                   fill = factor(variable,
                                 levels = c('Heat Flow',
                                            'Gradient',
                                            'Conductivity',
                                            'Elevation'))),
               color = 'black',
               alpha = 0.8,
               outlier.alpha = 1,
               show.legend = F) +
  labs(y = NULL, x = NULL, title = 'Lesser Antilles Summary') +
  scale_fill_viridis_d(option = 'D') +
  facet_wrap(~factor(variable,
                     levels = c('Heat Flow',
                                'Gradient',
                                'Conductivity',
                                'Elevation')),
             ncol = 1, scales = 'free_x') +
  theme_void() +
  theme(strip.text = element_text(size = 11, margin = margin(0, 0, 8, 0)),
        axis.text.x = element_text(),
        plot.margin = margin(32, 12, 32, 12),
        panel.spacing.y = unit(14, 'pt'),
        plot.title = element_text(size = 14,
                                  margin = margin(8, 0, 8, 0),
                                  hjust = 0.5),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA))

# Save Plot
ggsave('../figs/la/la-summary.png',
       device = 'png',
       dpi = 330,
       width = 7,
       height = 4.66,
       bg = 'transparent')

# Plot heat flow (map view)
shp.hf.pos.la %>%
  ggplot() +
  geom_sf(data = bbox_widen(st_bbox(shp.buf.la),
                            crs = proj4.robin.pacific,
                            c('left' = 0.1,
                              'right' = 0.1,
                              'top' = 0.1,
                              'bottom' = 0.1)),
          fill = 'cornflowerblue',
          alpha = 0.2,
          color = NA) +
  geom_sf(data = shp.con.la, color = 'white', alpha = 0.25, size = 0.25) +
  geom_sf(data = shp.volc.la, aes(shape = 'volcano'), color = 'deeppink4', alpha = 0.5) +
  geom_sf(data = shp.seg.la, fill = NA, color = 'black', size = 1.25) +
  geom_sf(data = shp.buf.la, fill = NA, color = 'black', size = 0.3, alpha = 0.8) +
  geom_sf(aes(color = round(hf), shape = pos)) +
  labs(title = 'Lesser Antilles Heat Flow', color = bquote(mWm^-2), shape = NULL) +
  scale_shape_manual(values = c('arc-side' = 15, 'outboard' = 19, 'volcano' = 2)) +
  scale_color_viridis_c(option = 'A') +
  theme_map(font_size = 11) +
  theme(
    plot.title = element_text(face = 'plain',
                              size = 18,
                              margin = margin(8, 0, 8, 0), hjust = 0.5),
    axis.text = element_text(),
    panel.border = element_blank(),
    panel.grid = element_line(size = 0.25, color = rgb(0.1, 0.1, 0.1, 0.5)),
    panel.background = element_blank(),
    panel.ontop = TRUE,
    plot.background = element_rect(fill = "transparent", color = NA)
  )

# Save Plot
ggsave('../figs/la/la-map.png',
       device = 'png',
       dpi = 330,
       width = 10,
       height = 10,
       bg = 'transparent')

# Plot all heat flow vs. distance to trench
shp.hf.pos.la %>%
  ggplot() +
  geom_point(aes(x = t.dist/1000, y = hf, shape = pos), color = 'black') +
  geom_point(data = shp.volc.la,
             aes(x = Tdis, y = 0, shape = 'volcano'), color = 'deeppink') +
  scale_shape_manual(values = c('arc-side' = 15, 'outboard' = 19, 'volcano' = 2)) +
  guides(shape = guide_legend(override.aes = list(color = 'black'))) +
  labs(x = bquote(Distance~~km),
       y = bquote(Heat~Flow~~mWm^-2),
       title = 'Lesser Antilles Heat Flow',
       shape = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA),
        legend.position = 'bottom')

# Save Plot
ggsave('../figs/la/la-hf-trench.png',
       device = 'png',
       dpi = 330,
       width = 5,
       height = 5,
       bg = 'transparent')

# Composite summary plot
k.la$k.plot + (k.la$v.plot / k.la$hist) +
  plot_layout(widths = c(1.5, 1)) +
  plot_annotation(title = 'Lesser Antilles Kriging Summary', tag_level = 'a') &
  theme(plot.title = element_text(size = 18,
                                  hjust = 0.5,
                                  margin = margin(8, 0, 8, 0)),
        legend.position = 'top',
        legend.justification = 'center',
        panel.spacing.y = unit(0, 'pt'),
        plot.tag = element_text(size = 14, face = 'plain'),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA))

# Save plot
ggsave('../figs/la/la.png',
       device = 'png',
       dpi = 330,
       width = 10,
       height = 7,
       bg = 'transparent')

# Summarise experimental variogram
cat('Number of point-pairs: Lesser Antilles', sep = '\n\n')
print(k.la$variogram %>%
  summarise(np.max = max(np),
            np.min = min(np),
            np.mean = round(mean(np)),
            np.median = median(np),
            .groups = 'keep'))

# Variogram model
cat('Variogram model: Lesser Antilles', sep = '\n\n')
print(k.la$v.model)

# Subsegment kriging

# Summarise
cat('Number of heat flow points: Lesser Antilles', sep = '\n\n')
print(k.la.subseg %>%
  purrr::map_df(~ .x$data) %>%
  group_by(subseg) %>%
  summarise(n = n(),
            .groups = 'keep'))

# Summarise experimental variograms
cat('Variogram point-pairs: Lesser Antilles', sep = '\n\n')
print(purrr::map_df(k.la.subseg,
              ~.x$variogram %>%
                summarise(np.max = max(np),
                          np.min = min(np),
                          np.mean = round(mean(np)),
                          np.median = median(np),
                          .groups = 'keep'),
              .id = 'subseg'))

# Summarise variogram models
cat('Variogram model summary: Lesser Antilles', sep = '\n\n')
print(purrr::map_df(k.la.subseg, ~.x$v.mod, .id = 'subseg') %>%
  mutate(segment = 'Lesser Antilles', .before = subseg))

# Composite summary plot (krige results)
purrr::map2(k.la.subseg, 1:length(k.la.subseg), ~{
  .x %>%
    purrr::pluck('k.plot') %>%
    wrap_plots() +
    plot_layout(widths = 1, heights = 1) +
    plot_annotation() &
    theme(plot.title = element_text(size = 16,
                                    hjust = 0.5,
                                    margin = margin(8, 0, 8, 0)),
          legend.position = 'top',
          legend.justification = 'center',
          plot.background = element_rect(fill = 'transparent', color = NA),
          panel.background = element_rect(fill = 'transparent', color = NA))
  # Save plot
  ggsave(paste0('../figs/la/la-subseg-krig', .y, '.png'),
         device = 'png',
         dpi = 330,
         width = 5,
         height = 5,
         bg = 'transparent')
})

# Save plot

# Composite summary plot (variograms)
purrr::map2(k.la.subseg, 1:length(k.la.subseg), ~{
  .x %>%
    purrr::pluck('v.plot') %>%
    wrap_plots() +
    plot_layout(widths = 1, heights = 1) +
    plot_annotation() &
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 16,
                                    hjust = 0.5,
                                    margin = margin(8, 0, 8, 0)),
          legend.position = 'top',
          plot.background = element_rect(fill = 'transparent', color = NA),
          panel.background = element_rect(fill = 'transparent', color = NA),
          legend.background = element_rect(fill = 'transparent', color = NA))
  # Save plot
  ggsave(paste0('../figs/la/la-subseg-variogram.png', .y, '.png'),
         device = 'png',
         dpi = 330,
         width = 5,
         height = 5,
         bg = 'transparent')
})
})
