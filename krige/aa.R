rm(list = ls())
suppressMessages(source('../functions.R'))
load('../data/hf.RData')

suppressWarnings({
# Compile and transform data__________________________________________________
# Alaska Aleutians

# Segment
shp.seg.aa <- shp.sa.segs.robin.pacific %>%
  filter(segment == 'Alaska Aleutians')

# Flat buffer
shp.buf.flat.aa <- shp.sa.segs.robin.pacific %>%
  filter(segment == 'Alaska Aleutians') %>%
  st_buffer(1000000, endCapStyle = 'FLAT')
  
# Round buffer
shp.buf.round.aa <- shp.sa.segs.robin.pacific %>% 
  filter(segment == 'Alaska Aleutians') %>% 
  st_buffer(100000)
  
# Contours
shp.con.aa <- shp.sa.countours.robin.pacific %>%
  filter(segment == 'Alaska Aleutians')

# Hf data
# Flat buffer
shp.hf.flat.aa <- shp.hf %>% st_join(shp.buf.flat.aa, left = F) %>%
  relocate(segment, .before = geometry)

# Round buffer
shp.hf.round.aa <- shp.hf %>% st_join(shp.buf.round.aa, left = F) %>%
  relocate(segment, .before = geometry)

# Hf data without simple feature geometry
d.hf.flat.aa <- shp.hf.flat.aa %>% st_set_geometry(NULL)
d.hf.round.aa <- shp.hf.round.aa %>% st_set_geometry(NULL)

# Determine hf position and distance relative to trench
shp.hf.flat.pos.aa <- shp.hf.aa %>%
  mutate(pos = pts_position(shp.hf.flat.aa, shp.seg.aa, 5e6, 5e6, 'updown', 'up'),
         t.dist = trench_distance(shp.hf.aa, pos, shp.seg.aa),
         .before = geometry)

# Volcanoes
shp.volc.aa <- rbind(
  shp.volc %>% st_intersection(shp.buf.aa),
  shp.volc.no.spread %>% st_intersection(shp.buf.aa))

# Split segment
shp.seg.aa.splt <- shp.seg.aa %>%
  splt(cut.prop = 18,
      buffer = F,
      buffer.dist = 1500000,
      cap.style = 'SQUARE') %>%
  mutate(segment = 'Alaska Aleutians', .before = geometry)

# Sample lines
shp.seg.aa.splt <- shp.seg.aa.splt[seq(1, nrow(shp.seg.aa.splt), 2),]

# Split buffer
shp.buf.aa.splt <- shp.seg.aa %>%
  st_zm() %>%
  splt(cut.prop = 18,
       buffer = T,
       buffer.dist = 1500000,
       cap.style = 'SQUARE') %>%
  mutate(segment = 'Alaska Aleutians', .before = geometry)

# Sample buffers
shp.buf.aa.splt <- shp.buf.aa.splt[seq(1, nrow(shp.buf.aa.splt), 2),]

# Quick look
st_coordinates(shp.hf.flat.aa) %>%
  as_tibble() %>%
  mutate(hf = shp.hf.flat.aa$hf) %>%
  ggplot() +
  geom_contour_fill(aes(x = X, y = Y, z = hf), kriging = T) +
  geom_sf(data = shp.seg.aa) +
  geom_sf(data = shp.hf.flat.aa, aes(color = hf)) +
  coord_sf(crs = proj4.robin.pacific) +
  scale_color_viridis_c(option = 'magma')

# Kriging ____________________________________________________________________
# Krige entire segment
cat('Kriging entire segment: Alaska Aleutians ...', sep = '\n\n')
k.aa <- suppressWarnings(
  krg(
    data = shp.hf.flat.pos.aa %>% select(-segment),
    lags = 10,
    lag.cutoff = 2,
    param = 'hf',
    buf = shp.buf.aa,
    seg = shp.seg.aa,
    contours = shp.con.aa,
    volc = shp.volc.aa,
    crs = proj4.robin.pacific,
    ngrid = 3e5,
    grid.method = 'regular',
    v.mod <- 'Sph',
    plot = T
  )
)

# Krige individual segments
cat('Kriging individual segments: Alaska Aleutians ...', sep = '\n\n')
k.aa.subseg <- purrr::map(unique(shp.buf.aa.splt$subseg) %>% sort(), ~{
  suppressWarnings(
    krg(
      data = shp.hf.aa %>% select(-segment),
      lags = 100,
      lag.cutoff = 2,
      param = 'hf',
      buf = shp.buf.aa.splt %>% filter(subseg == .x),
      seg = shp.seg.aa.splt %>% filter(subseg == .x),
      contours = shp.con.aa,
      volc = shp.volc.aa,
      crs = proj4.robin.pacific,
      ngrid = 3e5,
      grid.method = 'hexagonal',
      v.mod <- c('Sph', 'Nug', 'Gau', 'Mat', 'Lin', 'Cir', 'Per', 'Wav', 'Log'),
      plot = T
    )
  )
})

# Summarise and visualize data________________________________________________
# Print data
cat('All data: Alaska Aleutians', sep = '\n\n')
print(d.hf.aa %>% select(c(hf, grad, con, ele)) %>% arrange(hf),
      n = nrow(d.hf.aa))

# Summarise data
cat('Data summary: Alaska Aleutians', sep = '\n\n')
print(d.hf.aa %>%
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
d.hf.aa %>%
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
  labs(y = NULL, x = NULL, title = 'Alaska Aleutians Summary') +
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
ggsave('../figs/aa/aa-summary.png',
       device = 'png',
       dpi = 330,
       width = 7,
       height = 4.66,
       bg = 'transparent')

# Plot heat flow (map view)
shp.hf.pos.aa %>%
  ggplot() +
  geom_sf(data = bbox_widen(st_bbox(shp.buf.aa),
                            crs = proj4.robin.pacific,
                            c('left' = 0.1,
                              'right' = 0.1,
                              'top' = 0.1,
                              'bottom' = 0.1)),
          fill = 'cornflowerblue',
          alpha = 0.2,
          color = NA) +
  geom_sf(data = shp.con.aa, color = 'white', alpha = 0.25, size = 0.25) +
  geom_sf(data = shp.volc.aa, aes(shape = 'volcano'), color = 'deeppink4', alpha = 0.5) +
  geom_sf(data = shp.seg.aa, fill = NA, color = 'black', size = 1.25) +
  geom_sf(data = shp.buf.aa, fill = NA, color = 'black', size = 0.3, alpha = 0.8) +
  geom_sf(aes(color = round(hf), shape = pos)) +
  labs(title = 'Alaska Aleutians Heat Flow', color = bquote(mWm^-2), shape = NULL) +
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
ggsave('../figs/aa/aa-map.png',
       device = 'png',
       dpi = 330,
       width = 10,
       height = 10,
       bg = 'transparent')

# Plot all heat flow vs. distance to trench
shp.hf.pos.aa %>%
  ggplot() +
  geom_point(aes(x = t.dist/1000, y = hf, shape = pos), color = 'black') +
  geom_point(data = shp.volc.aa,
             aes(x = Tdis, y = 0, shape = 'volcano'), color = 'deeppink') +
  scale_shape_manual(values = c('arc-side' = 15, 'outboard' = 19, 'volcano' = 2)) +
  guides(shape = guide_legend(override.aes = list(color = 'black'))) +
  labs(x = bquote(Distance~~km),
       y = bquote(Heat~Flow~~mWm^-2),
       title = 'Alaska Aleutians Heat Flow',
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
ggsave('../figs/aa/aa-hf-trench.png',
       device = 'png',
       dpi = 330,
       width = 5,
       height = 5,
       bg = 'transparent')

# Composite summary plot
k.aa$k.plot + (k.aa$v.plot / k.aa$hist) +
  plot_layout(widths = c(1.5, 1)) +
  plot_annotation(title = 'Alaska Aleutians Kriging Summary', tag_level = 'a') &
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
ggsave('../figs/aa/aa.png',
       device = 'png',
       dpi = 330,
       width = 10,
       height = 7,
       bg = 'transparent')

# Summarise experimental variogram
cat('Number of point-pairs: Alaska Aleutians', sep = '\n\n')
print(k.aa$variogram %>%
  summarise(np.max = max(np),
            np.min = min(np),
            np.mean = round(mean(np)),
            np.median = median(np),
            .groups = 'keep'))

# Variogram model
cat('Variogram model: Alaska Aleutians', sep = '\n\n')
print(k.aa$v.model)

# Subsegment kriging

# Summarise
cat('Number of heat flow points: Alaska Aleutians', sep = '\n\n')
print(k.aa.subseg %>%
  purrr::map_df(~ .x$data) %>%
  group_by(subseg) %>%
  summarise(n = n(),
            .groups = 'keep'))

# Summarise experimental variograms
cat('Variogram point-pairs: Alaska Aleutians', sep = '\n\n')
print(purrr::map_df(k.aa.subseg,
              ~.x$variogram %>%
                summarise(np.max = max(np),
                          np.min = min(np),
                          np.mean = round(mean(np)),
                          np.median = median(np),
                          .groups = 'keep'),
              .id = 'subseg'))

# Summarise variogram models
cat('Variogram model summary: Alaska Aleutians', sep = '\n\n')
print(purrr::map_df(k.aa.subseg, ~.x$v.mod, .id = 'subseg') %>%
  mutate(segment = 'Alaska Aleutians', .before = subseg))

# Composite summary plot (krige results)
purrr::map2(k.aa.subseg, 1:length(k.aa.subseg), ~{
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
  ggsave(paste0('../figs/aa/aa-subseg-krig', .y, '.png'),
         device = 'png',
         dpi = 330,
         width = 5,
         height = 5,
         bg = 'transparent')
})

# Save plot

# Composite summary plot (variograms)
purrr::map2(k.aa.subseg, 1:length(k.aa.subseg), ~{
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
  ggsave(paste0('../figs/aa/aa-subseg-variogram.png', .y, '.png'),
         device = 'png',
         dpi = 330,
         width = 5,
         height = 5,
         bg = 'transparent')
})
})
