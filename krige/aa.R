suppressMessages(source('../functions.R'))
load('../data/hf.RData')

# Compile and transform data
cat('Compiling and transforming spatial data for Alaska Aleutians:
    draw segment
    draw buffer
    draw contours
    crop heat flow
    draw bounding box
    calculate distances from data points to trench
    get volcano positions
    \n')
# Alaska Aleutians

# Kriging options
lags <- 15
lag.cutoff <- 3
ngrid <- 1e3
grid.method <- 'hexagonal'
rotate.angle <- 10
v.mod <- c('Sph', 'Gau', 'Exp', 'Lin')
idp <- 2

cat('Kriging Alaska Aleutians heat flow...\nKriging options:',
    'lags:', lags,
    'lag cutoff (proportion of max distance):', paste0('1/', lag.cutoff),
    'grid size:', ngrid,
    'sampling method:', grid.method,
    'variogram model:', v.mod,
    'inverse weighting power:', idp,
    sep = '\n')

# Segment
shp.seg.aa <- shp.sa.segs.robin.pacific %>%
  filter(segment == 'Alaska Aleutians')

# Flat buffer
shp.buf.aa <- shp.sa.segs.robin.pacific %>%
  filter(segment == 'Alaska Aleutians') %>%
  st_buffer(1000000, endCapStyle = 'FLAT')

# Contours
shp.con.aa <- shp.sa.countours.robin.pacific %>%
  filter(segment == 'Alaska Aleutians')

# Heat flow
shp.hf.aa <- shp.hf %>% st_join(shp.buf.aa, left = F) %>%
  relocate(segment, .before = geometry)

# Volcanoes
suppressWarnings(
  shp.volc.aa <- rbind(
    shp.volc %>% st_intersection(shp.buf.aa),
    shp.volc.no.spread %>% st_intersection(shp.buf.aa)))

# Rotated projection
proj.rot <- rotate_proj(shp.hf.aa, rotate.angle)

# Viridis pallettes
v.color.pal <- scale_color_viridis_c(option = 'magma',
                                     limits = c(min(shp.hf.aa$hf), max(shp.hf.aa$hf)))
v.fill.pal <- scale_fill_viridis_c(option = 'magma',
                                   limits = c(min(shp.hf.aa$hf), max(shp.hf.aa$hf)))

# Kriging
# Krige entire segment
k <- model_variogram(
  data = shp.hf.aa,
  lags = lags,
  lag.cutoff = lag.cutoff,
  param = 'hf',
  ngrid = ngrid,
  grid.method = grid.method,
  grid.rotate = rotate.angle,
  v.mod = vgm(model = v.mod),
  krige = T,
  interp.power = 2
)
cat('\nSaving variogram to figs/aa/variogram.png\n')
ggsave(filename = '../figs/aa/variogram.png',
       plot = k$variogram.plot + labs(title = 'Alaska Aleutians Variogram'),
       device = 'png',
       type = 'cairo',
       width = 7,
       height = 4,
       bg = 'transparent')
cat('Drawing interpolation sampling grid\n')
p1 <- k$interp.results %>% 
  ggplot() +
  geom_sf(data = shp.seg.aa) +
  geom_sf(data = shp.hf.aa, aes(color = hf), size = 0.3) +
  geom_sf(data = shp.con.aa, color = 'black', alpha = 0.5) +
  geom_sf(data = shp.seg.aa, color = 'black') +
  geom_sf(data = shp.volc.aa, color = 'black', shape = 2) +
  geom_sf(size = 0.2, shape = 15) +
  labs(color = bquote(mWm^-2),
       shape = 'Position',
       title = 'Alaska Aleutians Interpolation Grid') +
  coord_sf(crs = proj.rot) +
  v.color.pal +
  scale_shape_manual(values = c('arc-side' = 17, 'outboard' = 15)) +
  guides(color = guide_colorbar(title.vjust = 1,
                                barwidth = unit(abs(st_bbox(shp.hf.aa)$xmax/5.5e5 -
                                                      st_bbox(shp.hf.aa)$xmin/5.5e5)/5,
                                                'in'))) +
  theme_map(font_size = 11) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color = rgb(0.5, 0.5, 0.5, 0.5)),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA),
        panel.ontop = T)
cat('Saving heat flow grid to figs/aa/hf_grid.png\n')
ggsave(filename = paste0('../figs/aa/hf_grid.png'),
       plot = p1,
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.hf.aa)$xmax/5.5e5 - st_bbox(shp.hf.aa)$xmin/5.5e5),
       height = abs(st_bbox(shp.hf.aa)$ymax/5.5e5 - st_bbox(shp.hf.aa)$ymin/5.5e5),
       bg = 'transparent')
# Draw interpolated results
cat('Drawing interpolated results\n')
p2 <- k$interp.results %>%
  st_coordinates() %>% 
  as_tibble() %>% 
  mutate(hf = k$interp.results$hf) %>%
  ggplot() +
  geom_contour_fill(aes(x = X, y = Y, z = hf)) +
  geom_sf(data = shp.hf.aa %>% st_transform(proj.rot), aes(color = hf), size = 0.3, show.legend = F) +
  geom_sf(data = shp.con.aa %>% st_transform(proj.rot), color = 'white', alpha = 0.5) +
  geom_sf(data = shp.seg.aa %>% st_transform(proj.rot), color = 'white') +
  geom_sf(data = shp.volc.aa %>% st_transform(proj.rot), color = 'white', shape = 2) +
  labs(fill = bquote(mWm^-2), title = 'Alaska Aleutians Interpolated Heat Flow') +
  v.color.pal +
  v.fill.pal +
  guides(fill = guide_colorbar(title.vjust = 1,
                               barwidth = unit(abs(st_bbox(shp.hf.aa)$xmax/5.5e5 -
                                                     st_bbox(shp.hf.aa)$xmin/5.5e5)/5,
                                               'in'))) +
  coord_sf(crs = proj.rot) +
  theme_map(font_size = 11) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color = rgb(0.5, 0.5, 0.5, 0.5)),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA),
        panel.ontop = T)
cat('Saving interpolated heat flow results to figs/aa/hf_interp.png\n')
ggsave(filename = paste0('../figs/aa/hf_interp.png'),
       plot = p2,
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.hf.aa)$xmax/5.5e5 - st_bbox(shp.hf.aa)$xmin/5.5e5),
       height = abs(st_bbox(shp.hf.aa)$ymax/5.5e5 - st_bbox(shp.hf.aa)$ymin/5.5e5),
       bg = 'transparent')
# Heat flow transects
cat('Drawing heat flow transects\n')
t <- k$interp.results %>%
  st_coordinates() %>%
  as_tibble() %>%
  mutate(hf = k$interp.results$hf) %>%
  group_by(X) %>%
  mutate(group = cur_group_id()) %>% 
  st_as_sf(coords = c(1,2), crs = proj.rot) %>% 
  st_transform(crs = proj4.wgs)
p3 <- t %>% ggplot() +
  geom_sf(aes(color = group), size = 0.2, shape = 15) +
  geom_sf(data = shp.hf.aa %>% st_transform(proj.rot), size = 0.3, show.legend = F) +
  geom_sf(data = shp.con.aa %>% st_transform(proj.rot), color = 'black', alpha = 0.5) +
  geom_sf(data = shp.seg.aa %>% st_transform(proj.rot), color = 'black') +
  geom_sf(data = shp.volc.aa %>% st_transform(proj.rot), color = 'black', shape = 2) +
  coord_sf(crs = proj.rot) +
  labs(color = 'Transect', title = 'Alaska Aleutians Transects') +
  scale_color_viridis_c(option = 'viridis') +
  guides(color = guide_colorbar(title.vjust = 1,
                                barwidth = unit(abs(st_bbox(shp.hf.aa)$xmax/5.5e5 -
                                                      st_bbox(shp.hf.aa)$xmin/5.5e5)/5,
                                                'in'))) +
  theme_map(font_size = 11) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color = rgb(0.5, 0.5, 0.5, 0.5)),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA),
        panel.ontop = T)
cat('Saving heat flow transect map to figs/aa/hf_trans_map.png\n')
ggsave(filename = paste0('../figs/aa/hf_trans_map.png'),
       plot = p3,
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.hf.aa)$xmax/5.5e5 - st_bbox(shp.hf.aa)$xmin/5.5e5),
       height = abs(st_bbox(shp.hf.aa)$ymax/5.5e5 - st_bbox(shp.hf.aa)$ymin/5.5e5),
       bg = 'transparent')
p4 <- t %>% 
  st_coordinates() %>% 
  as_tibble() %>% 
  mutate(hf = t$hf,
         group = t$group) %>% 
  group_by(group) %>% 
  ggplot() +
  geom_path(aes(x = Y, y = hf, group = group, color = group), size = 0.3, show.legend = F) +
  scale_color_viridis_c(option = 'viridis') +
  labs(x = 'Latitude',
       y = bquote(Heat~Flow~mWm^-2),
       color = 'Transect',
       title = 'Alaska Aleutians Transects') +
  theme_classic(base_size = 11) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color = rgb(0.5, 0.5, 0.5, 0.5)),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA),
        title = element_text(face = 'bold'))
cat('Saving heat flow transects to figs/aa/hf_trans.png\n')
ggsave(filename = paste0('../figs/aa/hf_trans.png'),
       plot = p4,
       device = 'png',
       type = 'cairo',
       width = 7,
       height = 4,
       bg = 'transparent')
# Composition
cat('Drawing composition\n')
p <- (p2 + p3) / ((k$variogram.plot +
                       labs(title = 'Alaska Aleutians Variogram') +
                       theme(title = element_text(face = 'bold'))) + p4) +
  plot_annotation(tag_levels = 'a',
                  theme = theme(plot.margin = margin(),
                                panel.background = element_rect(fill = 'transparent', color = NA),
                                plot.background = element_rect(fill = 'transparent', color = NA),
                                legend.background = element_rect(fill = 'transparent', color = NA))) +
  plot_layout(guides = 'collect',
              heights = c(2, 1)) &
  theme(legend.position = 'bottom',
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA))
cat('Saving composition to figs/aa/comp.png\n')
ggsave(filename = paste0('../figs/aa/comp.png'),
       plot = p,
       device = 'png',
       type = 'cairo',
       width = 13,
       height = 8,
       bg = 'transparent')