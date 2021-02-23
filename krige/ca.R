suppressMessages(source('../functions.R'))
load('../data/hf.RData')

# Compile and transform data
cat('Compiling and transforming spatial data for Central America:
    draw segment
    draw buffer
    draw contours
    crop heat flow
    draw bounding box
    calculate distances from data points to trench
    get volcano positions
    \n')
# Central America

# Kriging options
lags <- 12
lag.cutoff <- 2
ngrid <- 1e3
grid.method <- 'hexagonal'
rotate.angle <- -50
v.mod <- c('Sph', 'Exp', 'Gau', 'Lin')
idp <- 2

cat('Kriging Central America heat flow...\nKriging options:',
    'lags:', lags,
    'lag cutoff (proportion of max distance):', paste0('1/', lag.cutoff),
    'grid size:', ngrid,
    'sampling method:', grid.method,
    'variogram model:', v.mod,
    'inverse weighting power:', idp,
    sep = '\n')

# Segment
shp.seg.ca <- shp.sa.segs.robin.pacific %>%
  filter(segment == 'Central America')

# Flat buffer
shp.buf.ca <- shp.sa.segs.robin.pacific %>%
  filter(segment == 'Central America') %>%
  st_buffer(1000000, endCapStyle = 'FLAT')

# Contours
shp.con.ca <- shp.sa.countours.robin.pacific %>%
  filter(segment == 'Central America')

# Heat flow
shp.hf.ca <- shp.hf %>% st_join(shp.buf.ca, left = F) %>%
  relocate(segment, .before = geometry)

# Volcanoes
suppressWarnings(
  shp.volc.ca <- rbind(
    shp.volc %>% st_intersection(shp.buf.ca),
    shp.volc.no.spread %>% st_intersection(shp.buf.ca)))

# Rotated projection
proj.rot <- rotate_proj(shp.hf.ca, rotate.angle)

# Viridis pallettes
v.color.pal <- scale_color_viridis_c(option = 'magma',
                                     limits = c(min(shp.hf.ca$hf), max(shp.hf.ca$hf)))
v.fill.pal <- scale_fill_viridis_c(option = 'magma',
                                   limits = c(min(shp.hf.ca$hf), max(shp.hf.ca$hf)))

# Kriging
# Krige entire segment
k <- model_variogram(
  data = shp.hf.ca,
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
cat('\nSaving variogram to figs/ca/variogram.png\n')
ggsave(filename = '../figs/ca/variogram.png',
       plot = k$variogram.plot + labs(title = 'Central America Variogram'),
       device = 'png',
       type = 'cairo',
       width = 7,
       height = 4,
       bg = 'transparent')
cat('Drawing interpolation sampling grid\n')
p1 <- k$interp.results %>% 
  ggplot() +
  geom_sf(data = shp.seg.ca) +
  geom_sf(data = shp.hf.ca, aes(color = hf), size = 0.3) +
  geom_sf(data = shp.con.ca, color = 'black', alpha = 0.5) +
  geom_sf(data = shp.seg.ca, color = 'black') +
  geom_sf(data = shp.volc.ca, color = 'black', shape = 2) +
  geom_sf(size = 0.2, shape = 15) +
  labs(color = bquote(mWm^-2),
       shape = 'Position',
       title = 'Central America Interpolation Grid') +
  coord_sf(crs = proj.rot) +
  scale_color_viridis_c(option = 'magma') +
  guides(color = guide_colorbar(title.vjust = 1,
                                barwidth = unit(abs(st_bbox(shp.hf.ca)$xmax/5.5e5 -
                                                      st_bbox(shp.hf.ca)$xmin/5.5e5)/3,
                                                'in'))) +
  theme_map(font_size = 11) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color = rgb(0.5, 0.5, 0.5, 0.5)),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA),
        panel.ontop = T)
cat('Saving heat flow grid to figs/ca/hf_grid.png\n')
ggsave(filename = paste0('../figs/ca/hf_grid.png'),
       plot = p1,
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.hf.ca)$xmax/5.5e5 - st_bbox(shp.hf.ca)$xmin/5.5e5),
       height = abs(st_bbox(shp.hf.ca)$ymax/5.5e5 - st_bbox(shp.hf.ca)$ymin/5.5e5),
       bg = 'transparent')
# Draw interpolated results
cat('Drawing interpolated results\n')
p2 <- k$interp.results %>%
  st_coordinates() %>% 
  as_tibble() %>% 
  mutate(hf = k$interp.results$hf) %>%
  ggplot() +
  geom_contour_fill(aes(x = X, y = Y, z = hf)) +
  geom_sf(data = shp.hf.ca %>% st_transform(proj.rot), aes(color = hf), size = 0.3, show.legend = F) +
  geom_sf(data = shp.con.ca %>% st_transform(proj.rot), color = 'white', alpha = 0.5) +
  geom_sf(data = shp.seg.ca %>% st_transform(proj.rot), color = 'white') +
  geom_sf(data = shp.volc.ca %>% st_transform(proj.rot), color = 'white', shape = 2) +
  labs(fill = bquote(mWm^-2), title = 'Central America Interpolated Heat Flow') +
  scale_color_viridis_c(option = 'magma') +
  scale_fill_viridis_c(option = 'magma') +
  guides(fill = guide_colorbar(title.vjust = 1,
                               barwidth = unit(abs(st_bbox(shp.hf.ca)$xmax/5.5e5 -
                                                     st_bbox(shp.hf.ca)$xmin/5.5e5)/3,
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
cat('Saving interpolated heat flow results to figs/ca/hf_interp.png\n')
ggsave(filename = paste0('../figs/ca/hf_interp.png'),
       plot = p2,
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.hf.ca)$xmax/5.5e5 - st_bbox(shp.hf.ca)$xmin/5.5e5),
       height = abs(st_bbox(shp.hf.ca)$ymax/5.5e5 - st_bbox(shp.hf.ca)$ymin/5.5e5),
       bg = 'transparent')
# Heat flow transects
cat('Drawing heat flow transects\n')
t <- k$interp.results %>%
  st_coordinates() %>%
  as_tibble() %>%
  mutate(hf = k$interp.results$hf) %>%
  group_by(Y) %>%
  mutate(group = cur_group_id()) %>% 
  st_as_sf(coords = c(1,2), crs = proj.rot) %>% 
  st_transform(crs = proj4.wgs)
p3 <- t %>% ggplot() +
  geom_sf(aes(color = group), size = 0.2, shape = 15) +
  geom_sf(data = shp.hf.ca %>% st_transform(proj.rot), size = 0.3, show.legend = F) +
  geom_sf(data = shp.con.ca %>% st_transform(proj.rot), color = 'black', alpha = 0.5) +
  geom_sf(data = shp.seg.ca %>% st_transform(proj.rot), color = 'black') +
  geom_sf(data = shp.volc.ca %>% st_transform(proj.rot), color = 'black', shape = 2) +
  coord_sf(crs = proj.rot) +
  labs(color = 'Transect', title = 'Central America Transects') +
  scale_color_viridis_c(option = 'viridis') +
  guides(color = guide_colorbar(title.vjust = 1,
                                barwidth = unit(abs(st_bbox(shp.hf.ca)$xmax/5.5e5 -
                                                      st_bbox(shp.hf.ca)$xmin/5.5e5)/3,
                                                'in'))) +
  theme_map(font_size = 11) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color = rgb(0.5, 0.5, 0.5, 0.5)),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA),
        panel.ontop = T)
cat('Saving heat flow transect map to figs/ca/hf_trans_map.png\n')
ggsave(filename = paste0('../figs/ca/hf_trans_map.png'),
       plot = p3,
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.hf.ca)$xmax/5.5e5 - st_bbox(shp.hf.ca)$xmin/5.5e5),
       height = abs(st_bbox(shp.hf.ca)$ymax/5.5e5 - st_bbox(shp.hf.ca)$ymin/5.5e5),
       bg = 'transparent')
p4 <- t %>% 
  st_coordinates() %>% 
  as_tibble() %>% 
  mutate(hf = t$hf,
         group = t$group) %>% 
  group_by(group) %>% 
  ggplot() +
  geom_path(aes(x = X, y = hf, group = group, color = group), size = 0.3, show.legend = F) +
  scale_color_viridis_c(option = 'viridis') +
  labs(x = 'Longitude',
       y = bquote(Heat~Flow~mWm^-2),
       color = 'Transect',
       title = 'Central America Transects') +
  theme_classic(base_size = 11) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color = rgb(0.5, 0.5, 0.5, 0.5)),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA),
        title = element_text(face = 'bold'))
cat('Saving heat flow transects to figs/ca/hf_trans.png\n')
ggsave(filename = paste0('../figs/ca/hf_trans.png'),
       plot = p4,
       device = 'png',
       type = 'cairo',
       width = 7,
       height = 4,
       bg = 'transparent')
# Composition
cat('Drawing composition\n')
p <- (p2 + p3) / ((k$variogram.plot +
                       labs(title = 'Central America Variogram') +
                       theme(title = element_text(face = 'bold'))) + p4) +
  plot_annotation(tag_levels = 'a',
                  theme = theme(plot.margin = margin(),
                                panel.background = element_rect(fill = 'transparent', color = NA),
                                plot.background = element_rect(fill = 'transparent', color = NA),
                                legend.background = element_rect(fill = 'transparent', color = NA))) +
  plot_layout(guides = 'collect',
              heights = c(3, 1)) &
  theme(legend.position = 'bottom',
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA))
cat('Saving composition to figs/ca/comp.png\n')
ggsave(filename = paste0('../figs/ca/comp.png'),
       plot = p,
       device = 'png',
       type = 'cairo',
       width = 11,
       height = 9,
       bg = 'transparent')

